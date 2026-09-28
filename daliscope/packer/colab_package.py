import sqlite3
import numpy as np
import os
import pandas as pd
from collections import defaultdict
import sys

class ProteinDataRepo:
    def __init__(self, db_path: str, shard_dir: str):
        self.db_path   = db_path
        self.shard_dir = shard_dir
        self.cache_shard = None
        self.cache_name  = None

    def fetch_batch(self, dali_ids: list, ns_id: int) -> dict:
        if not dali_ids: return {}

        if not os.path.isfile(self.db_path):
            raise FileNotFoundError(f"Protein ledger not found: {self.db_path}")
        from pathlib import Path
        conn = sqlite3.connect(Path(self.db_path).resolve().as_uri() + "?mode=ro", uri=True)
        cur = conn.cursor()

        # 1. Create a temporary table to hold the IDs we want to fetch
        cur.execute("CREATE TEMPORARY TABLE target_ids (d_id TEXT)")

        # 2. Bulk insert the IDs into the temporary table
        # We use a list of tuples for executemany
        id_data = [(str(d),) for d in dali_ids]
        cur.executemany("INSERT INTO target_ids VALUES (?)", id_data)

        # 3. Perform a JOIN to fetch all metadata in one shot
        # This bypasses the SQL variable limit entirely
        query = """
            SELECT p.dali_id, p.n_res, p.umid, p.shard_path, p.cx, p.cy, p.cz 
            FROM proteins p
            JOIN target_ids t ON p.dali_id = t.d_id
            WHERE p.ns_id = ?
        """

        rows = cur.execute(query, (ns_id,)).fetchall()
        for row in rows[:10]: print(row)
        conn.close() # Temporary table is automatically dropped on close

        # --- Remaining logic is the same ---
        shard_map = defaultdict(list)
        for d_id, n_res, umid, s_name, cx, cy, cz in rows:
            shard_map[s_name].append((d_id, umid, cx, cy, cz))

        results = {}
        for s_name, entries in shard_map.items():
            path = os.path.join(self.shard_dir, s_name)
            if not os.path.exists(path):
                print(f"WARNING: Shard file missing: {path}. Skipping {len(entries)} proteins.")
                continue

            try:
                # mmap_mode='r' keeps memory usage low for huge monoliths
                shard = np.load(path, mmap_mode='r')
            except Exception as e:
                print(f"ERROR: Failed to load shard {s_name}: {e}")
                continue

            for d_id, umid, cx, cy, cz in entries:
                try:
                    results[d_id] = {
                        "umid": umid,
                        "coords": shard[f"{umid}_coords"],
                        "plddt": shard[f"{umid}_plddt"],
                        "seq": shard[f"{umid}_seq"],
                        "dssp": shard[f"{umid}_dssp"],
                        "compnd": shard[f"{umid}_compnd"],
                        "centroid": np.array([cx, cy, cz])
                        #"cx": cx,
                        #"cy": cy,
                        #"cz": cz
                    }
                except KeyError as e:
                    print(f"WARNING: UMID {umid} not found in shard {s_name}")

        return results

def get_umid(ns_id: int, dali_id: str) -> str:
    """Standardized UMID construction."""
    return f"{ns_id}.{dali_id}"

def create_colab_package(target_dali_ids: list,
                          ns_id:           int,
                          db_path:         str,
                          shard_dir:       str,
                          output_name:     str = "colab_pack") -> None:

    repo = ProteinDataRepo(db_path, shard_dir)
    # Using the temporary table JOIN method discussed previously to handle 10k+ IDs
    all_data = repo.fetch_batch(list(target_dali_ids), ns_id)

    if len(all_data) < len(target_dali_ids):
        missing = len(target_dali_ids) - len(all_data)
        print(f"PROCEEDING WITH CAUTION: {missing} proteins could not be retrieved.")

    # Monolithic buffers
    all_coords_list = []
    all_plddt_list = []
    all_seq_list = []
    all_dssp_list = []
    pack_ledger = []

    current_offset = 0
    n_ok = 0

    for d_id in target_dali_ids:
        data = all_data.get(d_id)
        if data is None:
            continue

        umid = get_umid(ns_id, d_id)

        # 1. Coordinates & Ground Truth length
        # Standardizing to float16 and uint8 saves massive amounts of space in the NPZ
        coords = data['coords'].astype(np.float16)
        n_res = coords.shape[0]
        plddt = data['plddt'].astype(np.uint8)[:n_res]
        centroid = data['centroid'].astype(np.float16)
        #cx = data['cx']
        #cy = data['cy']
        #cz = data['cz']

        # 2. Unified Sequence & DSSP Handling
        # We handle strings, bytes, or arrays and force-clip to n_res to ensure alignment
        def to_uint8_arr(raw_data, length):
            if isinstance(raw_data, (bytes, str)):
                b = raw_data.encode() if isinstance(raw_data, str) else raw_data
                return np.frombuffer(b, dtype=np.uint8)[:length]
            return np.array(raw_data).view(np.uint8)[:length]

        seq_arr = to_uint8_arr(data['seq'], n_res)
        dssp_arr = to_uint8_arr(data['dssp'], n_res)

        # Safety Check: Total Alignment (Should never fail with the clipping above)
        if not (n_res == len(seq_arr) == len(dssp_arr)):
            print(f"CRITICAL ERROR: {umid} length mismatch! Skipping.")
            continue

        # 3. Robust COMPND string reconstruction
        compnd_raw = data['compnd']
        if isinstance(compnd_raw, np.ndarray):
            if compnd_raw.dtype.kind in 'iu': # Numeric ASCII
                compnd_str = "".join(map(chr, compnd_raw)).strip()
            elif compnd_raw.ndim == 0: # 0-d object
                val = compnd_raw.item()
                compnd_str = val.decode('utf-8') if isinstance(val, bytes) else str(val)
            else:
                compnd_str = str(compnd_raw)
        else:
            compnd_str = str(compnd_raw)

        # 4. Accumulate
        all_coords_list.append(coords)
        all_plddt_list.append(plddt)
        all_seq_list.append(seq_arr)
        all_dssp_list.append(dssp_arr)

        pack_ledger.append({
            'umid': umid,
            'dali_id': d_id,
            'length': n_res,
            'start_idx': current_offset,
            'end_idx': current_offset + n_res,
            'compnd': compnd_str,
            'centroid': centroid,
            #'cx': cx,
            #'cy': cy,
            #'cz': cz
        })

        current_offset += n_res
        n_ok += 1

    if n_ok == 0:
        print("Error: No proteins were successfully processed.")
        return

    # 5. Final Concatenation
    full_coords = np.concatenate(all_coords_list, axis=0)
    full_plddt  = np.concatenate(all_plddt_list, axis=0)
    full_seq    = np.concatenate(all_seq_list, axis=0)
    full_dssp   = np.concatenate(all_dssp_list, axis=0)

    # 6. Save to NPZ (Compression is vital for Colab transfers)
    np.savez_compressed(f"{output_name}.npz",
                        coords=full_coords,
                        plddt=full_plddt,
                        seq=full_seq,
                        dssp=full_dssp)

    # 7. Save the Ledger
    pd.DataFrame(pack_ledger).to_csv(f"{output_name}_index.tsv", sep='\t', index=False)

    print(f"Successfully packed {n_ok} proteins ({full_coords.shape[0]} total residues).")
    print(f"Outputs: {output_name}.npz, {output_name}_index.tsv")

def main(argv=None):
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description="Build coordinate/sequence arrays from a protein ledger")
    parser.add_argument("namespace", type=int, choices=(1, 2, 3, 4))
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--shards", type=Path, required=True)
    parser.add_argument("--id-mapping", type=Path, default=Path("id_mapping.tsv"))
    parser.add_argument("--output", default="colab_pack")
    args = parser.parse_args(argv)
    df = pd.read_csv(args.id_mapping, sep="\t")
    create_colab_package(df["target_id"].drop_duplicates().tolist(), args.namespace,
                         str(args.database), str(args.shards), args.output)
    if not Path(f"{args.output}.npz").is_file():
        raise RuntimeError("No proteins were successfully packed")


if __name__ == "__main__":
    main()
