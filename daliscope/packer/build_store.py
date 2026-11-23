"""Build a local DaliScope structure store for the targets in one DALI search."""

import argparse
from contextlib import closing
import csv
import gzip
from pathlib import Path
import sqlite3
import tempfile

import gemmi
import numpy as np

from daliscope.packer.make_colab_pack import parse_namespace
from daliscope.packer.minimal_parser import read_alignments
from daliscope.packer.query_meta import extract_quoted


def read_dat(path):
    """Read the first DaliLite coordinate block and its sequence annotations."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().split()
        if len(header) < 4 or header[0] != ">>>>":
            raise ValueError("invalid DALI .dat header")
        length, n_sse = map(int, header[2:4])
        if length < 1 or n_sse < 0:
            raise ValueError("invalid DALI residue or secondary-structure count")
        for _ in range(n_sse):
            handle.readline()
        numbers = []
        while len(numbers) < 3 * length:
            line = handle.readline().rstrip("\r\n")
            if not line.strip():
                raise ValueError("truncated DALI coordinate block")
            numbers.extend(float(line[i:i + 8]) for i in range(0, len(line), 8)
                           if line[i:i + 8].strip())
        coords = np.asarray(numbers, dtype=np.float32).reshape(length, 3)
        annotations = {}
        for line in handle:
            parts = line.split(maxsplit=1)
            if len(parts) == 2 and parts[0] in ("-sequence", "-dssp", "-compnd"):
                annotations[parts[0]] = extract_quoted(line.strip())
    seq = annotations.get("-sequence", "").upper()
    dssp = annotations.get("-dssp", "")
    if len(seq) != length or len(dssp) != length or not np.isfinite(coords).all():
        raise ValueError("invalid coordinates or sequence/DSSP length in DALI .dat")
    return coords, seq, dssp, annotations.get("-compnd", "")


def model_atoms(path, chain_id):
    structure = gemmi.read_structure(str(path))
    if not len(structure):
        raise ValueError(f"no models in {path}")
    chain = structure[0].find_chain(chain_id)
    if chain is None:
        raise ValueError(f"chain {chain_id!r} missing in {path}")
    atoms = []
    for residue in chain:
        if gemmi.find_tabulated_residue(residue.name).is_amino_acid():
            atom = residue.find_atom("CA", "*")
            if atom and atom.element.name == "C":
                atoms.append(atom)
    return atoms


def target_arrays(row, base, namespace):
    dat = row.get("dat") or ""
    if not dat:
        raise ValueError("missing dat path in manifest")
    coords, seq, dssp, description = read_dat(base / dat)
    confidence = np.full(len(coords), 100, dtype=np.uint8)
    model = row.get("model")
    if namespace not in (1, 5) and not model:
        raise ValueError("the predicted model used for the search is required")
    if model:
        atoms = model_atoms(base / model, (row.get("chain") or "A").strip())
        xyz = np.array([list(a.pos) for a in atoms], dtype=np.float32)
        # DaliLite writes C-alpha coordinates to one decimal place. Comparing
        # them also catches a wrong chain or residue order of the same length.
        if xyz.shape != coords.shape or not np.allclose(xyz, coords, atol=0.051, rtol=0):
            raise ValueError("model and .dat C-alpha coordinates or residue order differ")
        coords = xyz
        if namespace != 1:
            scores = np.array([a.b_iso for a in atoms], dtype=np.float32)
            if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 100)):
                raise ValueError("predicted-model C-alpha confidence must be in [0, 100]")
            confidence = scores.astype(np.uint8)
    centroid = coords.mean(axis=0)
    arrays = {
        "coords": coords - centroid,
        "plddt": confidence,
        "seq": np.frombuffer(seq.encode("ascii"), dtype=np.uint8),
        "dssp": np.frombuffer(dssp.encode("ascii"), dtype=np.uint8),
        "compnd": np.array(description),
    }
    return centroid, arrays


def collect_targets(results, manifest):
    summary, _ = read_alignments(results)
    wanted = set(summary["target_id"]) if not summary.empty else set()
    if not wanted:
        raise ValueError("the search contains no target IDs")
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target_id", "dat"}.issubset(reader.fieldnames or ()):
            raise ValueError("the manifest requires target_id and dat columns")
        rows = {}
        for row in reader:
            target = (row["target_id"] or "").strip()
            if target in wanted:
                if target in rows:
                    raise ValueError(f"duplicate manifest entry for {target}")
                rows[target] = row
    missing = sorted(wanted.difference(rows))
    if missing:
        raise ValueError(f"missing {len(missing)} manifest entries, including {missing[0]}")
    return rows


def build_store(results, namespace, manifest, output):
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    rows = collect_targets(results, manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        stage = Path(temporary) / "store"
        shards = stage / "shards"
        shards.mkdir(parents=True)
        with closing(sqlite3.connect(stage / "master_ledger.db")) as conn:
            conn.execute("""CREATE TABLE proteins (
                dali_id TEXT, ns_id INTEGER, n_res INTEGER, umid TEXT, shard_path TEXT,
                cx REAL, cy REAL, cz REAL, PRIMARY KEY (ns_id, dali_id))""")
            batch = {}
            for i, (target, row) in enumerate(rows.items()):
                try:
                    centroid, arrays = target_arrays(row, manifest.parent, namespace)
                except (ValueError, OSError, RuntimeError) as error:
                    raise ValueError(f"{target}: {error}") from error
                key = f"p{i}"
                shard = f"shard-{i // 256:05d}.npz"
                batch.update({f"{key}_{field}": value for field, value in arrays.items()})
                conn.execute("INSERT INTO proteins VALUES (?,?,?,?,?,?,?,?)",
                             (target, namespace, len(arrays["seq"]), key, shard, *map(float, centroid)))
                if (i + 1) % 256 == 0 or i + 1 == len(rows):
                    np.savez_compressed(shards / shard, **batch)
                    batch.clear()
            conn.commit()
        if output.exists():
            raise FileExistsError(f"output already exists: {output}")
        stage.rename(output)
    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="DALI alignment TSV or text report")
    parser.add_argument("namespace", type=parse_namespace, help="same namespace used by daliscope-pack")
    parser.add_argument("manifest", type=Path, help="TSV mapping target_id to dat, model and chain")
    parser.add_argument("--output-dir", required=True, type=Path, help="new directory for the ledger and shards")
    args = parser.parse_args(argv)
    try:
        count = build_store(args.results, args.namespace, args.manifest, args.output_dir)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")
    print(f"Wrote {count} targets to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
