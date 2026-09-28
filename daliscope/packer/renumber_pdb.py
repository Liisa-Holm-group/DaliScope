import sys
import gzip

SOLVENT_BLACKLIST = {
    "HOH", "WAT", "TIP", "DOD", "SO4", "PO4", "EDT", "ACT",
    "PEG", "PGE", "PG4", "EDO", "GOL", "DMS", "MES", "TRS",
    "CL", "NA", "K", "MG", "CA", "2PA"
}
METAL_WHITELIST = {
    "ZN", "FE", "CU", "MN", "MG", "CA", "CO", "NI", "MO", "W"
}
BACKBONE = {"N", "CA", "C", "O"}

def is_relevant_hetatm(line):
    res_name = line[17:20].strip().upper()
    if res_name in SOLVENT_BLACKLIST and res_name not in METAL_WHITELIST:
        return False
    element  = line[76:78].strip().upper()
    if element in METAL_WHITELIST or res_name in METAL_WHITELIST:
        return True
    return len(res_name) >= 3

def extract_and_renumber_pdb(pdb_path, chain_id, output_path):
    open_func = gzip.open if pdb_path.endswith('.gz') else open

    if True:
        # --- Pass 1: collect ATOM residues to validate backbones ---
        backbone_seen = {}  # res_key -> set of backbone atoms found
        with open_func(pdb_path, 'rt') as f:
            for line in f:
                if line.startswith("ENDMDL"):
                    break
                if not line.startswith("ATOM  ") and not line.startswith("HETATM"): continue
		# alt. location must be blank of
                if line[16] not in [' ','A']: continue
                if (line.startswith("ATOM  ") or line.startswith("HETATM")) and line[21].strip() == chain_id:
                    key = (line[22:26].strip(), line[26].strip())
                    atom = line[12:16].strip()
                    if key not in backbone_seen:
                        backbone_seen[key] = set()
                    backbone_seen[key].add(atom)

        complete = {k for k, atoms in backbone_seen.items()
                    if BACKBONE.issubset(atoms)}

        # --- Pass 2: write lines in original order ---
        renumber = {}   # res_key -> new sequential number
        counter = 1

        with open_func(pdb_path, 'rt') as f, open(output_path, 'w') as out:
            for line in f:
                if line.startswith("ENDMDL"):
                    break

                if len(line) < 26:
                    continue

                if line[21].strip() != chain_id:
                    continue

                elif line.startswith("HETATM") and line[12] != ' ' and is_relevant_hetatm(line):
                        out.write(line)

                elif line.startswith("ATOM  ") or line.startswith("HETATM"):
                    key = (line[22:26].strip(), line[26].strip())
                    if key not in complete:
                        continue
                    if key not in renumber:
                        renumber[key] = counter
                        counter += 1
                    new_num = renumber[key]
                    out.write(line[:22] + f"{new_num:>4}" + " " + line[27:])

            out.write("TER\nEND\n")
        return True

#    except Exception:
    return False

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"USAGE: {sys.argv[0]} <pdbfile> <chainid> <outfile>")
        sys.exit()
    pdbfile = sys.argv[1]
    chainid = sys.argv[2]
    outfile = sys.argv[3]
    extract_and_renumber_pdb(pdbfile, chainid, outfile)
