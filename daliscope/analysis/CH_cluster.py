import numpy as np
import itertools
from collections import defaultdict

def parse_candidate_atoms(pdb_str, chain=None):
    """Extract Cys SG and His ND1/NE2 atoms from a PDB-format string."""
    atoms = []
    for line in pdb_str.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        record_chain = line[21]
        if chain and record_chain != chain:
            continue
        resname = line[17:20].strip()
        atom_name = line[12:16].strip()
        resnum = int(line[22:26])
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        if resname == "CYS" and atom_name == "SG":
            atoms.append(dict(resnum=resnum, resname="CYS", atom="SG",
                               chain=record_chain, coord=np.array([x, y, z])))
        elif resname == "HIS" and atom_name in ("ND1", "NE2"):
            atoms.append(dict(resnum=resnum, resname="HIS", atom=atom_name,
                               chain=record_chain, coord=np.array([x, y, z])))
    return atoms


def find_disulfides(atoms, bond_range=(1.9, 2.5)):
    """Cys SG-SG pairs within covalent bonding distance."""
    cys = [a for a in atoms if a["resname"] == "CYS"]
    bonds = []
    for a, b in itertools.combinations(cys, 2):
        d = np.linalg.norm(a["coord"] - b["coord"])
        if bond_range[0] <= d <= bond_range[1]:
            bonds.append({"res1": a["resnum"], "res2": b["resnum"],
                           "chain": a["chain"], "distance": round(d, 2)})
    return bonds


def cluster_zinc_candidates(atoms, disulfides, cluster_cutoff=4.6,
                             min_ligands=3, max_ligands=4):
    """
    Group Cys SG / His ND1/NE2 atoms into spatial clusters and score
    them as candidate metal-binding sites based on tetrahedral geometry.
    Cys residues already in a disulfide are excluded (their thiolate
    isn't free to coordinate a metal).
    """
    disulfide_res = {r for b in disulfides for r in (b["res1"], b["res2"])}
    candidates = [a for a in atoms if not (a["resname"] == "CYS" and a["resnum"] in disulfide_res)]
    n = len(candidates)
    coords = np.array([a["coord"] for a in candidates])

    # union-find clustering by distance
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i, j in itertools.combinations(range(n), 2):
        if np.linalg.norm(coords[i] - coords[j]) <= cluster_cutoff:
            union(i, j)

    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    sites = []
    for idxs in groups.values():
        # one atom per residue (His may contribute ND1 and NE2 both being
        # in range; keep whichever sits closer to the group's centroid)
        by_res = defaultdict(list)
        for i in idxs:
            a = candidates[i]
            by_res[(a["chain"], a["resnum"])].append(i)

        if len(by_res) < min_ligands:
            continue

        rough_centroid = coords[idxs].mean(axis=0)
        rep = [min(v, key=lambda i: np.linalg.norm(coords[i] - rough_centroid))
               for v in by_res.values()]

        # cap to the most tightly clustered max_ligands if an oversized group formed
        if len(rep) > max_ligands:
            centroid = coords[rep].mean(axis=0)
            rep = sorted(rep, key=lambda i: np.linalg.norm(coords[i] - centroid))[:max_ligands]

        rep_coords = coords[rep]
        centroid = rep_coords.mean(axis=0)
        d_centroid = np.linalg.norm(rep_coords - centroid, axis=1)
        pairwise = [np.linalg.norm(coords[a] - coords[b]) for a, b in itertools.combinations(rep, 2)]

        n_cys = sum(1 for i in rep if candidates[i]["resname"] == "CYS")
        n_his = sum(1 for i in rep if candidates[i]["resname"] == "HIS")

        mean_c, std_c = d_centroid.mean(), d_centroid.std()
        no_bonded_pairs = all(p > 3.0 for p in pairwise)  # excludes stray disulfide-range pairs

        if len(rep) in (3, 4) and 1.8 <= mean_c <= 3.0 and no_bonded_pairs:
            confidence = "high" if std_c < 0.5 else "medium"
        elif len(rep) == 2:
            confidence = "low (incomplete/partial site)"
        else:
            confidence = "low"

        sites.append({
            "residues": sorted({(candidates[i]["resname"], candidates[i]["resnum"],
                                  candidates[i]["chain"]) for i in rep}, key=lambda x: x[1]),
            "n_ligands": len(rep),
            "composition": f"Cys{n_cys}His{n_his}" if n_his else f"Cys{n_cys}",
            "centroid": tuple(round(c, 2) for c in centroid),
            "mean_dist_to_centroid": round(mean_c, 2),
            "std_dist_to_centroid": round(std_c, 2),
            "min_pairwise": round(min(pairwise), 2) if pairwise else None,
            "max_pairwise": round(max(pairwise), 2) if pairwise else None,
            "confidence": confidence,
        })

    return sorted(sites, key=lambda s: s["residues"][0][1])


def evaluate_cys_his_sites(pdb_str, chain=None, disulfide_range=(1.9, 2.5),
                            zn_cluster_cutoff=4.6, min_zn_ligands=3):
    """
    Scan a PDB structure string for:
      - disulfide bonds (Cys SG-SG ~2.05 A)
      - candidate metal-binding clusters (Cys SG / His ND1,NE2 converging
        on a shared point with tetrahedral-like spacing, ~3.6-5.0 A pairwise,
        ~2.0-2.8 A to centroid)

    Returns dict with 'disulfides' and 'zinc_sites' lists.
    Purely geometric heuristic on a metal-free AlphaFold model — treat
    'high' confidence as a strong candidate, not proof; cross-check against
    sequence motifs (CxxC spacing, zinc-finger consensus) and, ideally,
    AlphaFill or homology to a metal-bound template.
    """
    atoms = parse_candidate_atoms(pdb_str, chain=chain)
    disulfides = find_disulfides(atoms, bond_range=disulfide_range)
    zinc_sites = cluster_zinc_candidates(atoms, disulfides,
                                          cluster_cutoff=zn_cluster_cutoff,
                                          min_ligands=min_zn_ligands)
    return {"disulfides": disulfides, "zinc_sites": zinc_sites}


def print_report(results):
    print(f"Disulfide bonds found: {len(results['disulfides'])}")
    for b in results["disulfides"]:
        print(f"  Cys{b['res1']}-Cys{b['res2']} (chain {b['chain']}): {b['distance']} A")
    print(f"\nCandidate metal-binding clusters: {len(results['zinc_sites'])}")
    for s in results["zinc_sites"]:
        res_str = ", ".join(f"{r}{n}" for r, n, c in s["residues"])
        print(f"  [{s['confidence']}] {s['composition']} — {res_str}")
        print(f"      pairwise range: {s['min_pairwise']}-{s['max_pairwise']} A, "
              f"mean/std dist to centroid: {s['mean_dist_to_centroid']}/{s['std_dist_to_centroid']} A")
