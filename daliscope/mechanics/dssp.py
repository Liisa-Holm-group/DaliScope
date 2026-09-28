"""
DSSP-style beta sheet definition: bridges -> ladders -> sheets, at
residue resolution, using SABA's confirmed H-bond pairs directly.

Replaces the earlier strand-level chimera detection/splitting entirely.
Sheet membership is read off real bridge connectivity, never off
segment-level graph topology or fuzzy contiguity heuristics -- so
strand splitting becomes an EXACT rule (cut wherever two consecutive
residues' sheet labels differ), not something needing tunable overlap/
block-length parameters to approximate.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import networkx as nx


@dataclass(frozen=True)
class BetaBridge:
    i: int
    j: int
    kind: str  # 'parallel' or 'antiparallel'

def old_find_beta_bridges(sheet_pairs: Iterable[Any]) -> List[BetaBridge]:
    """
    sheet_pairs: raw residue-level SheetPair(i, j, kind) from
        saba.find_beta_pairs (already-confirmed H-bond pairs).

    A genuine BRIDGE requires two CONSECUTIVE confirmed H-bond pairs of
    the SAME type:
        parallel:     (i, j) and (i+1, j+1) both confirmed
        antiparallel: (i, j) and (i+1, j-1) both confirmed
    A single isolated H-bond pair with no consecutive partner does NOT
    count as a bridge -- same "joint confirmation, not a single isolated
    match" principle used for the original parallel/antiparallel core
    tests in saba.py.

    Returns every (i, j) pair that is confirmed as part of at least one
    bridge (both members of each valid consecutive pair are included).
    """
    print("#find_beta_bridges sheet_pairs", sheet_pairs)
    pair_lookup: Dict[Tuple[int, int], str] = {(sp.i, sp.j): sp.kind for sp in sheet_pairs}
    bridges: Dict[Tuple[int, int], str] = {}

    for sp in sheet_pairs:
        i, j, kind = sp.i, sp.j, sp.kind
        if kind == "parallel":
            step = (i + 1, j + 1)
        else:
            step = (i + 1, j - 1)

        if pair_lookup.get(step) == kind:
            bridges[(i, j)] = kind
            bridges[step] = kind

    return [BetaBridge(i, j, kind) for (i, j), kind in bridges.items()]


def old_build_ladders(bridges: List[BetaBridge]) -> List[List[BetaBridge]]:
    """
    Group bridges into maximal chains of CONSECUTIVE same-type bridges
    -- DSSP's "ladder". Not needed for sheet connectivity itself (bridge
    edges alone give that via connected components), but preserves
    per-stretch parallel/antiparallel character as useful metadata
    (e.g. for the sheet-topology fingerprint).
    """
    by_type: Dict[str, Dict[Tuple[int, int], BetaBridge]] = {
        "parallel": {}, "antiparallel": {}
    }
    for b in bridges:
        by_type[b.kind][(b.i, b.j)] = b

    visited = set()
    ladders: List[List[BetaBridge]] = []

    for kind, lookup in by_type.items():
        step = (1, 1) if kind == "parallel" else (1, -1)
        for (i, j), b in lookup.items():
            if (i, j) in visited:
                continue
            # walk backward to the start of this chain
            ci, cj = i, j
            while (ci - step[0], cj - step[1]) in lookup:
                ci, cj = ci - step[0], cj - step[1]
            # walk forward, collecting the whole chain
            chain = []
            while (ci, cj) in lookup:
                chain.append(lookup[(ci, cj)])
                visited.add((ci, cj))
                ci, cj = ci + step[0], cj + step[1]
            ladders.append(chain)

    return ladders


def old_compute_sheet_labels(ladders: List[List[BetaBridge]]) -> Dict[int, int]:
    """
    residue -> sheet_id, via connected components of LADDERS (not raw
    bridge edges directly) -- two ladders merge into one sheet if they
    share ANY residue, and every residue belonging to either merged
    ladder inherits the shared sheet label. This matches DSSP's actual
    "sheets are ladders connected by shared residues" definition.

    Connecting individual residues directly (skipping the ladder level)
    silently shatters ordinary sheets: in a normal 3-strand sheet,
    residue 20 (shared between the strand1-strand2 and strand2-strand3
    ladders) correctly bridges them, but residue 21 has no DIRECT edge
    to residue 20 at all (they're not H-bond partners of each other,
    just both members of strand2) -- so a residue-only graph produces
    one tiny disconnected component per "column" of the sheet instead
    of one unified sheet.
    """
    residue_to_ladders: Dict[int, set] = defaultdict(set)
    for lidx, ladder in enumerate(ladders):
        residues = set()
        for b in ladder:
            residues.add(b.i)
            residues.add(b.j)
        for r in residues:
            residue_to_ladders[r].add(lidx)

    G = nx.Graph()
    G.add_nodes_from(range(len(ladders)))
    for lidxs in residue_to_ladders.values():
        lidxs = list(lidxs)
        for k in range(1, len(lidxs)):
            G.add_edge(lidxs[0], lidxs[k])

    ladder_to_sheet: Dict[int, int] = {}
    for sheet_id, component in enumerate(nx.connected_components(G)):
        for lidx in component:
            ladder_to_sheet[lidx] = sheet_id

    labels: Dict[int, int] = {}
    for r, lidxs in residue_to_ladders.items():
        labels[r] = ladder_to_sheet[next(iter(lidxs))]
    return labels


def old_split_strands_by_sheet_membership(
    strands: List[Tuple[int, int]],
    sheet_pairs: Iterable[Any],
) -> Tuple[List[Tuple[int, int]], Dict[int, List[int]], Dict[int, int]]:
    """
    The direct replacement for find_chimeric_strands/split_chimeric_segments.
    Cuts a SABA strand wherever two CONSECUTIVE residues have DIFFERING,
    both-defined sheet labels -- an exact rule, no tunable overlap/
    block-length heuristics needed, because sheet_labels comes straight
    from real bridge connectivity rather than an aggregated proxy.

    Residues with NO confirmed bridge (sheet_labels doesn't cover them)
    do not force a cut on their own -- only a genuine disagreement
    between two labeled, differing sheets does.

    Returns:
        new_strands: refined (start, end) segment list
        old_to_new: old segment id -> list of new segment ids
        sheet_labels: residue -> sheet_id (for downstream sheet-collapse use)
    """
    bridges = find_beta_bridges(sheet_pairs)
    ladders = build_ladders(bridges)
    sheet_labels = compute_sheet_labels(ladders)

    new_strands: List[Tuple[int, int]] = []
    old_to_new: Dict[int, List[int]] = {}

    for seg_id, (start, end) in enumerate(strands):
        cut_points = []
        last_label: Optional[int] = sheet_labels.get(start)
        for r in range(start + 1, end + 1):
            label = sheet_labels.get(r)
            if last_label is not None and label is not None and label != last_label:
                cut_points.append(r)
            if label is not None:
                last_label = label

        if not cut_points:
            old_to_new[seg_id] = [len(new_strands)]
            new_strands.append((start, end))
            continue

        new_ids = []
        cursor = start
        for cp in cut_points:
            new_ids.append(len(new_strands))
            new_strands.append((cursor, cp - 1))
            cursor = cp
        new_ids.append(len(new_strands))
        new_strands.append((cursor, end))
        old_to_new[seg_id] = new_ids

    new_strands.sort()
    return new_strands, old_to_new, sheet_labels

"""DSSP-style beta sheet definition: bridges -> ladders -> sheets, at residue

resolution, using SABA's confirmed H-bond pairs directly.

Replaces the earlier strand-level chimera detection/splitting entirely.
"""

from dataclasses import dataclass
from typing import Any, Iterable, List, Tuple
import numpy as np


@dataclass(frozen=True)
class BetaBridge:
    i: int
    j: int
    kind: str  # 'parallel' or 'antiparallel'


def find_beta_bridges(sheet_pairs: Iterable[Any]) -> List[BetaBridge]:
    """Finds consecutive H-bond pairs of the same type and returns them as

    BetaBridge instances.
    """
    if len(sheet_pairs) == 0:
        return []

    # Extract arrays for vectorised pair matching
    i_arr = np.fromiter((sp.i for sp in sheet_pairs), dtype=np.int32)
    j_arr = np.fromiter((sp.j for sp in sheet_pairs), dtype=np.int32)
    is_p = np.fromiter((sp.kind == "parallel" for sp in sheet_pairs), dtype=bool)

    # Pack (i, j, kind) into unique 64-bit keys: bit 63 indicates kind
    # (0: antiparallel, 1: parallel)
    keys = (
        (i_arr.astype(np.int64) << 32)
        | (j_arr.astype(np.int64) & 0x7FFFFFFF)
        | (is_p.astype(np.int64) << 62)
    )
    lookup_set = set(keys)

    # Compute step candidates
    step_i = i_arr + 1
    step_j = np.where(is_p, j_arr + 1, j_arr - 1)

    step_keys = (
        (step_i.astype(np.int64) << 32)
        | (step_j.astype(np.int64) & 0x7FFFFFFF)
        | (is_p.astype(np.int64) << 62)
    )

    # Vectorised check if step key exists in original lookup
    valid_mask = np.fromiter((k in lookup_set for k in step_keys), dtype=bool)

    if not np.any(valid_mask):
        return []

    # Collect unique bridges
    bridge_keys = set()
    valid_indices = np.flatnonzero(valid_mask)
    for idx in valid_indices:
        bridge_keys.add(keys[idx])
        bridge_keys.add(step_keys[idx])

    # Unpack keys back to BetaBridge objects
    bridges = []
    for k in bridge_keys:
        kind = "parallel" if (k >> 62) & 1 else "antiparallel"
        i = int(k >> 32) & 0x3FFFFFFF
        j = int(k & 0x7FFFFFFF)
        bridges.append(BetaBridge(i, j, kind))

    return bridges


def build_ladders(bridges: List[BetaBridge]) -> List[List[BetaBridge]]:
    """Group bridges into maximal chains of CONSECUTIVE same-type bridges

    (DSSP's "ladder").
    """
    if not bridges:
        return []

    by_type = {"parallel": {}, "antiparallel": {}}
    for b in bridges:
        by_type[b.kind][(b.i, b.j)] = b

    visited = set()
    ladders: List[List[BetaBridge]] = []

    for kind, lookup in by_type.items():
        step = (1, 1) if kind == "parallel" else (1, -1)
        for (i, j), b in lookup.items():
            if (i, j) in visited:
                continue

            # Walk backward to chain start
            ci, cj = i, j
            while (ci - step[0], cj - step[1]) in lookup:
                ci, cj = ci - step[0], cj - step[1]

            # Walk forward to collect full chain
            chain = []
            while (ci, cj) in lookup:
                chain.append(lookup[(ci, cj)])
                visited.add((ci, cj))
                ci, cj = ci + step[0], cj + step[1]
            ladders.append(chain)

    return ladders


def compute_sheet_labels(
    ladders: List[List[BetaBridge]], max_res: int
) -> np.ndarray:
    """Computes a 1D array mapping residue index -> sheet_id.

    Returns -1 for residues not assigned to any beta sheet ladder.
    """
    labels = np.full(max_res, -1, dtype=np.int32)
    if not ladders:
        return labels

    # Build residue to ladder mapping
    res_to_ladders = {}
    for lidx, ladder in enumerate(ladders):
        for b in ladder:
            res_to_ladders.setdefault(b.i, set()).add(lidx)
            res_to_ladders.setdefault(b.j, set()).add(lidx)

    # Disjoint Set Union (Union-Find) to merge connected ladders
    parent = list(range(len(ladders)))

    def find(i):
        path = []
        while parent[i] != i:
            path.append(i)
            i = parent[i]
        for node in path:
            parent[node] = i
        return i

    def union(i, j):
        root_i, root_j = find(i), find(j)
        if root_i != root_j:
            parent[root_i] = root_j

    for lidxs in res_to_ladders.values():
        l_list = list(lidxs)
        for k in range(1, len(l_list)):
            union(l_list[0], l_list[k])

    # Map root ladder component IDs to contiguous sheet IDs (0, 1, 2...)
    component_map = {}
    sheet_counter = 0

    for r, lidxs in res_to_ladders.items():
        root = find(next(iter(lidxs)))
        if root not in component_map:
            component_map[root] = sheet_counter
            sheet_counter += 1
        labels[r] = component_map[root]

    return labels


def split_strands_by_sheet_membership(
    strands: np.ndarray,
    sheet_pairs: Iterable[Any],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cuts SABA strand segments wherever two CONSECUTIVE residues have

    differing defined sheet labels.

    Parameters
    ----------
    strands : np.ndarray
        (N, 2) array of (start, end) strand range indices.
    sheet_pairs : Iterable[Any]
        Raw residue-level SheetPair objects or array.

    Returns
    -------
    new_strands : np.ndarray
        (M, 2) integer array of refined strand ranges.
    old_to_new : np.ndarray
        1D object array where old_to_new[i] is a list of new strand indices.
    sheet_labels : np.ndarray
        1D integer array mapping residue_idx -> sheet_id (-1 for unlabeled).
    """
    bridges = find_beta_bridges(sheet_pairs)
    ladders = build_ladders(bridges)

    # Determine maximum residue index for spatial mapping array
    max_res = 0
    if len(strands) > 0:
        max_res = max(max_res, int(np.max(strands)) + 2)
    if len(sheet_pairs) > 0:
        max_res = max(
            max_res, max(max(sp.i, sp.j) for sp in sheet_pairs) + 2
        )

    sheet_labels = compute_sheet_labels(ladders, max_res)

    new_strands_list: List[Tuple[int, int]] = []
    old_to_new_list: List[List[int]] = []

    for seg_id, (start, end) in enumerate(strands):
        start, end = int(start), int(end)
        cut_points = []
        last_label = sheet_labels[start] if sheet_labels[start] != -1 else None

        for r in range(start + 1, end + 1):
            lbl = sheet_labels[r]
            curr_label = lbl if lbl != -1 else None

            if (
                last_label is not None
                and curr_label is not None
                and curr_label != last_label
            ):
                cut_points.append(r)
            if curr_label is not None:
                last_label = curr_label

        if not cut_points:
            old_to_new_list.append([len(new_strands_list)])
            new_strands_list.append((start, end))
            continue

        new_ids = []
        cursor = start
        for cp in cut_points:
            new_ids.append(len(new_strands_list))
            new_strands_list.append((cursor, cp - 1))
            cursor = cp
        new_ids.append(len(new_strands_list))
        new_strands_list.append((cursor, end))
        old_to_new_list.append(new_ids)

    # Convert results to NumPy arrays
    new_strands = (
        np.array(new_strands_list, dtype=np.int32)
        if new_strands_list
        else np.empty((0, 2), dtype=np.int32)
    )

    if len(new_strands) > 0:
        # Ensure new_strands are sorted by start index
        sort_order = np.argsort(new_strands[:, 0])
        new_strands = new_strands[sort_order]

    old_to_new = np.array(old_to_new_list, dtype=object)

    return new_strands, old_to_new, sheet_labels
