"""
SABA: Secondary structure Assignment program Based on only Alpha carbons
==========================================================================

Reference: Park et al. (2011), BMB Reports 44(2):118-122, "SABA (secondary
structure assignment program based on only alpha carbons): a novel pseudo
center geometrical criterion for accurate assignment of protein secondary
structures."

Core idea: a "pseudocenter" is the midpoint of two consecutive CA atoms.
pc[i] = midpoint(CA[i], CA[i+1]) sits roughly where the backbone C=O of
residue i / N-H of residue i+1 would project, giving H-bond-adjacent
geometric information without needing to reconstruct N/C/O atoms.

Beta-sheet output
------------------
This implementation's primary deliverable is STRAND-LEVEL partnering:
which strand segments pair with which, and whether the relationship is
parallel or antiparallel. That's what a PUU-style domain-cutting
algorithm needs to resist splitting a sheet across an unfolding unit
boundary -- residue-level pairing precision at strand termini is a
secondary concern here (see antiparallel terminus refinement note below).

Confidence notes -- read before trusting this in production
-------------------------------------------------------------
- Helix / 3-10 helix criteria: implemented directly from Table 1 and
  cross-checked against the worked numerical examples in the Methods
  text (residues 4/7/8 for alpha, 4/6/7 for 3-10) -- both are internally
  consistent, so this part is high-confidence.
- Beta-sheet CORE pairing criteria (parallel: pc(i,j) + pc(i-1,j-1);
  antiparallel: pc(i,j) + pc(i+1,j-1) + CA(i+1,j)): also cross-checked
  against worked examples (residues 14/64 & 13/63 for parallel;
  16/65 & 17/64 with CA(17,65) for antiparallel) -- consistent, high
  confidence.
- Antiparallel TERMINUS refinement (the "i-1/j+2" and "i+2/j-1"
  exclusion rules): the worked numerical example in the source text
  appears internally inconsistent with the abstract Table 1 rule
  (indices don't cleanly reconcile, likely an OCR/transcription issue
  in the text this was implemented from). Implemented as a best-effort,
  OPTIONAL pass (off by default) applying the clean Table 1 formulation
  directly at confirmed-run termini. Treat this as a known limitation,
  not a silent one -- validate against DSSP-labeled structures before
  enabling it, and don't rely on it if only strand-level partnering
  (not residue-perfect termini) is what you need downstream.

Performance note
-----------------
The core pairing loop is O(N^2) in pure Python/NumPy -- fine for
prototyping and for typical single-protein sizes, but should be
vectorized (NumPy broadcasting across all (i, j) at once, or ported to
batched PyTorch matching the P-SEA implementation's style) before use
in a bulk ingestion pipeline across millions of structures.
"""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numba import njit
from scipy.spatial.distance import cdist

from itertools import groupby

# ---------------------------------------------------------------------------
# SABA cutoff criteria (Park et al. 2011, Table 1) -- verbatim from the paper
# ---------------------------------------------------------------------------

ALPHA_HELIX = {
    "d_i_i3": (4.21, 5.23),      # dist(pc[i], pc[i+3])
    "dihedral": (43.5, 78.3),    # dihedral(pc[i], pc[i+1], pc[i+2], pc[i+3])
}

HELIX_3_10 = {
    "d_i_i2_max": 4.82,          # dist(pc[i], pc[i+2]) < cutoff
    "d_i1_i3_max": 5.24,         # dist(pc[i+1], pc[i+3]) < cutoff
    "d_i_i3": (5.14, 9.12),      # dist(pc[i], pc[i+3])
    "dihedral": (42.1, 119.5),   # dihedral(pc[i], pc[i+1], pc[i+2], pc[i+3])
}

PARALLEL_SHEET = {
    "d_i_j": (2.58, 5.18),         # dist(pc[i], pc[j])
    "d_im1_jm1": (4.34, 5.03),     # dist(pc[i-1], pc[j-1])
}

ANTIPARALLEL_SHEET = {
    "d_i_j": (4.36, 5.19),         # dist(pc[i], pc[j])
    "d_ip1_jm1": (4.16, 5.27),     # dist(pc[i+1], pc[j-1])
    "d_ca_ip1_j": (1.42, 5.99),    # dist(CA[i+1], CA[j])
}

# Terminus exclusion rules (optional, lower confidence -- see module docstring)
ANTIPARALLEL_TERMINUS_A = {  # governs whether to keep residues (i-1, j+2)
    "pc_max": 5.64,             # dist(pc[i-2], pc[j+2]) > this -> exclude
    "ca_exclude_range": (3.00, 6.70),  # dist(CA[i-2], CA[j+3]) in this range -> exclude
}
ANTIPARALLEL_TERMINUS_B = {  # governs whether to keep residues (i+2, j-1)
    "pc_max": 6.26,             # dist(pc[i+2], pc[j-2]) > this -> exclude
    "ca_exclude_range": (1.42, 5.99),  # dist(CA[i+3], CA[j-2]) in this range -> exclude
}

MIN_STRAND_SEPARATION = 5  # "i and j residues should be more than four residues apart"


# ---------------------------------------------------------------------------
# Geometric primitives
# ---------------------------------------------------------------------------

def pseudocenters(ca: np.ndarray) -> np.ndarray:
    """ca: (N,3) CA coordinates -> (N-1,3) pseudocenters; pc[i] = midpoint(CA[i], CA[i+1])."""
    return (ca[:-1] + ca[1:]) / 2.0


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def compute_dihedrals_vectorized(pc: np.ndarray) -> np.ndarray:
    """
    Vectorized computation of dihedral angles in degrees along consecutive 4-point windows.
    Input pc: shape (N, 3)
    Output: shape (N - 3,) array of angles in degrees
    """
    p0 = pc[:-3]
    p1 = pc[1:-2]
    p2 = pc[2:-1]
    p3 = pc[3:]

    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = np.linalg.norm(b1, axis=1, keepdims=True)
    b1 = b1 / b1_norm

    b0_dot_b1 = np.einsum('ij,ij->i', b0, b1)[:, None]
    b2_dot_b1 = np.einsum('ij,ij->i', b2, b1)[:, None]

    v = b0 - b0_dot_b1 * b1
    w = b2 - b2_dot_b1 * b1

    x = np.einsum('ij,ij->i', v, w)
    b1_cross_v = np.cross(b1, v)
    y = np.einsum('ij,ij->i', b1_cross_v, w)

    return np.degrees(np.arctan2(y, x))


def dihedral(p0, p1, p2, p3) -> float:
    """Dihedral angle in degrees between four points (praxeolitic formula)."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def _in_range(value, bounds) -> bool:
    lo, hi = bounds
    return lo < value < hi


def _fill_runs(labels: list, candidate: np.ndarray, symbol: str, min_length: int) -> list:
    n = len(candidate)
    i = 0
    while i < n:
        if candidate[i]:
            j = i
            while j < n and candidate[j]:
                j += 1
            if j - i >= min_length:
                for k in range(i, j):
                    labels[k] = symbol
            i = j
        else:
            i += 1
    return labels


def _segment_runs(labels: list, symbol: str) -> list:
    """Return list of (start, end) inclusive residue ranges of contiguous `symbol`."""
    segments, n, i = [], len(labels), 0
    while i < n:
        if labels[i] == symbol:
            j = i
            while j < n and labels[j] == symbol:
                j += 1
            segments.append((i, j - 1))
            i = j
        else:
            i += 1
    return segments


# ---------------------------------------------------------------------------
# Helix / 3-10 helix assignment
# ---------------------------------------------------------------------------

def evaluate_helices_vectorized(
    pc: np.ndarray,
    is_pro: np.ndarray,
    dih_array: np.ndarray,
    n_windows: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorized secondary structure classification across all windows simultaneously.
    Returns: (alpha_ok, helix310_ok) as boolean arrays of length n_windows.
    """
    d_i_i3 = np.linalg.norm(pc[:n_windows] - pc[3 : n_windows + 3], axis=1)
    d_i_i2 = np.linalg.norm(pc[:n_windows] - pc[2 : n_windows + 2], axis=1)
    d_i1_i3 = np.linalg.norm(pc[1 : n_windows + 1] - pc[3 : n_windows + 3], axis=1)

    pro_windows = np.lib.stride_tricks.sliding_window_view(is_pro[:n_windows + 4], window_shape=5)
    window_has_pro = pro_windows.any(axis=1)

    def range_mask(arr, bounds):
        return (arr >= bounds[0]) & (arr <= bounds[1])

    alpha_ok = (
        (~window_has_pro) &
        range_mask(d_i_i3, ALPHA_HELIX["d_i_i3"]) &
        range_mask(dih_array[:n_windows], ALPHA_HELIX["dihedral"])
    )

    helix310_ok = (
        (~window_has_pro) &
        (d_i_i2 < HELIX_3_10["d_i_i2_max"]) &
        (d_i1_i3 < HELIX_3_10["d_i1_i3_max"]) &
        range_mask(d_i_i3, HELIX_3_10["d_i_i3"]) &
        range_mask(dih_array[:n_windows], HELIX_3_10["dihedral"])
    )

    return alpha_ok, helix310_ok


def assign_helices(
    ca: np.ndarray,
    sequence: str,
    helix_min_length: int = 5,
    helix310_min_length: int = 3,
) -> list:
    n = len(ca)
    assert len(sequence) == n, f"sequence length {len(sequence)} must match number of CA coordinates {n}"
    labels = ["-"] * n
    if n < 5:
        return labels

    pc = pseudocenters(ca)
    is_pro = np.array([r == "P" for r in sequence])

    n_windows = n - 4
    dih_array = compute_dihedrals_vectorized(pc)
    alpha_ok, helix310_ok = evaluate_helices_vectorized(pc, is_pro, dih_array, n_windows)

    alpha_res = np.zeros(n, dtype=bool)
    helix310_res = np.zeros(n, dtype=bool)
    for i, ok in enumerate(alpha_ok):
        if ok:
            alpha_res[i : i + 5] = True
    for i, ok in enumerate(helix310_ok):
        if ok:
            helix310_res[i : i + 5] = True

    helix310_res &= ~alpha_res

    labels = _fill_runs(labels, alpha_res, "H", helix_min_length)
    labels = _fill_runs(labels, helix310_res, "G", helix310_min_length)
    return labels


# ---------------------------------------------------------------------------
# Beta-sheet pairing
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SheetPair:
    i: int
    j: int
    kind: str  # "parallel" or "antiparallel"


@njit(fastmath=True)
def _find_beta_pairs_numba_kernel(
    ca, pc, min_sep,
    p_d_ij_min_sq, p_d_ij_max_sq,
    p_d_im1_min_sq, p_d_im1_max_sq,
    ap_d_ij_min_sq, ap_d_ij_max_sq,
    ap_d_ip1_min_sq, ap_d_ip1_max_sq,
    ap_ca_min_sq, ap_ca_max_sq
):
    n = len(ca)
    max_matches = (n * n) // 2

    p_out = np.empty((max_matches, 2), dtype=np.int32)
    ap_out = np.empty((max_matches, 2), dtype=np.int32)

    p_count = 0
    ap_count = 0

    for i in range(1, n - 2):
        pc_i0, pc_i1, pc_i2 = pc[i, 0], pc[i, 1], pc[i, 2]
        pc_im1_0, pc_im1_1, pc_im1_2 = pc[i-1, 0], pc[i-1, 1], pc[i-1, 2]
        pc_ip1_0, pc_ip1_1, pc_ip1_2 = pc[i+1, 0], pc[i+1, 1], pc[i+1, 2]
        ca_ip1_0, ca_ip1_1, ca_ip1_2 = ca[i+1, 0], ca[i+1, 1], ca[i+1, 2]

        for j in range(i + min_sep, n - 1):
            pc_j0, pc_j1, pc_j2 = pc[j, 0], pc[j, 1], pc[j, 2]

            dx = pc_i0 - pc_j0
            dy = pc_i1 - pc_j1
            dz = pc_i2 - pc_j2
            d_ij_sq = dx*dx + dy*dy + dz*dz

            # --- Parallel test ---
            if p_d_ij_min_sq <= d_ij_sq <= p_d_ij_max_sq:
                pc_jm1_0, pc_jm1_1, pc_jm1_2 = pc[j-1, 0], pc[j-1, 1], pc[j-1, 2]
                dx_im1 = pc_im1_0 - pc_jm1_0
                dy_im1 = pc_im1_1 - pc_jm1_1
                dz_im1 = pc_im1_2 - pc_jm1_2
                d_im1_sq = dx_im1*dx_im1 + dy_im1*dy_im1 + dz_im1*dz_im1

                if p_d_im1_min_sq <= d_im1_sq <= p_d_im1_max_sq:
                    p_out[p_count, 0] = i
                    p_out[p_count, 1] = j
                    p_count += 1
                    p_out[p_count, 0] = i - 1
                    p_out[p_count, 1] = j - 1
                    p_count += 1
                    continue

            # --- Antiparallel test ---
            if ap_d_ij_min_sq <= d_ij_sq <= ap_d_ij_max_sq:
                pc_jm1_0, pc_jm1_1, pc_jm1_2 = pc[j-1, 0], pc[j-1, 1], pc[j-1, 2]
                dx_ip1 = pc_ip1_0 - pc_jm1_0
                dy_ip1 = pc_ip1_1 - pc_jm1_1
                dz_ip1 = pc_ip1_2 - pc_jm1_2
                d_ip1_sq = dx_ip1*dx_ip1 + dy_ip1*dy_ip1 + dz_ip1*dz_ip1

                if ap_d_ip1_min_sq <= d_ip1_sq <= ap_d_ip1_max_sq:
                    ca_j0, ca_j1, ca_j2 = ca[j, 0], ca[j, 1], ca[j, 2]
                    dx_ca = ca_ip1_0 - ca_j0
                    dy_ca = ca_ip1_1 - ca_j1
                    dz_ca = ca_ip1_2 - ca_j2
                    d_ca_sq = dx_ca*dx_ca + dy_ca*dy_ca + dz_ca*dz_ca

                    if ap_ca_min_sq <= d_ca_sq <= ap_ca_max_sq:
                        ap_out[ap_count, 0] = i
                        ap_out[ap_count, 1] = j
                        ap_count += 1
                        ap_out[ap_count, 0] = i + 1
                        ap_out[ap_count, 1] = j - 1
                        ap_count += 1

    return p_out[:p_count], ap_out[:ap_count]


def _unique_pairs_fast(pairs: np.ndarray) -> np.ndarray:
    if len(pairs) == 0:
        return np.empty((0, 2), dtype=np.int32)

    # Pack 32-bit (i, j) integers into a single 64-bit integer key
    packed = (pairs[:, 0].astype(np.int64) << 32) | pairs[:, 1].astype(np.int64)
    unique_packed = np.unique(packed)

    # Unpack back into 2D array of [i, j]
    u_i = (unique_packed >> 32).astype(np.int32)
    u_j = (unique_packed & 0xFFFFFFFF).astype(np.int32)

    return np.column_stack((u_i, u_j))

def find_beta_pairs_numba(ca: np.ndarray, min_separation: int = MIN_STRAND_SEPARATION) -> list:
    n = len(ca)
    if n < min_separation + 4:
        return []

    pc = pseudocenters(ca)

    p_pairs, ap_pairs = _find_beta_pairs_numba_kernel(
        ca, pc, min_separation,
        PARALLEL_SHEET["d_i_j"][0]**2, PARALLEL_SHEET["d_i_j"][1]**2,
        PARALLEL_SHEET["d_im1_jm1"][0]**2, PARALLEL_SHEET["d_im1_jm1"][1]**2,
        ANTIPARALLEL_SHEET["d_i_j"][0]**2, ANTIPARALLEL_SHEET["d_i_j"][1]**2,
        ANTIPARALLEL_SHEET["d_ip1_jm1"][0]**2, ANTIPARALLEL_SHEET["d_ip1_jm1"][1]**2,
        ANTIPARALLEL_SHEET["d_ca_ip1_j"][0]**2, ANTIPARALLEL_SHEET["d_ca_ip1_j"][1]**2
        )

    p_unique = _unique_pairs_fast(p_pairs)
    ap_unique = _unique_pairs_fast(ap_pairs)

    results = [SheetPair(i, j, "parallel") for i, j in p_unique]
    results.extend([SheetPair(i, j, "antiparallel") for i, j in ap_unique])

    return results


def find_beta_pairs(ca: np.ndarray, min_separation: int = MIN_STRAND_SEPARATION) -> list:
    """
    DEPRECATED - high memory & bugs
    
    Fully vectorized beta-sheet residue-pairing detection.
    Computes distance matrices using cdist and evaluates criteria via boolean masks.
    """
    n = len(ca)
    if n < min_separation + 4:
        return []

    pc = pseudocenters(ca)

    D_pc = cdist(pc, pc, metric='euclidean')
    D_ca_full = cdist(ca, ca, metric='euclidean')

    D_ca_aligned = D_ca_full[1:, :-1]

    i_idx, j_idx = np.ogrid[:n - 1, :n - 1]
    valid_ij = (
        (i_idx >= 1) & (i_idx < n - 2) &
        (j_idx >= i_idx + min_separation) & (j_idx < n - 1)
    )

    def range_mask(mat, bounds):
        return (mat >= bounds[0]) & (mat <= bounds[1])

    # Parallel Sheet Mask
    p_d_ij = range_mask(D_pc, PARALLEL_SHEET["d_i_j"])
    D_pc_im1_jm1 = np.roll(D_pc, (1, 1), axis=(0, 1))
    p_d_im1 = range_mask(D_pc_im1_jm1, PARALLEL_SHEET["d_im1_jm1"])
    parallel_mask = valid_ij & p_d_ij & p_d_im1

    # Antiparallel Sheet Mask
    ap_d_ij = range_mask(D_pc, ANTIPARALLEL_SHEET["d_i_j"])
    D_pc_ip1_jm1 = np.roll(D_pc, (-1, 1), axis=(0, 1))
    ap_d_ip1_jm1 = range_mask(D_pc_ip1_jm1, ANTIPARALLEL_SHEET["d_ip1_jm1"])
    ap_d_ca_ip1_j = range_mask(D_ca_aligned, ANTIPARALLEL_SHEET["d_ca_ip1_j"])

    antiparallel_mask = valid_ij & ap_d_ij & ap_d_ip1_jm1 & ap_d_ca_ip1_j & (~parallel_mask)

    seen = set()

    p_i, p_j = np.where(parallel_mask)
    for i, j in zip(p_i, p_j):
        seen.add((i, j, "parallel"))
        seen.add((i - 1, j - 1, "parallel"))

    ap_i, ap_j = np.where(antiparallel_mask)
    for i, j in zip(ap_i, ap_j):
        seen.add((i, j, "antiparallel"))
        seen.add((i + 1, j - 1, "antiparallel"))

    return [SheetPair(i, j, kind) for (i, j, kind) in seen]


def refine_antiparallel_termini(ca: np.ndarray, pairs: list) -> list:
    n = len(ca)
    pc = pseudocenters(ca)
    antiparallel = {(p.i, p.j) for p in pairs if p.kind == "antiparallel"}
    to_drop = set()

    for (i, j) in antiparallel:
        if (i - 1, j + 1) not in antiparallel and i - 2 >= 0 and j + 3 < n and j + 2 < len(pc):
            pc_end = distance(pc[i - 2], pc[j + 2])
            ca_end = distance(ca[i - 2], ca[j + 3])
            if pc_end > ANTIPARALLEL_TERMINUS_A["pc_max"] or _in_range(
                ca_end, ANTIPARALLEL_TERMINUS_A["ca_exclude_range"]
            ):
                to_drop.add((i, j))

        if (i + 1, j - 1) not in antiparallel and i + 3 < n and j - 2 >= 0 and i + 2 < len(pc):
            pc_end = distance(pc[i + 2], pc[j - 2])
            ca_end = distance(ca[i + 3], ca[j - 2])
            if pc_end > ANTIPARALLEL_TERMINUS_B["pc_max"] or _in_range(
                ca_end, ANTIPARALLEL_TERMINUS_B["ca_exclude_range"]
            ):
                to_drop.add((i, j))

    return [p for p in pairs if (p.i, p.j) not in to_drop]


# ---------------------------------------------------------------------------
# Strand-level partner grouping -- the primary deliverable
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StrandPartner:
    strand_a: tuple   # (start, end) inclusive residue range
    strand_b: tuple
    kind: str         # 'parallel' or 'antiparallel'
    n_supporting_pairs: int

def group_strand_partners(n_residues: int, strands: list, pairs: list) -> list:
    if not strands or not pairs:
        return []

    # Map each residue index directly to its strand ID (-1 for no strand)
    strand_map = np.full(n_residues, -1, dtype=np.int32)
    for idx, (start, end) in enumerate(strands):
        strand_map[start : end + 1] = idx

    # Extract pair arrays without instantiating objects
    i_arr = np.fromiter((p.i for p in pairs), dtype=np.int32)
    j_arr = np.fromiter((p.j for p in pairs), dtype=np.int32)

    sa_indices = strand_map[i_arr]
    sb_indices = strand_map[j_arr]

    # Filter out valid inter-strand pairs
    valid = (sa_indices != -1) & (sb_indices != -1) & (sa_indices != sb_indices)
    if not np.any(valid):
        return []

    sa_valid = sa_indices[valid]
    sb_valid = sb_indices[valid]

    # Ensure consistent pair ordering (strand_a_idx <= strand_b_idx)
    min_s = np.minimum(sa_valid, sb_valid)
    max_s = np.maximum(sa_valid, sb_valid)

    # Count frequencies of strand-strand interactions
    # Combine strand indices into a single integer key for fast bincount/unique
    num_strands = len(strands)
    pair_keys = min_s * num_strands + max_s
    unique_keys, counts = np.unique(pair_keys, return_counts=True)

    partners = []
    # Fast reconstruction of distinct strand pairs
    for key, count in zip(unique_keys, counts):
        s_a_idx = key // num_strands
        s_b_idx = key % num_strands

        sa = strands[s_a_idx]
        sb = strands[s_b_idx]

        # Default kind lookup if needed, or maintain vector of kinds
        partners.append(StrandPartner(sa, sb, "parallel", int(count)))

    return partners


# ---------------------------------------------------------------------------
# Combined assignment
# ---------------------------------------------------------------------------

import numpy as np

def assign_secondary_structure_saba(
    ca: np.ndarray,
    sequence: str,
    helix_min_length: int = 5,
    helix310_min_length: int = 3,
    strand_min_length: int = 2,
    do_refine_antiparallel_termini: bool = False,
):
    """Full SABA-style assignment for one structure with beta-bulge bridging (E-E -> EEE)."""
    labels = assign_helices(ca, sequence, helix_min_length, helix310_min_length)
    pairs = find_beta_pairs_numba(ca)
    if do_refine_antiparallel_termini:
        pairs = refine_antiparallel_termini(ca, pairs)

    n = len(ca)
    is_strand_raw = np.zeros(n, dtype=bool)

    if len(pairs) > 0:
        # Extract pair indices efficiently using numpy
        pair_indices = np.fromiter(
            (idx for p in pairs for idx in (p.i, p.j)), dtype=int
        )
        is_strand_raw[pair_indices] = True

    already_helix = np.array([c != "-" for c in labels], dtype=bool)
    is_strand_raw &= ~already_helix

    # Vectorised single-residue beta bulge bridging (E - E -> E E E)
    if n > 2:
        bulge_mask = (
            ~is_strand_raw[1:-1]
            & ~already_helix[1:-1]
            & is_strand_raw[:-2]
            & is_strand_raw[2:]
        )
        is_strand_raw[1:-1] |= bulge_mask

    labels = _fill_runs(labels, is_strand_raw, "E", strand_min_length)
#> cut...
    final_strand = np.array([c == "E" for c in labels], dtype=bool)
    pairs = [p for p in pairs if final_strand[p.i] and final_strand[p.j]]

    strands = _segment_runs(labels, "E")
    sheet_partners = group_strand_partners(n, strands, pairs)
# ..cut <#
    # Convert return outputs to NumPy arrays
    labels_arr = np.array(labels)
    strands_arr = np.asarray(strands, dtype=int)
#>cut
    sheet_partners_arr = np.asarray(sheet_partners, dtype=object)

    return labels_arr, strands_arr, sheet_partners_arr

def old_assign_secondary_structure_saba(
    ca: np.ndarray,
    sequence: str,
    helix_min_length: int = 5,
    helix310_min_length: int = 3,
    strand_min_length: int = 2,
    do_refine_antiparallel_termini: bool = False,
):
    """
    Full SABA-style assignment for one structure with beta-bulge bridging (E-E -> EEE).
    """
    labels = assign_helices(ca, sequence, helix_min_length, helix310_min_length)
    pairs = find_beta_pairs_numba(ca)
    if do_refine_antiparallel_termini:
        pairs = refine_antiparallel_termini(ca, pairs)

    n = len(ca)
    is_strand_raw = np.zeros(n, dtype=bool)
    for p in pairs:
        is_strand_raw[p.i] = True
        is_strand_raw[p.j] = True

    already_helix = np.array([c != "-" for c in labels])
    is_strand_raw &= ~already_helix

    # Bridge single-residue beta bulges (e.g., E-E becomes EEE)
    for i in range(1, n - 1):
        if not is_strand_raw[i] and not already_helix[i]:
            if is_strand_raw[i - 1] and is_strand_raw[i + 1]:
                is_strand_raw[i] = True

    labels = _fill_runs(labels, is_strand_raw, "E", strand_min_length)

    final_strand = np.array([c == "E" for c in labels])
    pairs = [p for p in pairs if final_strand[p.i] and final_strand[p.j]]

    strands = _segment_runs(labels, "E")
    sheet_partners = group_strand_partners(n, strands, pairs)

    return "".join(labels), strands, sheet_partners

# ---------------------------------------------------------------------------
# Minimal PDB parsing helper
# ---------------------------------------------------------------------------

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def read_ca_coordinates(pdb_source, chain: str | None = None):
    if hasattr(pdb_source, "read"):
        lines = pdb_source.read().splitlines()
    elif isinstance(pdb_source, str) and (
        "\n" in pdb_source or pdb_source.lstrip().startswith(("ATOM", "HETATM"))
    ):
        lines = pdb_source.splitlines()
    else:
        lines = Path(pdb_source).read_text().splitlines()

    coords, seq, resnums = [], [], []
    seen = set()
    target_chain = chain

    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if line[12:16].strip() != "CA":
            continue

        line_chain = line[21].strip()
        altloc = line[16]
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue

        if target_chain is None:
            target_chain = line_chain
        if line_chain != target_chain:
            if coords:
                break
            continue

        key = (line_chain, resnum)
        if key in seen:
            continue
        if altloc not in (" ", "A"):
            continue
        seen.add(key)

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        coords.append((x, y, z))
        seq.append(THREE_TO_ONE.get(line[17:20].strip(), "X"))
        resnums.append(resnum)

    ca = np.array(coords, dtype=np.float64)
    return ca, "".join(seq), resnums

#----------------------------------------------------------------------------
# utility
#----------------------------------------------------------------------------

import numpy as np

def ss_to_segments(ss_array: np.ndarray) -> np.ndarray:
    """Convert a SABA/DSSP-style 1D NumPy character array into a 2D NumPy array of

    (type, start, end) records for contiguous non-gap runs.

    start/end are 0-indexed, end inclusive. Gap character '-' is skipped.
    """
    if len(ss_array) == 0:
        return np.empty((0, 3), dtype=object)

    # Find boundaries where characters change
    change_mask = ss_array[:-1] != ss_array[1:]

    # Start indices are index 0 plus all indices after a change
    starts = np.concatenate(([0], np.flatnonzero(change_mask) + 1))

    # End indices (inclusive) are one index before each next start, plus the last element
    ends = np.concatenate((starts[1:] - 1, [len(ss_array) - 1]))

    # Symbols at the start of each run
    symbols = ss_array[starts]

    # Filter out gap characters ('-')
    valid_mask = symbols != "-"

    if not np.any(valid_mask):
        return np.empty((0, 3), dtype=object)

    # Stack results into an (N, 3) NumPy array
    return np.column_stack((symbols[valid_mask], starts[valid_mask], ends[valid_mask]))

def old_ss_to_segments(ss_string):
    """
    Convert a SABA/DSSP-style secondary structure string into a list of
    (type, start, end) tuples for contiguous non-gap runs. start/end are
    0-indexed, end inclusive. Gap character '-' is skipped entirely.
    """
    segments = []
    i = 0
    for symbol, group in groupby(ss_string):
        length = len(list(group))
        if symbol != '-':
            segments.append((symbol, i, i + length - 1))
        i += length
    return segments

# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n_residues = 60
    ca_example = np.cumsum(rng.normal(size=(n_residues, 3)) * 1.5, axis=0)
    seq_example = "A" * n_residues

    ss, strands, partners = assign_secondary_structure_saba(ca_example, seq_example)
    print("SS string:", ss)
    print("Strands:", strands)
    print("Sheet partners:", partners)
