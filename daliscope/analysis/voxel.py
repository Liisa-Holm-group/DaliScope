import numpy as np
from numba import njit

@njit
def _numba_lcs_core(seq_a, seq_b, min_length):
    m, n = len(seq_a), len(seq_b)
    if m == 0 or n == 0:
        return 0.0, np.empty(0, dtype=seq_a.dtype)

    # 1. Compute DP Matrix at C-speed
    dp = np.zeros((m + 1, n + 1), dtype=np.int32)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq_a[i - 1] == seq_b[j - 1]:
                dp[i, j] = dp[i - 1, j - 1] + 1
            else:
                dp[i, j] = max(dp[i - 1, j], dp[i, j - 1])

    lcs_len = dp[m, n]
    score = lcs_len / min(m, n) if min(m, n) > 0 else 0.0

    if lcs_len < min_length:
        return score, np.empty(0, dtype=seq_a.dtype)

    # 2. Fast Backtrace to recover cells (pre-allocating max possible size)
    covered_arr = np.empty(lcs_len, dtype=seq_a.dtype)
    idx = lcs_len - 1

    i, j = m, n
    while i > 0 and j > 0:
        if seq_a[i - 1] == seq_b[j - 1]:
            covered_arr[idx] = seq_a[i - 1]
            idx -= 1
            i -= 1
            j -= 1
        elif dp[i - 1, j] >= dp[i, j - 1]:
            i -= 1
        else:
            j -= 1

    return score, covered_arr

def lcs_score_and_coverage_fast(trace_a: list, trace_b: list, min_length: int = 5) -> tuple:
    # Convert incoming traces to flat NumPy integer arrays for Numba
    seq_a = np.array([c for c, _ in trace_a], dtype=np.int64)
    seq_b = np.array([c for c, _ in trace_b], dtype=np.int64)

    score, covered_arr = _numba_lcs_core(seq_a, seq_b, min_length)

    # Return a set to preserve the rest of your pipeline's downstream logic
    return score, set(covered_arr)

def lcs_score_and_coverage(trace_a: list, trace_b: list,
                           min_length: int = 5) -> tuple:
    """
    LCS with backtrace — returns both a normalized score AND the
    actual set of matched cells, which is what's needed to update
    the explained-weight pool correctly.
    """
    seq_a = [c for c, _ in trace_a]
    seq_b = [c for c, _ in trace_b]
    m, n = len(seq_a), len(seq_b)

    if m == 0 or n == 0:
        return 0.0, set()

    dp = np.zeros((m+1, n+1), dtype=np.int32)
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq_a[i-1] == seq_b[j-1]:
                dp[i,j] = dp[i-1,j-1] + 1
            else:
                dp[i,j] = max(dp[i-1,j], dp[i,j-1])

    lcs_len = dp[m, n]
    score   = lcs_len / min(m, n) if min(m, n) > 0 else 0.0

    if lcs_len < min_length:
        return score, set()   # too short to count as a real match

    # backtrace to recover the actual matched cells
    covered = set()
    i, j = m, n
    while i > 0 and j > 0:
        if seq_a[i-1] == seq_b[j-1]:
            covered.add(seq_a[i-1])
            i -= 1; j -= 1
        elif dp[i-1, j] >= dp[i, j-1]:
            i -= 1
        else:
            j -= 1

    return score, covered

import numpy as np
from collections import Counter


# ================================================================== #
# Grid discretization                                                 #
# ================================================================== #

def grid_trace(coords: np.ndarray,
               anchor_mask: np.ndarray,
               cell_size: float = 4.5) -> tuple:
    """
    Convert a target's CA trace (already superimposed onto the anchor
    frame) into a sequence of occupied grid cells, restricted to
    residues OUTSIDE the anchor domain.

    Parameters
    ----------
    coords      : (L, 3) float32 — full CA trace for this target,
                  in the shared anchor-superimposed frame
    anchor_mask : (L,) bool — True for residues belonging to the
                  anchor domain (these are excluded from the trace)
    cell_size   : Å — grid cell edge length

    Returns
    -------
    cell_sequence : list of int cell-hash values, IN CHAIN ORDER,
                    consecutive duplicate cells collapsed (a residue
                    sitting in the same cell as the previous one
                    contributes no new information to the trace)
    """
    non_anchor_coords = coords[~anchor_mask]

    if len(non_anchor_coords) == 0:
        return []

    cell_keys = np.floor(non_anchor_coords / cell_size).astype(np.int64)
    OFFSET = 1 << 16
    hashed = (
        (cell_keys[:, 0] + OFFSET) * (1 << 34) +
        (cell_keys[:, 1] + OFFSET) * (1 << 17) +
        (cell_keys[:, 2] + OFFSET)
    )

    # collapse consecutive duplicates — a stretch of residues in one
    # cell is one "visit" to that cell, not N separate visits
    cell_sequence = []
    for h in hashed:
        if not cell_sequence or cell_sequence[-1] != h:
            cell_sequence.append(int(h))

    return cell_sequence


def cell_atom_counts(coords: np.ndarray,
                     anchor_mask: np.ndarray,
                     cell_size: float = 4.5) -> Counter:
    """
    Atom (residue) count per occupied cell — used as the weight
    in the Round 1 weighted Jaccard.
    """
    non_anchor_coords = coords[~anchor_mask]
    if len(non_anchor_coords) == 0:
        return Counter()

    cell_keys = np.floor(non_anchor_coords / cell_size).astype(np.int64)
    OFFSET = 1 << 16
    hashed = (
        (cell_keys[:, 0] + OFFSET) * (1 << 34) +
        (cell_keys[:, 1] + OFFSET) * (1 << 17) +
        (cell_keys[:, 2] + OFFSET)
    )
    return Counter(hashed.tolist())


# ================================================================== #
# Round 1: weighted Jaccard                                          #
# ================================================================== #

def weighted_jaccard(counts_i: Counter, counts_j: Counter) -> float:
    """
    Weighted Jaccard similarity between two targets' cell-occupancy
    profiles, weighted by atom count per cell.

    sum(min(w_i, w_j) over shared cells) / sum(max(w_i, w_j) over all cells)

    This is the standard weighted Jaccard (Ruzicka similarity).
    """
    all_cells = set(counts_i) | set(counts_j)
    if not all_cells:
        return 0.0

    num = sum(min(counts_i.get(c, 0), counts_j.get(c, 0)) for c in all_cells)
    den = sum(max(counts_i.get(c, 0), counts_j.get(c, 0)) for c in all_cells)

    return num / den if den > 0 else 0.0


def round1_triage(all_counts: list, threshold: float = 0.2) -> list:
    """
    All-pairs weighted Jaccard, returning pairs above threshold
    for Round 2 refinement. O(N^2) in pair count but each comparison
    is cheap (Counter intersection, not coordinate-level work).

    For large N, this should be preceded by a cheaper coarse filter
    (e.g. comparing centroid-of-occupied-cells, or LSH on the cell
    sets) — flagged below as a scaling consideration.
    """
    N = len(all_counts)
    candidates = []
    for i in range(N):
        for j in range(i + 1, N):
            sim = weighted_jaccard(all_counts[i], all_counts[j])
            if sim >= threshold:
                candidates.append((i, j, sim))
    return candidates


# ================================================================== #
# Round 2: LCS refinement on ordered cell sequences                  #
# ================================================================== #

def lcs_length(seq_a: list, seq_b: list) -> int:
    """
    Standard LCS length via DP. O(|A| * |B|).

    For grid trace sequences this is fast — cell sequences are much
    shorter than residue counts (consecutive duplicates collapsed),
    typically tens of entries for a 50-100 residue auxiliary domain.
    """
    m, n = len(seq_a), len(seq_b)
    if m == 0 or n == 0:
        return 0

    # rolling 1D DP — O(min(m,n)) memory
    if m < n:
        seq_a, seq_b = seq_b, seq_a
        m, n = n, m

    prev = [0] * (n + 1)
    for i in range(1, m + 1):
        curr = [0] * (n + 1)
        for j in range(1, n + 1):
            if seq_a[i-1] == seq_b[j-1]:
                curr[j] = prev[j-1] + 1
            else:
                curr[j] = max(prev[j], curr[j-1])
        prev = curr
    return prev[n]


def lcs_similarity(seq_a: list, seq_b: list) -> float:
    """
    Normalised LCS — fraction of the shorter sequence's cells that
    participate in the longest common ordered subsequence.
    1.0 = one trace's cell sequence is fully embedded in the other's,
    in matching order.
    """
    if not seq_a or not seq_b:
        return 0.0
    L = lcs_length(seq_a, seq_b)
    return L / min(len(seq_a), len(seq_b))


def round2_refine(candidates: list,
                  cell_sequences: list) -> list:
    """
    Refine Round 1 candidates with LCS on ordered cell sequences.
    Returns (i, j, jaccard_sim, lcs_sim) for all candidates,
    sorted by lcs_sim descending.
    """
    refined = []
    for i, j, jsim in candidates:
        lsim = lcs_similarity(cell_sequences[i], cell_sequences[j])
        refined.append((i, j, jsim, lsim))

    refined.sort(key=lambda x: -x[3])
    return refined

def compute_voxel_weights(target_voxel_sets: list) -> dict:
    """
    voxel -> log(occupancy) weight, where occupancy = number of
    targets whose voxel set contains this voxel.
    """
    from collections import Counter
    counts = Counter()
    for vset in target_voxel_sets:
        counts.update(vset)
    return {v: np.log1p(c) for v, c in counts.items()}   # log1p avoids log(0)/log(1)=0 edge cases

def greedy_max_coverage_fast(target_voxel_sets: list,
                             voxel_weights: dict,
                             k: int,
                             target_labels: list = None,
                             min_gain: float = 0.0) -> list:
    """
    Same algorithm, with lazy-greedy (priority queue) speedup.
    Marginal gains only decrease as `covered` grows (submodularity),
    so we can use a max-heap with lazy re-evaluation — avoids
    recomputing gain for targets that clearly can't be the best.
    """
    import heapq

    N = len(target_voxel_sets)
    if target_labels is None:
        target_labels = list(range(N))

    covered      = set()
    total_weight = sum(voxel_weights.values())
    selected     = []

    # initial gains — full voxel set weight for each target
    heap = []
    for i in range(N):
        gain = sum(voxel_weights.get(v, 0.0) for v in target_voxel_sets[i])
        heapq.heappush(heap, (-gain, i, 0))   # 0 = "freshness" counter

    cumulative_gain = 0.0
    iteration = 0

    while heap and len(selected) < k:
        neg_gain, i, stamp = heapq.heappop(heap)

        if stamp < iteration:
            # stale — recompute actual current gain (submodularity:
            # true gain can only be <= the stored stale value)
            new_voxels = target_voxel_sets[i] - covered
            actual_gain = sum(voxel_weights.get(v, 0.0) for v in new_voxels)
            heapq.heappush(heap, (-actual_gain, i, iteration))
            continue

        gain = -neg_gain
        if gain <= min_gain:
            print(f"Stopping at {len(selected)} representatives — "
                  f"top remaining gain {gain:.3f} below threshold")
            break

        new_voxels = target_voxel_sets[i] - covered
        covered |= new_voxels
        cumulative_gain += gain
        iteration += 1

        selected.append({
            'idx':                  i,
            'label':                target_labels[i],
            'gain':                 gain,
            'cumulative_weight':    cumulative_gain,
            'cumulative_fraction':  cumulative_gain / total_weight,
        })
        print(f"  [{len(selected)}/{k}] selected {target_labels[i]}: "
              f"gain={gain:.2f}, "
              f"cumulative coverage={100*cumulative_gain/total_weight:.1f}%")

    return selected

def representatives_via_coverage(cell_sequences: list,
                                 target_labels: list,
                                 k: int) -> list:
    """
    Reframe the grid-trace cell sets (not the sequences — the SETS,
    for coverage purposes) as input to weighted max coverage,
    exactly as developed earlier in this conversation.
    """
    voxel_sets = [set(seq) for seq in cell_sequences]
    weights    = compute_voxel_weights(voxel_sets)   # log(occupancy)
    return greedy_max_coverage_fast(
        voxel_sets, weights, k=k, target_labels=target_labels
    )

import numpy as np


def build_anchor_mask(
        project,
        df,
        anchor_q_start: int,
        anchor_q_end: int):
    """
    Returns
    -------
    dict {target_id: (anchor_mask, full_coords)}

        anchor_mask : (n_res,) bool

        full_coords : (n_res,3) float32
            coordinates transformed using
            df[['rotation','translation']]
    """

    segs = project.segments

    meta = (
        df[
            [
                'alignment_id',
                'target_id',
                'start_idx',
                'target_length',
                'R',
                't'
            ]
        ]
        .drop_duplicates(subset='alignment_id')
    )

    result = {}

    for _, row in meta.iterrows():

        aln_id = row['alignment_id']
        target_id = row['target_id']

        start_idx = int(row['start_idx'])
        end_idx = start_idx + int(row['target_length'])

        coords = project.coords[start_idx:end_idx]

        R = np.asarray(row['R'], dtype=np.float32)
        t = np.asarray(row['t'], dtype=np.float32)

        # x' = xR + t
        # R.T is correct for Dali, which uses Fortran indexing
        full_coords = coords @ R.T + t

        n_res = len(full_coords)

        anchor_mask = np.zeros(n_res, dtype=bool)

        aln_segs = segs[segs['alignment_id'] == aln_id]

        for _, seg in aln_segs.iterrows():

            qs = int(seg['q_start'])
            ss = int(seg['s_start'])
            l = int(seg['length'])

            q_lo = qs
            q_hi = qs + l

            overlap_lo = max(q_lo, anchor_q_start)
            overlap_hi = min(q_hi, anchor_q_end)

            if overlap_hi <= overlap_lo:
                continue

            offset_in_seg = overlap_lo - qs
            n_overlap = overlap_hi - overlap_lo

            s_lo = ss + offset_in_seg
            s_hi = s_lo + n_overlap

            s_lo = max(0, s_lo)
            s_hi = min(n_res, s_hi)

            if s_hi > s_lo:
                anchor_mask[s_lo:s_hi] = True

        result[target_id] = (
            anchor_mask,
            full_coords.astype(np.float32)
        )

    return result

def grid_trace_with_dwell(coords: np.ndarray,
                          anchor_mask: np.ndarray,
                          cell_size: float = 4.5) -> list:
    """
    Run-length-encoded cell trace: each entry is (cell_hash, dwell_count).
    Preserves dwelling time as an explicit attribute rather than
    either discarding it (collapsed) or conflating it with revisits
    (uncollapsed).

    A long dwell (e.g. 8 consecutive residues in one cell) suggests
    a helix curling tightly in place; rapid cell-to-cell movement
    with dwell=1 suggests an extended strand.
    """
    non_anchor_coords = coords[~anchor_mask]
    if len(non_anchor_coords) == 0:
        return []

    cell_keys = np.floor(non_anchor_coords / cell_size).astype(np.int64)
    OFFSET = 1 << 16
    hashed = (
        (cell_keys[:, 0] + OFFSET) * (1 << 34) +
        (cell_keys[:, 1] + OFFSET) * (1 << 17) +
        (cell_keys[:, 2] + OFFSET)
    )

    trace = []
    for h in hashed:
        if trace and trace[-1][0] == h:
            trace[-1] = (h, trace[-1][1] + 1)
        else:
            trace.append((h, 1))
    return trace

def compute_cell_occupancy(all_coords: list,
                           all_anchor_masks: list,
                           cell_size: float = 4.5) -> Counter:
    """
    First pass over ALL targets: count how many DISTINCT TARGETS
    (not atoms) visit each cell. This is occupancy in the sense of
    "how many targets agree this region exists", which is the
    right denominator for filtering — NOT atom count, which would
    be dominated by how many residues happen to sit there.

    Parameters
    ----------
    all_coords       : list of (n_res_i, 3) arrays, one per target
    all_anchor_masks : list of (n_res_i,) bool arrays, one per target
    cell_size        : Å
    """
    target_occupancy = Counter()

    for coords, mask in zip(all_coords, all_anchor_masks):
        non_anchor = coords[~mask]
        if len(non_anchor) == 0:
            continue

        cell_keys = np.floor(non_anchor / cell_size).astype(np.int64)
        OFFSET = 1 << 16
        hashed = (
            (cell_keys[:, 0] + OFFSET) * (1 << 34) +
            (cell_keys[:, 1] + OFFSET) * (1 << 17) +
            (cell_keys[:, 2] + OFFSET)
        )

        # count each cell ONCE per target, regardless of how many
        # atoms from this target sit there — that's the whole point
        unique_cells_this_target = set(hashed.tolist())
        for c in unique_cells_this_target:
            target_occupancy[c] += 1

    return target_occupancy


def filter_low_occupancy_cells(target_occupancy: Counter,
                               n_targets: int,
                               min_fraction: float = 0.05,
                               min_count: int = None) -> set:
    """
    Decide which cells survive — visited by at least min_fraction
    of all targets (or an absolute min_count, whichever you prefer
    to reason about).

    A hairball cell visited by only 1-2 targets out of thousands
    gets dropped. A cell that's part of a genuinely recurrent
    auxiliary domain, visited by even a modest fraction of targets,
    survives.

    Returns the set of cell hashes to KEEP.
    """
    if min_count is None:
        min_count = max(2, int(np.ceil(min_fraction * n_targets)))

    kept = {c for c, n in target_occupancy.items() if n >= min_count}

    print(f"Cell occupancy filter: {len(target_occupancy)} total cells, "
          f"{len(kept)} kept (visited by >= {min_count}/{n_targets} targets), "
          f"{len(target_occupancy) - len(kept)} dropped as hairball/noise")

    return kept

def grid_trace_filtered(coords: np.ndarray,
                        anchor_mask: np.ndarray,
                        kept_cells: set,
                        cell_size: float = 4.5) -> list:
    """
    Same as grid_trace, but residues landing in a low-occupancy
    (hairball) cell are dropped from the trace entirely — they
    don't even appear as a gap-causing entry, they're simply absent,
    as if that part of the chain weren't there for this analysis.
    """
    non_anchor_coords = coords[~anchor_mask]
    if len(non_anchor_coords) == 0:
        return []

    cell_keys = np.floor(non_anchor_coords / cell_size).astype(np.int64)
    OFFSET = 1 << 16
    hashed = (
        (cell_keys[:, 0] + OFFSET) * (1 << 34) +
        (cell_keys[:, 1] + OFFSET) * (1 << 17) +
        (cell_keys[:, 2] + OFFSET)
    )

    trace = []
    for h in hashed:
        h = int(h)
        if h not in kept_cells:
            continue   # hairball/noise cell — skip silently
        if trace and trace[-1][0] == h:
            trace[-1] = (h, trace[-1][1] + 1)
        else:
            trace.append((h, 1))
    return trace

def build_filtered_traces(all_coords: list,
                          all_anchor_masks: list,
                          target_labels: list,
                          cell_size: float = 4.5,
                          min_fraction: float = 0.05) -> dict:
    """
    Two-pass construction:
      Pass 1: occupancy census across ALL targets, determine
              which cells are real (recurrent) vs. hairball (rare)
      Pass 2: build each target's trace using only kept cells
    """
    n_targets = len(all_coords)

    print(f"Pass 1: computing cell occupancy across {n_targets} targets...")
    occupancy = compute_cell_occupancy(all_coords, all_anchor_masks, cell_size)
    kept_cells = filter_low_occupancy_cells(
        occupancy, n_targets, min_fraction=min_fraction
    )

    print(f"Pass 2: building filtered traces...")
    traces = {}
    for label, coords, mask in zip(target_labels, all_coords, all_anchor_masks):
        traces[label] = grid_trace_filtered(coords, mask, kept_cells, cell_size)

    print('Done')

    return traces

def weighted_jaccard_matrix_exact(all_counts: list) -> np.ndarray:
    """
    Exact weighted Jaccard (Ruzicka similarity) for all pairs,
    vectorized via sparse matrices — no Python-level pair loop.

    weighted_jaccard(i,j) = sum(min(w_i,w_j)) / sum(max(w_i,w_j))

    min/max don't factor into a single matrix product the way dot
    products do, but they CAN be computed efficiently using the
    identity:
        sum(min(a,b)) = sum(a) + sum(b) - sum(max(a,b))
    and max can be handled via a different sparse trick, OR more
    simply: compute it densely if N and the cell vocabulary are
    small enough, which is the practical path here.
    """
    from scipy.sparse import csr_matrix

    all_cells = sorted(set().union(*(c.keys() for c in all_counts)))
    cell_idx  = {c: i for i, c in enumerate(all_cells)}
    n_cells   = len(all_cells)
    N         = len(all_counts)

    print(f"{N} targets, {n_cells} distinct cells")

    # dense matrix — check memory first
    mem_gb = N * n_cells * 4 / 1e9
    print(f"Dense matrix would need ~{mem_gb:.2f} GB")

    rows, cols, vals = [], [], []
    for i, counts in enumerate(all_counts):
        for cell, w in counts.items():
            rows.append(i); cols.append(cell_idx[cell]); vals.append(w)

    W = csr_matrix((vals, (rows, cols)), shape=(N, n_cells),
                   dtype=np.float32).toarray()   # dense from here

    # exact min/max weighted Jaccard via broadcasting, chunked over i
    # to avoid materializing the full (N, N, n_cells) tensor at once
    sim = np.zeros((N, N), dtype=np.float32)
    chunk = 10

    for i0 in range(0, N, chunk):
        i1 = min(i0 + chunk, N)
        Wi = W[i0:i1, None, :]      # (chunk, 1, n_cells)
        Wj = W[None, :, :]          # (1, N, n_cells)

        mins = np.minimum(Wi, Wj).sum(axis=2)   # (chunk, N)
        maxs = np.maximum(Wi, Wj).sum(axis=2)   # (chunk, N)

        with np.errstate(divide='ignore', invalid='ignore'):
            sim[i0:i1] = np.where(maxs > 0, mins / maxs, 0.0)

        print(f"  {i1}/{N} rows done")

    return sim

def trace_to_weighted_set(trace: list) -> Counter:
    """
    Convert a run-length trace [(cell, dwell), ...] into a weighted
    occupancy Counter {cell: total_dwell} for Round 1 Jaccard.
    Revisits (same cell appearing non-adjacently) are summed.
    """
    counts = Counter()
    for cell, dwell in trace:
        counts[cell] += dwell
    return counts
