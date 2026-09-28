"""Graph construction for infomap"""

import numpy as np
from typing import Any, List, Tuple, Dict
from scipy import sparse
from collections import defaultdict, deque


"""
Tandem repeat detection from the sequentially-ordered segment adjacency
matrix.

Signal (per discussion): a genuine tandem-repeat stack shows a CONSISTENT,
LINEARLY-OCCUPIED band at a small diagonal offset p (the repeat period,
in segments) -- e.g. for period 3, (i,i+3), (i+1,i+4), (i+2,i+5), ... are
ALL present, not just occasionally -- while diagonals beyond a fixed
high-diagonal cutoff (k>3, per prior discussion) stay sparse. This is
what distinguishes an open, linear repeat stack (armadillo, TPR, LRR)
from both a closed ring (which would show a strong contact at the
wrap-around offset connecting first and last unit) and an ordinary
compact globular fold (which typically has real, scattered long-range
3D contacts that a simple periodicity check would NOT explain away).
"""

def _diagonal(binary: np.ndarray, k: int) -> np.ndarray:
    """Boolean array: entry i = whether (i, i+k) is a contact, i in [0, n-k)."""
    n = binary.shape[0]
    idx = np.arange(n - k)
    return binary[idx, idx + k]


def _high_diagonal_occupancy(binary: np.ndarray, start: int, end_exclusive: int, cutoff: int) -> float:
    """Mean occupancy over all diagonal offsets k > cutoff, restricted to
    the [start, end_exclusive) segment range."""
    hits, total = 0, 0
    n_region = end_exclusive - start
    for k in range(cutoff + 1, n_region):
        idx = np.arange(start, end_exclusive - k)
        if len(idx) == 0:
            continue
        hits += int(binary[idx, idx + k].sum())
        total += len(idx)
    return hits / total if total > 0 else 0.0


def detect_tandem_repeats(
    segment_adj: np.ndarray,
    min_period: int = 2,
    max_period: int = 3,
    high_diagonal_cutoff: int = 3,
    occupancy_threshold: float = 0.8,
    sparsity_threshold: float = 0.15,
    min_repeat_units: int = 3,
    weight_threshold: float = 0.0,
) -> List[Dict]:
    """
    segment_adj: (n_segments, n_segments) symmetric contact/weight matrix,
        segments in SEQUENTIAL chain order (essential -- diagonal offsets
        only mean "distance along the chain, in repeat units" if the
        matrix is sequence-ordered).
    min_period: smallest period tested. Defaults to 2, NOT 1 -- period=1
        just means consecutive segments touch each other, which is
        ordinary chain connectivity, not a meaningful repeat unit, and
        would otherwise trivially out-compete genuine candidate periods
        whenever local packing happens to also be dense.
    max_period: largest period tested (segments per repeat unit). Kept
        <= high_diagonal_cutoff by design, so a period's own defining
        diagonal is never itself misclassified as "should be sparse".
    high_diagonal_cutoff: diagonals k > this must stay sparse (the fixed
        k>3 boundary from prior discussion).
    occupancy_threshold: minimum fraction of (i, i+p) pairs that must be
        present within a candidate run for it to count as a genuine
        LINEAR periodic band, not just scattered low-offset contacts.
    sparsity_threshold: maximum allowed occupancy at high diagonals
        (k > high_diagonal_cutoff) within a candidate region.
    min_repeat_units: minimum number of repeat units required (run
        length / period) before a candidate is accepted as real.
    weight_threshold: a matrix entry counts as "contact present" if its
        weight exceeds this.

    Returns: list of candidate repeat regions, sorted best-first, with
    overlapping candidates resolved by greedy non-overlapping selection.
        [{'start', 'end', 'period', 'n_units',
          'low_occupancy', 'high_occupancy', 'score'}, ...]
    (start/end are inclusive segment indices.)
    """
    n = segment_adj.shape[0]
    binary = segment_adj > weight_threshold
    candidates: List[Dict] = []

    for p in range(min_period, max_period + 1):
        diag = _diagonal(binary, p)
        i = 0
        while i < len(diag):
            if not diag[i]:
                i += 1
                continue

            # Greedily extend a run from i while cumulative occupancy
            # stays >= occupancy_threshold (tolerates the occasional
            # missed contact real repeats show, rather than requiring
            # a perfectly unbroken run).
            j = i
            hits, total = 1, 1
            while j + 1 < len(diag):
                hits_next = hits + (1 if diag[j + 1] else 0)
                total_next = total + 1
                if hits_next / total_next >= occupancy_threshold:
                    j += 1
                    hits, total = hits_next, total_next
                else:
                    break

            start, end = i, j + p  # inclusive segment range
            n_units = (end - start) / p

            if n_units >= min_repeat_units:
                low_occ = hits / total
                high_occ = _high_diagonal_occupancy(binary, start, end + 1, high_diagonal_cutoff)
                if high_occ <= sparsity_threshold:
                    candidates.append({
                        "start": start, "end": end, "period": p,
                        "n_units": n_units,
                        "low_occupancy": low_occ,
                        "high_occupancy": high_occ,
                        "score": low_occ - high_occ,
                    })

            i = j + 1

    candidates.sort(key=lambda c: c["score"], reverse=True)
    selected: List[Dict] = []
    covered = np.zeros(n, dtype=bool)
    for c in candidates:
        seg_range = slice(c["start"], c["end"] + 1)
        if covered[seg_range].any():
            continue
        selected.append(c)
        covered[seg_range] = True

    return selected

def graph_triplets_to_segment_adj(graph_triplets, n_segments=None):
    """
    Convert graph_triplets = [(u, v, w), ...] into a dense symmetric
    (n_segments, n_segments) matrix, segment-index-ordered -- required
    for detect_tandem_repeats, which reads diagonal offsets as
    "distance along the chain in repeat units".

    n_segments: if None, inferred as max(u, v) + 1 across all triplets.
    Pass explicitly if the LAST segment could be isolated (e.g. no
    surviving edges) -- inference would silently undercount and drop
    it from the matrix entirely, rather than error.
    """
    if n_segments is None:
        n_segments = max(max(u, v) for u, v, _w in graph_triplets) + 1

    segment_adj = np.zeros((n_segments, n_segments))
    for u, v, w in graph_triplets:
      if u != v: # excl
        segment_adj[u, v] = w
        segment_adj[v, u] = w  # harmless no-op when u == v (self-loops)

    return segment_adj

"""
Repeat-cluster adapter + combined contraction, reusing create_cluster_graph_gaussian
directly rather than writing a parallel contraction path.

Sequencing (sheets first, matches the reasoning that H-bond evidence is
direct/unambiguous while periodicity is a derived statistical pattern
best read off the already-cleaned graph):

  1. identify_beta_sheets_and_bifurcations() -> sheet_clusters
     (NOTE: as currently written, this does NOT catch the common
     degree-2 chimeric-strand-bridging-two-sheets case -- see message
     text. That needs a separate, upstream fix before pairings/
     StrandPartner aggregation, not something this module addresses.)

  2. Build the PLAIN, segment-id-ordered graph (trivial singleton
     "clusters", i.e. no sheet collapsing) via create_cluster_graph_gaussian,
     restricted to segments NOT already claimed by a sheet cluster --
     sheet cluster IDs are not sequence-ordered, so repeat detection
     cannot run on the post-sheet-collapse graph directly.

  3. detect_tandem_repeats() on that plain, ordered graph -> repeat
     candidates, in terms of raw segment_id ranges.

  4. Convert repeat candidates into the same List[List[(seg_id, range)]]
     shape sheet clusters already use (this module's job).

  5. ONE combined call to create_cluster_graph_gaussian with
     sheet_clusters + repeat_clusters + remaining singleton segments --
     no second, parallel contraction mechanism needed.
"""

def repeat_candidates_to_clusters(
    repeat_candidates: List[Dict],
    segments: List[Tuple[str, int, int]],
) -> List[List[Tuple[int, Tuple[int, int]]]]:
    """
    Convert detect_tandem_repeats() output into the SAME
    List[List[(seg_id, (res_start, res_end))]] shape
    identify_beta_sheets_and_bifurcations() already produces, so both
    can be concatenated and passed to create_cluster_graph_gaussian
    unchanged.
 
    repeat_candidates: output of detect_tandem_repeats -- start/end are
        SEGMENT indices (inclusive) into `segments`, NOT residue indices.
    segments: the same segments list used everywhere else in the
        pipeline, (sse_type, res_start, res_end) per entry.
    """
    clusters = []
    for cand in repeat_candidates:
        member_seg_ids = range(cand["start"], cand["end"] + 1)
        cluster = [(seg_id, (segments[seg_id][1], segments[seg_id][2]))
                   for seg_id in member_seg_ids]
        clusters.append(cluster)
    return clusters

import numpy as np

def claimed_segment_ids(clusters) -> np.ndarray:
    """Extracts all unique segment IDs already claimed by prior clustering passes."""
    if not clusters:
        return np.empty(0, dtype=np.int32)

    # Extract values if dictionary, or use input sequence
    cluster_iterable = clusters.values() if isinstance(clusters, dict) else clusters

    # Extract IDs regardless of whether elements are scalars or tuples
    flat_ids = []
    for cluster in cluster_iterable:
        for item in cluster:
            flat_ids.append(
                item[0] if isinstance(item, (tuple, list)) else item
            )

    return (
        np.unique(np.array(flat_ids, dtype=np.int32))
        if flat_ids
        else np.empty(0, dtype=np.int32)
    )

def old_claimed_segment_ids(clusters) -> np.ndarray:
    """Extracts all unique segment IDs already claimed by prior clustering passes.

    Parameters
    ----------
    clusters : Dict[int, List[int] | np.ndarray] or np.ndarray (object array)
        The sheet clusters mapping cluster IDs to lists/arrays of segment IDs.

    Returns
    -------
    np.ndarray
        1D NumPy integer array of unique claimed segment IDs.
    """
    if not clusters:
        return np.empty(0, dtype=np.int32)

    # Handle dictionary of cluster_id -> list/array of segment IDs
    if isinstance(clusters, dict):
        if not clusters:
            return np.empty(0, dtype=np.int32)
        # Flatten all segment ID lists into a single vectorised sequence
        all_ids = [
            seg_id
            for seg_list in clusters.values()
            for seg_id in (
                seg_list.tolist()
                if isinstance(seg_list, np.ndarray)
                else seg_list
            )
        ]
        return np.unique(np.array(all_ids, dtype=np.int32))

    # Handle 1D NumPy object array of arrays/lists
    if isinstance(clusters, np.ndarray):
        if len(clusters) == 0:
            return np.empty(0, dtype=np.int32)
        all_ids = np.concatenate(
            [
                np.asarray(item, dtype=np.int32)
                for item in clusters
                if len(item) > 0
            ]
        )
        return np.unique(all_ids)

    return np.empty(0, dtype=np.int32)


def old_claimed_segment_ids(clusters: List[List[Tuple[int, Tuple[int, int]]]]) -> set:
    """Segment ids already grouped by a prior clustering pass (e.g. sheets)."""
    claimed = set()
    print("#claimed_segment_ids clusters", clusters)
    for cluster in clusters.keys():
        for seg_id in clusters[cluster]:
            claimed.add(seg_id)

#    for cluster in clusters:
#        for seg_id, _range in cluster:
#            claimed.add(seg_id)
    return claimed


def build_full_ordered_clusters(
    segments: Any,
) -> List[List[Tuple[int, Tuple[int, int]]]]:
    """Trivial singleton clustering over EVERY segment, in segment-id order

    -- guarantees cluster_id == seg_id exactly (no gaps, no compaction),
    which is what create_cluster_graph_gaussian needs to preserve
    sequential ordering for detect_tandem_repeats to interpret diagonal
    offsets correctly.

    DELIBERATELY does not exclude sheet-claimed segments -- excluding
    entries here causes create_cluster_graph_gaussian's own fallback
    logic (which assigns any un-clustered segment a new TRAILING
    singleton ID) to silently scramble node index vs. seg_id wherever
    an excluded segment sits anywhere but the very end of the sequence.
    Filter candidates AFTER detect_tandem_repeats runs instead (see
    filter_repeat_candidates below), never before.
    """
    if isinstance(segments, np.ndarray):
        # Extract start and end columns as integers directly from numpy array
        starts = segments[:, 1].astype(int)
        ends = segments[:, 2].astype(int)
        return [[(seg_id, (int(starts[seg_id]), int(ends[seg_id])))] for seg_id in range(len(segments))]

    return [[(seg_id, (int(seg[1]), int(seg[2])))] for seg_id, seg in enumerate(segments)]

def old_build_full_ordered_clusters(
    segments: List[Tuple[str, int, int]],
) -> List[List[Tuple[int, Tuple[int, int]]]]:
    """
    Trivial singleton clustering over EVERY segment, in segment-id order
    -- guarantees cluster_id == seg_id exactly (no gaps, no compaction),
    which is what create_cluster_graph_gaussian needs to preserve
    sequential ordering for detect_tandem_repeats to interpret diagonal
    offsets correctly.
 
    DELIBERATELY does not exclude sheet-claimed segments -- excluding
    entries here causes create_cluster_graph_gaussian's own fallback
    logic (which assigns any un-clustered segment a new TRAILING
    singleton ID) to silently scramble node index vs. seg_id wherever
    an excluded segment sits anywhere but the very end of the sequence.
    Filter candidates AFTER detect_tandem_repeats runs instead (see
    filter_repeat_candidates below), never before.
    """
    return [[(seg_id, (seg[1], seg[2]))] for seg_id, seg in enumerate(segments)]


def filter_repeat_candidates(
    repeat_candidates: List[Dict],
    sheet_claimed: set,
) -> List[Dict]:
    """
    Reject any repeat candidate whose segment range overlaps a segment
    already claimed by sheet detection, even partially -- sheet
    membership is direct H-bond evidence and takes priority; a partial
    overlap would be ambiguous, not something to silently resolve.
    """
    kept = []
    for cand in repeat_candidates:
        span = range(cand["start"], cand["end"] + 1)
        if any(seg_id in sheet_claimed for seg_id in span):
            continue
        kept.append(cand)
    return kept

from typing import Dict


def assert_clusters_disjoint(*cluster_lists) -> None:
    """Defensive check: ensures no segment ID appears in multiple cluster lists."""
    seen: Dict[int, int] = {}

    for list_idx, cluster_list in enumerate(cluster_lists):
        clusters = (
            cluster_list.values()
            if isinstance(cluster_list, dict)
            else cluster_list
        )

        for cluster in clusters:
            for item in cluster:
                seg_id = int(
                    item[0] if isinstance(item, (tuple, list)) else item
                )

                if seg_id in seen:
                    raise ValueError(
                        f"Segment {seg_id} claimed by both cluster list "
                        f"{seen[seg_id]} and cluster list {list_idx} -- "
                        f"would silently corrupt segment_to_cluster mapping."
                    )
                seen[seg_id] = list_idx

def old_assert_clusters_disjoint(
    *cluster_lists: List[List[Tuple[int, Tuple[int, int]]]],
) -> None:
    """
    Defensive check: no segment id may appear in more than one cluster
    across the combined lists. create_cluster_graph_gaussian's
    segment_to_cluster dict silently OVERWRITES on a duplicate key with
    no warning -- a segment could get silently reassigned from its
    (authoritative) sheet cluster to a (less certain) repeat cluster
    with zero indication anything went wrong. Call this before the
    final combined create_cluster_graph_gaussian call, every time.
    """
    seen: Dict[int, int] = {}
    print("# cluster_lists", cluster_lists)
    for list_idx, cluster_list in enumerate(cluster_lists):
        # Handle dict vs list/tuple for cluster_list
        clusters = cluster_list.values() if isinstance(cluster_list, dict) else cluster_list

        for cluster in clusters:
            # Handle if cluster is directly a list/set of IDs or a list of tuples like (seg_id, _range)
            for item in cluster:
                seg_id = item[0] if isinstance(item, (tuple, list)) else item

                if seg_id in seen:
                    raise ValueError(
                        f"segment {seg_id} claimed by both cluster list "
                        f"{seen[seg_id]} and cluster list {list_idx} -- "
                        f"would silently corrupt segment_to_cluster"
                    )
                seen[seg_id] = list_idx

from typing import Any, List, Tuple, Union, Dict, Set

from typing import Any, List, Set


def merge_cluster_lists(
    *cluster_lists: Any,
    n_segments: int,
) -> List[Any]:
    """Concatenates all cluster lists and appends remaining unclaimed segment IDs

    as singletons up to n_segments.
    """
    merged: List[Any] = []
    claimed: Set[int] = set()

    for cluster_list in cluster_lists:
        clusters = (
            cluster_list.values()
            if isinstance(cluster_list, dict)
            else cluster_list
        )

        for cluster in clusters:
            merged.append(cluster)
            for item in cluster:
                seg_id = int(
                    item[0] if isinstance(item, (tuple, list)) else item
                )
                claimed.add(seg_id)

    # Append unclaimed segment IDs as singletons
#    for seg_id in range(n_segments):
#        if seg_id not in claimed:
#            merged.append([(seg_id, (0, 0))])
#            claimed.add(seg_id)

    return merged

def old_merge_cluster_lists(
    *cluster_lists: List[List[Tuple[int, Tuple[int, int]]]],
    n_segments: int,
) -> List[List[Tuple[int, Tuple[int, int]]]]:
    """
    Concatenate sheet clusters + repeat clusters + anything still
    unclaimed, as singletons -- the combined list create_cluster_graph_gaussian
    should get for the ONE final contraction call.
    """
    print('#merge_cluster_lists', cluster_lists)
    merged = []
    claimed = set()
    for cluster_list in cluster_lists:
        for cluster in cluster_list:
            merged.append(cluster)
            for seg_id, _range in cluster:
                claimed.add(seg_id)
    return merged

# ==============================================================================
# 2. Beta-Sheet & Bifurcation Identification (with Spoke Pruning)
# ==============================================================================

def identify_beta_sheets_and_bifurcations(secondary_struct_tuple, segments):
    """
    Identifies beta sheets and bifurcated strands. 
    Prunes degree-1 spokes connected to degree > 2 hubs before clustering.
    Returns clusters as lists of (segment_id, (start, end)).
    """
    _, strand_ranges, pairings = secondary_struct_tuple

    # Map strand range to 0-based index in `segments`
    range_to_seg_id = {}
    for seg_id, seg in enumerate(segments):
        range_to_seg_id[(seg[1], seg[2])] = seg_id

    for idx, strand in enumerate(strand_ranges):
        if strand not in range_to_seg_id:
            range_to_seg_id[strand] = idx

    # Build adjacency list
    adj = defaultdict(set)
    for partner in pairings:
        u, v = partner.strand_a, partner.strand_b
        adj[u].add(v)
        adj[v].add(u)

    for strand in strand_ranges:
        if strand not in adj:
            adj[strand] = set()

    # Graph Cleaning: Iteratively prune degree-1 spokes connected to hubs (deg > 2)
    modified = True
    while modified:
        modified = False
        edges_to_remove = []
        for u in list(adj.keys()):
            if len(adj[u]) == 1:
                v = next(iter(adj[u]))
                if len(adj[v]) > 2:
                    edges_to_remove.append((u, v))

        if edges_to_remove:
            for u, v in edges_to_remove:
                if u in adj[v]:
                    adj[u].remove(v)
                    adj[v].remove(u)
                    modified = True

    # Classify degrees post-pruning
    bifurcated_strands = set()
    normal_strands = set()

    for strand, neighbors in adj.items():
        if len(neighbors) > 2:
            bifurcated_strands.add(strand)
        else:
            normal_strands.add(strand)

    clusters = []

    # Bifurcated strands -> Singletons
    for strand in sorted(bifurcated_strands, key=lambda s: range_to_seg_id[s]):
        seg_id = range_to_seg_id[strand]
        clusters.append([(seg_id, strand)])

    # Connected components (BFS) for non-bifurcated strands
    visited = set()
    for strand in sorted(normal_strands, key=lambda s: range_to_seg_id[s]):
        if strand not in visited:
            component = []
            queue = deque([strand])
            visited.add(strand)

            while queue:
                curr = queue.popleft()
                component.append((range_to_seg_id[curr], curr))

                for neighbor in adj[curr]:
                    if neighbor in normal_strands and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

            clusters.append(sorted(component, key=lambda x: x[0]))

    return clusters


# ==============================================================================
# 3. Cluster-Level Gaussian Graph Construction (with Self-Loops)
# ==============================================================================

from typing import Any, List, Tuple
import numpy as np
from scipy import sparse


def create_cluster_graph_gaussian(
    n_residues: int,
    adj: sparse.coo_matrix,
    segments: Any,
    clusters: Any,
    sigma: float = 3.0,
    diag: int = 1,
) -> List[Tuple[int, int, float]]:
    """Converts residue distance COO matrix to SSE cluster graph with Gaussian weights.

    Includes (diag=0) / Excludes (diag=1, default) intra-cluster self-loops (u, u).
    """
    n_segs = len(segments) if segments is not None else 0
    if n_segs == 0 or not clusters or adj.nnz == 0:
        return []

    # 1. Normalize cluster iteration (handles dict, list, or numpy object array)
    if isinstance(clusters, dict):
        cluster_iterable = list(clusters.values())
    else:
        cluster_iterable = list(clusters)

    n_clusters_input = len(cluster_iterable)

    # Array mapping segment_id -> cluster_id, initialized to -1
    seg_to_cluster = np.full(n_segs, -1, dtype=np.int32)

    # 2. Vectorised mapping of segment IDs to cluster IDs
    for cluster_id, cluster in enumerate(cluster_iterable):
        for item in cluster:
            seg_id = int(item[0] if isinstance(item, (tuple, list)) else item)
            if 0 <= seg_id < n_segs:
                seg_to_cluster[seg_id] = cluster_id

    # Assign unclaimed non-strand segments sequential singleton cluster IDs
    unclaimed = np.flatnonzero(seg_to_cluster == -1)
    if len(unclaimed) > 0:
        seg_to_cluster[unclaimed] = np.arange(
            n_clusters_input, n_clusters_input + len(unclaimed), dtype=np.int32
        )

    num_clusters = n_clusters_input + len(unclaimed)

    # 3. Vectorised residue_cluster_map construction
    residue_cluster_map = np.full(n_residues, -1, dtype=np.int32)

    if isinstance(segments, np.ndarray):
        starts = segments[:, 1].astype(np.int32)
        ends = segments[:, 2].astype(np.int32)
    else:
        starts = np.fromiter((int(s[1]) for s in segments), dtype=np.int32)
        ends = np.fromiter((int(s[2]) for s in segments), dtype=np.int32)

    # Slice residues for each segment into the residue_cluster_map
    for seg_id in range(n_segs):
        c_id = seg_to_cluster[seg_id]
        s_res, e_res = starts[seg_id], ends[seg_id]
        residue_cluster_map[s_res : e_res + 1] = c_id

    # 4. COO Upper Triangle Slicing & Masking
    adj_upper = sparse.triu(adj, k=diag).tocoo()

    u = residue_cluster_map[adj_upper.row]
    v = residue_cluster_map[adj_upper.col]

    # Valid mask ensuring both residues belong to clusters and distance > 0
    valid_mask = (u != -1) & (v != -1) & (adj_upper.data > 0)
    if not np.any(valid_mask):
        return []

    u_valid = u[valid_mask]
    v_valid = v[valid_mask]
    dists = adj_upper.data[valid_mask]

    # Gaussian weighting
    weights = np.exp(-0.5 * ((dists / sigma) ** 2))

    # 5. Fast Array Aggregation using 64-bit bincount keys
    min_node = np.minimum(u_valid, v_valid).astype(np.int64)
    max_node = np.maximum(u_valid, v_valid).astype(np.int64)

    keys = min_node * num_clusters + max_node
    summed_weights = np.bincount(
        keys, weights=weights, minlength=num_clusters * num_clusters
    )

    nonzero_keys = np.flatnonzero(summed_weights)
    if len(nonzero_keys) == 0:
        return []

    cluster_u = (nonzero_keys // num_clusters).astype(np.int32)
    cluster_v = (nonzero_keys % num_clusters).astype(np.int32)
    cluster_w = summed_weights[nonzero_keys]

    return list(
        zip(cluster_u.tolist(), cluster_v.tolist(), cluster_w.tolist())
    )

def old_create_cluster_graph_gaussian(
    n_residues: int,
    adj: sparse.coo_matrix,
    segments: List[Tuple[str, int, int]],
    clusters: List[List[Tuple[int, Tuple[int, int]]]],
    sigma: float = 3.0,
    diag: int = 1,
) -> List[Tuple[int, int, float]]:
    """
    Converts residue distance COO matrix to SSE cluster graph with Gaussian weights.
    Includes (diag=0) / Excludes (diag=1, default) intra-cluster self-loops (u, u).
    """
    if not segments or not clusters or adj.nnz == 0:
        return []

    # Map segments to cluster IDs
    print("#clusters",clusters)
    segment_to_cluster = {}
    for cluster_id, cluster in enumerate(clusters):
        for item in cluster:
            # Extract seg_id whether item is a tuple like (0, (0, 0)) or a plain integer
            seg_id = item[0] if isinstance(item, (tuple, list)) else item
            segment_to_cluster[int(seg_id)] = cluster_id

    # Non-strand segments (helices/turns/coils) get unique singleton cluster IDs
    current_cluster_id = len(clusters)
    for seg_id in range(len(segments)):
        if seg_id not in segment_to_cluster:
            segment_to_cluster[seg_id] = current_cluster_id
            current_cluster_id += 1

    num_clusters = current_cluster_id

    # Map residue index -> cluster ID
    residue_cluster_map = np.full(n_residues, -1, dtype=np.int32)
    for seg_id, seg in enumerate(segments):
        cluster_id = segment_to_cluster[seg_id]
        start_res, end_res = seg[1], seg[2]
        residue_cluster_map[start_res : end_res + 1] = cluster_id

    # Upper triangle including diagonal (k=0)
    adj_upper = sparse.triu(adj, k=diag).tocoo()

    u = residue_cluster_map[adj_upper.row]
    v = residue_cluster_map[adj_upper.col]

    # Valid mask allowing intra-cluster (u == v)
    valid_mask = (u != -1) & (v != -1) & (adj_upper.data > 0)
    if not np.any(valid_mask):
        return []

    u_valid = u[valid_mask]
    v_valid = v[valid_mask]
    dists = adj_upper.data[valid_mask]

    weights = np.exp(-(dists ** 2) / (2.0 * (sigma ** 2)))

    # Aggregate weights per cluster pair
    min_node = np.minimum(u_valid, v_valid)
    max_node = np.maximum(u_valid, v_valid)

    keys = min_node * num_clusters + max_node
    summed_weights = np.bincount(keys, weights=weights, minlength=num_clusters * num_clusters)

    nonzero_keys = np.nonzero(summed_weights)[0]
    cluster_u = (nonzero_keys // num_clusters).astype(int)
    cluster_v = (nonzero_keys % num_clusters).astype(int)
    cluster_w = summed_weights[nonzero_keys]

    return list(zip(cluster_u.tolist(), cluster_v.tolist(), cluster_w.tolist()))
