from infomap import Infomap
import numpy as np
from typing import List, Tuple, Dict
from scipy.spatial import cKDTree

"""
SASA-proxy agglomerative refinement of Infomap's final residue assignments.

Replaces the min-density background model entirely -- that rule was
tuned to fix one specific asymmetric failure (a dense beta-core
inflating expectations against a looser helical domain) and generalized
poorly to the symmetric case (a beta sheet split into two fragments that
are EACH individually dense, so min(density_A, density_B) stays high on
both sides and the real interface has to clear an artificially high bar).

This version uses ONE Delaunay-derived signal throughout, for both the
merge-candidate prefilter and the actual scoring -- no separate 8A-cutoff
contact matrix living in parallel with the real edge weights, which was
the concrete bug in the previous version (adj_csr was computed and used
ONLY for the nnz>0 prefilter, then discarded in favor of a cruder,
inconsistent contact definition for the actual scoring).

SASA proxy: bounded soft occupancy from Delaunay edge lengths (Cα-only),
following the half-sphere-exposure precedent (Hamelryck 2005) that this
class of proxy is a real, validated technique, not a compromise --
though note HSE itself uses DIRECTIONAL NEIGHBOR COUNT, not mean
distance alone; this implementation follows the mean-edge-length
formulation as specified, with degree/count-based refinement flagged
as a known, easy extension point if empirical validation shows it's
needed (see docstring note in compute_residue_occupancy).
"""

from typing import Dict, List, Tuple, Union

import numpy as np
from scipy.sparse import coo_matrix, spmatrix
from scipy.spatial.distance import cdist

# OPTIMIZE TREE PARTITION

"""
Dynamic-programming tree partitioning of Infomap's hierarchy.

Replaces the global-depth sweep with a proper optimal-substructure DP,
allowing DIFFERENT effective cut depths on different branches -- exactly
what Infomap's own varying path lengths were already hinting at (a
branch whose path terminates early is a branch Infomap itself judged
didn't need further splitting).

Correctness: multi-way normalized cut is additive over the chosen
partition -- Ncut(partition) = sum over clusters C of cut(C,rest)/assoc(C)
-- which is exactly the property a tree DP needs (optimal substructure +
additive cost). For every node v in the hierarchy tree:

    DP(v) = min(
        cost_as_one_cluster(v),
        sum(DP(child) for child in children(v)) + alpha
    )

`alpha` is a fixed complexity penalty paid once per split decision
(same spirit as CART cost-complexity pruning) -- a split is only taken
if the Ncut improvement it buys exceeds alpha, replacing the earlier
fixed accept/reject threshold with a genuine cost-complexity tradeoff
baked into the optimization itself, not a post-hoc filter.
"""

from typing import Dict, List, Tuple, Any

import numpy as np
import scipy.sparse as sp


def build_tree(hierarchy: Dict[Any, Tuple[int, ...]]):
    """
    Build the hierarchy tree from path-prefix relationships.
    Returns:
        children: dict, path-prefix -> list of immediate child prefixes
        leaves_under: dict, path-prefix -> set of original node ids
                      whose path has this prefix (i.e. all leaves in
                      that subtree)
        root: the top-level prefix, ()
    """
    all_prefixes = set()
    leaves_under: Dict[Tuple[int, ...], set] = {}

    for node, path in hierarchy.items():
        for depth in range(len(path) + 1):
            prefix = path[:depth]
            all_prefixes.add(prefix)
            leaves_under.setdefault(prefix, set()).add(node)

    children: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = {p: [] for p in all_prefixes}
    for p in all_prefixes:
        if len(p) == 0:
            continue
        parent = p[:-1]
        if p not in children[parent]:
            children[parent].append(p)

    return children, leaves_under, ()



def contribution_as_one_cluster(leaf_set, adj_matrix, node_index, m, degree) -> float:
    """
    Newman modularity contribution of treating leaf_set as ONE cluster:
    internal_weight (double-counted, full symmetric submatrix) minus the
    null-model expectation. Unlike Ncut, the whole-graph-as-one-cluster
    case sits at a NEUTRAL baseline (0), not an unbeatable minimum --
    real splits can genuinely IMPROVE on this, spurious ones make it worse.
    """
    idx = np.array([node_index[n] for n in leaf_set])
    mask = np.zeros(adj_matrix.shape[0], dtype=bool)
    mask[idx] = True
    internal = adj_matrix[mask][:, mask].sum()  # double-counted, matches standard formula
    deg_sum = degree[mask].sum()
    return float(internal - (deg_sum ** 2) / (2.0 * m))


def _edge_list_to_matrix(edges, nodes):
    """Convert undirected (source, target, weight) records to a symmetric matrix."""
    order = sorted(nodes)
    lookup = {node: index for index, node in enumerate(order)}
    matrix = np.zeros((len(order), len(order)), dtype=float)
    for source, target, weight in edges:
        i, j = lookup[source], lookup[target]
        matrix[i, j] += float(weight)
        if i != j:
            matrix[j, i] += float(weight)
    return matrix, order


def dp_optimal_partition(
    hierarchy: Dict[Any, Tuple[int, ...]],
    adj,
    alpha: float = 0.0,
    min_cluster_size: int = 2,
):
    """
    Returns:
        partition: dict node_id -> chosen cluster label (a path prefix,
            possibly at DIFFERENT depths for different branches)
        total_modularity: the optimal DP value at the root (in modularity
            units -- higher is better; 0 = no better than random null model)
        choices: dict path-prefix -> "leaf" | "kept" | "split"

    alpha: optional flat penalty subtracted per split decision, for extra
        parsimony beyond what the null-model correction already provides
        on its own (modularity already penalizes weak/small clusters
        implicitly; alpha=0 is a reasonable starting default, unlike the
        Ncut-based version where a nonzero threshold was load-bearing).
    """
    node_order = sorted(hierarchy.keys())
    if not (hasattr(adj, "shape") or sp.issparse(adj)):
        # Accept an undirected weighted edge list without an external module.
        adj, matrix_node_order = _edge_list_to_matrix(adj, hierarchy.keys())
        assert matrix_node_order == node_order

    if sp.issparse(adj):
        adj = adj.toarray()

    node_index = {n: i for i, n in enumerate(node_order)}
    children, leaves_under, root = build_tree(hierarchy)

    m = 0.5 * adj.sum()
    degree = adj.sum(axis=1)

    all_prefixes = sorted(leaves_under.keys(), key=len, reverse=True)

    dp_value: Dict[Tuple[int, ...], float] = {}
    choice: Dict[Tuple[int, ...], str] = {}

    for p in all_prefixes:
        leaf_set = leaves_under[p]
        kid_list = children.get(p, [])

        contrib_here = contribution_as_one_cluster(leaf_set, adj, node_index, m, degree)

        if not kid_list or len(leaf_set) < min_cluster_size:
            dp_value[p] = contrib_here
            choice[p] = "leaf"
            continue

        split_value = sum(dp_value[c] for c in kid_list) - alpha

        if split_value > contrib_here:
            dp_value[p] = split_value
            choice[p] = "split"
        else:
            dp_value[p] = contrib_here
            choice[p] = "kept"

        print(p, choice[p], split_value, contrib_here)

    chosen_clusters: List[Tuple[int, ...]] = []

    def backtrack(p):
        if choice[p] in ("leaf", "kept"):
            chosen_clusters.append(p)
        else:
            for c in children[p]:
                backtrack(c)

    backtrack(root)

    partition = {}
    for cluster_prefix in chosen_clusters:
        for node in leaves_under[cluster_prefix]:
            partition[node] = cluster_prefix

    return partition, dp_value[root] / (2.0 * m), choice

def map_nodes_to_selected_clusters(hierarchy: dict, choices: dict) -> dict:
    """
    Maps each node to its deepest ancestor marked as 'kept'.
    
    Parameters:
    -----------
    hierarchy : dict
        Mapping from node index to path tuple, e.g., {25: (2, 2, 1), ...}
    choices : dict
        Mapping from path tuple to decision string ('kept' or 'leaf').
        
    Returns:
    --------
    dict
        Mapping from node index to chosen cluster tuple.
    """
    node_to_cluster = {}

    for node, path in hierarchy.items():
        matched_cluster = None

        # Check prefixes bottom-up (from full path down to shortest)
        # e.g., for (2, 2, 1): checks (2, 2, 1), then (2, 2), then (2,)
        for i in range(len(path), 0, -1):
            prefix = path[:i]
            if choices.get(prefix) == 'kept':
                matched_cluster = prefix
                break  # Found the deepest 'kept' decision

        # Fallback if no ancestor was marked 'kept'
        if matched_cluster is None:
            matched_cluster = path

        node_to_cluster[node] = matched_cluster

    return node_to_cluster

from typing import Any, Dict, List, Tuple
import numpy as np


def get_core_communities_from_pruned_tree(
    hierarchy: Dict[int, Tuple[int, ...]],
    choices: Dict[Tuple[int, ...], str],
    segments: Any,
    clusters: Any,
) -> Dict[int, np.ndarray]:
    """Extracts pruned tree cluster assignments, converts path tuples into

    numeric cluster indices (0, 1, 2...), and expands SSE segments into
    individual residue indices for py3Dmol visualization.

    Parameters
    ----------
    hierarchy : Dict[int, Tuple[int, ...]]
        Mapping from Infomap graph node_id (segment/sheet cluster ID) to path
        tuple.
    choices : Dict[Tuple[int, ...], str]
        Optimization choices dictionary mapping path tuples to 'kept' or
        'leaf'.
    segments : Any (np.ndarray of shape (N, 3) or List[Tuple[str, int, int]])
        List/array of SSE segments as (sse_type, start_residue, end_residue).
    clusters : Any (Dict[int, np.ndarray], List[np.ndarray], or np.ndarray)
        Sheet/strand clusters mapping cluster ID to segment IDs.

    Returns
    -------
    Dict[int, np.ndarray]
        Mapping from numeric cluster ID (0, 1, 2...) to 1D int32 array of
        residue indices.
    """
    n_segs = len(segments) if segments is not None else 0
    if n_segs == 0 or not hierarchy:
        return {}

    # 1. Map graph node IDs (Infomap cluster IDs) to SSE segment IDs
    cluster_to_segments: Dict[int, List[int]] = {}

    # Handle dictionary of cluster_id -> array/list of segment IDs
    if isinstance(clusters, dict):
        cluster_iterable = list(clusters.values())
    else:
        cluster_iterable = list(clusters)

    for cluster_id, cluster in enumerate(cluster_iterable):
        seg_ids = []
        for item in cluster:
            # Safely handle both scalar seg_id (np.int64/int) and tuple (seg_id, range)
            seg_id = int(item[0] if isinstance(item, (tuple, list)) else item)
            seg_ids.append(seg_id)
        cluster_to_segments[cluster_id] = seg_ids

    # Non-strand singleton clusters (helices, turns, coils)
    current_cluster_id = len(cluster_iterable)
    claimed_segs = {
        seg_id
        for segs in cluster_to_segments.values()
        for seg_id in segs
    }

    for seg_id in range(n_segs):
        if seg_id not in claimed_segs:
            cluster_to_segments[current_cluster_id] = [seg_id]
            current_cluster_id += 1

    # 2. Derive node_id to selected cluster path mapping bottom-up
    node_to_selected_cluster: Dict[int, Tuple[int, ...]] = {}
    for node_id, path in hierarchy.items():
        matched_cluster = None
        # Bottom-up search for the deepest 'kept' ancestor
        for i in range(len(path), 0, -1):
            prefix = path[:i]
            if choices.get(prefix) == "kept":
                matched_cluster = prefix
                break

        if matched_cluster is None:
            matched_cluster = path

        node_to_selected_cluster[node_id] = matched_cluster

    # 3. Create a stable mapping from cluster path tuples to integer IDs (0, 1, 2...)
    unique_selected_clusters = sorted(
        list(set(node_to_selected_cluster.values()))
    )
    cluster_tuple_to_int_id = {
        tuple_path: int_id
        for int_id, tuple_path in enumerate(unique_selected_clusters)
    }

    # Extract start and end boundaries for vectorised residue expansion
    if isinstance(segments, np.ndarray):
        starts = segments[:, 1].astype(np.int32)
        ends = segments[:, 2].astype(np.int32)
    else:
        starts = np.fromiter((int(s[1]) for s in segments), dtype=np.int32)
        ends = np.fromiter((int(s[2]) for s in segments), dtype=np.int32)

    # 4. Aggregate residue indices under numeric cluster IDs
    raw_core_communities: Dict[int, List[int]] = {}

    for node_id, chosen_cluster_tuple in node_to_selected_cluster.items():
        numeric_cluster_id = cluster_tuple_to_int_id[chosen_cluster_tuple]

        if numeric_cluster_id not in raw_core_communities:
            raw_core_communities[numeric_cluster_id] = []

        seg_ids = cluster_to_segments.get(node_id, [])
        for seg_id in seg_ids:
            if 0 <= seg_id < n_segs:
                start_res, end_res = starts[seg_id], ends[seg_id]
                raw_core_communities[numeric_cluster_id].extend(
                    range(start_res, end_res + 1)
                )

    # Convert residue index lists to sorted, unique 1D NumPy int32 arrays
    opt_core_communities: Dict[int, np.ndarray] = {}
    for cid, res_list in raw_core_communities.items():
        if res_list:
            opt_core_communities[cid] = np.unique(
                np.array(res_list, dtype=np.int32)
            )
        else:
            opt_core_communities[cid] = np.empty(0, dtype=np.int32)

    return opt_core_communities


def patch_unstructured_residues(
    final_residue_assignments: Union[List[int], np.ndarray],
    opt_domains: Union[List[int], np.ndarray]
) -> List[int]:
    """
    Patches opt_domains by carrying over -1 assignments (unstructured/unassigned residues)
    from final_residue_assignments, and re-indexes remaining clusters cleanly.
    """
    patched = list(opt_domains)

    for i, orig_idx in enumerate(final_residue_assignments):
        if orig_idx == -1:
            patched[i] = -1

    unique_clusters = sorted(set(c for c in patched if c != -1))
    cluster_mapping = {old_id: new_id for new_id, old_id in enumerate(unique_clusters)}

    final_patched = [
        cluster_mapping[c] if c != -1 else -1
        for c in patched
    ]

    return final_patched

# ---------------------------------------------------------------------------
# Delaunay adjacency: build once, use everywhere (prefilter AND scoring)
# ---------------------------------------------------------------------------

def build_delaunay_neighbors(
    adj: Union[spmatrix, List[Tuple[Tuple[int, int], float]]],
    n_residues: int,
) -> List[set]:
    """
    adj: scipy.sparse matrix (coo_matrix, or anything with .tocoo()) --
    Delaunay-derived weighted adjacency, OR (for backwards compat) a
    [((u, v), w), ...] edge list.

    Returns: list of length n_residues, each entry a set of Delaunay-
    connected neighbor residue indices. Ignores edge weights/`.data`
    deliberately -- real Euclidean distances are recomputed from
    ca_coords for the occupancy calculation, rather than assuming a
    particular meaning for the stored weight (raw distance vs.
    Gaussian-kernel-transformed vs. other -- given this pipeline's
    adjacency is built via a Gaussian-weighted construction elsewhere,
    .data here is very likely NOT a raw distance, which makes
    recomputing from ca_coords the correct choice, not just a cautious one).

    Symmetric assumed; works whether the matrix stores each edge once
    (upper triangle) or twice (full symmetric) -- sets make this
    idempotent either way.
    """
    neighbors = [set() for _ in range(n_residues)]

    if isinstance(adj, spmatrix) or hasattr(adj, "tocoo"):
        coo = adj.tocoo()
        for u, v in zip(coo.row, coo.col):
            u, v = int(u), int(v)
            if u == v:
                continue  # skip self-loops if any are present
            neighbors[u].add(v)
            neighbors[v].add(u)
    else:
        for (u, v), _w in adj:
            neighbors[u].add(v)
            neighbors[v].add(u)

    return neighbors


# ---------------------------------------------------------------------------
# SASA proxy from Delaunay edge lengths
# ---------------------------------------------------------------------------

def _capped_neighbor_distances(
    i: int,
    neighbors: set,
    ca_coords: np.ndarray,
    d_max: float,
    min_seq_dist: int,
) -> np.ndarray:
    """
    Distances from residue i to its Delaunay neighbors, EXCLUDING
    trivially-close sequence neighbors (peptide-bond adjacency, which
    is always short regardless of true 3D burial) and CAPPING anything
    beyond d_max (convex-hull tessellation artifacts -- long edges
    represent empty solvent space, not real packing).
    """
    valid = [j for j in neighbors if abs(j - i) >= min_seq_dist]
    if not valid:
        return np.array([])
    d = np.linalg.norm(ca_coords[valid] - ca_coords[i], axis=1)
    return np.minimum(d, d_max)


def compute_residue_occupancy(
    ca_coords: np.ndarray,
    neighbors: List[set],
    d_core: float = 4.5,
    d_max: float = 10.0,
    min_seq_dist: int = 4,
) -> np.ndarray:
    """
    Per-residue occupancy phi_i in [0, 1]: 1 = fully buried (core-like
    packing), 0 = fully exposed. SASA_i = 1 - phi_i.

    NOTE: uses mean neighbor distance only, matching the originally
    specified formula. Real HSE-style measures also use DIRECTIONAL
    NEIGHBOR COUNT (how many neighbors, not just how close on average)
    -- two residues with identical mean distance but very different
    neighbor counts are NOT equally buried in reality. If validation
    against real structures shows this matters, the fix is a simple
    combination of degree and mean distance, not a different objective.
    """
    n = len(ca_coords)
    phi = np.zeros(n)
    for i in range(n):
        d = _capped_neighbor_distances(i, neighbors[i], ca_coords, d_max, min_seq_dist)
        if len(d) == 0:
            phi[i] = 0.0  # no qualifying neighbors -- treat as fully exposed
            continue
        d_mean = d.mean()
        phi[i] = np.clip((d_max - d_mean) / (d_max - d_core), 0.0, 1.0)
    return phi


def _burial_contribution(d_ij: float, d_core: float, d_max: float) -> float:
    """Same clamp/linear transform as residue occupancy, applied to one
    Delaunay edge: short edge (near d_core) -> strong burial contribution
    (~1); long edge (near/beyond d_max) -> negligible (~0)."""
    return float(np.clip((d_max - d_ij) / (d_max - d_core), 0.0, 1.0))


def compute_interface_burial(
    res_a: List[int],
    res_b: List[int],
    ca_coords: np.ndarray,
    neighbors: List[set],
    d_core: float = 4.5,
    d_max: float = 10.0,
    min_seq_dist: int = 4,
) -> Tuple[float, float]:
    """
    "Buried surface area upon merger" proxy for candidate domains A, B --
    computed directly from the EXISTING Delaunay cross-edges between
    them, no re-tessellation needed. Sequence-adjacent cross-edges are
    excluded for the same reason as within-residue occupancy: a
    peptide-bond boundary between two domains is always short
    regardless of real 3D packing, and would otherwise inflate the
    score for every sequence-adjacent domain pair trivially.

    Returns: (raw_score, normalized_score). Normalization uses the
    geometric mean of the two candidate sizes (sqrt(len_a * len_b)) so
    a large domain absorbing a small fragment isn't penalized purely by
    size mismatch -- WORTH CALIBRATING against real validation data,
    this is the single most likely parameter to need adjustment.
    """
    b_set = set(res_b)
    raw = 0.0
    for i in res_a:
        for j in neighbors[i]:
            if j in b_set and abs(i - j) >= min_seq_dist:
                d_ij = np.linalg.norm(ca_coords[i] - ca_coords[j])
                raw += _burial_contribution(d_ij, d_core, d_max)

    denom = np.sqrt(len(res_a) * len(res_b))
    normalized = raw / denom if denom > 0 else 0.0
    return float(raw), float(normalized)


# ---------------------------------------------------------------------------
# Agglomerative refinement loop
# ---------------------------------------------------------------------------

def agglomerate_infomap_clusters(
    ca_coords: np.ndarray,
    initial_partition: Union[List[int], Dict[int, List[int]], np.ndarray],
    adj: Union[coo_matrix, spmatrix],
    d_core: float = 4.5,
    d_max: float = 10.0,
    min_seq_dist: int = 4,
    burial_threshold: float = 0.5,
) -> Tuple[List[int], float]:
    """
    Refine Infomap's final residue assignments by greedily merging the
    pair of domains with the strongest SASA-proxy interface burial,
    stopping once no remaining candidate pair clears `burial_threshold`
    (normalized score). No density optimization anywhere in this path.

    burial_threshold: WORTH CALIBRATING empirically -- this replaces
    both the old gamma and min-density-floor parameters with a single,
    more directly interpretable knob (how strong must the real
    3D-proximity interface be, in occupancy units, to justify merging).
    """
    n_residues = len(ca_coords)
    neighbors = build_delaunay_neighbors(adj, n_residues)

    if isinstance(initial_partition, (list, np.ndarray)):
        partition_arr = np.array(initial_partition)
        domains = {
            int(c): np.where(partition_arr == c)[0].tolist()
            for c in np.unique(partition_arr)
        }
    else:
        domains = {c: list(res) for c, res in initial_partition.items()}

    while len(domains) > 1:
        best_pair = None
        best_score = burial_threshold  # only accept merges clearing this bar

        domain_ids = list(domains.keys())
        for i in range(len(domain_ids)):
            id_a = domain_ids[i]
            res_a = domains[id_a]
            for j in range(i + 1, len(domain_ids)):
                id_b = domain_ids[j]
                res_b = domains[id_b]

                # SAME Delaunay adjacency used for prefilter and scoring --
                # no separate contact definition living in parallel.
                has_any_edge = any(
                    n in set(res_b) for i_res in res_a for n in neighbors[i_res]
                )
                if not has_any_edge:
                    continue

                _raw, normalized = compute_interface_burial(
                    res_a, res_b, ca_coords, neighbors, d_core, d_max, min_seq_dist
                )
                if normalized > best_score:
                    best_score = normalized
                    best_pair = (id_a, id_b)

        if best_pair is None:
            break

        id_a, id_b = best_pair
        print("# merge", best_pair, best_score)
        domains[id_a] = sorted(domains[id_a] + domains[id_b])
        del domains[id_b]

    final_domains = [-1] * n_residues
    for new_id, (_, res_list) in enumerate(domains.items()):
        for res_idx in res_list:
            final_domains[res_idx] = new_id

    quality = compute_sasa_domain_score(domains, ca_coords, neighbors, d_core, d_max, min_seq_dist)
    print("# partition quality score", round(quality,4))
    return final_domains, round(quality, 4)


def compute_sasa_domain_score(
    domains: Dict[int, List[int]],
    ca_coords: np.ndarray,
    neighbors: List[set],
    d_core: float = 4.5,
    d_max: float = 10.0,
    min_seq_dist: int = 4,
) -> float:
    """
    Reporting/diagnostic score for a final partition: mean per-domain
    SASA-proxy per residue (lower = more compact/globular, consistent
    with the ~0.2-0.4 compact vs ~0.7-1.0 extended range from the
    original design notes). Purely descriptive -- NOT used inside the
    greedy merge loop's stopping decision.
    """
    if not domains:
        return 0.0
    scores = []
    for res_list in domains.values():
        if not res_list:
            continue
        phi = compute_residue_occupancy(ca_coords, neighbors, d_core, d_max, min_seq_dist)
        # only average over this domain's own residues
        domain_sasa = (1.0 - phi[res_list]).mean()
        scores.append(domain_sasa)
    return float(np.mean(scores)) if scores else 0.0


def run_infomap_hierarchy(im, graph_triplets: List[Tuple[int, int, float]]) -> Dict[int, Tuple[int, ...]]:
    """Runs Infomap and extracts the multi-level cluster path per node.
    Instantiate infomap outside with either
        im = Infomap("--two-level --silent") or
        im = Infomap()
    """
    for u, v, w in graph_triplets:
        im.add_link(int(u), int(v), float(w))

    im.run()

    hierarchy = {}
    for node in im.tree:
        if node.is_leaf:
            hierarchy[node.node_id] = node.path[:-1]
    return hierarchy

# ==============================================================================
# 5. Infomap Core Community Extractor
# ==============================================================================

from typing import Any, Dict, List
import numpy as np


def get_core_communities_from_infomap(
    im: Any,
    segments: Any,
    clusters: Any,
) -> Dict[int, np.ndarray]:
    """Extracts module assignments from Infomap and maps cluster IDs back

    through SSE segments to individual residue indices.

    Parameters
    ----------
    im : Infomap
        Infomap instance after running optimization (containing .nodes).
    segments : Any (np.ndarray of shape (N, 3) or List[Tuple[str, int, int]])
        List/array of SSE segments as (sse_type, start_residue, end_residue).
    clusters : Any (Dict[int, np.ndarray], List[np.ndarray], or np.ndarray)
        Sheet/strand clusters mapping cluster ID to segment IDs.

    Returns
    -------
    Dict[int, np.ndarray]
        Mapping from 0-indexed module ID to 1D int32 array of residue indices.
    """
    n_segs = len(segments) if segments is not None else 0
    if n_segs == 0 or im is None:
        return {}

    # 1. Normalize cluster iteration (handles dict, list, or numpy object array)
    if isinstance(clusters, dict):
        cluster_iterable = list(clusters.values())
    else:
        cluster_iterable = list(clusters)

    cluster_to_segments: Dict[int, List[int]] = {}

    for cluster_id, cluster in enumerate(cluster_iterable):
        seg_ids = []
        for item in cluster:
            # Safely handle both scalar seg_id (np.int64/int) and tuple (seg_id, range)
            seg_id = int(item[0] if isinstance(item, (tuple, list)) else item)
            seg_ids.append(seg_id)
        cluster_to_segments[cluster_id] = seg_ids

    # 2. Process non-strand singleton clusters (helices, turns, coils)
    current_cluster_id = len(cluster_iterable)
    claimed_segs = {
        seg_id
        for segs in cluster_to_segments.values()
        for seg_id in segs
    }

    for seg_id in range(n_segs):
        if seg_id not in claimed_segs:
            cluster_to_segments[current_cluster_id] = [seg_id]
            current_cluster_id += 1

    # 3. Extract start and end boundaries for vectorised residue expansion
    if isinstance(segments, np.ndarray):
        starts = segments[:, 1].astype(np.int32)
        ends = segments[:, 2].astype(np.int32)
    else:
        starts = np.fromiter((int(s[1]) for s in segments), dtype=np.int32)
        ends = np.fromiter((int(s[2]) for s in segments), dtype=np.int32)

    # 4. Map Infomap modules to residue index lists
    raw_core_communities: Dict[int, List[int]] = {}

    for node in im.nodes:
        cluster_id = node.node_id
        module_id = node.module_id - 1  # 0-indexed module ID

        if module_id not in raw_core_communities:
            raw_core_communities[module_id] = []

        seg_ids = cluster_to_segments.get(cluster_id, [])
        for seg_id in seg_ids:
            if 0 <= seg_id < n_segs:
                start_res, end_res = starts[seg_id], ends[seg_id]
                raw_core_communities[module_id].extend(
                    range(start_res, end_res + 1)
                )

    # Convert residue index lists to sorted, unique 1D NumPy int32 arrays
    core_communities: Dict[int, np.ndarray] = {}
    for mod_id, res_list in raw_core_communities.items():
        if res_list:
            core_communities[mod_id] = np.unique(
                np.array(res_list, dtype=np.int32)
            )
        else:
            core_communities[mod_id] = np.empty(0, dtype=np.int32)

    return core_communities


def assign_non_core_residues(
    xyz: np.ndarray,
    core_communities: Dict[int, List[int]],
    dist_threshold: float = 16.0
) -> np.ndarray:
    """
    Assigns non-CORE residues to the closest community containing a CORE residue.
    
    Args:
        xyz: (N, 3) float array of residue coordinates.
        core_communities: Dict mapping community_id -> list of residue indices.
                          e.g., {0: [10, 11, 12], 1: [45, 46, 47]}
        dist_threshold: Maximum allowed distance (Å) to assign to a community.

    Returns:
        assignments: (N,) int array where assignments[i] is community_id, 
                     or -1 if unassigned/too far.
    """
    n_residues = len(xyz)
    assignments = np.full(n_residues, -1, dtype=np.int32)

    # 1. Flatten CORE residues and build a mapping to their community IDs
    core_indices = []
    core_comm_ids = []

    for comm_id, res_indices in core_communities.items():
      try:
        # Set existing core assignments directly
        assignments[res_indices] = comm_id

        core_indices.extend(res_indices)
        core_comm_ids.extend([comm_id] * len(res_indices))
      except:
        continue

    if not core_indices:
        return assignments

    core_indices = np.array(core_indices, dtype=np.int32)
    core_comm_ids = np.array(core_comm_ids, dtype=np.int32)

    # 2. Identify non-CORE residues
    non_core_mask = assignments == -1
    non_core_indices = np.where(non_core_mask)[0]

    if len(non_core_indices) == 0:
        return assignments

    # 3. Build a SINGLE cKDTree of all CORE residue coordinates
    core_pts = xyz[core_indices]
    tree = cKDTree(core_pts)

    # 4. Single-pass query for all non-CORE points at once
    # distance_upper_bound stops tree search early for points > threshold
    non_core_pts = xyz[non_core_indices]
    distances, nearest_core_local_indices = tree.query(
        non_core_pts,
        k=1,
        distance_upper_bound=dist_threshold
    )

    # 5. Map nearest CORE points back to their community IDs
    # Points farther than distance_upper_bound return index == len(core_pts)
    valid_mask = nearest_core_local_indices < len(core_pts)
    valid_non_core_idx = non_core_indices[valid_mask]
    valid_core_local_idx = nearest_core_local_indices[valid_mask]

    # Assign community IDs to valid non-CORE residues
    assignments[valid_non_core_idx] = core_comm_ids[valid_core_local_idx]

    return assignments

#----------------------------------------------------------------
# Repair final partitioning: close gaps, reassign short spans
#----------------------------------------------------------------

import numpy as np
from scipy.spatial.distance import cdist
from typing import List, Tuple, Dict, Set

# ----------------------------------------------------------------
# Repair final partitioning: close gaps, reassign short spans
# ----------------------------------------------------------------

def residue_assignments_to_domain_ranges(
    final_domains: np.ndarray,
) -> np.ndarray:
    """Converts 1D array of per-residue domain assignments into 2D array of range spans.

    Parameters
    ----------
    final_domains : np.ndarray
        1D integer array of domain IDs for each residue (-1 or negative for unassigned).

    Returns
    -------
    np.ndarray
        2D int32 array of shape (M, 3), where each row is:
        [domain_id, start_residue_1based, end_residue_1based].
    """
    final_domains = np.asarray(final_domains, dtype=np.int32)
    valid_mask = final_domains >= 0
    if not np.any(valid_mask):
        return np.empty((0, 3), dtype=np.int32)

    unique_domains = np.unique(final_domains[valid_mask])
    range_list = []

    for dom_id in unique_domains:
        res_indices = np.where(final_domains == dom_id)[0]
        if res_indices.size == 0:
            continue

        # Detect breaks in continuous sequence
        breaks = np.where(np.diff(res_indices) != 1)[0]
        starts = np.insert(res_indices[breaks + 1], 0, res_indices[0])
        ends = np.append(res_indices[breaks], res_indices[-1])

        # Convert to 1-based inclusive indices
        for s, e in zip(starts, ends):
            range_list.append([dom_id, s + 1, e + 1])

    if not range_list:
        return np.empty((0, 3), dtype=np.int32)

    return np.array(range_list, dtype=np.int32)


from typing import Tuple
import numpy as np


def parse_sheets(
    sheets_data: np.ndarray, n_residues: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Parses sheet structure into vectorized residue-cluster lookups.

    Handles homogeneous 2D/3D NumPy arrays as well as ragged/nested Python structures
    (e.g. Dict[int, List], List[List[Tuple]], or object arrays).

    Parameters
    ----------
    sheets_data : Any
        Sheet clusters data structure. Supported formats:
        - Dict[int, List[Tuple[int, int]]] mapping cluster_id to list of (start, end)
        - List[List[Tuple[int, Tuple[int, int]]]] or List[List[Tuple[int, int]]]
        - 2D/3D np.ndarray
    n_residues : int
        Total number of residues in the protein structure.

    Returns
    -------
    res_to_cluster : np.ndarray
        1D int32 array of length `n_residues`, where each index holds the cluster ID or -1.
    cluster_bounds : np.ndarray
        2D int32 array of shape (C, 3) where rows are [cluster_id, min_res_idx, max_res_idx + 1].
    """
    res_to_cluster = np.full(n_residues, -1, dtype=np.int32)
    if sheets_data is None or len(sheets_data) == 0:
        return res_to_cluster, np.empty((0, 3), dtype=np.int32)

    flat_records = []

    # 1. Parse dictionary structure: {cluster_id: [(start, end), ...]} or {cluster_id: [seg_ids...]}
    if isinstance(sheets_data, dict):
        for cid, cluster_items in sheets_data.items():
            for item in cluster_items:
                if isinstance(item, (tuple, list, np.ndarray)):
                    if len(item) == 2 and isinstance(item[1], (tuple, list)):
                        # Format: (strand_idx, (start, end))
                        s, e = item[1][0], item[1][1]
                    else:
                        s, e = item[0], item[1]
                    flat_records.append([int(cid), int(s), int(e)])
                else:
                    # Single segment/residue ID
                    flat_records.append([int(cid), int(item), int(item)])

    # 2. Parse nested lists / object arrays / ragged lists
    else:
        # Check if it's already a clean 2D numeric NumPy array
        if (
            isinstance(sheets_data, np.ndarray)
            and sheets_data.dtype != object
            and sheets_data.ndim == 2
        ):
            if sheets_data.shape[1] == 4:
                # [cluster_id, strand_id, start, end]
                for cid, _, s, e in sheets_data:
                    flat_records.append([int(cid), int(s), int(e)])
            elif sheets_data.shape[1] == 3:
                # [cluster_id, start, end]
                for cid, s, e in sheets_data:
                    flat_records.append([int(cid), int(s), int(e)])
        else:
            # Fallback for ragged/nested list or object array:
            # e.g., List[List[(strand_idx, (start, end))]] or List[List[(start, end)]]
            cluster_iterable = (
                sheets_data.items()
                if hasattr(sheets_data, "items")
                else enumerate(sheets_data)
            )
            for cid, cluster_list in cluster_iterable:
                if cluster_list is None:
                    continue
                for item in cluster_list:
                    if isinstance(item, (tuple, list, np.ndarray)):
                        if len(item) == 2 and isinstance(
                            item[1], (tuple, list, np.ndarray)
                        ):
                            # Format: (strand_idx, (start, end))
                            s, e = item[1][0], item[1][1]
                        elif len(item) >= 2:
                            # Format: (start, end) or (cluster_id, start, end)
                            s, e = item[-2], item[-1]
                        else:
                            s = e = item[0]
                        flat_records.append([int(cid), int(s), int(e)])
                    else:
                        flat_records.append([int(cid), int(item), int(item)])

    if not flat_records:
        return res_to_cluster, np.empty((0, 3), dtype=np.int32)

    # Convert flattened list into a clean 2D int32 NumPy array
    sheets_arr = np.array(flat_records, dtype=np.int32)

    cluster_ids = sheets_arr[:, 0]
    starts = sheets_arr[:, 1]
    ends = sheets_arr[:, 2]

    # Convert 1-indexed range boundaries to 0-indexed half-open intervals
    s_indices = np.maximum(0, starts - 1)
    e_indices = np.minimum(n_residues, ends)

    for cid, s, e in zip(cluster_ids, s_indices, e_indices):
        res_to_cluster[s:e] = cid

    unique_clusters = np.unique(cluster_ids[cluster_ids >= 0])
    cluster_bounds_list = []

    for cid in unique_clusters:
        res_mask = res_to_cluster == cid
        res_idxs = np.where(res_mask)[0]
        if res_idxs.size > 0:
            cluster_bounds_list.append([cid, res_idxs[0], res_idxs[-1] + 1])

    if not cluster_bounds_list:
        cluster_bounds = np.empty((0, 3), dtype=np.int32)
    else:
        cluster_bounds = np.array(cluster_bounds_list, dtype=np.int32)

    return res_to_cluster, cluster_bounds

def repair_domain_assignments_3d(
    coords: np.ndarray,
    domain_ranges: np.ndarray,
    sheets_data: np.ndarray,
    contact_cutoff: float = 8.0,
    w_seq: float = 3.0,
    w_contact: float = 1.0,
    min_core_len: int = 15,
    max_iter: int = 15,
) -> Tuple[np.ndarray, np.ndarray]:
    """Repairs domain assignments while preserving protected beta-sheets and repeat elements.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (N, 3) with residue 3D coordinates.
    domain_ranges : np.ndarray
        2D array of shape (M, 3) with rows [domain_id, start_1based, end_1based],
        or 1D array of domain assignments of length N.
    sheets_data : np.ndarray
        2D array of shape (K, 4) or (K, 3) containing strand ranges per sheet cluster.
    contact_cutoff : float
        Distance cutoff (in Angstroms) for defining 3D contacts.
    w_seq : float
        Weight penalty/support factor for sequence proximity.
    w_contact : float
        Weight factor for 3D spatial contact energy.
    min_core_len : int
        Minimum length threshold to lock anchor core spans.
    max_iter : int
        Maximum optimization iterations.

    Returns
    -------
    labels : np.ndarray
        1D int32 array of shape (N,) containing final domain assignments per residue.
    domain_ranges_out : np.ndarray
        2D int32 array of shape (M, 3) containing final domain range spans:
        [domain_id, start_1based, end_1based].
    """
    n_residues = len(coords)
    labels = np.full(n_residues, -1, dtype=np.int32)

    # 1. Parse initial domain labels
    domain_ranges = np.asarray(domain_ranges, dtype=np.int32)
    if domain_ranges.ndim == 1 and len(domain_ranges) == n_residues:
        labels = domain_ranges.copy()
    elif domain_ranges.ndim == 2 and domain_ranges.shape[1] == 3:
        for dom_id, start, end in domain_ranges:
            s_idx = max(0, start - 1)
            e_idx = min(n_residues, end)
            labels[s_idx:e_idx] = dom_id

    unique_domains = np.unique(labels[labels >= 0])
    num_domains = (
        int(np.max(unique_domains)) + 1 if unique_domains.size > 0 else 0
    )
    if num_domains == 0:
        return labels, np.empty((0, 3), dtype=np.int32)

    # 2. Parse sheet elements
    res_to_cluster, _ = parse_sheets(sheets_data, n_residues)
    unique_sheet_clusters = np.unique(res_to_cluster[res_to_cluster >= 0])

    # Unify initial cluster labels by majority vote
    for cid in unique_sheet_clusters:
        c_mask = res_to_cluster == cid
        c_labels = labels[c_mask & (labels != -1)]
        if c_labels.size > 0:
            counts = np.bincount(c_labels)
            maj_domain = np.argmax(counts)
            labels[c_mask] = maj_domain

    # 3. Build 3D Contact Matrix
    dist_mat = cdist(coords, coords)
    contact_mask = (dist_mat <= contact_cutoff) & (dist_mat > 0)
    contact_weights = np.zeros_like(dist_mat)
    contact_weights[contact_mask] = np.exp(-dist_mat[contact_mask] / 4.0)

    # 4. Identify Anchors (Spans >= min_core_len)
    is_core = np.zeros(n_residues, dtype=bool)
    if n_residues > 0:
        label_changes = np.where(labels[:-1] != labels[1:])[0] + 1
        split_starts = np.insert(label_changes, 0, 0)
        split_ends = np.append(label_changes, n_residues)

        for s, e in zip(split_starts, split_ends):
            lbl = labels[s]
            if (e - s) >= min_core_len and lbl != -1:
                is_core[s:e] = True

    # 5. Optimization Loop
    for iteration in range(max_iter):
        changes = 0

        # --- PASS A: Reassign Non-Sheet Residues Individually ---
        new_labels = labels.copy()
        unprotected_mask = (res_to_cluster == -1) & ~(
            is_core & (iteration < 3)
        )
        unprotected_indices = np.where(unprotected_mask)[0]

        for r in unprotected_indices:
            current_lbl = labels[r]
            domain_scores = np.zeros(num_domains, dtype=np.float64)

            # Sequence context
            seq_min = max(0, r - 2)
            seq_max = min(n_residues, r + 3)
            for j in range(seq_min, seq_max):
                if j != r and labels[j] != -1:
                    domain_scores[labels[j]] += w_seq / abs(r - j)

            # 3D Contact context
            contacts = np.where(contact_mask[r])[0]
            for j in contacts:
                if j != r and labels[j] != -1:
                    domain_scores[labels[j]] += (
                        w_contact * contact_weights[r, j]
                    )

            best_domain = int(np.argmax(domain_scores))
            best_score = domain_scores[best_domain]
            current_score = (
                domain_scores[current_lbl] if current_lbl != -1 else -1.0
            )

            if best_domain != current_lbl and best_score > current_score + 0.1:
                print("# reassign",r, current_lbl, ' -> ', best_domain)
                new_labels[r] = best_domain
                changes += 1

        labels = new_labels

        # --- PASS B: Atomic Reassignment of Sheet Clusters ---
        for cid in unique_sheet_clusters:
            cluster_mask = res_to_cluster == cid
            cluster_res = np.where(cluster_mask)[0]
            if cluster_res.size == 0:
                continue

            curr_cluster_lbl = labels[cluster_res[0]]

            # Skip if locked core
            if np.any(is_core[cluster_mask]) and iteration < 3:
                continue

            cluster_scores = np.zeros(num_domains, dtype=np.float64)

            # Aggregate external support across ALL residues in the sheet cluster
            for r in cluster_res:
                contacts = np.where(contact_mask[r])[0]
                for j in contacts:
                    if not cluster_mask[j] and labels[j] != -1:
                        cluster_scores[labels[j]] += (
                            w_contact * contact_weights[r, j]
                        )

                if (
                    r - 1 >= 0
                    and not cluster_mask[r - 1]
                    and labels[r - 1] != -1
                ):
                    cluster_scores[labels[r - 1]] += w_seq
                if (
                    r + 1 < n_residues
                    and not cluster_mask[r + 1]
                    and labels[r + 1] != -1
                ):
                    cluster_scores[labels[r + 1]] += w_seq

            best_cluster_domain = int(np.argmax(cluster_scores))

            if (
                best_cluster_domain != curr_cluster_lbl
                and cluster_scores[best_cluster_domain]
                > cluster_scores[curr_cluster_lbl] + 0.5
            ):
                print("# cluster reassign",curr_cluster_lbl, " -> ", best_cluster_domain) 
                labels[cluster_mask] = best_cluster_domain
                changes += 1

        # --- PASS C: Fragment-Level Reassignment (Unprotected Loop Chunks) ---
        label_changes = np.where(labels[:-1] != labels[1:])[0] + 1
        split_starts = np.insert(label_changes, 0, 0)
        split_ends = np.append(label_changes, n_residues)

        for start, end in zip(split_starts, split_ends):
            span_len = end - start
            lbl = labels[start]

            if span_len < min_core_len and not np.all(is_core[start:end]):
                if np.all(res_to_cluster[start:end] == -1):
                    frag_indices = np.arange(start, end)
                    frag_scores = np.zeros(num_domains, dtype=np.float64)

                    for r in frag_indices:
                        contacts = np.where(contact_mask[r])[0]
                        external_contacts = contacts[
                            (contacts < start) | (contacts >= end)
                        ]
                        for j in external_contacts:
                            if labels[j] != -1:
                                frag_scores[labels[j]] += (
                                    w_contact * contact_weights[r, j]
                                )

                    left_lbl = labels[start - 1] if start > 0 else -1
                    right_lbl = labels[end] if end < n_residues else -1

                    if left_lbl != -1:
                        frag_scores[left_lbl] += w_seq * 2.0
                    if right_lbl != -1:
                        frag_scores[right_lbl] += w_seq * 2.0

                    best_frag_domain = int(np.argmax(frag_scores))
                    if (
                        best_frag_domain != lbl
                        and frag_scores[best_frag_domain]
                        > frag_scores[lbl] + 0.5
                    ):
                        print("# loop reassign",start,end,' -> ', best_frag_domain)
                        labels[start:end] = best_frag_domain
                        changes += 1

        if changes == 0:
            break

    # 6. Convert assignments to range array format
    domain_ranges_out = residue_assignments_to_domain_ranges(labels)

    return labels, domain_ranges_out
