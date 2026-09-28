"""
structal_similarity.py

STRUCTAL similarity pipeline for a DaliScope project view.

(0) Select up to max_targets rows from project.views[view_name].
(1) Apply each target's precomputed (R, t) -- the domain-clipped refined
    initial superimposition already computed upstream -- to its full CA
    coordinate slice. No Kabsch refit here; R/t are taken as given.
(2) For every target pair, generate candidate residue-residue matches via
    a vectorized KDTree radius join over the FULL superimposed structures
    (not just the query-defined core, so shared structure beyond that
    core can extend the apparent rigid core), then resolve them into a
    single sequential (non-crossing) alignment via a numba-JIT dynamic
    program -- this is what rules out the biologically-invalid
    non-sequential matches that plain (mutual) nearest-neighbor matching
    can produce.
(3) Build a tree and a matching heatmap from the similarity matrix.
    Two options, both driven by the same distance
        D(i,j) = sim(i,i) - sim(i,j) + sim(j,j) - sim(j,i)
              = sim(i,i) + sim(j,j) - 2*sim(i,j)   (sim symmetric)
    which is exact/non-negative for this score since a DP chain can
    match at most min(length_i, length_j) residues:

    (a) Neighbor-Joining (recommended) -- build_nj_tree + plot_nj_tree +
        plot_similarity_heatmap_ordered. Produces a proper additive tree:
        each leaf gets its OWN branch length toward a merge (large,
        long-self-similarity proteins get long branches; small ones stay
        close in), rather than forcing every leaf to the same depth.
        Note: NJ trees are inherently unrooted; plot_nj_tree applies
        midpoint rooting purely for a sensible layout, and NJ topology
        can legitimately differ from UPGMA's when self-similarities vary
        a lot (it isn't forced to be ultrametric).
    (b) UPGMA (scipy linkage) -- plot_similarity_heatmap +
        plot_similarity_dendrogram. Standard cladogram: every leaf ends
        at the same depth, so a leaf-pair's branch length is split
        equally regardless of each protein's own size.

Usage (Neighbor-Joining, recommended)
--------------------------------------
    target_ids, raw, norm, lengths = compute_similarity_matrix(
        project, view_name="my_view", max_targets=1000, d0=5.0,
    )
    labels = make_labels(project.views["my_view"], target_ids)

    tree = build_nj_tree(target_ids, raw)
    fig_tree, order = plot_nj_tree(tree, labels=labels)
    fig_heat = plot_similarity_heatmap_ordered(target_ids, raw, order, labels=labels)

Usage (UPGMA)
-------------
    fig_heat, ordered_ids, Z = plot_similarity_heatmap(target_ids, raw, labels=labels)
    fig_dend = plot_similarity_dendrogram(target_ids, Z=Z, labels=labels)
"""

import ast
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

try:
    import numba
    _HAVE_NUMBA = True
except ImportError:
    _HAVE_NUMBA = False
    import warnings
    warnings.warn(
        "numba not found -- the DP equivalence-resolution step will fall "
        "back to a much slower pure-Python loop. `pip install numba` "
        "(or `pip install numba --break-system-packages` in Colab) for "
        "the vectorized/JIT path this module is designed around."
    )

try:
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    import Bio.Phylo as Phylo
    _HAVE_BIOPYTHON = True
except ImportError:
    _HAVE_BIOPYTHON = False

try:
    from pyroaring import BitMap
    _HAVE_PYROARING = True
except ImportError:
    _HAVE_PYROARING = False


# =======================================================================
# (2a) Fenwick max-tree + DP chain resolution (numba-JIT)
#
# b values are used directly as Fenwick positions (no rank-compression
# needed): b is literally the KDTree atom index into structure j, which
# is already a dense 0..len_j-1 integer range.
# =======================================================================

if _HAVE_NUMBA:
    @numba.njit(cache=True)
    def _fenwick_update_max(tree, i, value):
        i += 1
        n = len(tree) - 1
        while i <= n:
            if value > tree[i]:
                tree[i] = value
            i += i & (-i)

    @numba.njit(cache=True)
    def _fenwick_query_max(tree, i):
        if i < 0:
            return 0.0
        i += 1
        res = 0.0
        while i > 0:
            if tree[i] > res:
                res = tree[i]
            i -= i & (-i)
        return res

    @numba.njit(cache=True)
    def _dp_chain_numba(a_sorted, b_sorted, score_sorted, n_b):
        """a_sorted, b_sorted, score_sorted already sorted by (a asc, b asc).
        Returns the max total score of a strictly-(a,b)-increasing chain --
        i.e. a valid, non-crossing, sequential alignment."""
        n = len(a_sorted)
        fenwick = np.zeros(n_b + 1, dtype=np.float64)
        global_best = 0.0
        idx = 0
        while idx < n:
            cur_a = a_sorted[idx]
            start = idx
            while idx < n and a_sorted[idx] == cur_a:
                idx += 1
            end = idx
            vals = np.empty(end - start, dtype=np.float64)
            for k in range(start, end):
                r = b_sorted[k]
                prefix = _fenwick_query_max(fenwick, r - 1)
                val = score_sorted[k] + prefix
                vals[k - start] = val
                if val > global_best:
                    global_best = val
            for k in range(start, end):
                r = b_sorted[k]
                _fenwick_update_max(fenwick, r, vals[k - start])
        return global_best


class _FenwickMax:
    """Pure-Python fallback, only used if numba is unavailable."""
    def __init__(self, n):
        self.n = n
        self.tree = [0.0] * (n + 1)

    def update_max(self, i, value):
        i += 1
        while i <= self.n:
            if value > self.tree[i]:
                self.tree[i] = value
            i += i & (-i)

    def query_max(self, i):
        if i < 0:
            return 0.0
        i += 1
        res = 0.0
        while i > 0:
            if self.tree[i] > res:
                res = self.tree[i]
            i -= i & (-i)
        return res


def _dp_chain_score_python_fallback(candidates):
    candidates = list(candidates)
    if not candidates:
        return 0.0
    b_vals = sorted(set(b for _, b, _ in candidates))
    b_rank = {b: r for r, b in enumerate(b_vals)}
    fenwick = _FenwickMax(len(b_vals))
    candidates.sort(key=lambda t: (t[0], t[1]))
    global_best = 0.0
    idx = 0
    n = len(candidates)
    while idx < n:
        cur_a = candidates[idx][0]
        batch = []
        while idx < n and candidates[idx][0] == cur_a:
            batch.append(candidates[idx])
            idx += 1
        batch_updates = []
        for (a, b, score) in batch:
            r = b_rank[b]
            prefix = fenwick.query_max(r - 1)
            val = score + prefix
            batch_updates.append((r, val))
            if val > global_best:
                global_best = val
        for (r, val) in batch_updates:
            fenwick.update_max(r, val)
    return global_best


# =======================================================================
# (1) Apply precomputed transforms
# =======================================================================

def _coerce_array(value):
    """R/t cells might already be np.ndarray, or a list, or (if the
    dataframe was round-tripped through CSV) a stringified list."""
    if isinstance(value, np.ndarray):
        return value
    if isinstance(value, (list, tuple)):
        return np.asarray(value, dtype=float)
    if isinstance(value, str):
        return np.asarray(ast.literal_eval(value), dtype=float)
    raise TypeError(f"Unrecognized R/t cell type: {type(value)}")


def apply_transforms(project, df, target_col="target_id",
                      r_col="R", t_col="t", verbose=True):
    """
    Returns
    -------
    dict: target_id -> (N_t, 3) float32 array, CA coords in the common
    (query-anchored) frame.

    Convention assumed: x_query_frame = x_target @ R.T + t   (row-vector
    convention). If the sanity-check RMSD printed below looks large
    (many Angstroms) rather than typical DALI-superposition RMSD
    (usually < 3-4 A over the aligned core), the stored convention is
    probably transposed/column-vector -- see the note at the bottom.
    """
    meta = project.master_metadata.set_index(target_col)
    coords_dict = {}
    check_rmsds = []

    segs_by_alignment = None
    if hasattr(project, "segments") and "alignment_id" in getattr(project.segments, "columns", []):
        segs_by_alignment = project.segments.groupby("alignment_id")

    for _, row in df.iterrows():
        target_id = row[target_col]
        if target_id not in meta.index:
            continue
        meta_row = meta.loc[target_id]
        start, end = int(meta_row["start_idx"]), int(meta_row["end_idx"])
        local_coords = np.asarray(project.coords[start:end], dtype=np.float64)

        R = _coerce_array(row[r_col])
        t = _coerce_array(row[t_col])
        transformed = local_coords @ R.T + t
        coords_dict[target_id] = transformed.astype(np.float32)

        # Optional sanity check: RMSD over the target's own aligned pairs
        # (from project.segments), transformed vs query_ca_coords.
        if segs_by_alignment is not None and "alignment_id" in meta_row.index:
            aln_id = meta_row["alignment_id"]
            if aln_id in segs_by_alignment.groups:
                segs = segs_by_alignment.get_group(aln_id)
                q_idx, t_idx = [], []
                for _, seg in segs.iterrows():
                    L = int(seg["length"])
                    q_idx.extend(range(int(seg["q_start"]), int(seg["q_start"]) + L))
                    t_idx.extend(range(int(seg["s_start"]), int(seg["s_start"]) + L))
                if q_idx:
                    q_pts = project.query_ca_coords[q_idx]
                    t_pts = transformed[t_idx]
                    rmsd = np.sqrt(np.mean(np.sum((q_pts - t_pts) ** 2, axis=1)))
                    check_rmsds.append(rmsd)

    if verbose and check_rmsds:
        check_rmsds = np.array(check_rmsds)
        print(f"[apply_transforms] sanity check over {len(check_rmsds)} targets: "
              f"median aligned-core RMSD = {np.median(check_rmsds):.2f} A "
              f"(range {check_rmsds.min():.2f}-{check_rmsds.max():.2f}).")
        print("  Typical refined-superposition RMSD should be a few A. "
              "If this is instead tens/hundreds of A, the R/t convention "
              "is probably transposed -- try `local_coords @ R + t` or "
              "`(R @ local_coords.T).T + t` instead.")

    return coords_dict


# =======================================================================
# (2) Pairwise STRUCTAL score via KDTree candidates + DP resolution
# =======================================================================

def build_trees(coords_dict):
    return {tid: cKDTree(coords) for tid, coords in coords_dict.items()}

    try:
        import numba
        _HAVE_NUMBA = True
    except ImportError:
        _HAVE_NUMBA = False
        print("numba not installed, running without number")
        # Dummy decorator that leaves functions uncompiled if numba is absent
        class numba:
            @staticmethod
            def njit(*args, **kwargs):
                return lambda fn: fn

if _HAVE_NUMBA:
    @numba.njit(cache=True)
    def _fenwick_update_max_idx(tree_val, tree_idx, i, value, src_idx):
        i += 1
        n = len(tree_val) - 1
        while i <= n:
            if value > tree_val[i]:
                tree_val[i] = value
                tree_idx[i] = src_idx
            i += i & (-i)

    @numba.njit(cache=True)
    def _fenwick_query_max_idx(tree_val, tree_idx, i):
        if i < 0:
            return 0.0, -1
        i += 1
        res = 0.0
        res_idx = -1
        while i > 0:
            if tree_val[i] > res:
                res = tree_val[i]
                res_idx = tree_idx[i]
            i -= i & (-i)
        return res, res_idx

    @numba.njit(cache=True)
    def _dp_chain_numba_backtrack(a_sorted, b_sorted, score_sorted, n_b):
        """Same DP as _dp_chain_numba, but also returns enough to
        reconstruct the actual matched (a,b) chain: the best-scoring
        chain's endpoint index and a parent pointer per candidate."""
        n = len(a_sorted)
        fenwick_val = np.zeros(n_b + 1, dtype=np.float64)
        fenwick_idx = np.full(n_b + 1, -1, dtype=np.int64)
        parent = np.full(n, -1, dtype=np.int64)
        dp_val = np.zeros(n, dtype=np.float64)

        global_best = 0.0
        global_best_idx = -1
        idx = 0
        while idx < n:
            cur_a = a_sorted[idx]
            start = idx
            while idx < n and a_sorted[idx] == cur_a:
                idx += 1
            end = idx
            for k in range(start, end):
                r = b_sorted[k]
                pv, pidx = _fenwick_query_max_idx(fenwick_val, fenwick_idx, r - 1)
                val = score_sorted[k] + pv
                dp_val[k] = val
                parent[k] = pidx
                if val > global_best:
                    global_best = val
                    global_best_idx = k
            for k in range(start, end):
                r = b_sorted[k]
                _fenwick_update_max_idx(fenwick_val, fenwick_idx, r, dp_val[k], k)
        return global_best, global_best_idx, parent


def pairwise_structal_alignment(coords_i, tree_i, coords_j, tree_j,
                                 d0=5.0, radius=None):
    """Like pairwise_structal_score, but returns the actual matched
    (a,b) index-pair chain (backtracked from the DP) instead of just the
    scalar score. Needed whenever you care WHICH residues matched, not
    just how well they matched -- e.g. recovering how far a rigid core
    extends beyond a fixed reference alignment, or building coverage
    inputs for a representative-selection algorithm.

    Returns (score, a_indices, b_indices) -- a_indices/b_indices are the
    matched positions in structure i / j respectively, in matching order
    (a_indices[k] is aligned to b_indices[k]), strictly increasing in
    both (a valid, non-crossing chain)."""
    if radius is None:
        radius = 2.0 * d0

    sdm = tree_i.sparse_distance_matrix(tree_j, max_distance=radius,
                                         output_type="coo_matrix")
    if sdm.nnz == 0:
        return 0.0, np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)

    a = sdm.row.astype(np.int64)
    b = sdm.col.astype(np.int64)
    score = 1.0 / (1.0 + (sdm.data / d0) ** 2)
    order = np.lexsort((b, a))
    a_s, b_s, s_s = a[order], b[order], score[order]
    n_b = len(coords_j)

    if not _HAVE_NUMBA:
        raise ImportError(
            "pairwise_structal_alignment requires numba (the pure-Python "
            "DP fallback doesn't implement backtracking). `pip install numba`."
        )

    best_score, best_idx, parent = _dp_chain_numba_backtrack(a_s, b_s, s_s, n_b)
    chain_a, chain_b = [], []
    cur = int(best_idx)
    while cur != -1:
        chain_a.append(a_s[cur])
        chain_b.append(b_s[cur])
        cur = int(parent[cur])
    chain_a.reverse()
    chain_b.reverse()
    return float(best_score), np.array(chain_a, dtype=np.int64), np.array(chain_b, dtype=np.int64)


def pairwise_structal_score(coords_i, tree_i, coords_j, tree_j,
                             d0=5.0, radius=None):
    """STRUCTAL score for one pair, resolving equivalences via DP over
    candidates found within `radius` (default 2*d0).

    Fully vectorized candidate generation: cKDTree.sparse_distance_matrix
    computes every (a,b) pair within radius AND its distance in one C-level
    call (no per-atom Python loop), then the score transform and sort are
    plain numpy array ops. Equivalence resolution runs in the numba-JIT
    DP above when available (falls back to a pure-Python DP otherwise,
    with a warning issued at import time)."""
    if radius is None:
        radius = 2.0 * d0

    sdm = tree_i.sparse_distance_matrix(tree_j, max_distance=radius,
                                         output_type="coo_matrix")
    if sdm.nnz == 0:
        return 0.0

    a = sdm.row.astype(np.int64)
    b = sdm.col.astype(np.int64)
    score = 1.0 / (1.0 + (sdm.data / d0) ** 2)

    if _HAVE_NUMBA:
        order = np.lexsort((b, a))
        n_b = len(coords_j)
        return float(_dp_chain_numba(a[order], b[order], score[order], n_b))
    else:
        candidates = list(zip(a.tolist(), b.tolist(), score.tolist()))
        return _dp_chain_score_python_fallback(candidates)


# =======================================================================
# Full matrix computation
# =======================================================================

def compute_similarity_matrix(project, view_name, max_targets=1000,
                               d0=5.0, radius_factor=2.0,
                               target_col="target_id", verbose=True):
    df = project.views[view_name].head(max_targets).copy()

    coords_dict = apply_transforms(project, df, target_col=target_col,
                                    verbose=verbose)
    target_ids = [tid for tid in df[target_col] if tid in coords_dict]
    n = len(target_ids)

    trees = build_trees(coords_dict)
    lengths = np.array([len(coords_dict[tid]) for tid in target_ids])

    raw = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        raw[i, i] = lengths[i]  # trivial self-alignment score = length

    radius = radius_factor * d0
    t0 = time.time()
    n_pairs = n * (n - 1) // 2
    done = 0
    for i in range(n):
        ti = target_ids[i]
        for j in range(i + 1, n):
            tj = target_ids[j]
            s = pairwise_structal_score(
                coords_dict[ti], trees[ti], coords_dict[tj], trees[tj],
                d0=d0, radius=radius,
            )
            raw[i, j] = raw[j, i] = s
            done += 1
        if verbose and n_pairs > 0 and (i % max(1, n // 20) == 0):
            elapsed = time.time() - t0
            frac = done / n_pairs if n_pairs else 1.0
            eta = elapsed / frac - elapsed if frac > 0 else float("nan")
            print(f"[compute_similarity_matrix] {done}/{n_pairs} pairs "
                  f"({frac:.0%}), elapsed {elapsed:.0f}s, ETA {eta:.0f}s")

    min_len = np.minimum.outer(lengths, lengths).astype(np.float64)
    norm = raw / min_len
    np.fill_diagonal(norm, 1.0)

    return target_ids, raw, norm, lengths


def similarity_to_distance(raw):
    """
    d(i,j) = sim(i,i) - sim(i,j) + sim(j,j) - sim(j,i)

    With a symmetric sim matrix this is sim(i,i) + sim(j,j) - 2*sim(i,j) --
    nodes sit at a "depth" equal to their own self-similarity (raw score
    against themselves, i.e. their length), and the distance between two
    nodes is how far apart their branches are above the point where they
    join. d(i,i) = 0 by construction.

    Guaranteed non-negative for this particular STRUCTAL/DP score: since
    the DP chain matches each residue at most once, score(i,j) <=
    min(length_i, length_j) <= min(sim(i,i), sim(j,j)) always holds, which
    gives d(i,j) >= |sim(i,i) - sim(j,j)| >= 0.
    """
    diag = np.diag(raw)
    d = diag[:, None] + diag[None, :] - raw - raw.T
    np.fill_diagonal(d, 0.0)
    return d


# =======================================================================
# (3) Heatmap + dendrogram: static seaborn/matplotlib versions
# =======================================================================

def _leaf_order(raw, linkage_method="average"):
    """Shared clustering step: returns (Z, leaf_order) so the heatmap and
    the dendrogram, drawn as separate static figures, stay consistent
    with each other."""
    dist = similarity_to_distance(raw)
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method=linkage_method)
    order = dendrogram(Z, no_plot=True)["leaves"]
    return Z, order


def make_labels(df, target_ids, target_col="target_id",
                 desc_col="description", max_desc_len=40):
    """Build 'target_id  description' labels, in target_ids order.
    Missing/absent description falls back to target_id alone."""
    label_map = {}
    for _, row in df.iterrows():
        tid = row[target_col]
        desc = str(row[desc_col]) if desc_col in df.columns and pd.notna(row.get(desc_col)) else ""
        if max_desc_len and len(desc) > max_desc_len:
            desc = desc[: max_desc_len - 1] + "\u2026"
        label_map[tid] = f"{tid}  {desc}".rstrip()
    return [label_map.get(tid, tid) for tid in target_ids]


def plot_similarity_heatmap(target_ids, raw, labels=None,
                             linkage_method="average",
                             title="STRUCTAL similarity",
                             cmap="viridis", figsize=None):
    """Static seaborn heatmap, square cells, rows/cols in dendrogram-
    leaf order. Returns (fig, ordered_target_ids, Z) -- pass that same Z
    into plot_similarity_dendrogram(..., Z=Z) to GUARANTEE the two
    figures share one tree, rather than each recomputing linkage from
    whatever matrix happens to be passed (which silently breaks if the
    two calls are given different matrices, e.g. one clipped/normalized
    and the other not)."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    n = len(target_ids)
    Z, order = _leaf_order(raw, linkage_method)
    ordered_raw = raw[np.ix_(order, order)]
    if labels is None:
        labels = target_ids
    ordered_labels = [labels[i] for i in order]

    if figsize is None:
        side = max(6.0, 0.28 * n)
        figsize = (side, side)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(ordered_raw, xticklabels=ordered_labels, yticklabels=ordered_labels,
                square=True, cmap=cmap, cbar_kws={"label": "STRUCTAL score"}, ax=ax)
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    ax.tick_params(axis="y", labelrotation=0, labelsize=7)
    plt.tight_layout()
    return fig, [target_ids[i] for i in order], Z


def plot_similarity_dendrogram(target_ids, raw=None, labels=None,
                                linkage_method="average",
                                title="STRUCTAL clustering", figsize=None,
                                Z=None):
    """Static matplotlib dendrogram (left-oriented so labels stay
    horizontal and readable).

    Pass `Z` (as returned by plot_similarity_heatmap) to guarantee this
    draws the SAME tree as the heatmap, regardless of what `raw` you
    happen to pass here -- if Z is given, raw is ignored for linkage
    purposes entirely. If Z is not given, raw is required and linkage is
    (re)computed from it -- only safe if you are certain it is the exact
    same matrix used for the corresponding heatmap call."""
    import matplotlib.pyplot as plt

    if Z is None:
        if raw is None:
            raise ValueError("Must supply either Z (from plot_similarity_heatmap) or raw.")
        Z, _ = _leaf_order(raw, linkage_method)

    n = Z.shape[0] + 1
    if labels is None:
        labels = target_ids
    if figsize is None:
        figsize = (max(6.0, 0.15 * n), max(4.0, 0.22 * n))

    fig, ax = plt.subplots(figsize=figsize)
    dendrogram(Z, labels=labels, ax=ax, orientation="left", leaf_font_size=7)
    # scipy's orientation="left" places leaf index 0 at the BOTTOM and the
    # last leaf at the TOP -- the opposite of plot_similarity_heatmap's row
    # order (row 0 at the top, standard imshow/seaborn convention). Without
    # this, the two figures show genuinely different leaf sequences, not
    # just a cosmetic flip -- invert so index 0 is at the top in both.
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel("distance: sim(i,i) + sim(j,j) - 2*sim(i,j)")
    plt.tight_layout()
    return fig



# =======================================================================
# (3b) Neighbor-Joining tree: proper asymmetric branch lengths
# =======================================================================

def build_nj_tree(target_ids, raw, root="midpoint"):
    """
    Build a Neighbor-Joining tree from the same distance
        D(i,j) = sim(i,i) + sim(j,j) - 2*sim(i,j)
    used elsewhere in this module. Unlike UPGMA, NJ does not force every
    leaf to the same depth: each leaf gets its own branch length, so
    e.g. two directly-joined leaves i,j split their shared distance as
    (sim(i,i)-sim(i,j)) to leaf i and (sim(j,j)-sim(i,j)) to leaf j --
    exactly matching sim(i,i)+sim(j,j)-2sim(i,j) = D(i,j), so this is a
    genuine re-derivation of the SAME distance, not a new computation.

    NJ trees are inherently unrooted; `root="midpoint"` applies midpoint
    rooting purely to get a sensible left-to-right layout, with no
    implied biological "ancestor" at that point. Pass root=None to skip.

    Tip names are kept as target_ids (stable, guaranteed-unique keys);
    use `labels=` in plot_nj_tree for display text instead.

    NJ topology can legitimately differ from the UPGMA/average-linkage
    topology used by plot_similarity_heatmap/plot_similarity_dendrogram,
    especially when self-similarities (diag(raw), e.g. protein length)
    vary a lot -- NJ isn't forced to produce an ultrametric tree.
    """
    if not _HAVE_BIOPYTHON:
        raise ImportError(
            "biopython is required for Neighbor-Joining trees. "
            "`pip install biopython` (or `pip install biopython "
            "--break-system-packages` in Colab)."
        )

    diag = np.diag(raw)
    D = diag[:, None] + diag[None, :] - raw - raw.T
    np.fill_diagonal(D, 0.0)
    D = np.maximum(D, 0.0)  # guard tiny floating-point negatives

    n = len(target_ids)
    lower_tri = [list(D[i, : i + 1]) for i in range(n)]
    dm = DistanceMatrix(names=list(target_ids), matrix=lower_tri)

    tree = DistanceTreeConstructor().nj(dm)
    if root == "midpoint":
        tree.root_at_midpoint()
    return tree


def plot_nj_tree(tree, target_ids, labels=None, title="STRUCTAL NJ tree", figsize=None):
    """
    Draw the NJ tree with Bio.Phylo, with target_id (+ optional
    description) tip labels.

    `target_ids` and `labels` (if a list, not a dict) must be aligned to
    each other, same convention as elsewhere in this module -- NOT
    aligned to tree.get_terminals()'s own traversal order, which is a
    different, tree-structure-dependent sequence.

    Returns (fig, order) where `order` is the ACTUAL rendered top-to-
    bottom leaf sequence -- extracted from the drawn text objects'
    positions, not assumed -- so it can be passed directly into
    plot_similarity_heatmap_ordered to guarantee the heatmap matches
    this exact tree, the same way Z-sharing does for the UPGMA path.
    """
    if not _HAVE_BIOPYTHON:
        raise ImportError(
            "biopython is required for Neighbor-Joining trees. "
            "`pip install biopython`."
        )
    import matplotlib.pyplot as plt

    if labels is None:
        label_map = {tid: tid for tid in target_ids}
    elif isinstance(labels, dict):
        label_map = labels
    else:
        label_map = dict(zip(target_ids, labels))

    def label_func(clade):
        if clade.is_terminal():
            return label_map.get(clade.name, clade.name)
        return None  # don't clutter the plot with internal "InnerN" labels

    n = len(target_ids)
    if figsize is None:
        figsize = (10.0, max(4.0, 0.28 * n))

    fig, ax = plt.subplots(figsize=figsize)
    Phylo.draw(tree, axes=ax, do_show=False, label_func=label_func)
    ax.set_title(title)

    # Bio.Phylo draws leaf labels as free-floating ax.text() annotations,
    # NOT tick labels -- layout tools like tight_layout() have no idea
    # they exist and won't reserve room for them, so long labels silently
    # run off the edge of the fixed-size canvas.
    #
    # Naively growing fig.set_size_inches() does NOT fix this on its own:
    # the axes' position is stored as a FRACTION of figure width, so
    # widening the figure re-stretches the axes proportionally too,
    # dragging the tree's rightmost tip outward right along with the
    # added space -- the gap never actually closes. Instead, pin the
    # axes to its current ABSOLUTE size/position and let the added
    # figure width become pure right-margin space for the labels to
    # extend into, which the data-to-pixel mapping is untouched by.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    orig_w, orig_h = fig.get_size_inches()
    bbox = ax.get_position()  # figure-fraction (x0, y0, width, height)

    max_right_px = max(
        (t.get_window_extent(renderer=renderer).x1 for t in ax.texts),
        default=orig_w * fig.dpi,
    )
    overflow_px = max_right_px - orig_w * fig.dpi
    if overflow_px > 0:
        extra_inches = overflow_px / fig.dpi + 0.2
        new_w = orig_w + extra_inches
        fig.set_size_inches(new_w, orig_h)
        # rescale the axes fractional position so its ABSOLUTE size/
        # placement (in inches) is unchanged -- all new width becomes
        # blank margin to the right of the existing drawing
        scale = orig_w / new_w
        ax.set_position([bbox.x0 * scale, bbox.y0, bbox.width * scale, bbox.height])

    # Extract the VERIFIED rendered order: match drawn text against the
    # known display-label strings rather than assuming a traversal order.
    display_to_tid = {label_map[tid]: tid for tid in target_ids}
    found = []
    for t in ax.texts:
        txt = t.get_text().strip()
        if txt in display_to_tid:
            found.append((t.get_position()[1], display_to_tid[txt]))

    ylim = ax.get_ylim()
    # smaller y = visually higher UNLESS the axis is non-inverted
    found.sort(key=lambda p: p[0], reverse=(ylim[0] < ylim[1]))
    order = [tid for _, tid in found]

    return fig, order


def plot_similarity_heatmap_ordered(target_ids, raw, order, labels=None,
                                     title="STRUCTAL similarity",
                                     cmap="viridis", figsize=None):
    """Static seaborn heatmap using an EXTERNALLY supplied leaf order
    (e.g. from plot_nj_tree) rather than computing its own clustering --
    guarantees this heatmap matches whatever tree `order` came from,
    regardless of what matrix was used to build that tree."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    id_to_idx = {tid: i for i, tid in enumerate(target_ids)}
    idx_order = [id_to_idx[tid] for tid in order]
    n = len(order)

    ordered_raw = raw[np.ix_(idx_order, idx_order)]
    if labels is None:
        labels = target_ids
    label_map = dict(zip(target_ids, labels))
    ordered_labels = [label_map[tid] for tid in order]

    if figsize is None:
        side = max(6.0, 0.28 * n)
        figsize = (side, side)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(ordered_raw, xticklabels=ordered_labels, yticklabels=ordered_labels,
                square=True, cmap=cmap, cbar_kws={"label": "STRUCTAL score"}, ax=ax)
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    ax.tick_params(axis="y", labelrotation=0, labelsize=7)
    plt.tight_layout()
    return fig


# =======================================================================
# (4) Export for iTOL
# =======================================================================

def export_newick(tree, path):
    """
    Write a Bio.Phylo tree (e.g. from build_nj_tree) to Newick format,
    for import into iTOL (https://itol.embl.de) or any other tree viewer.

    Tip names are kept as target_ids -- stable, Newick-safe identifiers
    (no commas/parens/colons to escape). For full display labels (e.g.
    target_id + description), use export_itol_labels to write a separate
    companion annotation file rather than embedding free text in the
    Newick file itself -- this is the standard iTOL workflow, and avoids
    Newick special-character escaping issues entirely.
    """
    if not _HAVE_BIOPYTHON:
        raise ImportError("biopython is required. `pip install biopython`.")
    Phylo.write(tree, path, "newick")


def export_itol_labels(target_ids, labels, path):
    """
    Write an iTOL 'LABELS' annotation file (see
    https://itol.embl.de/help/labels_template.txt) mapping each
    target_id to a display label (e.g. target_id + description).

    Usage in iTOL: upload the Newick file from export_newick, then drag
    and drop this file onto the tree in the browser to rename the leaves
    -- no need to re-upload anything, and the underlying tree/branch
    lengths are untouched.

    `labels` can be a list aligned to target_ids (same convention as
    elsewhere in this module) or a dict keyed by target_id.
    """
    label_map = labels if isinstance(labels, dict) else dict(zip(target_ids, labels))
    with open(path, "w") as f:
        f.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for tid in target_ids:
            label = str(label_map.get(tid, tid)).replace("\t", " ").replace("\n", " ")
            f.write(f"{tid}\t{label}\n")


# =======================================================================
# (5) Coverage alignments for greedy max-coverage representative selection
# =======================================================================

def compute_all_pairs_coverage_alignments(project, view_name, d0=5.0,
                                           radius_factor=2.0, max_targets=1000,
                                           verbose=True, checkpoint_path=None,
                                           checkpoint_every=0.1):
    """
    All-vs-all target-target STRUCTAL-DP rescoring -- the SAME all-pairs
    loop as compute_similarity_matrix -- but instead of discarding the
    matched-residue chain down to a scalar score, accumulate for each
    target the GLOBAL residue indices (into project.coords) that got
    DP-matched in ANY of its pairwise alignments against every other
    target.

    This is the shared-coordinate-space coverage input a greedy
    max-coverage selector needs: column i's bitmap = union, over every
    other target j, of the specific residues of BOTH i and j that
    structurally aligned when i and j were DP-compared. Two different
    targets' bitmaps can and should overlap whenever they both align
    well to the same region of some third target -- that's exactly what
    lets greedy max-coverage pick a small, complementary representative
    set instead of just the largest structures.

    Memory: accumulates into a pyroaring BitMap per target INCREMENTALLY,
    pair by pair, rather than collecting raw numpy arrays from all
    ~N^2/2 pairs and deduplicating only at the end -- that earlier
    approach holds every pair's matched indices in memory simultaneously
    (hundreds of millions of int64 values at N=1000), which is what
    silently OOM-kills the kernel with no Python traceback. Roaring
    bitmaps also compress near-contiguous matched runs (the common case
    here) far more compactly than a plain array ever could, so peak
    memory stays bounded by N compressed bitmaps, not by pair count.

    checkpoint_path : optional str/Path. If given, the partial alignments
    dict is pickled to this path every `checkpoint_every` fraction of
    pairs completed (default every 10%), so a crash doesn't lose all
    prior work -- reload with pickle.load and pass to
    FastGreedyMaxCoverage directly, or re-run to pick up where you'd
    left off (this function does not auto-resume; it re-checkpoints from
    scratch each call, but at least the LAST checkpoint survives a crash).

    Returns
    -------
    dict: target_id -> 1D np.ndarray[uint32] of GLOBAL residue indices
    (into project.coords) matched by this target across all pairwise DP
    alignments. Feed directly into FastGreedyMaxCoverage(target_alignments).
    """
    if not _HAVE_PYROARING:
        raise ImportError(
            "pyroaring is required for memory-safe coverage accumulation. "
            "`pip install pyroaring` (or `pip install pyroaring "
            "--break-system-packages` in Colab)."
        )

    df = project.views[view_name].head(max_targets).copy()
    coords_dict = apply_transforms(project, df, verbose=verbose)
    target_ids = list(coords_dict.keys())
    n = len(target_ids)

    meta = project.master_metadata.set_index("target_id")
    start_offset = {tid: int(meta.loc[tid, "start_idx"]) for tid in target_ids}

    trees = build_trees(coords_dict)
    covered = {tid: BitMap() for tid in target_ids}

    radius = radius_factor * d0
    t0 = time.time()
    n_pairs = n * (n - 1) // 2
    done = 0
    next_checkpoint_frac = checkpoint_every

    for ii in range(n):
        ti = target_ids[ii]
        off_i = start_offset[ti]
        for jj in range(ii + 1, n):
            tj = target_ids[jj]
            off_j = start_offset[tj]
            _, a_idx, b_idx = pairwise_structal_alignment(
                coords_dict[ti], trees[ti], coords_dict[tj], trees[tj],
                d0=d0, radius=radius,
            )
            if len(a_idx):
                # CREDIT THE OTHER SIDE to each column: column ti gets tj's
                # matched positions (in tj's own global range), and column
                # tj gets ti's matched positions (in ti's own global range).
                # This is what makes representative bitmaps span multiple
                # OTHER targets' residues and actually overlap with each
                # other -- crediting each column with its own matched
                # residues (what this said before) confines every target's
                # bitmap to its own disjoint coordinate slice, so no two
                # different targets' bitmaps could ever intersect at all.
                covered[ti].update((b_idx + off_j).astype(np.uint32))
                covered[tj].update((a_idx + off_i).astype(np.uint32))
            done += 1

        frac = done / n_pairs if n_pairs else 1.0
        if verbose and n_pairs > 0 and (ii % max(1, n // 20) == 0):
            elapsed = time.time() - t0
            eta = elapsed / frac - elapsed if frac > 0 else float("nan")
            print(f"[compute_all_pairs_coverage_alignments] {done}/{n_pairs} "
                  f"({frac:.0%}), elapsed {elapsed:.0f}s, ETA {eta:.0f}s")

        if checkpoint_path is not None and frac >= next_checkpoint_frac:
            import pickle
            partial = {tid: np.array(bm, dtype=np.uint32) for tid, bm in covered.items()}
            with open(checkpoint_path, "wb") as f:
                pickle.dump(partial, f)
            if verbose:
                print(f"[compute_all_pairs_coverage_alignments] checkpoint written "
                      f"at {frac:.0%} -> {checkpoint_path}")
            next_checkpoint_frac += checkpoint_every

    return {tid: np.array(bm, dtype=np.uint32) for tid, bm in covered.items()}


def assign_to_representatives_by_score(project, view_name, selected,
                                        target_ids=None, d0=5.0,
                                        radius_factor=2.0, max_targets=1000,
                                        normalize=False, verbose=True):
    """
    Assign every target to whichever of the `selected` k representatives
    it has the highest raw STRUCTAL score against -- a distance-weighted
    structural-fit measure (near-perfect matches contribute close to
    1.0 per residue, marginal matches near the DP's radius cutoff
    contribute much less), generally more informative for cluster
    affinity than bitsum's plain matched-residue COUNT.

    Only computes N x k pairwise comparisons (each target vs each
    representative), not the full N x N similarity matrix -- much
    cheaper than compute_similarity_matrix for this purpose.

    Representatives self-assign with no special-casing needed (unlike
    the bitsum version): score(r,r) = length_r exactly, and any DP chain
    against a DIFFERENT structure can match at most
    min(length_r, length_other) <= length_r residues each contributing
    at most 1.0, so no other representative can ever out-score r against
    itself.

    normalize : bool
        If True, divide each score by min(own_length, rep_length) before
        comparing (bounds scores to roughly [0,1], fraction of the
        smaller structure explained) -- use this if you're worried a
        very large representative could win purely by size rather than
        genuine fit quality. Off by default (raw score, as requested).

    Returns
    -------
    pandas.DataFrame: target_id, representative, score (score actually
    used for the assignment decision -- raw or normalized per
    `normalize`), raw_score, own_length.
    """
    df = project.views[view_name].head(max_targets).copy()
    coords_dict = apply_transforms(project, df, verbose=verbose)
    if target_ids is None:
        target_ids = list(coords_dict.keys())

    trees = build_trees(coords_dict)
    radius = radius_factor * d0

    rows = []
    t0 = time.time()
    n_total = len(target_ids)
    for i, tid in enumerate(target_ids):
        own_length = len(coords_dict[tid])
        best_rep, best_score, best_raw = None, -1.0, 0.0
        for r in selected:
            raw_score = pairwise_structal_score(
                coords_dict[tid], trees[tid], coords_dict[r], trees[r],
                d0=d0, radius=radius,
            )
            score = raw_score / min(own_length, len(coords_dict[r])) if normalize else raw_score
            if score > best_score:
                best_score = score
                best_raw = raw_score
                best_rep = r
        rows.append(dict(target_id=tid, representative=best_rep,
                          score=best_score, raw_score=best_raw,
                          own_length=own_length))
        if verbose and n_total > 20 and i % max(1, n_total // 10) == 0:
            print(f"[assign_to_representatives_by_score] {i}/{n_total}, "
                  f"elapsed {time.time()-t0:.1f}s")

    return pd.DataFrame(rows)


def assign_to_representatives(project, alignments, selected):
    """
    Assign every target to whichever of the `selected` representatives
    (from FastGreedyMaxCoverage.solve()) explains the most of ITS OWN
    residues -- bitsum = how many of target t's own global residue
    positions appear in representative r's coverage bitmap.

    Important: under the corrected compute_all_pairs_coverage_alignments
    semantics, alignments[t] does NOT contain t's own residues at all --
    it contains OTHER targets' residues that t explains. So this can't
    be computed as alignments[t] & alignments[r] (that would compare two
    bitmaps that, for a non-representative t, live in largely disjoint
    territory from t's own identity). Instead, for each non-representative
    target t we build its own full [start_idx, end_idx) range as a
    BitMap and intersect that against each representative's coverage
    bitmap -- i.e. "how many of t's positions did r explain, across all
    pairwise comparisons." Representatives are assigned to themselves
    directly (see below), not via this bitsum comparison -- a
    representative's own bitmap structurally never contains its own
    range (credit always goes to the other side of every comparison it
    appears in), so it could never win a bitsum contest against itself.

    Parameters
    ----------
    project : object with .master_metadata (target_id, start_idx, end_idx)
    alignments : dict target_id -> np.ndarray[uint32] or pyroaring.BitMap
        From compute_all_pairs_coverage_alignments.
    selected : list of target_id
        The k representatives chosen by FastGreedyMaxCoverage.solve().

    Returns
    -------
    pandas.DataFrame with columns: target_id, representative, bitsum,
    own_length (target's own residue count, for context/normalization).
    representative is None (bitsum 0) for a target none of the
    representatives explain any of.
    """
    if not _HAVE_PYROARING:
        raise ImportError("pyroaring is required. `pip install pyroaring`.")

    def to_bitmap(arr):
        return arr if isinstance(arr, BitMap) else BitMap(np.asarray(arr, dtype=np.uint32))

    meta = project.master_metadata.set_index("target_id")
    rep_bitmaps = {r: to_bitmap(alignments[r]) for r in selected}
    selected_set = set(selected)

    rows = []
    for tid in alignments.keys():
        start, end = int(meta.loc[tid, "start_idx"]), int(meta.loc[tid, "end_idx"])
        own_length = end - start

        if tid in selected_set:
            # A representative's own bitmap structurally never contains its
            # OWN range (credit always goes to the OTHER side of every
            # comparison it appears in), so it can never win its own
            # bitsum comparison against itself -- without this, a
            # representative could end up assigned to some OTHER
            # representative it happens to align well with (or unassigned
            # entirely), rather than anchoring its own cluster as intended.
            rows.append(dict(target_id=tid, representative=tid,
                              bitsum=own_length, own_length=own_length))
            continue

        own_range = BitMap(range(start, end))
        best_rep, best_overlap = None, -1
        for r, rep_bm in rep_bitmaps.items():
            overlap = len(rep_bm & own_range)
            if overlap > best_overlap:
                best_overlap = overlap
                best_rep = r
        if best_overlap <= 0:
            best_rep = None
        rows.append(dict(target_id=tid, representative=best_rep,
                          bitsum=max(best_overlap, 0), own_length=own_length))

    return pd.DataFrame(rows)


def print_clusters_with_annotations(assignment_df, parent_df, target_col="target_id",
                                     score_col="bitsum",
                                     annotation_cols=("description", "pfam", "clan")):
    """
    Print each representative's cluster members, annotated by merging
    target_id against parent_df's annotation columns.

    assignment_df : output of assign_to_representatives (score_col
        defaults to "bitsum") or assign_to_representatives_by_score
        (pass score_col="score" for that one).
    parent_df : full (redundant) dataframe carrying target_col plus
        whichever of annotation_cols it has -- missing columns are
        silently skipped rather than raising, since parent_df may not
        carry all three depending on what's been annotated so far.
    """
    annot_cols = [c for c in annotation_cols if c in parent_df.columns]
    missing = [c for c in annotation_cols if c not in parent_df.columns]
    if missing:
        print(f"(note: parent_df is missing columns {missing}, skipping those)")

    merged = assignment_df.merge(parent_df[[target_col] + annot_cols],
                                  on=target_col, how="left")
    is_float = pd.api.types.is_float_dtype(merged[score_col])

    # largest clusters first; UNASSIGNED (representative is NaN) goes last
    cluster_sizes = merged.groupby("representative", dropna=False).size()
    order = [r for r in cluster_sizes.sort_values(ascending=False).index if pd.notna(r)]
    if merged["representative"].isna().any():
        order.append(None)

    for rep in order:
        group = merged[merged["representative"] == rep] if pd.notna(rep) else merged[merged["representative"].isna()]
        label = rep if rep is not None else "UNASSIGNED (no overlap with any representative)"
        print(f"\n=== Cluster: {label}  (n={len(group)}) ===")
        for _, row in group.sort_values(score_col, ascending=False).iterrows():
            score_str = f"{row[score_col]:8.2f}" if is_float else f"{row[score_col]:6d}"
            ann = "  ".join(f"{c}={row[c]}" for c in annot_cols if pd.notna(row.get(c)))
            marker = " <- representative" if row[target_col] == rep else ""
            print(f"  {row[target_col]:15s} {score_col}={score_str}  {ann}{marker}")


def compute_query_coverage_alignments(project, view_name, d0=5.0,
                                       radius_factor=2.0, max_targets=1000,
                                       verbose=True):
    """
    For each target (up to max_targets), re-run the STRUCTAL-DP
    equivalence search between the QUERY and that target's FULL
    superimposed structure -- not just the originally-fixed
    query-defined segments -- and return the set of QUERY residue
    indices each target's DP alignment actually matches.

    Why query positions, not target positions: a greedy max-coverage
    selector needs a residue-ID space that's SHARED across targets so
    bitmaps can meaningfully overlap. Each target's own local coordinate
    indices (into project.coords) are disjoint from every other target's
    by construction, so "coverage" over that space would just degenerate
    into picking the largest structures. Query positions are the one
    coordinate system every target is already superimposed into, so
    they're the natural shared universe here.

    Why re-run the DP (rather than reuse project.segments): this
    recovers any EXTENDED rigid core -- structure the target shares with
    the query beyond what the original fixed alignment captured -- the
    same "how far does the core extend" question from the similarity-
    matrix analysis, now applied query-vs-target instead of target-vs-
    target.

    Returns
    -------
    dict: target_id -> 1D np.ndarray[uint32] of sorted unique QUERY
    residue indices matched by that target. Feed directly into
    FastGreedyMaxCoverage(target_alignments).
    """
    df = project.views[view_name].head(max_targets).copy()
    coords_dict = apply_transforms(project, df, verbose=verbose)

    query_coords = np.asarray(project.query_ca_coords, dtype=np.float64)
    query_tree = cKDTree(query_coords)
    radius = radius_factor * d0

    alignments = {}
    t0 = time.time()
    for i, (tid, coords) in enumerate(coords_dict.items()):
        tree_t = cKDTree(coords)
        _, q_idx, _ = pairwise_structal_alignment(
            query_coords, query_tree, coords, tree_t, d0=d0, radius=radius,
        )
        alignments[tid] = np.unique(q_idx).astype(np.uint32)
        if verbose and len(coords_dict) > 20 and i % max(1, len(coords_dict) // 10) == 0:
            print(f"[compute_query_coverage_alignments] {i}/{len(coords_dict)}, "
                  f"elapsed {time.time()-t0:.1f}s")

    return alignments


# =======================================================================
# (6) Cluster-vs-function agreement: bin-count-invariant metrics
# =======================================================================

def evaluate_cluster_functional_agreement(cluster_df, parent_df, id_col,
                                           cluster_col, functional_col,
                                           target_col="target_id"):
    """
    Quantify how well ANY clustering (infomap modules, STRUCTAL nearest-
    representative assignment, anything with an id + cluster-label
    column) agrees with a functional annotation column -- using metrics
    that are NOT biased by how many clusters or functional categories
    exist. This matters because eyeballing a printed module list can't
    tell "genuinely more functionally pure" apart from "just has fewer
    bins to scatter across" -- a clustering with 14 groups will always
    look tidier at a glance than one with 40, independent of true purity.

    NMI (normalized mutual information) and ARI (adjusted Rand index)
    are both corrected for this: NMI is normalized against the entropy
    of each labeling, and ARI is corrected for chance agreement, so
    neither one rewards having fewer clusters for its own sake.

    Parameters
    ----------
    cluster_df : DataFrame with an id column and a cluster-label column
        (e.g. your modules_df, or the output of
        assign_to_representatives_by_score with cluster_col="representative").
    parent_df : DataFrame with target_col + functional_col
        (e.g. "pfam" or "clan" -- pick ONE; results will differ since
        clan is coarser than pfam, and a single structural module can
        legitimately span multiple Pfam families within one clan).
    id_col : the id column name in cluster_df (may differ from target_col,
        e.g. infomap's own "node_name" or similar).

    Returns
    -------
    (summary: dict, purity_df: DataFrame)
    summary has n_items, n_clusters, n_functional_categories, NMI, ARI,
    homogeneity, completeness, v_measure, weighted_purity.
    purity_df has one row per cluster: n members, majority functional
    label, and purity (fraction matching that majority label).
    """
    from sklearn.metrics import (
        normalized_mutual_info_score, adjusted_rand_score,
        homogeneity_completeness_v_measure,
    )

    merged = cluster_df[[id_col, cluster_col]].rename(columns={id_col: target_col})
    merged = merged.merge(parent_df[[target_col, functional_col]], on=target_col, how="inner")
    merged = merged.dropna(subset=[cluster_col, functional_col])

    labels_true = merged[functional_col].astype(str)
    labels_pred = merged[cluster_col].astype(str)

    nmi = normalized_mutual_info_score(labels_true, labels_pred)
    ari = adjusted_rand_score(labels_true, labels_pred)
    homogeneity, completeness, v_measure = homogeneity_completeness_v_measure(labels_true, labels_pred)

    purity_rows = []
    for cluster_id, group in merged.groupby(cluster_col):
        counts = group[functional_col].value_counts()
        majority_label = counts.index[0]
        purity = counts.iloc[0] / len(group)
        purity_rows.append(dict(cluster=cluster_id, n=len(group),
                                 majority_label=majority_label, purity=purity))
    purity_df = pd.DataFrame(purity_rows).sort_values("n", ascending=False)
    weighted_purity = sum(r["n"] * r["purity"] for r in purity_rows) / len(merged)

    summary = dict(
        n_items=len(merged),
        n_clusters=merged[cluster_col].nunique(),
        n_functional_categories=merged[functional_col].nunique(),
        NMI=nmi, ARI=ari,
        homogeneity=homogeneity, completeness=completeness, v_measure=v_measure,
        weighted_purity=weighted_purity,
    )
    return summary, purity_df


def permutation_test_cluster_agreement(cluster_df, parent_df, id_col,
                                        cluster_col, functional_col,
                                        target_col="target_id", metric="NMI",
                                        n_permutations=1000, random_state=0):
    """
    Test whether cluster_col's agreement with functional_col is
    significantly better than chance, by comparing the observed NMI/ARI
    against a null distribution built by randomly permuting the
    functional labels many times.

    This is the robust complement to evaluate_cluster_functional_agreement:
    individual metrics from that function can each be trivially gamed by
    a degenerate clustering (completeness=1.0 with everyone in one
    cluster; homogeneity=purity=1.0 with every item its own cluster --
    see module notes).

    IMPORTANT LIMITATION (verified empirically, not just reasoned about):
    label-permutation only reorders which item holds which existing
    label value -- it can NEVER change how many groups exist or how big
    they are, on EITHER side. So this test can only ask "does this exact
    assignment beat a random reshuffle within the SAME group-size
    structure" -- it cannot validly judge whether the number of clusters
    itself is reasonable. For an all-singletons clustering specifically,
    the null distribution has EXACTLY ZERO variance (confirmed by direct
    comparison: shuffling either side gives identical, zero-spread null
    distributions when one side is a bijection), so its p-value reflects
    this mathematical degeneracy, not real (non-)significance. A warning
    is raised when n_clusters is 1 or equals n_items for exactly this
    reason -- always sanity-check n_clusters directly (e.g. via
    evaluate_cluster_functional_agreement's summary) rather than relying
    on this test alone to catch degenerate solutions.

    Returns
    -------
    dict: observed, null_mean, null_std, z_score,
    p_value (one-sided: fraction of null permutations >= observed,
    with +1 smoothing), n_permutations.
    """
    from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
    metric_fn = {"NMI": normalized_mutual_info_score, "ARI": adjusted_rand_score}[metric]

    merged = cluster_df[[id_col, cluster_col]].rename(columns={id_col: target_col})
    merged = merged.merge(parent_df[[target_col, functional_col]], on=target_col, how="inner")
    merged = merged.dropna(subset=[cluster_col, functional_col])

    labels_true = merged[functional_col].astype(str).to_numpy()
    labels_pred = merged[cluster_col].astype(str).to_numpy()

    observed = metric_fn(labels_true, labels_pred)

    n_items = len(merged)
    n_clusters = merged[cluster_col].nunique()
    if n_clusters == 1 or n_clusters == n_items:
        import warnings
        warnings.warn(
            f"n_clusters={n_clusters} out of {n_items} items -- this is a "
            "degenerate partition (either everyone in one cluster, or "
            "every item its own cluster). Label-permutation tests CANNOT "
            "validly flag this: permuting either side only reorders which "
            "item holds which existing label, it never changes group "
            "count/sizes, and for an all-singletons clustering the null "
            "distribution has EXACTLY ZERO variance (verified empirically "
            "-- permuting never perturbs the metric when one side is a "
            "bijection), so p/z here reflect this mathematical degeneracy, "
            "not genuine (non-)significance. Do not trust this result; "
            "inspect n_clusters directly instead."
        )

    rng = np.random.default_rng(random_state)
    null_scores = np.empty(n_permutations)
    for i in range(n_permutations):
        null_scores[i] = metric_fn(rng.permutation(labels_true), labels_pred)

    null_mean, null_std = float(null_scores.mean()), float(null_scores.std())
    z = (observed - null_mean) / null_std if null_std > 0 else float("inf")
    p_value = (int(np.sum(null_scores >= observed)) + 1) / (n_permutations + 1)

    return dict(observed=float(observed), null_mean=null_mean, null_std=null_std,
                z_score=float(z), p_value=p_value, n_permutations=n_permutations)


# =======================================================================
# (7) Cross-tabulate max-coverage representatives against infomap modules
# =======================================================================

def cross_tabulate_representatives_modules(project, alignments, selected,
                                            modules_df, id_col, module_col,
                                            target_col="target_id"):
    """
    Reconcile global max-coverage representatives with infomap modules
    as a REPORT, not a constraint: representatives are chosen globally
    and unconstrained (so cross-module structural redundancy -- e.g.
    several functionally-distinct modules sharing real structural
    overlap -- can actually surface), then this cross-tabulates where
    they landed and how far their explanatory reach extends.

    Returns three DataFrames:

    rep_module_df : one row per selected representative
        target_id, module -- which module each global rep happens to
        fall in.

    module_summary_df : one row per module
        module, n_members, total_residues, n_reps_in_module,
        coverage_fraction (residues explained by ANY selected rep),
        endogenous_fraction (explained specifically by a rep that is
        ITSELF a member of this module),
        exogenous_fraction (explained only by a rep from a DIFFERENT
        module -- this is the cross-module-redundancy signal; high
        exogenous_fraction with low n_reps_in_module means this module
        is structurally well-covered by representatives infomap placed
        elsewhere, a candidate for NOT needing a targeted top-up).

    module_x_module_df : long-format explanation matrix
        explaining_module, explained_module, residues_explained --
        for every representative, which modules does ITS coverage reach
        into. explaining_module == explained_module is the endogenous
        diagonal; off-diagonal entries are the cross-module signal
        (e.g. several dioxygenase modules mutually explaining each
        other's residues despite infomap having split them by function).
    """
    if not _HAVE_PYROARING:
        raise ImportError("pyroaring is required. `pip install pyroaring`.")

    def to_bitmap(arr):
        return arr if isinstance(arr, BitMap) else BitMap(np.asarray(arr, dtype=np.uint32))

    meta = project.master_metadata.set_index("target_id")
    module_map = modules_df.set_index(id_col)[module_col].to_dict()

    # each target's own full residue range, and which module it's in
    own_ranges = {}
    target_module = {}
    for tid in alignments.keys():
        start, end = int(meta.loc[tid, "start_idx"]), int(meta.loc[tid, "end_idx"])
        own_ranges[tid] = BitMap(range(start, end))
        target_module[tid] = module_map.get(tid)

    modules = sorted({m for m in target_module.values() if pd.notna(m)})

    # module -> union of its members' own ranges (total residue mass)
    module_range = {m: BitMap() for m in modules}
    module_members = {m: [] for m in modules}
    for tid, m in target_module.items():
        if pd.notna(m):
            module_range[m] |= own_ranges[tid]
            module_members[m].append(tid)

    # --- rep_module_df ---
    rep_module_rows = [dict(target_id=r, module=target_module.get(r)) for r in selected]
    rep_module_df = pd.DataFrame(rep_module_rows)

    # reps grouped by their own module
    reps_by_module = {}
    for r in selected:
        m = target_module.get(r)
        reps_by_module.setdefault(m, []).append(r)

    rep_bitmaps = {r: to_bitmap(alignments[r]) for r in selected}
    global_covered = BitMap()
    for bm in rep_bitmaps.values():
        global_covered |= bm

    # --- module_summary_df ---
    summary_rows = []
    for m in modules:
        total = module_range[m]
        n_total = len(total)
        covered = total & global_covered
        n_covered = len(covered)

        endogenous = BitMap()
        for r in reps_by_module.get(m, []):
            endogenous |= rep_bitmaps[r]
        endog_covered = total & endogenous
        n_endog = len(endog_covered)
        n_exog = n_covered - n_endog

        summary_rows.append(dict(
            module=m,
            n_members=len(module_members[m]),
            total_residues=n_total,
            n_reps_in_module=len(reps_by_module.get(m, [])),
            coverage_fraction=(n_covered / n_total) if n_total else 0.0,
            endogenous_fraction=(n_endog / n_total) if n_total else 0.0,
            exogenous_fraction=(n_exog / n_total) if n_total else 0.0,
        ))
    module_summary_df = pd.DataFrame(summary_rows).sort_values("total_residues", ascending=False)

    # --- module_x_module_df ---
    mxm_rows = []
    for r in selected:
        explaining_module = target_module.get(r)
        bm = rep_bitmaps[r]
        for explained_module in modules:
            overlap = bm & module_range[explained_module]
            if len(overlap) > 0:
                mxm_rows.append(dict(
                    representative=r,
                    explaining_module=explaining_module,
                    explained_module=explained_module,
                    residues_explained=len(overlap),
                ))
    mxm_long = pd.DataFrame(mxm_rows)
    if len(mxm_long):
        module_x_module_df = (mxm_long.groupby(["explaining_module", "explained_module"])
                               ["residues_explained"].sum().reset_index()
                               .sort_values("residues_explained", ascending=False))
    else:
        module_x_module_df = pd.DataFrame(columns=["explaining_module", "explained_module", "residues_explained"])

    return rep_module_df, module_summary_df, module_x_module_df


# =======================================================================
# (8) Module cohesion: within- vs between-module STRUCTAL score sampling
# =======================================================================

def evaluate_module_cohesion(project, view_name, modules_df, target_col="target_id",
                              module_col="module", n_within=1000, n_between=1000,
                              d0=5.0, radius_factor=2.0, modules=None,
                              random_state=0, verbose=True):
    """
    Per module: sample random WITHIN-module pairs (both members of this
    module) and random BETWEEN-module pairs (one member of this module,
    one from the POOL of everyone in every other module), and compute
    fresh STRUCTAL-DP scores for each sampled pair.

    This is a genuine out-of-sample check, not circular, even if the
    module assignment itself came from a STRUCTAL-derived graph: a
    module typically has far more members than any single node's k1
    neighbor list, so most randomly-sampled within-module pairs were
    NEVER directly compared during graph construction. Testing whether
    THOSE pairs still score highly is a real test of whether community
    membership predicts structural similarity beyond what was used to
    assign it.

    Sampling handles small modules correctly: if a module's total
    possible within-pairs (n_members choose 2) is below n_within, ALL
    pairs are used (no sampling, no replacement) rather than erroring
    or duplicating. Between-pairs are drawn i.i.d. (member, non-member)
    -- the non-member pool is large enough in practice that duplicate
    pairs are negligible for a descriptive/distributional evaluation
    like this.

    Cost warning: this is expensive at full scale. 10,000+10,000 pairs
    across ~243 modules is ~4.9M STRUCTAL-DP calls -- likely hours. Use
    `modules=` to restrict to a specific subset (recommended default
    workflow: evaluate the modules you're actually interested in, not
    all of them blindly), and/or lower n_within/n_between from the
    defaults here.

    Parameters
    ----------
    project, view_name : as elsewhere -- used to fetch R/t and apply
        transforms for every target appearing in modules_df.
    modules_df : DataFrame with target_col, module_col (e.g. from
        infomap's node/module output).
    n_within, n_between : max pairs sampled per module (per-category;
        actual within-count is capped at C(n_members,2) automatically).
    modules : optional list of module ids to evaluate. None = all
        modules in modules_df (see cost warning above).

    Returns
    -------
    (summary_df, raw_scores) where summary_df has one row per module
    (module, n_members, n_within_sampled, n_between_sampled,
    within_mean/median/std, between_mean/median/std, separation =
    within_median - between_median), and raw_scores is
    {module: {"within": array, "between": array}} for plotting
    (e.g. histograms/violins) beyond the summary statistics.
    """
    rng = np.random.default_rng(random_state)

    df = project.views[view_name]
    df = df[df[target_col].isin(modules_df[target_col])].copy()
    coords_dict = apply_transforms(project, df, target_col=target_col, verbose=verbose)
    trees = build_trees(coords_dict)
    radius = radius_factor * d0

    module_map = modules_df.set_index(target_col)[module_col].to_dict()
    all_targets = [t for t in coords_dict.keys() if t in module_map]
    members_by_module = {}
    for t in all_targets:
        members_by_module.setdefault(module_map[t], []).append(t)

    eval_modules = modules if modules is not None else sorted(members_by_module.keys())

    if verbose:
        total_est = len(eval_modules) * (n_within + n_between)
        print(f"[evaluate_module_cohesion] evaluating {len(eval_modules)} modules, "
              f"up to {n_within}+{n_between} pairs each (~{total_est:,} STRUCTAL calls "
              f"worst case -- actual may be less for small modules)")

    def score_pair(ti, tj):
        return pairwise_structal_score(coords_dict[ti], trees[ti],
                                        coords_dict[tj], trees[tj],
                                        d0=d0, radius=radius)

    summary_rows = []
    raw_scores = {}

    for mi, m in enumerate(eval_modules):
        members = members_by_module.get(m, [])
        non_members = [t for t in all_targets if module_map[t] != m]
        n_members = len(members)

        # --- within-module pairs ---
        max_within_pairs = n_members * (n_members - 1) // 2
        if max_within_pairs <= n_within:
            within_pairs = [(members[a], members[b])
                             for a in range(n_members) for b in range(a + 1, n_members)]
        else:
            within_pairs = set()
            while len(within_pairs) < n_within:
                a, b = rng.choice(n_members, size=2, replace=False)
                within_pairs.add((members[a], members[b]) if a < b else (members[b], members[a]))
            within_pairs = list(within_pairs)

        within_scores = np.array([score_pair(ti, tj) for ti, tj in within_pairs])

        # --- between-module pairs (pooled across ALL other modules) ---
        n_between_actual = min(n_between, n_members * len(non_members)) if non_members else 0
        between_scores = np.empty(n_between_actual)
        for k in range(n_between_actual):
            ti = members[rng.integers(n_members)]
            tj = non_members[rng.integers(len(non_members))]
            between_scores[k] = score_pair(ti, tj)

        summary_rows.append(dict(
            module=m, n_members=n_members,
            n_within_sampled=len(within_scores), n_between_sampled=len(between_scores),
            within_mean=float(within_scores.mean()) if len(within_scores) else np.nan,
            within_median=float(np.median(within_scores)) if len(within_scores) else np.nan,
            within_std=float(within_scores.std()) if len(within_scores) else np.nan,
            between_mean=float(between_scores.mean()) if len(between_scores) else np.nan,
            between_median=float(np.median(between_scores)) if len(between_scores) else np.nan,
            between_std=float(between_scores.std()) if len(between_scores) else np.nan,
            separation=(float(np.median(within_scores) - np.median(between_scores))
                        if len(within_scores) and len(between_scores) else np.nan),
        ))
        raw_scores[m] = {"within": within_scores, "between": between_scores}

        if verbose and (mi % max(1, len(eval_modules) // 20) == 0):
            print(f"[evaluate_module_cohesion] {mi+1}/{len(eval_modules)} modules done")

    summary_df = pd.DataFrame(summary_rows).sort_values("separation", ascending=False)
    return summary_df, raw_scores
