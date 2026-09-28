"""
structural_fingerprint_patch.py

Three concrete changes to structural_fingerprint.py, each verified
independently against the original before being combined:

1. norms computation (in build_structural_fingerprint_index): Python loop
   over n_targets -> vectorized segment-sum via cumsum. VERIFIED bit-
   identical to the original (including the empty-segment/zero-feature
   target edge case) on synthetic data.

2. build_knn_graph's symmetrization: Python loop over up to N*k unique
   undirected edges (picking the max-scoring direction per duplicate
   pair) -> fully vectorized via sort + first-occurrence, same technique
   used for pair-deduplication elsewhere in this project. VERIFIED
   identical edge sets against the original, and benchmarked: 5.5x
   speedup in isolation at N=10,000/k=200 scale (2M directed edges,
   4.46s -> 0.81s). Neither change increases memory -- both operate on
   the SAME arrays the original already builds.

3. topk's outer per-target loop -> parallelized via ThreadPoolExecutor.
   Each target's computation only READS shared index/raw arrays and
   writes to its own disjoint output row, so this is safe without any
   locking. VERIFIED correctness (identical neighbors/scores to serial).
   NOT verified for actual speedup: this dev sandbox has exactly 1 CPU
   core, so no thread pool can show real parallelism here regardless of
   correctness. The numpy operations used (searchsorted, bincount,
   argpartition, boolean masking on large arrays) release the GIL during
   their C-level execution, which is the theoretical basis for expecting
   a real speedup on a genuine multi-core machine -- but BENCHMARK THIS
   YOURSELF on your actual environment before relying on it; I cannot
   confirm the speedup claim from here.

   Threading (not multiprocessing) is used deliberately: it has ZERO
   memory-duplication risk, since threads share the process's address
   space natively. multiprocessing.ProcessPoolExecutor is an alternative
   if threading doesn't give enough speedup on your machine (e.g. if too
   much time is spent in Python-level orchestration between numpy calls,
   which stays GIL-bound) -- but it comes with a real memory landmine:
   it MUST use the "fork" start method (Linux default) so the large
   shared index arrays are copy-on-write rather than pickled/duplicated
   per worker. If your environment's multiprocessing default is "spawn"
   (macOS, Windows, or newer Python defaults on some Linux configs),
   each worker process would get its OWN COPY of index_hashes/
   index_targets/index_weights/raw_hashes/raw_weights -- multiplying
   memory by worker count, which directly fights the memory constraint
   this whole module was designed around. Explicitly set the start
   method (multiprocessing.get_context("fork")) rather than relying on
   the platform default if you go this route.

Usage
-----
    from structural_fingerprint import build_structural_fingerprint_index
    from structural_fingerprint_patch import topk_threaded, build_knn_graph_fast

    idx = build_structural_fingerprint_index(...)  # unchanged
    neighbors, scores = topk_threaded(idx, k=200, max_workers=8)  # benchmark max_workers yourself
    source, target, score = build_knn_graph_fast(neighbors, scores)
"""

import numpy as np
from concurrent.futures import ThreadPoolExecutor


# =======================================================================
# 1. Vectorized norms (drop-in replacement for the Python loop inside
#    build_structural_fingerprint_index)
# =======================================================================

def compute_norms_vectorized(target_offsets: np.ndarray,
                              raw_weights: np.ndarray) -> np.ndarray:
    """
    Per-target L2 norm of raw_weights, vectorized via a cumsum-based
    segment sum (correctly handles zero-feature/empty-segment targets,
    verified against the original Python-loop version).
    """
    sq = raw_weights.astype(np.float64) ** 2
    cum = np.concatenate(([0.0], np.cumsum(sq)))
    seg_sums = cum[target_offsets[1:]] - cum[target_offsets[:-1]]
    return np.sqrt(seg_sums).astype(np.float32)


# =======================================================================
# 2. Vectorized build_knn_graph symmetrization
# =======================================================================

def build_knn_graph_fast(neighbors: np.ndarray, scores: np.ndarray,
                          min_score: float = 0.0, symmetrize: bool = True):
    """
    Drop-in replacement for build_knn_graph, taking topk's raw
    (neighbors, scores) arrays directly (skips the intermediate
    candidate_pairs() call, same data either way).

    Symmetrization is vectorized: sort directed edges by
    (lo, hi, -score) so the best-scoring copy of each undirected pair
    comes first, then keep just the first occurrence per pair -- no
    per-group Python loop. VERIFIED identical edge sets to the
    original; 5.5x faster in isolation at N=10000/k=200.
    """
    n, k = neighbors.shape
    source = np.repeat(np.arange(n, dtype=np.int32), k)
    target = neighbors.ravel()
    score = scores.ravel()
    keep = (target >= 0) & (score > min_score) if min_score > 0 else (target >= 0)
    source, target, score = source[keep], target[keep], score[keep]

    if not symmetrize:
        return source, target, score

    lo = np.minimum(source, target)
    hi = np.maximum(source, target)

    if lo.size == 0:
        return (np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float32))

    order = np.lexsort((-score, hi, lo))
    lo_s, hi_s, score_s = lo[order], hi[order], score[order]

    is_new_pair = np.empty(lo_s.size, dtype=bool)
    is_new_pair[0] = True
    is_new_pair[1:] = (lo_s[1:] != lo_s[:-1]) | (hi_s[1:] != hi_s[:-1])

    out_source = lo_s[is_new_pair].astype(np.int32, copy=False)
    out_target = hi_s[is_new_pair].astype(np.int32, copy=False)
    out_score = score_s[is_new_pair].astype(np.float32, copy=False)
    return out_source, out_target, out_score


# =======================================================================
# 3. Parallelized topk (threaded)
# =======================================================================

def _topk_one_target(idx, target_idx: int, k: int, min_score: float, n: int):
    """Single-target body, identical logic to the original topk()'s loop
    body -- extracted so it can be dispatched to a thread pool. Reads
    only shared, immutable arrays on `idx`; writes nothing shared."""
    start = int(idx.target_offsets[target_idx])
    end = int(idx.target_offsets[target_idx + 1])
    source_hashes = idx.raw_hashes[start:end]
    source_weights = idx.raw_weights[start:end]
    keep = source_weights > idx.weight_epsilon
    empty = (np.full(k, -1, dtype=np.int32), np.zeros(k, dtype=np.float32))
    if not np.any(keep):
        return target_idx, *empty
    source_hashes = source_hashes[keep]
    source_weights = source_weights[keep]

    query_hashes = (source_hashes[:, None] + idx.kernel_delta_hashes[None, :]).ravel()
    query_weights = (source_weights[:, None] * idx.kernel_weights[None, :]).ravel()

    n_index = idx.index_hashes.size
    left = np.searchsorted(idx.index_hashes, query_hashes, side="left")
    hit = (left < n_index) & (idx.index_hashes[np.minimum(left, n_index - 1)] == query_hashes)
    if not np.any(hit):
        return target_idx, *empty

    hit_left = left[hit]
    hit_hashes = query_hashes[hit]
    hit_source_weights = query_weights[hit]
    right = np.searchsorted(idx.index_hashes, hit_hashes, side="right")
    multiplicity = right - hit_left
    total_hits = int(multiplicity.sum())
    if total_hits == 0:
        return target_idx, *empty

    repeated_left = np.repeat(hit_left, multiplicity)
    repeated_starts = np.repeat(np.cumsum(multiplicity) - multiplicity, multiplicity)
    local_offsets = np.arange(total_hits, dtype=np.int64) - repeated_starts
    index_positions = repeated_left + local_offsets
    candidate_targets = idx.index_targets[index_positions]
    candidate_weights = idx.index_weights[index_positions]
    source_contrib = np.repeat(hit_source_weights, multiplicity)
    contributions = (source_contrib * candidate_weights).astype(np.float32, copy=False)

    candidate_scores = np.bincount(candidate_targets, weights=contributions,
                                    minlength=n).astype(np.float32, copy=False)
    candidate_scores[target_idx] = 0.0

    denom = idx.norms[target_idx] * idx.norms
    valid = (candidate_scores > min_score) & (denom > 0)
    if not np.any(valid):
        return target_idx, *empty
    candidate_scores[valid] /= denom[valid]
    candidate_scores[~valid] = 0.0

    n_valid = int(np.count_nonzero(candidate_scores > min_score))
    if n_valid == 0:
        return target_idx, *empty

    take = min(k, n_valid)
    valid_indices = np.flatnonzero(candidate_scores > min_score)
    if valid_indices.size > take:
        partial = np.argpartition(candidate_scores[valid_indices], -take)[-take:]
        selected = valid_indices[partial]
    else:
        selected = valid_indices
    order = np.argsort(candidate_scores[selected])[::-1]
    selected = selected[order]

    m = selected.size
    nb = np.full(k, -1, dtype=np.int32)
    sc = np.zeros(k, dtype=np.float32)
    nb[:m] = selected.astype(np.int32, copy=False)
    sc[:m] = candidate_scores[selected]
    return target_idx, nb, sc


def topk_threaded(idx, k: int = 200, min_score: float = 0.0, max_workers: int = 8):
    """
    Threaded drop-in for StructuralFingerprintIndex.topk(). Correctness
    VERIFIED identical to the serial version. Speedup NOT verified in
    this environment (1 CPU core available) -- benchmark max_workers on
    your actual machine; the numpy ops used here release the GIL during
    C-level execution, which is why this should parallelize in principle
    on a real multi-core machine, but "should in principle" is not the
    same as "measured," and I couldn't measure it here.

    Zero memory-duplication risk (unlike multiprocessing): threads share
    the process's address space, no fork/spawn/pickling involved.
    """
    n = idx.n_targets
    if n < 2:
        return (np.empty((n, 0), dtype=np.int32), np.empty((n, 0), dtype=np.float32))
    k = min(int(k), n - 1)
    if k <= 0:
        return (np.empty((n, 0), dtype=np.int32), np.empty((n, 0), dtype=np.float32))

    neighbors = np.full((n, k), -1, dtype=np.int32)
    scores_out = np.zeros((n, k), dtype=np.float32)
    if idx.index_hashes.size == 0:
        return neighbors, scores_out

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for target_idx, nb, sc in ex.map(
            lambda t: _topk_one_target(idx, t, k, min_score, n), range(n)
        ):
            neighbors[target_idx] = nb
            scores_out[target_idx] = sc

    return neighbors, scores_out
