from typing import Tuple, Optional
import numpy as np
from scipy.spatial import Delaunay, QhullError
from scipy import sparse
from numba import njit

@njit(fastmath=True)
def _extract_unique_delaunay_keys(simplices: np.ndarray) -> np.ndarray:
    """
    Directly extracts unique 64-bit packed edges (lo << 32 | hi) from tetrahedra
    using a Numba typed set. Avoids np.unique and sorting overhead entirely.
    """
    unique_keys = set()
    n_tets = simplices.shape[0]

    # 6 edges per tetrahedron: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    pairs = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32)

    for i in range(n_tets):
        tet = simplices[i]
        for p in range(6):
            u = tet[pairs[p, 0]]
            v = tet[pairs[p, 1]]
            if u < v:
                lo, hi = np.uint64(u), np.uint64(v)
            else:
                lo, hi = np.uint64(v), np.uint64(u)

            key = (lo << 32) | hi
            unique_keys.add(key)

    out = np.empty(len(unique_keys), dtype=np.uint64)
    idx = 0
    for k in unique_keys:
        out[idx] = k
        idx += 1
    return out


def delaunay_edges(
    points: np.ndarray,
    return_sparse: bool = True,
    dist_cutoff: Optional[float] = None,
    normalize: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Optional[sparse.coo_matrix]]:
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must be an (N, 3) array")

    pts = np.asarray(points, dtype=float)
    if normalize:
        mu = pts.mean(axis=0)
        sigma = pts.std(axis=0)
        sigma[sigma == 0.0] = 1.0
        pts_norm = (pts - mu) / sigma
    else:
        pts_norm = pts

    # 1. Run Qhull Delaunay Tessellation with fast jiggle flag
    try:
        tri = Delaunay(pts_norm, qhull_options="QJ")
    except QhullError as e:
        raise RuntimeError(
            "Delaunay (Qhull) failed. Data may be coplanar or degenerate: " + str(e)
        )

    simplices = tri.simplices  # Shape: (n_tets, 4)
    if simplices.size == 0:
        edges = np.empty((0, 2), dtype=np.int32)
        distances = np.empty((0,), dtype=float)
        adj = sparse.coo_matrix(([], ([], [])), shape=(pts.shape[0], pts.shape[0])) if return_sparse else None
        return edges, distances, adj

    # 2. Extract unique 64-bit keys directly via Numba JIT
    keys = _extract_unique_delaunay_keys(simplices)

    # 3. Unpack keys directly into (M, 2) edge array
    edges = np.empty((len(keys), 2), dtype=np.int32)
    edges[:, 0] = (keys >> 32).astype(np.int32)
    edges[:, 1] = (keys & 0xFFFFFFFF).astype(np.int32)

    # 4. Lazy Distance Calculation
    need_distances = return_sparse or (dist_cutoff is not None)

    if need_distances:
        diffs = pts[edges[:, 0]] - pts[edges[:, 1]]
        sq_dists = (diffs * diffs).sum(axis=1)

        if dist_cutoff is not None:
            dist_mask = sq_dists <= (float(dist_cutoff) ** 2)
            edges = edges[dist_mask]
            distances = np.sqrt(sq_dists[dist_mask])
        else:
            distances = np.sqrt(sq_dists)
    else:
        distances = np.empty((0,), dtype=float)

    # 5. Sparse Matrix Creation
    adj = None
    if return_sparse:
        n = pts.shape[0]
        m = len(edges)

        rows = np.empty(2 * m, dtype=np.int32)
        cols = np.empty(2 * m, dtype=np.int32)
        data = np.empty(2 * m, dtype=distances.dtype)

        rows[:m] = edges[:, 0]
        rows[m:] = edges[:, 1]

        cols[:m] = edges[:, 1]
        cols[m:] = edges[:, 0]

        data[:m] = distances
        data[m:] = distances

        adj = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))

    return edges, distances, adj
