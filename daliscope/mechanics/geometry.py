# ------------------------------------------------------------------ #
# Pure function — no self, fully testable in isolation               #
# ------------------------------------------------------------------ #

import numpy as np

def _parse_query_ca_coords(query_pdb_str) -> np.ndarray:
    """
    Extract CA coordinates from self.query_pdb_str.
    Returns (query_length, 3) float32 array in original (non-centroid-subtracted) frame.
    """
    if query_pdb_str is None:
        raise ValueError("query_pdb_str is None — Query.pdb not found in archive")

    coords = []
    for line in query_pdb_str.splitlines():
        #if line[:4] == "ATOM" and line[12:16].strip() == "CA":
        # alt location index blank or A
        if line[12:16].strip() == "CA" and line[16] in [' ','A']:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append((x, y, z))

    ca = np.array(coords, dtype=np.float32)
    if len(ca) == 0:
        raise ValueError("No CA atoms found in Query.pdb")
    return ca


def svd_superimpose(P: np.ndarray, Q: np.ndarray):
    """
    Superimpose P onto Q using the Kabsch/SVD algorithm.

    Parameters
    ----------
    P : (n, 3) float32  — subject (target domain) CA coordinates
    Q : (n, 3) float32  — reference (query) CA coordinates

    Returns
    -------
    R    : (3, 3) rotation matrix
    t    : (3,)   translation vector
    rmsd : float  RMSD after superimposition

    Such that:  Q ≈ (R @ P.T).T + t
    Requires n >= 3.
    """
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"P and Q must both be (n,3), got {P.shape} and {Q.shape}")
    if len(P) < 3:
        raise ValueError(f"Need at least 3 aligned points, got {len(P)}")

    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)

    H        = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)

    # correct for improper rotation (reflection)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, float(d)])
    R = (Vt.T @ D @ U.T).astype(np.float32)
    t = (Q.mean(axis=0) - R @ P.mean(axis=0)).astype(np.float32)

    diff = (R @ P.T).T + t - Q
    rmsd = float(np.sqrt((diff ** 2).sum(axis=1).mean()))

    return R, t, rmsd
