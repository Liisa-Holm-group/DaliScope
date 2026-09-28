# daliscope/analysis/structural_fingerprint.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Voxel convention
# ---------------------------------------------------------------------------

CELL_SIZE_DEFAULT = 4.5

# Must match daliscope.analysis.voxel
OFFSET = 1 << 16
SHIFT_Y = 17
SHIFT_X = 34
CELL_MASK = (1 << 17) - 1

MAX_CELL_INDEX = OFFSET - 1
MIN_CELL_INDEX = -OFFSET


def _hash_cells(cell_xyz: np.ndarray) -> np.ndarray:
    """
    Hash integer voxel coordinates using the same convention as voxel.py.

    Existing convention:

        (x + OFFSET) * 2^34
      + (y + OFFSET) * 2^17
      + (z + OFFSET)

    Returns int64 hashes. Protein-scale coordinates are comfortably within
    the signed-int64 range with the existing 17-bit-per-axis representation.
    """
    cell_xyz = np.asarray(cell_xyz)

    if cell_xyz.ndim != 2 or cell_xyz.shape[1] != 3:
        raise ValueError("cell_xyz must have shape (N, 3)")

    cell_xyz = cell_xyz.astype(np.int64, copy=False)

    if (
        np.any(cell_xyz < MIN_CELL_INDEX)
        or np.any(cell_xyz > MAX_CELL_INDEX)
    ):
        raise ValueError(
            "Voxel coordinate outside the existing 17-bit hash range: "
            f"[{MIN_CELL_INDEX}, {MAX_CELL_INDEX}]"
        )

    x = cell_xyz[:, 0] + OFFSET
    y = cell_xyz[:, 1] + OFFSET
    z = cell_xyz[:, 2] + OFFSET

    return (
        x * (1 << SHIFT_X)
        + y * (1 << SHIFT_Y)
        + z
    ).astype(np.int64, copy=False)


def _hash_delta(dx: int, dy: int, dz: int) -> int:
    """Hash difference corresponding to a voxel displacement."""
    return (
        dx * (1 << SHIFT_X)
        + dy * (1 << SHIFT_Y)
        + dz
    )


# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------

def make_cross7_kernel(
    neighbor_weight: float = 0.25,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Seven-point local kernel:

        center:       1
        +/- x/y/z:    neighbor_weight

    The weights are normalized to sum to one.

    This is the recommended kernel for the ~10M-feature regime because it
    requires only seven inverted-index probes per occupied voxel.
    """
    if neighbor_weight < 0:
        raise ValueError("neighbor_weight must be >= 0")

    offsets = np.asarray(
        [
            [0, 0, 0],
            [1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 1],
            [0, 0, -1],
        ],
        dtype=np.int8,
    )

    weights = np.asarray(
        [
            1.0,
            neighbor_weight,
            neighbor_weight,
            neighbor_weight,
            neighbor_weight,
            neighbor_weight,
            neighbor_weight,
        ],
        dtype=np.float32,
    )

    weights /= weights.sum()

    return offsets, weights


def make_gaussian27_kernel(
    cell_size: float = CELL_SIZE_DEFAULT,
    sigma: Optional[float] = None,
    radius_voxels: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Truncated isotropic Gaussian kernel.

    This produces up to 27 offsets for radius_voxels=1.

    It is useful for benchmarking, but the cross7 kernel is substantially
    cheaper for very large datasets.
    """
    if radius_voxels < 0:
        raise ValueError("radius_voxels must be >= 0")

    if sigma is None:
        sigma = cell_size / 2.0

    if sigma <= 0:
        raise ValueError("sigma must be > 0")

    offsets: List[Tuple[int, int, int]] = []
    weights: List[float] = []

    for dx in range(-radius_voxels, radius_voxels + 1):
        for dy in range(-radius_voxels, radius_voxels + 1):
            for dz in range(-radius_voxels, radius_voxels + 1):
                r = cell_size * np.sqrt(
                    dx * dx + dy * dy + dz * dz
                )

                w = np.exp(-0.5 * (r / sigma) ** 2)

                offsets.append((dx, dy, dz))
                weights.append(float(w))

    offsets_arr = np.asarray(offsets, dtype=np.int8)
    weights_arr = np.asarray(weights, dtype=np.float32)
    weights_arr /= weights_arr.sum()

    return offsets_arr, weights_arr


def make_kernel(
    kind: str = "cross7",
    *,
    cell_size: float = CELL_SIZE_DEFAULT,
    neighbor_weight: float = 0.25,
    sigma: Optional[float] = None,
    radius_voxels: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Construct the kernel used for soft voxel matching.

    Parameters
    ----------
    kind
        "cross7" or "gaussian".
    """
    kind = kind.lower()

    if kind == "cross7":
        return make_cross7_kernel(neighbor_weight)

    if kind in {"gaussian", "gaussian27"}:
        return make_gaussian27_kernel(
            cell_size=cell_size,
            sigma=sigma,
            radius_voxels=radius_voxels,
        )

    raise ValueError(
        f"Unknown kernel {kind!r}; expected 'cross7' or 'gaussian'."
    )


# ---------------------------------------------------------------------------
# IDF
# ---------------------------------------------------------------------------

def capped_idf(
    occupancy: np.ndarray,
    n_targets: int,
    *,
    alpha: float = 1.0,
    cap: float = 2.0,
) -> np.ndarray:
    """
    Smoothed, capped inverse-document-frequency.

        IDF(v) = min(
            log((N + alpha) / (n_v + alpha)),
            cap
        )

    occupancy is the number of distinct targets containing voxel v.

    Because each target contributes at most one entry per voxel, occupancy is
    simply the number of entries for that voxel in the global inverted index.
    """
    if n_targets <= 0:
        raise ValueError("n_targets must be > 0")

    if alpha <= 0:
        raise ValueError("alpha must be > 0")

    if cap <= 0:
        raise ValueError("cap must be > 0")

    occupancy = np.asarray(occupancy, dtype=np.float32)

    result = np.log(
        (float(n_targets) + alpha)
        / (occupancy + alpha)
    )

    np.minimum(result, cap, out=result)

    return result.astype(np.float32, copy=False)


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class StructuralFingerprintIndex:
    """
    Memory-conscious structural fingerprint index.

    Raw representation
    ------------------
    raw_hashes/raw_weights are grouped by target. Each occupied voxel occurs
    exactly once per target.

    Inverted representation
    ------------------------
    index_hashes/index_targets/index_weights are sorted by voxel hash and are
    used for fast spatial candidate generation.

    No globally kernel-expanded fingerprint is stored.
    """

    target_ids: List[str]

    # Target -> raw feature ranges.
    target_offsets: np.ndarray

    # Raw per-target voxel representation.
    raw_hashes: np.ndarray
    raw_weights: np.ndarray

    # Inverted voxel index.
    index_hashes: np.ndarray
    index_targets: np.ndarray
    index_weights: np.ndarray

    # Per-target L2 norm of the weighted raw fingerprint.
    norms: np.ndarray

    # Kernel.
    kernel_offsets: np.ndarray
    kernel_weights: np.ndarray
    kernel_delta_hashes: np.ndarray

    # Parameters.
    cell_size: float
    radial_r0: Optional[float]
    query_com: Optional[np.ndarray]

    # Diagnostics.
    idf_alpha: float
    idf_cap: float
    weight_epsilon: float

    # Optional diagnostic vocabulary.
    vocabulary_hashes: Optional[np.ndarray] = None
    vocabulary_occupancy: Optional[np.ndarray] = None
    vocabulary_idf: Optional[np.ndarray] = None

    @property
    def n_targets(self) -> int:
        return len(self.target_ids)

    @property
    def n_raw_features(self) -> int:
        return self.raw_hashes.size

    @property
    def n_index_features(self) -> int:
        return self.index_hashes.size

    def target_index(self, target_id: str) -> int:
        return self.target_ids.index(target_id)

    def target_range(self, target_idx: int) -> Tuple[int, int]:
        return (
            int(self.target_offsets[target_idx]),
            int(self.target_offsets[target_idx + 1]),
        )

    def topk(
        self,
        k: int = 200,
        *,
        min_score: float = 0.0,
        return_scores: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Retrieve top-k structural candidates for every target.

        Returns
        -------
        neighbors
            int32 array, shape (N, k), containing target indices.
            -1 indicates no candidate.

        scores
            float32 array, shape (N, k), containing normalized soft
            voxel-overlap scores.

        Notes
        -----
        This is deliberately a candidate-generation score, not an
        approximation of STRUCTAL. STRUCTAL should still rerank the
        resulting candidates.
        """
        n = self.n_targets

        if n < 2:
            return (
                np.empty((n, 0), dtype=np.int32),
                np.empty((n, 0), dtype=np.float32),
            )

        k = min(int(k), n - 1)

        if k <= 0:
            return (
                np.empty((n, 0), dtype=np.int32),
                np.empty((n, 0), dtype=np.float32),
            )

        neighbors = np.full(
            (n, k),
            -1,
            dtype=np.int32,
        )

        scores_out = np.zeros(
            (n, k),
            dtype=np.float32,
        )

        n_index = self.index_hashes.size

        if n_index == 0:
            return neighbors, scores_out

        # ------------------------------------------------------------------
        # Process one target at a time.
        #
        # This is intentional: it keeps temporary candidate arrays O(K * M_t)
        # rather than O(N_targets * K * M_t).
        # ------------------------------------------------------------------

        for target_idx in range(n):
            start = int(self.target_offsets[target_idx])
            end = int(self.target_offsets[target_idx + 1])

            source_hashes = self.raw_hashes[start:end]
            source_weights = self.raw_weights[start:end]

            keep = source_weights > self.weight_epsilon

            if not np.any(keep):
                continue

            source_hashes = source_hashes[keep]
            source_weights = source_weights[keep]

            # Generate spatially tolerant voxel queries.
            #
            # Shape before flattening:
            #     n_source_voxels x n_kernel_offsets
            #
            query_hashes = (
                source_hashes[:, None]
                + self.kernel_delta_hashes[None, :]
            ).ravel()

            query_weights = (
                source_weights[:, None]
                * self.kernel_weights[None, :]
            ).ravel()

            # Find matching voxel ranges in the inverted index.
            left = np.searchsorted(
                self.index_hashes,
                query_hashes,
                side="left",
            )

            hit = (
                (left < n_index)
                & (self.index_hashes[
                    np.minimum(left, n_index - 1)
                ] == query_hashes)
            )

            if not np.any(hit):
                continue

            hit_left = left[hit]
            hit_hashes = query_hashes[hit]
            hit_source_weights = query_weights[hit]

            right = np.searchsorted(
                self.index_hashes,
                hit_hashes,
                side="right",
            )

            multiplicity = right - hit_left

            total_hits = int(multiplicity.sum())

            if total_hits == 0:
                continue

            # Expand matched voxel ranges.
            #
            # For the usual case where a voxel is present in only a few
            # targets, this remains very small.
            repeated_left = np.repeat(
                hit_left,
                multiplicity,
            )

            repeated_starts = np.repeat(
                np.cumsum(multiplicity) - multiplicity,
                multiplicity,
            )

            local_offsets = (
                np.arange(total_hits, dtype=np.int64)
                - repeated_starts
            )

            index_positions = (
                repeated_left + local_offsets
            )

            candidate_targets = self.index_targets[
                index_positions
            ]

            candidate_weights = self.index_weights[
                index_positions
            ]

            source_contrib = np.repeat(
                hit_source_weights,
                multiplicity,
            )

            contributions = (
                source_contrib
                * candidate_weights
            ).astype(np.float32, copy=False)

            # Accumulate candidate scores.
            candidate_scores = np.bincount(
                candidate_targets,
                weights=contributions,
                minlength=n,
            ).astype(np.float32, copy=False)

            # Self is not a candidate.
            candidate_scores[target_idx] = 0.0

            # Normalize by raw fingerprint norms.
            #
            # This is deliberately cheap. The score is only used for
            # candidate generation; STRUCTAL performs the final comparison.
            denom = self.norms[target_idx] * self.norms

            valid = (
                (candidate_scores > min_score)
                & (denom > 0)
            )

            if not np.any(valid):
                continue

            candidate_scores[valid] /= denom[valid]
            candidate_scores[~valid] = 0.0

            n_valid = int(np.count_nonzero(candidate_scores > min_score))

            if n_valid == 0:
                continue

            take = min(k, n_valid)

            # Select top-k without sorting all N targets.
            valid_indices = np.flatnonzero(
                candidate_scores > min_score
            )

            if valid_indices.size > take:
                partial = np.argpartition(
                    candidate_scores[valid_indices],
                    -take,
                )[-take:]

                selected = valid_indices[partial]
            else:
                selected = valid_indices

            order = np.argsort(
                candidate_scores[selected]
            )[::-1]

            selected = selected[order]

            m = selected.size

            neighbors[target_idx, :m] = selected.astype(
                np.int32,
                copy=False,
            )

            scores_out[target_idx, :m] = candidate_scores[
                selected
            ]

        if return_scores:
            return neighbors, scores_out

        return neighbors, np.empty((0, 0), dtype=np.float32)

    def candidate_pairs(
        self,
        k: int = 200,
        *,
        min_score: float = 0.0,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return directed k-NN candidate edges.

        Returns
        -------
        source
            int32 source target indices.

        target
            int32 candidate target indices.

        score
            float32 fingerprint scores.
        """
        neighbors, scores = self.topk(
            k=k,
            min_score=min_score,
            return_scores=True,
        )

        n, kk = neighbors.shape

        source = np.repeat(
            np.arange(n, dtype=np.int32),
            kk,
        )

        target = neighbors.ravel()
        score = scores.ravel()

        keep = target >= 0

        return (
            source[keep],
            target[keep],
            score[keep],
        )


# ---------------------------------------------------------------------------
# Per-target voxelization
# ---------------------------------------------------------------------------

def _voxelize_target(
    coords: np.ndarray,
    anchor_mask: Optional[np.ndarray],
    cell_size: float,
    query_com: Optional[np.ndarray],
    radial_r0: Optional[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Voxelize one target.

    Returns
    -------
    hashes
        Unique voxel hashes.

    dwell
        Number of residues/CA atoms in each voxel.

    radial
        Radial damping factor for each voxel.
    """
    coords = np.asarray(coords, dtype=np.float32)

    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(
            "Coordinates must have shape (N, 3)"
        )

    if anchor_mask is not None:
        anchor_mask = np.asarray(anchor_mask, dtype=bool)

        if anchor_mask.shape != (coords.shape[0],):
            raise ValueError(
                "anchor_mask must have one boolean per coordinate"
            )

        coords = coords[~anchor_mask]

    if coords.shape[0] == 0:
        return (
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.uint16),
            np.empty(0, dtype=np.float32),
        )

    # Same voxelization convention as voxel.py.
    cell_xyz = np.floor(
        coords / cell_size
    ).astype(np.int64)

    hashes = _hash_cells(cell_xyz)

    # np.unique gives one entry per occupied voxel. Since each target is
    # processed independently, every returned voxel is automatically a
    # binary document for that target.
    unique_hashes, first, dwell = np.unique(
        hashes,
        return_index=True,
        return_counts=True,
    )

    unique_hashes = unique_hashes.astype(
        np.int64,
        copy=False,
    )

    # Counts are <= number of residues. uint16 is sufficient for ordinary
    # protein domains and halves memory relative to int64.
    if dwell.max(initial=0) <= np.iinfo(np.uint16).max:
        dwell = dwell.astype(np.uint16, copy=False)
    else:
        dwell = dwell.astype(np.uint32, copy=False)

    if radial_r0 is None:
        radial = np.ones(
            unique_hashes.size,
            dtype=np.float32,
        )
    else:
        if query_com is None:
            raise ValueError(
                "query_com is required when radial_r0 is not None"
            )

        # The first occurrence of each unique hash provides its integer
        # voxel coordinate.
        unique_cell_xyz = cell_xyz[first]

        voxel_centers = (
            unique_cell_xyz.astype(np.float32) + 0.5
        ) * cell_size

        delta = voxel_centers - query_com[None, :]

        radius = np.sqrt(
            np.sum(delta * delta, axis=1)
        )

        radial = (
            1.0
            / (1.0 + (radius / radial_r0) ** 2)
        ).astype(np.float32)

    return unique_hashes, dwell, radial


# ---------------------------------------------------------------------------
# Index construction
# ---------------------------------------------------------------------------

def build_structural_fingerprint_index(
    target_labels,
    all_coords,
    all_anchor_masks,
    *,
    query_com=None,
    cell_size=4.5,
    radial_r0=20.0,
    idf_alpha=1.0,
    idf_cap=2.0,
    kernel="cross7",
    kernel_neighbor_weight=0.25,
    kernel_sigma=None,
    kernel_radius_voxels=1,
    weight_epsilon=1e-7,
    keep_diagnostics=False,
) -> StructuralFingerprintIndex:
    """
    Build the scalable structural fingerprint index.

    Parameters
    ----------
    coords_dict
        Mapping:

            target_id -> (N_i, 3) float coordinates

        Coordinates MUST already be transformed into the common DALI/query
        frame.

    anchor_masks
        Optional mapping:

            target_id -> boolean array

        True entries are excluded from the fingerprint. This should normally
        contain the DALI-aligned query/anchor domain.

        If None, all coordinates are treated as non-anchor coordinates.

    query_com
        Query-domain center of mass in the same common coordinate frame.

    cell_size
        Voxel size in Å. 4.5 Å matches voxel.py.

    radial_r0
        Radial damping length in Å.

            D(r) = 1 / (1 + (r/r0)^2)

        Set to None to disable radial damping.

    idf_alpha
        Additive smoothing in the IDF formula.

    idf_cap
        Maximum IDF value.

    kernel
        "cross7" is recommended for ~10K targets.
        "gaussian" produces a 27-point kernel by default.

    kernel_neighbor_weight
        Neighbor weight for cross7.

    kernel_sigma
        Gaussian sigma in Å.

    kernel_radius_voxels
        Radius of Gaussian kernel in voxels.

    weight_epsilon
        Features below this weighted threshold are omitted from the
        inverted index.

    keep_diagnostics
        If True, retain the global voxel vocabulary, occupancy, and IDF.
        These arrays can be large at 10M features, so False is recommended
        unless diagnostics are needed.

    Notes
    -----
    The implementation intentionally does NOT construct:

        target x voxel

    dense matrices, nor does it construct globally kernel-expanded sparse
    matrices.

    At 10M occupied voxels, the raw representation is roughly:

        hashes       80 MB
        dwell        20-40 MB
        weights      40 MB

    plus the inverted index and temporary sorting arrays.
    """

    if not target_labels:
        raise ValueError("target_labels is empty")

    if len(all_coords) != len(target_labels):
        raise ValueError(
            "target_labels and all_coords must have the same length"
        )

    if len(all_anchor_masks) != len(target_labels):
        raise ValueError(
            "target_labels and all_anchor_masks must have the same length"
        )

    target_ids = list(target_labels)
    n_targets = len(target_ids)

    if query_com is not None:
        query_com_arr = np.asarray(
            query_com,
            dtype=np.float32,
        )

        if query_com_arr.shape != (3,):
            raise ValueError(
                "query_com must have shape (3,)"
            )
    else:
        query_com_arr = None

    if radial_r0 is not None and radial_r0 <= 0:
        raise ValueError(
            "radial_r0 must be > 0 or None"
        )

    # ---------------------------------------------------------------
    # 1. Per-target voxelization
    # ---------------------------------------------------------------

    hash_chunks: List[np.ndarray] = []
    dwell_chunks: List[np.ndarray] = []
    radial_chunks: List[np.ndarray] = []

    offsets = np.zeros(
        n_targets + 1,
        dtype=np.int64,
    )

    for target_idx in range(n_targets):
        coords = all_coords[target_idx]
        anchor_mask = all_anchor_masks[target_idx]

        hashes, dwell, radial = _voxelize_target(
            coords=coords,
            anchor_mask=anchor_mask,
            cell_size=cell_size,
            query_com=query_com_arr,
            radial_r0=radial_r0,
        )

        hash_chunks.append(hashes)
        dwell_chunks.append(dwell)
        radial_chunks.append(radial)

        offsets[target_idx + 1] = (
            offsets[target_idx] + hashes.size
        )

    raw_hashes = (
        np.concatenate(hash_chunks)
        if hash_chunks
        else np.empty(0, dtype=np.int64)
    )

    raw_dwell = (
        np.concatenate(dwell_chunks)
        if dwell_chunks
        else np.empty(0, dtype=np.uint16)
    )

    raw_radial = (
        np.concatenate(radial_chunks)
        if radial_chunks
        else np.empty(0, dtype=np.float32)
    )

    # Free the many small arrays as soon as possible.
    del hash_chunks
    del dwell_chunks
    del radial_chunks

    n_features = raw_hashes.size

    if n_features == 0:
        raise ValueError(
            "No non-anchor voxels were found"
        )

    # ---------------------------------------------------------------
    # 2. Global inverted voxel ordering
    #
    # Every target has at most one entry for a given voxel, so after
    # per-target np.unique(), the number of entries in a voxel group is
    # exactly its document frequency n_v.
    # ---------------------------------------------------------------

    order = np.argsort(
        raw_hashes,
        kind="mergesort",
    )

    sorted_hashes = raw_hashes[order]

    voxel_start = np.empty(
        n_features,
        dtype=bool,
    )

    voxel_start[0] = True
    voxel_start[1:] = (
        sorted_hashes[1:]
        != sorted_hashes[:-1]
    )

    voxel_starts = np.flatnonzero(voxel_start)

    voxel_ends = np.empty_like(voxel_starts)
    voxel_ends[:-1] = voxel_starts[1:]
    voxel_ends[-1] = n_features

    occupancy = (
        voxel_ends - voxel_starts
    ).astype(np.int32)

    vocabulary_hashes = sorted_hashes[
        voxel_starts
    ].copy()

    # ---------------------------------------------------------------
    # 3. IDF
    # ---------------------------------------------------------------

    vocabulary_idf = capped_idf(
        occupancy,
        n_targets=n_targets,
        alpha=idf_alpha,
        cap=idf_cap,
    )

    # Expand one IDF value per raw voxel entry in sorted order.
    sorted_idf = np.repeat(
        vocabulary_idf,
        occupancy,
    ).astype(np.float32, copy=False)

    # Return IDF values to target-grouped raw order.
    raw_idf = np.empty(
        n_features,
        dtype=np.float32,
    )

    raw_idf[order] = sorted_idf

    # ---------------------------------------------------------------
    # 4. Final per-voxel fingerprint weight
    #
    #      dwell * IDF * radial damping
    # ---------------------------------------------------------------

    raw_weights = (
        raw_dwell.astype(np.float32)
        * raw_idf
        * raw_radial
    ).astype(np.float32)

    # Raw feature norm used for inexpensive normalization during candidate
    # retrieval.
    norms = np.zeros(
        n_targets,
        dtype=np.float32,
    )

    for target_idx in range(n_targets):
        start = int(offsets[target_idx])
        end = int(offsets[target_idx + 1])

        w = raw_weights[start:end]

        norms[target_idx] = np.sqrt(
            np.sum(w * w, dtype=np.float64)
        )

    # ---------------------------------------------------------------
    # 5. Build sparse inverted index
    #
    # We retain only features with non-negligible weight. This naturally
    # removes voxels whose IDF is exactly zero, i.e. voxels occupied by every
    # target.
    # ---------------------------------------------------------------

    sorted_weights = raw_weights[order]

    # Target index for each sorted entry. Since raw arrays are target-grouped,
    # generate this without storing a second target array during voxelization.
    if n_targets <= np.iinfo(np.uint16).max:
        target_dtype = np.uint16
    else:
        target_dtype = np.uint32

    raw_targets = np.repeat(
        np.arange(n_targets, dtype=target_dtype),
        np.diff(offsets),
    )

    sorted_targets = raw_targets[order]

    keep = sorted_weights > weight_epsilon

    index_hashes = sorted_hashes[keep].copy()
    index_targets = sorted_targets[keep].copy()
    index_weights = sorted_weights[keep].copy()

    # ---------------------------------------------------------------
    # 6. Kernel
    # ---------------------------------------------------------------

    kernel_offsets, kernel_weights = make_kernel(
        kernel,
        cell_size=cell_size,
        neighbor_weight=kernel_neighbor_weight,
        sigma=kernel_sigma,
        radius_voxels=kernel_radius_voxels,
    )

    kernel_delta_hashes = np.asarray(
        [
            _hash_delta(
                int(dx),
                int(dy),
                int(dz),
            )
            for dx, dy, dz in kernel_offsets
        ],
        dtype=np.int64,
    )

    # ---------------------------------------------------------------
    # 7. Optional diagnostics
    # ---------------------------------------------------------------

    if keep_diagnostics:
        diag_hashes = vocabulary_hashes
        diag_occupancy = occupancy
        diag_idf = vocabulary_idf
    else:
        diag_hashes = None
        diag_occupancy = None
        diag_idf = None

    return StructuralFingerprintIndex(
        target_ids=target_ids,
        target_offsets=offsets,
        raw_hashes=raw_hashes,
        raw_weights=raw_weights,
        index_hashes=index_hashes,
        index_targets=index_targets,
        index_weights=index_weights,
        norms=norms,
        kernel_offsets=kernel_offsets,
        kernel_weights=kernel_weights,
        kernel_delta_hashes=kernel_delta_hashes,
        cell_size=float(cell_size),
        radial_r0=radial_r0,
        query_com=query_com_arr,
        idf_alpha=float(idf_alpha),
        idf_cap=float(idf_cap),
        weight_epsilon=float(weight_epsilon),
        vocabulary_hashes=diag_hashes,
        vocabulary_occupancy=diag_occupancy,
        vocabulary_idf=diag_idf,
    )


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

def build_query_mask(
    n_residues: int,
    domain_ranges: list[tuple[int, int]],
) -> np.ndarray:
    """
    Build a boolean mask from inclusive residue ranges.

    Example
    -------
    domain_ranges (1-based) = [(38, 166)]

    marks residues 38..166 inclusive as True.
    """
    mask = np.zeros(n_residues, dtype=bool)

    for start, end in domain_ranges:
        if start > end:
            raise ValueError(
                f"Invalid domain range: ({start}, {end})"
            )

        if start < 1 or end > n_residues:
            raise ValueError(
                f"Domain range ({start}, {end}) is outside "
                f"1..{n_residues}"
            )

        mask[start-1:end] = True

    return mask


def compute_query_com(
    query_coords: np.ndarray,
    query_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Compute a geometric center of mass from query coordinates.

    Parameters
    ----------
    query_coords
        Query coordinates in the same coordinate frame as the targets.

    query_mask
        Optional boolean mask selecting the query/anchor domain.
    """
    query_coords = np.asarray(
        query_coords,
        dtype=np.float32,
    )

    if query_coords.ndim != 2 or query_coords.shape[1] != 3:
        raise ValueError(
            "query_coords must have shape (N, 3)"
        )

    if query_mask is not None:
        query_mask = np.asarray(
            query_mask,
            dtype=bool,
        )

        if query_mask.shape != (query_coords.shape[0],):
            raise ValueError(
                "query_mask must have one boolean per coordinate"
            )

        query_coords = query_coords[query_mask]

    if query_coords.shape[0] == 0:
        raise ValueError(
            "No coordinates available for query COM"
        )

    return query_coords.mean(
        axis=0,
        dtype=np.float32,
    )


def build_knn_graph(
    fingerprint_index: StructuralFingerprintIndex,
    *,
    k: int = 200,
    min_score: float = 0.0,
    symmetrize: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the sparse candidate graph.

    Parameters
    ----------
    k
        Number of directed candidates per target.

    symmetrize
        If True, retain the union of directed k-NN edges.

    Returns
    -------
    source, target, score
        Sparse graph edge arrays.
    """
    source, target, score = fingerprint_index.candidate_pairs(
        k=k,
        min_score=min_score,
    )

    if not symmetrize:
        return source, target, score

    # Make an undirected union of the directed k-NN graph.
    lo = np.minimum(source, target)
    hi = np.maximum(source, target)

    pairs = np.column_stack((lo, hi))

    order = np.lexsort(
        (
            pairs[:, 1],
            pairs[:, 0],
        )
    )

    pairs = pairs[order]
    ordered_scores = score[order]

    if pairs.shape[0] == 0:
        return (
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.float32),
        )

    unique = np.ones(
        pairs.shape[0],
        dtype=bool,
    )

    unique[1:] = np.any(
        pairs[1:] != pairs[:-1],
        axis=1,
    )

    starts = np.flatnonzero(unique)

    # If both directions exist, retain the larger score.
    ends = np.empty_like(starts)
    ends[:-1] = starts[1:]
    ends[-1] = pairs.shape[0]

    n_edges = starts.size

    out_source = np.empty(
        n_edges,
        dtype=np.int32,
    )
    out_target = np.empty(
        n_edges,
        dtype=np.int32,
    )
    out_score = np.empty(
        n_edges,
        dtype=np.float32,
    )

    for i, (start, end) in enumerate(
        zip(starts, ends)
    ):
        block_scores = ordered_scores[start:end]
        best = start + int(
            np.argmax(block_scores)
        )

        out_source[i] = pairs[best, 0]
        out_target[i] = pairs[best, 1]
        out_score[i] = ordered_scores[best]

    return (
        out_source,
        out_target,
        out_score,
    )


# ---------------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # coords_dict:
    #
    #     {
    #         target_id: transformed_CA_coordinates,
    #         ...
    #     }
    #
    # anchor_masks:
    #
    #     {
    #         target_id: boolean_mask,
    #         ...
    #     }
    #
    # query_com:
    #
    #     np.array([x, y, z], dtype=np.float32)

    pass
