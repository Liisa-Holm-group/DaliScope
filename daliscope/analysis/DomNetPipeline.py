import daliscope
import daliscope.core.project
import daliscope.mechanics.delaunay
import daliscope.mechanics.saba
import daliscope.mechanics.segmentation
import daliscope.mechanics.graph
import daliscope.mechanics.partitioning
import daliscope.mechanics.dssp
import daliscope.optics.py3Dmol_viewer

from infomap import Infomap
import numpy as np
from collections import defaultdict

import math
from collections import defaultdict
from typing import Dict, Optional, Tuple, Union
import numpy as np


def contract_graph(
    edges: np.ndarray,
    cluster_of: Dict[int, int],
    gaussian_sigma: Optional[float] = None,
) -> np.ndarray:
    """
    ONE generic contraction primitive, reused for BOTH residue->segment
    and segment->final-contracted-graph steps.

    edges: 2D NumPy array of shape (N, 3) or structured array containing (u, v, w).
    cluster_of: node id -> cluster id. Missing keys or cluster id < 0
        EXCLUDE that node from the contracted graph entirely.
    gaussian_sigma: if given, apply exp(-w^2 / (2*sigma^2)) to each edge
        weight before summing. If None, sum raw weights as-is.

    Returns: NumPy array of shape (M, 3) where columns are (cluster_u, cluster_v, summed_weight),
             or structured NumPy array if input edges was structured.
    """
    if len(edges) == 0:
        if edges.dtype.names:
            return np.empty(0, dtype=[("src", int), ("dst", int), ("weight", float)])
        return np.empty((0, 3), dtype=float)

    # Check if edges is a structured array or a 2D ndarray
    is_structured = False #edges.dtype.names is not None
    if is_structured:
        src = edges["src"]
        dst = edges["dst"]
        weights = edges["weight"]
    else:
        src = edges[:, 0]
        dst = edges[:, 1]
        weights = edges[:, 2]

    # Pre-apply Gaussian scaling vectorized if sigma is provided
    if gaussian_sigma is not None:
        weights = np.exp(-(weights ** 2) / (2.0 * gaussian_sigma ** 2))

    agg: Dict[Tuple[int, int], float] = defaultdict(float)

    # Aggregate edge weights across cluster boundaries
    for u, v, w in zip(src, dst, weights):
        cu = cluster_of.get(int(u), -1)
        cv = cluster_of.get(int(v), -1)

        # Exclude unmapped nodes (-1) and self-loops (cu == cv)
        if cu < 0 or cv < 0 or cu == cv:
            continue

        key = (cu, cv) if cu < cv else (cv, cu)
        agg[key] += float(w)

    if not agg:
        if is_structured:
            return np.empty(0, dtype=[("src", int), ("dst", int), ("weight", float)])
        return np.empty((0, 3), dtype=float)

    # Reconstruct into NumPy array
    result_data = [(cu, cv, w) for (cu, cv), w in agg.items()]

    if is_structured:
        dtype = [("src", int), ("dst", int), ("weight", float)]
        return np.array(result_data, dtype=dtype)

    return np.array(result_data, dtype=float)

# ==============================================================================
# DomNet: Protein Domain Inference Pipeline
# ==============================================================================

class DomNetPipeline:
    """
    Reusable DomNet protein domain inference pipeline.

    The pipeline is instantiated once with algorithmic configuration and can
    then be run repeatedly on different proteins or coordinate sources.

    Input coordinates are intentionally NOT managed by this class. The
    caller supplies C-alpha coordinates and the corresponding sequence to
    ``run()``. Coordinates may come from ``project``, a PDB/mmCIF parser,
    or any other compatible source.

    Pipeline structure
    ------------------

        DomNetPipeline
        │
        ├── configuration
        │
        ├── run(xyz, sequence)
        │
        ├── secondary structure
        │   ├── ss_str
        │   └── segments
        │
        ├── structural organization
        │   └── sheets
        │
        ├── community detection
        │   ├── hierarchy
        │   ├── core_communities
        │   └── opt_core_communities
        │
        └── domain inference
            ├── domains
            ├── residue_assignments
            ├── domain_strings
            ├── repaired_strings
            └── score

    Typical usage
    --------------

        pipeline = DomNetPipeline(
            im=im3,
            sigma=24.0,
            core_dist_threshold=10.0,
            opt_dist_threshold=12.0,
            alpha=0.0,
            min_cluster_size=2,
        )

        xyz, sequence, dssp = project.retrieve_coords_sequ_dssp(
            protein_id=152
        )

        pipeline.run(xyz, sequence)

        #print(pipeline.domain_strings)

        daliscope.optics.py3Dmol_viewer.visualize_communities_ca_sticks(
            pipeline.xyz,
            pipeline.residue_assignments,
            pipeline.segments,
        )

    The same pipeline instance can then be reused:

        for protein_id in protein_ids:
            xyz, sequence, dssp = project.retrieve_coords_sequ_dssp(
                protein_id=protein_id
            )

            pipeline.run(xyz, sequence)

            #print(pipeline.domain_strings)

    Public input state
    ------------------
    xyz : ndarray
        C-alpha coordinates for the protein from the most recent run.
    sequence : str
        Protein sequence from the most recent run.
    n_residues : int
        Number of residues in the current protein.

    Public structural results
    -------------------------
    ss_str : str
        Secondary-structure assignment string.
    segments : list
        SSE segment tuples.
    sheets : list
        Final beta-sheet / strand structural clusters.

    Public community results
    ------------------------
    hierarchy : object
        Infomap hierarchy.
    core_communities : object
        Initial CORE communities.
    opt_core_communities : object
        CORE communities after hierarchy optimization.

    Public domain results
    --------------------
    domains : object
        Final domain representation.
    residue_assignments : object
        Final residue-to-domain assignments.
    domain_strings : object
        Domain assignments represented as residue ranges.
    repaired_strings : object
        Domain ranges after 3D repair.
    score : float
        Final domain agglomeration score.

    Configuration
    -------------
    sigma : float
        Gaussian weighting parameter for the cluster graph.
    core_dist_threshold : float
        Distance threshold for initial non-core residue assignment.
    opt_dist_threshold : float
        Distance threshold after hierarchy optimization.
    alpha : float
        DP optimization regularization parameter.
    min_cluster_size : int
        Minimum cluster size for DP optimization.

    Notes
    -----
    Intermediate computational state is deliberately kept private.

    The pipeline does not retrieve protein data and has no dependency on the
    ``project`` data-management layer beyond the example usage above. This
    keeps coordinate acquisition separate from structural domain inference.
    """

    def __init__(
        self,
        im_params="",
        sigma=4.0,
        core_dist_threshold=10.0,
        opt_dist_threshold=12.0,
        alpha=0.0,
        min_cluster_size=2,
    ):
        # ------------------------------------------------------------------
        # Algorithm configuration
        # ------------------------------------------------------------------

        self.im = Infomap(im_params)
        # Infomap()
        # Infomap(" --num-trials 10 --two-level --silent ")
        # increased walk time counter-acts small sigma
        # Infomap(" --num-trials 10 --two-level --silent --markov-time 1.5 ")

        self.sigma = sigma
        self.core_dist_threshold = core_dist_threshold
        self.opt_dist_threshold = opt_dist_threshold
        self.alpha = alpha
        self.min_cluster_size = min_cluster_size

        # ------------------------------------------------------------------
        # Current input
        #
        # These are populated by run() and replaced on every subsequent run.
        # ------------------------------------------------------------------

        self.xyz = None
        self.sequence = None
        self.n_residues = None

        # ------------------------------------------------------------------
        # Public scientific results
        # ------------------------------------------------------------------

        self.ss_str = None
        self.segments = None
        self.sheets = None

        self.hierarchy = None
        self.core_communities = None
        self.opt_core_communities = None

        self.domains = None
        self.residue_assignments = None

        self.domain_strings = None
        self.repaired_strings = None

        self.score = None

        # ------------------------------------------------------------------
        # Private implementation state
        # ------------------------------------------------------------------

        self._delaunay_edges = None
        self._delaunay_distances = None
        self._distance_adj = None

        self._strand_ranges = None
        self._sheet_partners = None
        self._saba_pairs = None

        self._sheet_claimed = None
        self._full_clusters = None
        self._full_graph_triplets = None
        self._full_segment_adj = None

        self._repeat_candidates_raw = None
        self._repeat_candidates = None
        self._repeat_clusters = None

        self._final_clusters = None
        self._graph_triplets = None

        self._final_residue_assignments = None

        self._partition = None
        self._total_modularity = None
        self._choices = None
        self._node_to_cluster = None

        self._opt_residue_assignments = None

    # ==========================================================================
    # Public API
    # ==========================================================================

    def run(self, xyz, sequence):
        """
        Run DomNet on a protein.

        Parameters
        ----------
        xyz : ndarray, shape (N, 3)
            C-alpha coordinates.
        sequence : str
            Protein amino-acid sequence.

        Returns
        -------
        DomNetPipeline_build_structural_clusters
            This pipeline instance.

        """

        self._reset_run_state()

        self.xyz = xyz
        self.sequence = sequence
        self.n_residues = len(sequence)

        self._build_geometry()
        self._assign_secondary_structure()
        print('1', self.ss_str)
        print('1 _sheet_partners', self._sheet_partners)
        self._build_structural_clusters()
        print('2 _graph_triplets')
        for u,v,w in self._graph_triplets: print(u,v,w)
        self._run_infomap()
        print('3 hierarchy:',self.hierarchy)
        self._optimize_partition()
        print('4 optimized partition',self._opt_residue_assignments)
        self._agglomerate_domains()
        print('5 after agglomeration',self.domains)
        self._repair_domains()
        print('6 after repairs',self.domain_strings)

        return self

    # ==========================================================================
    # State management
    # ==========================================================================

    def _reset_run_state(self):
        """
        Clear all protein-specific state before a new run.

        Configuration is deliberately preserved.
        """

        self.xyz = None
        self.sequence = None
        self.n_residues = None

        # Public results
        self.ss_str = None
        self.segments = None
        self.sheets = None

        self.hierarchy = None
        self.core_communities = None
        self.opt_core_communities = None

        self.domains = None
        self.residue_assignments = None

        self.domain_strings = None
        self.repaired_strings = None

        self.score = None

        # Private intermediates
        self._delaunay_edges = None
        self._delaunay_distances = None
        self._distance_adj = None

        self._strand_ranges = None
        self._sheet_partners = None
        self._saba_pairs = None

        self._sheet_claimed = None
        self._full_clusters = None
        self._full_graph_triplets = None
        self._full_segment_adj = None

        self._repeat_candidates_raw = None
        self._repeat_candidates = None
        self._repeat_clusters = None

        self._final_clusters = None
        self._graph_triplets = None

        self._final_residue_assignments = None

        self._partition = None
        self._total_modularity = None
        self._choices = None
        self._node_to_cluster = None

        self._opt_residue_assignments = None

    # ==========================================================================
    # Read-only diagnostic properties
    # ==========================================================================

    @property
    def distance_adj(self):
        """Residue-distance adjacency matrix for the current protein."""
        return self._distance_adj

    @property
    def graph_triplets(self):
        """Final graph as ``(u, v, weight)`` triplets."""
        return self._graph_triplets

    @property
    def partition(self):
        """Optimal partition returned by DP optimization."""
        return self._partition

    @property
    def total_modularity(self):
        """Total modularity associated with the selected partition."""
        return self._total_modularity

    @property
    def choices(self):
        """Selected hierarchy choices from DP optimization."""
        return self._choices

    @property
    def node_to_cluster(self):
        """Mapping from hierarchy nodes to selected clusters."""
        return self._node_to_cluster

    @property
    def final_residue_assignments(self):
        """
        Residue assignments before hierarchy optimization.

        Primarily useful for diagnostics. The preferred final result is
        ``residue_assignments``.
        """
        return self._final_residue_assignments

    @property
    def opt_residue_assignments(self):
        """
        Residue assignments after optimized CORE assignment and before
        final domain repair.
        """
        return self._opt_residue_assignments

    # ==========================================================================
    # Pipeline stages
    # ==========================================================================

    def _build_geometry(self):
        """Construct the residue-distance representation."""

        (
            self._delaunay_edges,
            self._delaunay_distances,
            self._distance_adj, # sparse square matrix
        ) = daliscope.mechanics.delaunay.delaunay_edges(
            self.xyz,
            return_sparse=True,
            dist_cutoff=20.0,
            normalize=False,
        )

    # --------------------------------------------------------------------------

    def _assign_secondary_structure(self):
        """
        Assign secondary structure and split chimeric strands using SABA
        hydrogen-bond topology.
        """

        (
            self.ss_str,           # np.array
            self._strand_ranges,   # np.array
            self._sheet_partners,  # np.array
        ) = daliscope.mechanics.saba.assign_secondary_structure_saba(
            ca=self.xyz,
            sequence=self.sequence,
            helix_min_length=5,
            helix310_min_length=3,
            strand_min_length=2,
            do_refine_antiparallel_termini=False,
        )
        #print(self.ss_str, self._strand_ranges, self._sheet_partners)

        self.segments = daliscope.mechanics.saba.ss_to_segments(
            self.ss_str
        )
        #print(self.ss_str)

        # Re-compute SABA pairs for H-bond strand-cutting logic.
        self._saba_pairs = daliscope.mechanics.saba.find_beta_pairs_numba(
            self.xyz
        )
        #print(self._saba_pairs)

        (
            self.segments,
            self._strand_ranges,
            self._sheet_partners,
            self.ss_str,
        ) = (
            daliscope.mechanics.segmentation
            .split_strands_by_saba_hbond_transitions(
                segments=self.segments,
                strand_ranges=self._strand_ranges,
                sheet_partners=self._sheet_partners,
                pairs=self._saba_pairs,
                ss_arr=self.ss_str,
                min_strand_len=2,
            )
        )

    # --------------------------------------------------------------------------

    def _build_structural_clusters(self):
        """Identify structural clusters and construct the contracted graph for Infomap."""
        secondary_struct_tuple = (
            self.ss_str,
            self._strand_ranges,
            self._sheet_partners,
        )

        saba_strands = self._strand_ranges

        # 1. Detect beta pairs and split strands by sheet membership
        raw_sheet_pairs = daliscope.mechanics.saba.find_beta_pairs_numba(self.xyz)
        strands, old_to_new, cluster_sheets = (
            daliscope.mechanics.dssp.split_strands_by_sheet_membership(
                saba_strands, raw_sheet_pairs
            )
        )

        # 2. Extract sheet clusters: maps sheet_id -> np.ndarray of segment indices
        unique_clusters = np.unique(cluster_sheets[cluster_sheets >= 0])
        sheet_clusters = {
            int(cid): np.flatnonzero(cluster_sheets == cid)
            for cid in unique_clusters
        }

        # Track segments claimed by sheets
        self._sheet_claimed = daliscope.mechanics.graph.claimed_segment_ids(
            sheet_clusters
        )

        # 3. BUILD SHEET-CONTRACTED CLUSTERS (Sheets = Super-Nodes, Helices/Loops = Single Nodes)
        # Merging sheet_clusters with all un-claimed single segments
        self._full_clusters = daliscope.mechanics.graph.merge_cluster_lists(
            sheet_clusters,
            {},  # No repeat clusters yet
            n_segments=len(self.segments),
        )

        # 4. Sheet-contracted graph used for tandem-repeat detection
        self._full_graph_triplets = (
            daliscope.mechanics.graph.create_cluster_graph_gaussian(
                self.n_residues,
                self._distance_adj,
                self.segments,
                self._full_clusters,  # Uses sheet-contracted super-nodes!
                sigma=self.sigma,
            )
        )

        self._full_segment_adj = (
            daliscope.mechanics.graph.graph_triplets_to_segment_adj(
                self._full_graph_triplets,
                n_segments=len(self._full_clusters), # dynamic node count
            )
        )

        # ------------------------------------------------------------------
        # Tandem Repeats Detection
        # ------------------------------------------------------------------

        self._repeat_candidates_raw = (
            daliscope.mechanics.graph.detect_tandem_repeats(self._full_segment_adj)
        )
        self._repeat_candidates = (
            daliscope.mechanics.graph.filter_repeat_candidates(
                self._repeat_candidates_raw,
                self._sheet_claimed,
            )
        )
        self._repeat_clusters = (
            daliscope.mechanics.graph.repeat_candidates_to_clusters(
                self._repeat_candidates,
                self.segments,
            )
        )

        # Defensive check
        daliscope.mechanics.graph.assert_clusters_disjoint(
            sheet_clusters,
            self._repeat_clusters,
        )

        # ------------------------------------------------------------------
        # Final Contracted Structural Clusters (Sheets + Repeats + Helices)
        # ------------------------------------------------------------------

        self._final_clusters = daliscope.mechanics.graph.merge_cluster_lists(
            sheet_clusters,
            self._repeat_clusters,
            n_segments=len(self.segments),
        )

        _raw_triplets = daliscope.mechanics.graph.create_cluster_graph_gaussian(
            self.n_residues,
            self._distance_adj,
            self.segments,
            self._final_clusters,
            sigma=self.sigma,
        )

        # Filter out self-edges
        dtype = [("src", int), ("dst", int), ("weight", float)]
        triplets_arr = np.array(_raw_triplets, dtype=dtype)
        self._graph_triplets = triplets_arr[triplets_arr["src"] != triplets_arr["dst"]]

        self.sheets = self._final_clusters

    def buggy_build_structural_clusters(self):
        """Identify structural clusters and construct the final graph used by Infomap."""
        secondary_struct_tuple = (
            self.ss_str,
            self._strand_ranges,
            self._sheet_partners,
        )

        saba_strands = self._strand_ranges
        #print("#saba_strands", saba_strands)

        raw_sheet_pairs = daliscope.mechanics.saba.find_beta_pairs_numba(self.xyz)
        #print("#raw_sheet_pairs", raw_sheet_pairs)

        strands, old_to_new, cluster_sheets = (
            daliscope.mechanics.dssp.split_strands_by_sheet_membership(
                saba_strands, raw_sheet_pairs
            )
        )
        print("#cluster_sheets", cluster_sheets) # nres array

        # Map each cluster_id to a 1D NumPy array of segment indices (skipping -1)
        unique_clusters = np.unique(cluster_sheets[cluster_sheets >= 0])
        sheet_clusters = {
            int(cid): np.flatnonzero(cluster_sheets == cid)
            for cid in unique_clusters
        }
        print("#sheet_clusters", sheet_clusters) # Dict: { cluster-id: residue-array }

        self._sheet_claimed = daliscope.mechanics.graph.claimed_segment_ids(
            sheet_clusters
        )

        # Build complete ordered cluster representation
        #print("# segments", self.segments)
        self._full_clusters = (
            daliscope.mechanics.graph.build_full_ordered_clusters(self.segments)
        )

        # Graph used for tandem-repeat detection
        #print("# _full_clusters", self._full_clusters)
        self._full_graph_triplets = (
            daliscope.mechanics.graph.create_cluster_graph_gaussian(
                self.n_residues,
                self._distance_adj,
                self.segments,
                self._full_clusters,
                sigma=self.sigma,
            )
        )

        #print("# _full_graph_triplets", self._full_graph_triplets)
        self._full_segment_adj = (
            daliscope.mechanics.graph.graph_triplets_to_segment_adj(
                self._full_graph_triplets,
                n_segments=len(self.segments),
            )
        )

        # ------------------------------------------------------------------
        # Tandem repeats
        # ------------------------------------------------------------------

        self._repeat_candidates_raw = (
            daliscope.mechanics.graph.detect_tandem_repeats(self._full_segment_adj)
        )
        #print("# _repeat_candidates_raw", self._repeat_candidates_raw)
        self._repeat_candidates = (
            daliscope.mechanics.graph.filter_repeat_candidates(
                self._repeat_candidates_raw,
                self._sheet_claimed,
            )
        )
        #print("# repeat_candidates", self._repeat_candidates)
        self._repeat_clusters = (
            daliscope.mechanics.graph.repeat_candidates_to_clusters(
                self._repeat_candidates,
                self.segments,
            )
        )
        #print("# _repeat_clusters", self._repeat_clusters)
        # Defensive check: sheet and repeat clusters must remain disjoint
        daliscope.mechanics.graph.assert_clusters_disjoint(
            sheet_clusters,
            self._repeat_clusters,
        )

        # ------------------------------------------------------------------
        # Final structural clusters
        # ------------------------------------------------------------------

        self._final_clusters = daliscope.mechanics.graph.merge_cluster_lists(
            sheet_clusters,
            self._repeat_clusters,
            n_segments=len(self.segments),
        )

        _graph_triplets = daliscope.mechanics.graph.create_cluster_graph_gaussian(
                self.n_residues,
                self._distance_adj,
                self.segments,
                self._final_clusters,
                sigma=self.sigma,
        )
        # no self edges
        # Create a boolean mask where src != dst
        dtype = [("src", int), ("dst", int), ("weight", float)]
        triplets_arr = np.array(_graph_triplets, dtype=dtype)
        # Filter out self-edges
        mask = triplets_arr["src"] != triplets_arr["dst"]
        self._graph_triplets = triplets_arr[mask]
        #print("# _graph_triplets", self._graph_triplets)
        self.sheets = self._final_clusters

    def old_build_structural_clusters(self):
        """
        Identify structural clusters and construct the final graph used
        by Infomap.
        """

        secondary_struct_tuple = (
            self.ss_str,
            self._strand_ranges,
            self._sheet_partners,
        )

        saba_strands = self._strand_ranges
        #print('#saba_strands', saba_strands)
        raw_sheet_pairs = daliscope.mechanics.saba.find_beta_pairs_numba(self.xyz)
        #print('#raw_sheet_pairs', raw_sheet_pairs)
        strands, old_to_new, cluster_sheets = daliscope.mechanics.dssp.split_strands_by_sheet_membership(saba_strands, raw_sheet_pairs)
        sheet_clusters = {} # {cluster_id: [seg_id_1, seg_id_2, ...]}
        #print("#cluster_sheets", cluster_sheets)

        # Extract valid unique cluster IDs (excluding unassigned -1)
        unique_clusters = np.unique(cluster_sheets[cluster_sheets >= 0])

        # Map each cluster_id to a 1D NumPy array of segment indices
        sheet_clusters = {
            int(cid): np.flatnonzero(cluster_sheets == cid) for cid in unique_clusters
        }

#        for seg_id, cluster_id in cluster_sheets.items():
#            sheet_clusters.setdefault(cluster_id, []).append(seg_id)
        #print('#sheet_clusters', sheet_clusters)

        self._sheet_claimed = (
            daliscope.mechanics.graph.claimed_segment_ids(
                sheet_clusters
            )
        )

        # Build complete ordered cluster representation.
        #print("# segments",self.segments)
        self._full_clusters = (
            daliscope.mechanics.graph.build_full_ordered_clusters(
                self.segments
            )
        )

        # Graph used for tandem-repeat detection.
        self._full_graph_triplets = (
            daliscope.mechanics.graph.create_cluster_graph_gaussian(
                self.n_residues,
                self._distance_adj,
                self.segments,
                self._full_clusters,
                sigma=self.sigma,
            )
        )

        self._full_segment_adj = (
            daliscope.mechanics.graph.graph_triplets_to_segment_adj(
                self._full_graph_triplets,
                n_segments=len(self.segments),
            )
        )

        # ------------------------------------------------------------------
        # Tandem repeats
        # ------------------------------------------------------------------

        self._repeat_candidates_raw = (
            daliscope.mechanics.graph.detect_tandem_repeats(
                self._full_segment_adj
            )
        )

        self._repeat_candidates = (
            daliscope.mechanics.graph.filter_repeat_candidates(
                self._repeat_candidates_raw,
                self._sheet_claimed,
            )
        )

        self._repeat_clusters = (
            daliscope.mechanics.graph.repeat_candidates_to_clusters(
                self._repeat_candidates,
                self.segments,
            )
        )

        # Sheet and repeat clusters must remain disjoint.
        daliscope.mechanics.graph.assert_clusters_disjoint(
            sheet_clusters,
            self._repeat_clusters,
        )

        # ------------------------------------------------------------------
        # Final structural clusters
        # ------------------------------------------------------------------

        self._final_clusters = (
            daliscope.mechanics.graph.merge_cluster_lists(
                sheet_clusters,
                self._repeat_clusters,
                n_segments=len(self.segments),
            )
        )

        self._graph_triplets = (
            daliscope.mechanics.graph.create_cluster_graph_gaussian(
                self.n_residues,
                self._distance_adj,
                self.segments,
                self._final_clusters,
                sigma=self.sigma,
            )
        )

        self.sheets = self._final_clusters

    # --------------------------------------------------------------------------

    def _run_infomap(self):
        """Run Infomap and construct initial CORE assignments."""

        self.hierarchy = (
            daliscope.mechanics.partitioning.run_infomap_hierarchy(
                self.im,
                self._graph_triplets,
            )
        )
        #print("# hierarchy", self.hierarchy)
        self.core_communities = (
            daliscope.mechanics.partitioning
            .get_core_communities_from_infomap(
                self.im,
                self.segments,
                self.sheets,
            )
        )
        #print("# core_communities", self.core_communities)
        self._final_residue_assignments = (
            daliscope.mechanics.partitioning.assign_non_core_residues(
                xyz=self.xyz,
                core_communities=self.core_communities,
                dist_threshold=self.core_dist_threshold,
            )
        )
        #print("# _final_residue_assignments", self._final_residue_assignments)

    # --------------------------------------------------------------------------

    def _optimize_partition(self):
        """Optimize the Infomap hierarchy using dynamic programming."""

        (
            self._partition,
            self._total_modularity,
            self._choices,
        ) = daliscope.mechanics.partitioning.dp_optimal_partition(
            self.hierarchy,
            self._distance_adj,
            alpha=self.alpha,
            min_cluster_size=self.min_cluster_size,
        )

        self._node_to_cluster = (
            daliscope.mechanics.partitioning
            .map_nodes_to_selected_clusters(
                self.hierarchy,
                self._choices,
            )
        )

        self.opt_core_communities = (
            daliscope.mechanics.partitioning
            .get_core_communities_from_pruned_tree(
                self.hierarchy,
                self._choices,
                self.segments,
                self.sheets,
            )
        )

        self._opt_residue_assignments = (
            daliscope.mechanics.partitioning.assign_non_core_residues(
                xyz=self.xyz,
                core_communities=self.opt_core_communities,
                dist_threshold=self.opt_dist_threshold,
            )
        )

    # --------------------------------------------------------------------------

    def _agglomerate_domains(self):
        """
        Agglomerate optimized Infomap communities into candidate domains.

        NOTE
        ----
        ``patch_unstructured_residues()`` is currently known to have
        incorrect/incomplete behavior. Its existing integration is retained
        here until its implementation is corrected.
        """

        self.domains, self.score = (
            daliscope.mechanics.partitioning
            .agglomerate_infomap_clusters(
                self.xyz,
                self._opt_residue_assignments,
                self._distance_adj,
            )
        )

        # TODO:
        # Correct patch_unstructured_residues() separately. Its current
        # implementation does not appear to use its arguments meaningfully.
        self.domains = (
            daliscope.mechanics.partitioning
            .patch_unstructured_residues(
                self._final_residue_assignments,
                self.domains,
            )
        )

    # --------------------------------------------------------------------------

    def _repair_domains(self):
        """Convert domains to residue ranges and perform final 3D repair."""

        #print("domains", self.domains)
        self.domain_strings = (
            daliscope.mechanics.partitioning
            .residue_assignments_to_domain_ranges(
                self.domains
            )
        )

        #print("domain_strings",self.domain_strings)
        #print("sheets", self.sheets)
        (
            self.residue_assignments,
            self.repaired_strings,
        ) = (
            daliscope.mechanics.partitioning
            .repair_domain_assignments_3d(
                self.xyz,
                self.domain_strings,
                self.sheets,
            )
        )

def convert_domains(domain_strings):
    # Group domain boundaries by chain/protein index (first column)
    chains = defaultdict(list)
    for row in domain_strings:
        chain_id, start, end = row
        chains[chain_id].append(f"{start}-{end}")

    # Combine each chain's boundaries with '_' and join chains with ', '
    return ", ".join("_".join(ranges) for chain_id, ranges in sorted(chains.items()))
