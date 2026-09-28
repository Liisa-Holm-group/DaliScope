# daliscope/core/project.py
"""
DaliScope Project — single central object replacing the awkward
DaliLakehouse + DaliScopeRegistry split.

Architecture
------------
The Project owns:
  - Pack (immutable bronze data: coords, sequences, segments, query)
  - Views (named materialized DataFrames)
  - Models (named analytical artifacts derived from views)
  - Plots (PlotSpec recipes — not figure objects)

Usage
-----
    project = Project.load_pack("path/to/pack.tar.gz")
    full_view = next(iter(project.views))
    project.views[full_view]                     # default full-length view
    project.add_view("trusted", trusted_df,      # register a view
                     parent=full_view,
                     function="threshold",
                     parameters={"z_score_min": 8})
    project.save("my_analysis/")
    project = Project.open("my_analysis/")       # restore across sessions
"""

import io
import json
import os
import pickle
import shutil
import tarfile
from dataclasses import dataclass, field, asdict, InitVar
from datetime import datetime
from pathlib import Path
from typing import Tuple, Any, Optional
from IPython.display import display

import numpy as np
import pandas as pd

from ..mechanics.metrics import best_pfam_per_alignment, best_pfam_clan_per_target, build_pileups
from ..mechanics.geometry import svd_superimpose, _parse_query_ca_coords
from ..mechanics.metrics import encode_dssp_optimized, nearest_neighbor_reordering_fast
from ..analysis.data_preview import show_heading
from ..optics.domain_cartoons import plot_domain_architecture_plotly, plot_domain_architecture_matplotlib, plot_domain_architecture
from ..optics.msa import plot_msa

from .. import __version__

class ProjectViz:
    """Visualization entry points for a Project. Owns display/renderer
    concerns and provides convenience plotting methods.
    """

    def __init__(self, project):
        self._project = project

    def architecture(
        self,
        view_or_df,
        **kwargs,
    ):
        """Dispatches to pure plot_domain_architecture function."""
        if isinstance(view_or_df, str):
            view_name = view_or_df
            df = self._project.views[view_name]
            title = f"{view_name} (n={len(df)})"
        elif isinstance(view_or_df, pd.DataFrame):
            df = view_or_df
            title = kwargs.pop("title", f"Custom DataFrame (n={len(df)})")
        else:
            raise TypeError(f"Expected str or pandas.DataFrame, got {type(view_or_df).__name__}")

        return plot_domain_architecture(df=df, title=title, **kwargs)

    def old_architecture(
        self,
        view_name,
        max_targets=10,
        pfam_col="pfam_domains",
        renderer="plotly",
        show_legend=True,
        group_by_architecture=True,
        representative_sort="z_score",
        sort_col=None,
        ascending=False,
        **kwargs,
    ):
        show_heading(view_name)
        df = self._project.views[view_name]

        if renderer == "plotly":
            fig = plot_domain_architecture_plotly(
                df=df,
                pfam_col=pfam_col,
                max_targets=max_targets,
                group_by_architecture=group_by_architecture,
                representative_sort=representative_sort,
                sort_col=sort_col,
                ascending=ascending,
                show_legend=show_legend,
                **kwargs,
            )
            display(fig)

        elif renderer == "matplotlib":
            fig = plot_domain_architecture_matplotlib(
                df=df,
                pfam_col=pfam_col,
                max_targets=max_targets,
                group_by_architecture=group_by_architecture,
                representative_sort=representative_sort,
                sort_col=sort_col,
                ascending=ascending,
                **kwargs,
            )
            # display(fig)

        else:
            raise ValueError(
                f"Unknown renderer: {renderer!r}"
            )

        #return fig

    def plot_msa(self, view_name, left, right, plot_type="logo", data_col="sequ_pileup"):
        df = self._project.views[view_name]
        show_heading( f"{view_name} ( n = {len(df)} ) " )
        return plot_msa(df, left, right, plot_type=plot_type, data_col=data_col)

# ================================================================== #
# Provenance — attached to every persistent object                   #
# ================================================================== #

@dataclass
class Provenance:
    created:           str         = field(default_factory=lambda: datetime.now().isoformat())
    function:          str         = ""
    parameters:        dict        = field(default_factory=dict)
    parent_objects:    list        = field(default_factory=list)
    daliscope_version: str         = __version__

    # InitVar means this argument is passed to __init__ but not saved as an attribute
    df: InitVar[pd.DataFrame] = None

    # The actual integer attribute stored on the dataclass
    N: int = 0

    def __post_init__(self, df):
        if df is not None:
            self.N = len(df)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Provenance":
        return cls(**d)


# ================================================================== #
# PlotSpec — recipe for a figure, not the figure itself              #
# ================================================================== #

@dataclass
class PlotSpec:
    """
    Lightweight recipe describing how to reproduce a plot.
    The figure object is never stored — always reproducible from
    renderer + source + params.
    """
    name:        str
    renderer:    str          # fully qualified function name e.g. "daliscope.optics.plots.scatter_plot"
    source:      str          # view name
    params:      dict         = field(default_factory=dict)
    caption:     str          = ""
    notes:       str          = ""      # user's scientific observations
    provenance:  Provenance   = field(default_factory=Provenance)
    png_path:    Optional[str] = None   # set after export_figure()

    def to_dict(self) -> dict:
        d = asdict(self)
        d["provenance"] = self.provenance.to_dict()
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "PlotSpec":
        d = d.copy()
        d["provenance"] = Provenance.from_dict(d.get("provenance", {}))
        return cls(**d)


# ================================================================== #
# Project                                                            #
# ================================================================== #

from typing import Optional, Any
import numpy as np
import pandas as pd


class Project:
    """
    Single central object for a DaliScope analysis session.

    Attributes
    ----------
    pack_path    : str — original tar.gz path
    views        : dict[str, pd.DataFrame]
    models       : dict[str, Any]
    plots        : dict[str, PlotSpec]
    provenance   : dict[str, Provenance] — one per view/model

    Pack data (set by _load_bronze / _build_silver)
    -------------------------------------------------
    coords           : (N_total, 3) float32
    target_seqs      : (N_total,)   int32
    target_dssp      : (N_total,)   int32
    query_pdb_str    : str
    query_ca_coords  : (L, 3) float32
    query_length     : int
    query_sequence   : str
    query_ascii      : (L,)   int32
    segments         : pd.DataFrame  — 0-based normalised
    master_metadata  : pd.DataFrame  — full silver join
    """

    # ---------------------------------------------------------------- #
    # Required files inside the tar                                    #
    # ---------------------------------------------------------------- #
    _REQUIRED_TSV = {
        "query_meta",
        "summary",
        "id_mapping",
        "colab_pack_index",
        "segments",
    }
    _REQUIRED_NPZ = {"coords", "seq", "dssp"}

    # ---------------------------------------------------------------- #
    # Construction                                                     #
    # ---------------------------------------------------------------- #

    def __init__(self):
        # pack data — set by _load_bronze / _build_silver
        self.pack_path: Optional[str] = None
        self.coords: Optional[np.ndarray] = None
        self.target_seqs: Optional[np.ndarray] = None
        self.target_dssp: Optional[np.ndarray] = None
        self.query_pdb_str: Optional[str] = None
        self.query_ca_coords: Optional[np.ndarray] = None
        self.query_length: Optional[int] = None
        self.query_sequence: Optional[str] = None
        self.query_ascii: Optional[np.ndarray] = None
        self.segments: Optional[pd.DataFrame] = None
        self.master_metadata: Optional[pd.DataFrame] = None

        # raw bronze tables (kept for diagnostics)
        self._query_meta: Optional[pd.DataFrame] = None
        self._summary: Optional[pd.DataFrame] = None
        self._id_mapping: Optional[pd.DataFrame] = None
        self._pack_index: Optional[pd.DataFrame] = None
        self._raw_segments: Optional[pd.DataFrame] = None
        self._pfam: Optional[pd.DataFrame] = None
        self._pfam_names: Optional[pd.DataFrame] = None
        self._ligands: Optional[pd.DataFrame] = None
        self._ligand_names: Optional[pd.DataFrame] = None

        # named assets
        self.views: dict[str, pd.DataFrame] = {}
        self.models: dict[str, Any] = {}
        self.plots: dict[str, PlotSpec] = {}
        self.provenance: dict[str, Provenance] = {}

        # internal flags
        self._bronze_loaded = False
        self._silver_built = False

        # for plotting views, e.g. project.viz.architecture(view_name)
        self.viz = ProjectViz(self)

    # ---------------------------------------------------------------- #
    # Query Metadata Properties & Convenience Methods                 #
    # ---------------------------------------------------------------- #

# ---------------------------------------------------------------- #
    # Query Metadata Properties & Convenience Methods                 #
    # ---------------------------------------------------------------- #

    @property
    def query_identifier(self) -> str:
        """Returns the query identifier (mapped from 'id')."""
        if self._query_meta is not None and not self._query_meta.empty:
            if "id" in self._query_meta.columns:
                return str(self._query_meta.iloc[0]["id"])
        return "N/A"

    @property
    def query_description(self) -> str:
        """Returns the query protein description (mapped from 'compnd')."""
        if self._query_meta is not None and not self._query_meta.empty:
            if "compnd" in self._query_meta.columns:
                return str(self._query_meta.iloc[0]["compnd"])
        return "N/A"

    @property
    def puu_domains(self) -> Optional[Any]:
        """Returns the PUU domain definitions (mapped from 'domain_string')."""
        if self._query_meta is not None and not self._query_meta.empty:
            if "domain_string" in self._query_meta.columns:
                return self._query_meta.iloc[0]["domain_string"]
        return None

    def show_query_meta(self) -> None:
        """Displays formatted metadata for the project query structure."""
        print(
            f"""
Identifier:  {self.query_identifier}
Length:      {self.query_length}
Description: {self.query_description}
PUU domains: {self.puu_domains}
        """.strip()
        )

    # ---------------------------------------------------------------- #
    # Entry Point                                                      #
    # ---------------------------------------------------------------- #

    @classmethod
    def load_pack(cls, tar_path: str, verbose: bool = True) -> "Project":
        """
        Main entry point. Load a colab-pack and build all silver data.
        Automatically registers a 'dom_0-<query_length>' view.

        Parameters
        ----------
        tar_path : path to .tar.gz colab-pack
        verbose  : print progress messages

        Returns
        -------
        Project with pack loaded and 'full' view registered.
        """
        project = cls()
        project.pack_path = str(Path(tar_path).resolve())
        project._load_bronze(tar_path)
        project._build_silver()
        project._register_full_view()
        return project

    # ---------------------------------------------------------------- #
    # Bronze — pure IO                                                  #
    # ---------------------------------------------------------------- #

    def _load_bronze(self, tar_path: str) -> None:
        raw_tsvs  = {}
        npz_data  = {}
        found_npz = False

        with tarfile.open(tar_path, "r:gz") as tar:
            for member in tar.getmembers():
                filename = member.name.split("/")[-1]
                if not filename:
                    continue
                fobj = tar.extractfile(member)
                if fobj is None:
                    continue

                if filename == "colab_pack.npz":
                    with np.load(io.BytesIO(fobj.read())) as data:
                        missing = self._REQUIRED_NPZ - set(data.files)
                        if missing:
                            raise KeyError(
                                f"colab_pack.npz missing arrays: {missing}"
                            )
                        npz_data["coords"] = data["coords"].astype(np.float32)
                        npz_data["seq"]    = data["seq"].astype(np.int32)
                        npz_data["dssp"]   = data["dssp"].astype(np.int32)
                    found_npz = True

                elif filename == "Query.pdb":
                    self.query_pdb_str   = fobj.read().decode("utf-8")
                    self.query_ca_coords = _parse_query_ca_coords(
                        self.query_pdb_str
                    )

                elif filename.endswith(".tsv"):
                    key = filename.replace(".tsv", "")
                    raw_tsvs[key] = pd.read_csv(fobj, sep="\t")

        if not found_npz:
            raise FileNotFoundError("colab_pack.npz not found in archive")
        missing_tsvs = self._REQUIRED_TSV - set(raw_tsvs)
        if missing_tsvs:
            raise FileNotFoundError(
                f"Required TSV files missing: {missing_tsvs}"
            )

        self.coords       = npz_data["coords"]
        self.target_seqs  = npz_data["seq"]
        self.target_dssp  = npz_data["dssp"]

        self._query_meta   = raw_tsvs["query_meta"]
        self._summary      = raw_tsvs["summary"]
        self._id_mapping   = raw_tsvs["id_mapping"]
        self._pack_index   = raw_tsvs["colab_pack_index"]
        self._raw_segments = raw_tsvs["segments"]
        self._pfam         = raw_tsvs.get("pfam", None)
        self._pfam_names   = raw_tsvs.get("mini_pfam_names", None)
        self._ligands         = raw_tsvs.get("mini_ligand", None)
        self._ligand_names   = raw_tsvs.get("mini_ligand_names", None)

        # ---------------------------------------------------------------- #
        # Pre-aggregate Pfam/Clan tuples at the Bronze level               #
        # ---------------------------------------------------------------- #
        if self._pfam is not None:
            # Create a target-level map of the string tuples
            pfam_summary = (
                self._pfam.groupby('target_name', sort=False)
                .agg(
                    pfam_tuples=('query_accession', lambda x: tuple(dict.fromkeys(x.dropna()))),
                    clan_tuples=('clan', lambda x: tuple(dict.fromkeys(x.dropna())))
                )
                .reset_index()
            )

            # Map them directly onto self._summary (or self._id_mapping)
            # depending on which key your silver builder uses for alignment joins.
            join_key = 'target_name' if 'target_name' in self._id_mapping.columns else 'target_id'

            self._id_mapping = self._id_mapping.merge(
                pfam_summary,
                left_on=join_key,
                right_on='target_name', 
                how='left'
            )

            # Clean up the extra merge key if they didn't match names exactly
            if join_key != 'target_name' and 'target_name' in self._id_mapping.columns:
                self._id_mapping.drop(columns=['target_name'], inplace=True)


        if not self.query_pdb_str:
            raise FileNotFoundError("Query.pdb not found in archive")
        if self._pfam is None or self._pfam.empty:
            self._pfam = None
            self._id_mapping["pfam_tuples"] = None
            self._id_mapping["clan_tuples"] = None

        self._bronze_loaded = True
        print(
            f"Pack loaded: {len(self.coords)} coord rows, "
            f"{len(self._raw_segments)} segments, "
            f"{len(self._summary)} alignments, "
            f"Query.pdb: {'yes' if self.query_pdb_str else 'not found'}"
        )

    def _seriate(self):
        pileups = self.df["dssp_pileup"].to_list()
        n_rows = len(pileups)
        n_pos = len(pileups[0])

        # convert to fixed-width uint8 array for Numba
        pileup_arr = np.zeros((n_rows, n_pos), dtype=np.uint8)
        for i, row in enumerate(pileups):
            pileup_arr[i, :] = np.frombuffer(row.encode("ascii"), dtype=np.uint8)

        bitmat = encode_dssp_optimized(pileup_arr, n_rows, n_pos)
        order = nearest_neighbor_reordering_fast(bitmat)
        return order

    # ---------------------------------------------------------------- #
    # Silver — joins and derived fields                                 #
    # ---------------------------------------------------------------- #

    def _build_silver(self) -> None:
        if not self._bronze_loaded:
            raise RuntimeError("Call _load_bronze() first")

        self.query_length   = int(self._query_meta["length"].iloc[0])
        self.query_sequence = str(self._query_meta["sequence"].iloc[0])
        self.query_ascii    = np.array(
            [ord(c) for c in self.query_sequence], dtype=np.int32
        )

        if len(self.query_ascii) != self.query_length:
            raise ValueError(
                f"query_meta length={self.query_length} but "
                f"sequence has {len(self.query_ascii)} characters"
            )

        if len(self.query_ca_coords) != self.query_length:
            raise ValueError("Query.pdb C-alpha count does not match query_meta length")

        # normalise segments: 1-based (DALI) → 0-based half-open
        segs            = self._raw_segments.copy()
        segs["q_start"] = segs["q_start"] - 1
        segs["s_start"] = segs["s_start"] - 1
        self.segments   = segs

        # master metadata join
        try:
            self.master_metadata = (
                self._summary[["alignment_id", "z_score", "protein_id"]]
                .merge(
                    self._id_mapping[["protein_id", "target_id", "pfam_tuples", "clan_tuples"]],
                    on="protein_id", how="inner"
                )
                .merge(
                    self._pack_index[
                        ["dali_id", "length", "start_idx", "end_idx", "compnd"]
                    ].rename(columns={
                        "dali_id": "target_id",
                        "compnd":  "description",
                        "length":  "target_length",
                    }),
                    on="target_id", how="inner"
                )
                .reset_index(drop=True)
            )
        except KeyError as e:
            raise KeyError(f"Column missing during silver join: {e}") from e

        # bounds validation
        max_offset = (
            self.master_metadata["start_idx"] +
            self.master_metadata["target_length"]
        ).max()
        if max_offset > len(self.target_seqs):
            raise ValueError(
                f"start_idx + target_length exceeds flat array "
                f"({max_offset} > {len(self.target_seqs)})"
            )

        n_dropped = len(self._summary) - len(self.master_metadata)

        self._compute_superimpositions()
        self._compute_pfam()
        self._build_pfam_domains()
        self._build_pileups()

        self._silver_built = True
        print(
            f"Silver built: {len(self.master_metadata)} alignments "
            f"({n_dropped} dropped), query length={self.query_length}"
        )

        if False:
            print("self._summary:", len(self._summary))
            print("self._id_mapping:", len(self._id_mapping))
            print("self._pack_index:", len(self._pack_index))
            print("self.master_metadata:", len(self.master_metadata))

    def _compute_superimpositions(self) -> None:
        rmsd_list, R_list, t_list = [], [], []
        seg_groups = self.segments.groupby("alignment_id")

        for _, row in self.master_metadata.iterrows():
            aln_id    = row["alignment_id"]
            start_idx = int(row["start_idx"])

            if aln_id not in seg_groups.groups:
                rmsd_list.append(np.nan); R_list.append(None); t_list.append(None)
                continue

            segs = seg_groups.get_group(aln_id)
            q_indices = np.concatenate([
                np.arange(qs, qs + l)
                for qs, l in zip(segs["q_start"].values, segs["length"].values)
            ]).astype(np.int64)
            s_indices = np.concatenate([
                np.arange(start_idx + ss, start_idx + ss + l)
                for ss, l in zip(segs["s_start"].values, segs["length"].values)
            ]).astype(np.int64)

            valid = (
                (q_indices >= 0) & (q_indices < self.query_length) &
                (s_indices >= 0) & (s_indices < len(self.coords))
            )
            q_indices = q_indices[valid]
            s_indices = s_indices[valid]

            if len(q_indices) < 3:
                rmsd_list.append(np.nan); R_list.append(None); t_list.append(None)
                continue

            Q = self.query_ca_coords[q_indices]
            P = self.coords[s_indices]

            try:
                R, t, rmsd = svd_superimpose(P, Q)
                rmsd_list.append(rmsd); R_list.append(R); t_list.append(t)
            except (np.linalg.LinAlgError, ValueError):
                rmsd_list.append(np.nan); R_list.append(None); t_list.append(None)

        self.master_metadata["rmsd"] = rmsd_list
        self.master_metadata["R"]    = R_list
        self.master_metadata["t"]    = t_list
        n_ok = self.master_metadata["rmsd"].notna().sum()
        print(f"  Superimpositions: {n_ok}/{len(self.master_metadata)} succeeded")

    def _compute_pfam_clan(self) -> None:
        if self._pfam is None:
            self.master_metadata["pfam"] = "Unassigned"
            self.master_metadata["clan"] = "Unassigned"
            print("  Pfam/clan: skipped (not in archive)")
            return

        # DataFrame indexed by target_id with columns ["pfam", "clan"]
        best = best_pfam_clan_per_target(self._pfam)

        # Join onto master_metadata via target_id
        self.master_metadata = self.master_metadata.join(
            best,
            on="target_id",
        )

        # Replace missing assignments
        self.master_metadata["pfam"] = self.master_metadata["pfam"].fillna("Unassigned")
        self.master_metadata["clan"] = self.master_metadata["clan"].fillna("Unassigned")

        n_pfam = (self.master_metadata["pfam"] != "Unassigned").sum()
        n_clan = (self.master_metadata["clan"] != "Unassigned").sum()

        print(
            f"  Pfam: {n_pfam}/{len(self.master_metadata)} assigned; "
            f"Clan: {n_clan}/{len(self.master_metadata)} assigned"
        )

    def _compute_pfam(self) -> None:
        if self._pfam is None:
            self.master_metadata["pfam"] = None
            print("  Pfam: skipped (not in archive)")
            return

        best_pfam_df = best_pfam_per_alignment(
            pfam_df         = self._pfam,
            segments_df     = self.segments,
            master_metadata = self.master_metadata,
        )

        self.master_metadata = self.master_metadata.merge(
            best_pfam_df, on='alignment_id', how='left'
        )

        n = self.master_metadata['pfam'].notna().sum()
        print(f"  Pfam: {n}/{len(self.master_metadata)} alignments assigned")


    def _build_pfam_domains(self) -> None:
        """
        Build ordered domain-architecture lists per target, at both
        Pfam-family and clan resolution. Merges 'pfam_domains' and
        'clan_domains' into master_metadata.
        """
        if self._pfam is None:
            self.master_metadata['pfam_domains'] = None
            self.master_metadata['clan_domains'] = None
            print("  Pfam domains: skipped (no pfam TSV in archive)")
            return

        from ..mechanics.metrics import build_pfam_domains

        pfam_df = build_pfam_domains(self._pfam, self._pfam_names, level='pfam')
        clan_df = build_pfam_domains(self._pfam, self._pfam_names, level='clan')

        self.master_metadata = (
            self.master_metadata
            .merge(pfam_df, on="target_id", how="left")
            .merge(clan_df, on="target_id", how="left")
        )

        n_pfam = self.master_metadata['pfam_domains'].notna().sum()
        n_clan = self.master_metadata['clan_domains'].notna().sum()
        print(f"  Pfam domains: {n_pfam}/{len(self.master_metadata)} targets")
        print(f"  Clan domains: {n_clan}/{len(self.master_metadata)} targets")


    def _build_pileups(self) -> None:
        pileup_df = build_pileups(
            segments        = self.segments,
            master_metadata = self.master_metadata,
            target_seqs     = self.target_seqs,
            target_dssp     = self.target_dssp,
            query_length    = self.query_length,
        )
        self.master_metadata = self.master_metadata.merge(
            pileup_df, on="alignment_id", how="left"
        )
        n_seq  = self.master_metadata["sequ_pileup"].notna().sum()
        n_dssp = self.master_metadata["dssp_pileup"].notna().sum() \
                 if "dssp_pileup" in self.master_metadata.columns else 0
        print(f"  Pileups: {n_seq} sequ, {n_dssp} dssp")

    # ---------------------------------------------------------------- #
    # Default full view                                                 #
    # ---------------------------------------------------------------- #

    def _register_full_view(self) -> None:
        from ..mechanics.metrics import register_domain_view
        pack_stem = Path(self.pack_path).stem.replace(".tar", "")
        view_name = f"dom_0-{self.query_length}" # #f"full_{pack_stem}"
        project = self
        domain_ranges = [(0, self.query_length)]
        register_domain_view( project, domain_ranges, view_name)
        print(f"  Default view '{view_name}' registered")


    # ---------------------------------------------------------------- #
    # View management                                                   #
    # ---------------------------------------------------------------- #

    def add_view(self, name: str, df: pd.DataFrame,
                 parent: str = None,
                 function: str = "",
                 parameters: dict = None) -> None:
        """
        Register a named view with provenance.

        Parameters
        ----------
        name       : unique view name
        df         : materialized DataFrame
        parent     : name of parent view this was derived from
        function   : name of function that created this view
        parameters : parameters passed to that function
        """
        if name in self.views:
            self.list_inventory()
            print(
                f"View '{name}' already exists — "
                f"choose a different name or delete it first: "
                f"del project.views['{name}']"
            )
            return False
        else:
            self.views[name] = df.copy()
            self.provenance[name] = Provenance(
                function        = function,
                parameters      = parameters or {},
                N               = len(df),
                parent_objects  = [parent] if parent else [],
            )
            return True

    def add_subset(self, name: str,
                   mask: pd.Series,
                   parent: str,
                   function: str = "row_filter",
                   parameters: dict = None) -> pd.DataFrame:
        """
        Register a row-filtered subset of an existing view.
        No recomputation — same columns as parent.

        mask is a Boolean Data Series, e.g.
        
        mask = parent_df["some_column"] > threshold
        mask = (parent_df["z_score"] > 5) & (parent_df["alignment_length"] > 50)

        Returns the subset DataFrame.
        """
        parent_df = self.views[parent]
        subset_df = parent_df[mask].reset_index(drop=True)
        _ = self.add_view(name, subset_df,
                      parent=parent,
                      function=function,
                      parameters=parameters or {})
        if _: print(f"  Subset '{name}': {len(subset_df)}/{len(parent_df)} rows")
        return subset_df

    def add_custom_subset(
        self,
        query_str: str,
        subset_name: str = "CUSTOM_SUBSET",
        source_view: str = 'full',
    ) -> pd.Series:
        """Evaluates a query string against a source view and registers the resulting subset.

        Parameters
        ----------
        query_str : str
            Pandas query expression (e.g., "z_score > 10 and rmsd < 5").
        subset_name : str, default "CUSTOM_SUBSET"
            Name under which to register the new view in `project.views`.
        source_view : str, default "full"
            Name of the existing view to query against.

        Returns
        -------
        pd.Series
            Boolean mask corresponding to the target hits in `source_view`.
        """
        if source_view not in self.views:
            source_view = self.list_inventory().iloc[0]['name']
            print(f"Source view '{source_view}' not found in project.views, using full data set.")

        source_df = self.views[source_view]

        # 1. Do-nothing default for empty query strings
        if not query_str or not query_str.strip():
            mask = pd.Series(False, index=source_df.index)
            print(f"No query string provided. Created empty subset '{subset_name}'.")
            return mask

        # 2. Evaluate boolean condition safely
        try:
            mask = source_df.eval(query_str)
            if not isinstance(mask, pd.Series) or mask.dtype != bool:
                raise ValueError(
                    "Query expression must evaluate to a boolean condition."
                )
        except Exception as e:
            print(f"Error evaluating query '{query_str}': {e}")
            mask = pd.Series(False, index=source_df.index)
            return mask

        # 3. Log statistics
        selected_count = mask.sum()
        print(
            f"Subset '{subset_name}': selected {selected_count:,} / {len(source_df):,} target hits."
        )

        # 4. Register new view internally
        self.add_subset(
            subset_name,
            mask,
            parent=source_view,
            function="add_custom_subset",
            parameters={"query_str": query_str, "source_view": source_view},
        )

        return mask

    # ---------------------------------------------------------------- #
    # Model management                                                  #
    # ---------------------------------------------------------------- #

    def add_model(self, name: str, artifact: Any,
                  parent: str = None,
                  function: str = "",
                  parameters: dict = None) -> None:
        """Register a named analytical artifact with provenance."""
        self.models[name] = artifact
        self.provenance[name] = Provenance(
            function        = function,
            parameters      = parameters or {},
            parent_objects  = [parent] if parent else [],
        )

    # ---------------------------------------------------------------- #
    # Plot management                                                   #
    # ---------------------------------------------------------------- #

    def bookmark_plot(self, name: str,
                      renderer: str,
                      source: str,
                      params: dict,
                      caption: str = "",
                      notes: str = "") -> PlotSpec:
        """
        Bookmark a plot as a recipe. Does not store the figure.

        Parameters
        ----------
        name     : unique plot name
        renderer : fully qualified renderer function name
        source   : view name used as data source
        params   : keyword arguments passed to renderer
        caption  : short title for the report
        notes    : user's scientific observations (multi-line text)
        """
        spec = PlotSpec(
            name       = name,
            renderer   = renderer,
            source     = source,
            params     = params,
            caption    = caption,
            notes      = notes,
            provenance = Provenance(
                function       = "bookmark_plot",
                parameters     = params,
                parent_objects = [source],
            )
        )
        self.plots[name] = spec
        return spec

    def export_figure(self, name: str,
                      session_dir: str,
                      fmt: str = "png",
                      dpi: int = 150) -> Path:
        """
        Render a bookmarked plot and save to disk.
        Requires the renderer function to be importable.
        """
        import importlib
        spec     = self.plots[name]
        out_dir  = Path(session_dir) / "figures"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{_safe_name(name)}.{fmt}"

        # resolve renderer
        module_path, func_name = spec.renderer.rsplit(".", 1)
        module   = importlib.import_module(module_path)
        renderer = getattr(module, func_name)

        # render
        df  = self.views[spec.source]
        fig = renderer(df, **spec.params)

        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        spec.png_path = str(out_path)
        print(f"  Exported '{name}' → {out_path}")
        return out_path

    # ---------------------------------------------------------------- #
    # Inventory                                                         #
    # ---------------------------------------------------------------- #

    def list_inventory(self) -> pd.DataFrame:
        """Summary of all registered assets."""
        rows = []
        for name, df in self.views.items():
            prov = self.provenance.get(name, Provenance())
            rows.append({
                "name":     name,
                "type":     "view",
                "rows":     len(df),
                "parent":   prov.parent_objects[0] if prov.parent_objects else "",
                "function": prov.function,
                "params":   prov.parameters,
                "created":  prov.created[:16],
            })
        for name, model in self.models.items():
            prov = self.provenance.get(name, Provenance())
            rows.append({
                "name":     name,
                "type":     "model",
                "rows":     "",
                "parent":   prov.parent_objects[0] if prov.parent_objects else "",
                "function": prov.function,
                "created":  prov.created[:16],
            })
        for name, spec in self.plots.items():
            rows.append({
                "name":     name,
                "type":     "plot",
                "rows":     spec.source,
                "parent":   spec.source,
                "function": spec.renderer.split(".")[-1],
                "created":  spec.provenance.created[:16],
            })
        return pd.DataFrame(rows)

    # ---------------------------------------------------------------- #
    # Persistence                                                       #
    # ---------------------------------------------------------------- #

    def save(self, session_dir: str, overwrite: bool = False) -> Path:
        """
        Save the project to a session directory.

        Layout
        ------
        session_dir/
            manifest.json
            views/
                {name}.parquet
                {name}_arrays.pkl   (for numpy object columns)
            models/
                {name}.pkl
            figures/
                (populated by export_figure())
        """
        out = Path(session_dir)
        if out.exists() and not overwrite:
            raise FileExistsError(
                f"{out} exists — pass overwrite=True"
            )
        if out.exists() and overwrite:
            shutil.rmtree(out)

        (out / "views").mkdir(parents=True)
        (out / "models").mkdir(parents=True)
        (out / "figures").mkdir(parents=True)

        manifest = {
            "saved_at":          datetime.now().isoformat(),
            "daliscope_version": __version__,
            "pack_path":         self.pack_path,
            "views":             {},
            "models":            {},
            "plots":             {},
        }

        # --- views ---
        for name, df in self.views.items():
            safe = _safe_name(name)
            array_cols  = [c for c in df.columns
                           if df[c].dtype == object
                           and _col_contains_arrays(df[c])]
            simple_cols = [c for c in df.columns if c not in array_cols]

            parquet_path = out / "views" / f"{safe}.parquet"
            df[simple_cols].to_parquet(parquet_path, index=True)

            arrays_path = None
            if array_cols:
                arrays_path = out / "views" / f"{safe}_arrays.pkl"
                with open(arrays_path, "wb") as f:
                    pickle.dump(
                        {c: df[c].tolist() for c in array_cols}, f,
                        protocol=pickle.HIGHEST_PROTOCOL
                    )

            prov = self.provenance.get(name, Provenance())
            manifest["views"][name] = {
                "parquet":    parquet_path.name,
                "arrays":     arrays_path.name if arrays_path else None,
                "array_cols": array_cols,
                "row_count":  len(df),
                "columns":    list(df.columns),
                "provenance": prov.to_dict(),
            }

        # --- models ---
        for name, artifact in self.models.items():
            safe = _safe_name(name)
            pkl_path = out / "models" / f"{safe}.pkl"
            with open(pkl_path, "wb") as f:
                pickle.dump(artifact, f, protocol=pickle.HIGHEST_PROTOCOL)
            prov = self.provenance.get(name, Provenance())
            manifest["models"][name] = {
                "pkl":        pkl_path.name,
                "provenance": prov.to_dict(),
            }

        # --- plots ---
        for name, spec in self.plots.items():
            manifest["plots"][name] = spec.to_dict()

        with open(out / "manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)

        print(
            f"Project saved → {out}\n"
            f"  {len(self.views)} views, "
            f"{len(self.models)} models, "
            f"{len(self.plots)} plots"
        )
        return out

    @classmethod
    def open(cls, session_dir: str) -> "Project":
        """
        Restore a project from a session directory.
        Pack data (coords, sequences) must be reloaded from the
        original pack_path recorded in the manifest.

        Parameters
        ----------
        session_dir : directory written by save()

        Returns
        -------
        Project with views, models, and plots restored.
        Pack is reloaded from its original path.
        """
        out = Path(session_dir)
        manifest_path = out / "manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(f"No manifest.json in {out}")

        with open(manifest_path) as f:
            manifest = json.load(f)

        saved_ver = manifest.get("daliscope_version", "0.0")
        if saved_ver != __version__:
            print(f"WARNING: saved with version {saved_ver}, "
                  f"current is {__version__}")

        project = cls()

        # reload pack from original path
        pack_path = manifest.get("pack_path")
        project.pack_path = pack_path
        if pack_path and Path(pack_path).exists():
            print(f"Reloading pack from {pack_path} ...")
            project._load_bronze(pack_path)
            project._build_silver()
        else:
            print(f"WARNING: pack not found at {pack_path!r} — "
                  f"pack data unavailable. Views are still accessible.")

        # restore views
        for name, meta in manifest["views"].items():
            parquet_path = out / "views" / meta["parquet"]
            df           = pd.read_parquet(parquet_path)
            if meta.get("arrays"):
                arrays_path = out / "views" / meta["arrays"]
                with open(arrays_path, "rb") as f:
                    array_data = pickle.load(f)
                for col, values in array_data.items():
                    df[col] = values
            df = df[meta.get("columns", list(df.columns))]
            project.views[name]      = df
            project.provenance[name] = Provenance.from_dict(
                meta.get("provenance", {})
            )

        # restore models
        for name, meta in manifest["models"].items():
            pkl_path = out / "models" / meta["pkl"]
            with open(pkl_path, "rb") as f:
                project.models[name] = pickle.load(f)
            project.provenance[name] = Provenance.from_dict(
                meta.get("provenance", {})
            )

        # restore plot specs
        for name, d in manifest["plots"].items():
            project.plots[name] = PlotSpec.from_dict(d)

        print(
            f"Project opened ← {out} "
            f"(saved {manifest['saved_at'][:16]})\n"
            f"  {len(project.views)} views, "
            f"{len(project.models)} models, "
            f"{len(project.plots)} plots"
        )
        return project

    def retrieve_coords_sequ_dssp(self, protein_id=None, target_id=None):
        """Get CA coordinates, sequence and dssp strings from project
        Either protein_id or target_id required.
        """
        # default: empty result
        xyz = np.empty((0, 3), dtype=float)
        sequ = ""
        dssp = ""
        if target_id is not None:
            match_col = 'target_id'
            match_item = target_id
        elif protein_id is not None:
            match_col = 'protein_id'
            match_item = protein_id
        else:
            print('# ERROR: protein_id or target_id must be defined.')
            return xyz, sequ, dssp
        try:
            start_idx, end_idx = self.master_metadata.loc[self.master_metadata[match_col] == match_item, ['start_idx', 'end_idx']].iloc[0]
        except:
            print(f" ERROR in retrieve_coords_sequ_dssp: {match_item}")
            return xyz, sequ, dssp

        sequ = self.target_seqs[start_idx:end_idx].astype(np.uint8).tobytes().decode('ascii')
        dssp = self.target_dssp[start_idx:end_idx].astype(np.uint8).tobytes().decode('ascii').replace('L','-') # for comparison
        xyz = self.coords[start_idx:end_idx,:]
        return xyz, sequ, dssp

    # ---------------------------------------------------------------- #
    # Diagnostics                                                       #
    # ---------------------------------------------------------------- #

    def status(self) -> None:
        print(f"DaliScope Project — {self.pack_path}")
        print(f"  Bronze loaded : {self._bronze_loaded}")
        print(f"  Silver built  : {self._silver_built}")
        if self._bronze_loaded:
            print(f"  coords        : {self.coords.shape}  "
                  f"dtype={self.coords.dtype}")
            print(f"  target_seqs   : {self.target_seqs.shape}")
            print(f"  Query.pdb     : "
                  f"{'yes' if self.query_pdb_str else 'not found'}")
        if self._silver_built:
            print(f"  query_length  : {self.query_length}")
            print(f"  master_meta   : {len(self.master_metadata)} rows, "
                  f"cols: {list(self.master_metadata.columns)}")
        print(f"  views         : {list(self.views.keys())}")
        print(f"  models        : {list(self.models.keys())}")
        print(f"  plots         : {list(self.plots.keys())}")


# ================================================================== #
# Module helpers                                                      #
# ================================================================== #

def _safe_name(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in name)


def _col_contains_arrays(col: pd.Series) -> bool:
    return any(isinstance(value, (np.ndarray, list, tuple, dict)) for value in col.dropna())


def initialize_project(pack_path: str, current_project: Optional[Project] = None,
                       verbose=True) -> Tuple[Optional[Project], Optional[pd.DataFrame], Optional[str]]:
    """
    Validates the path to a colab-pack, handles auto-saving/warnings for existing 
    projects, loads the new pack, and extracts the default view.

    Returns
    -------
    tuple (Project or None, DataFrame or None, view name or None)
        The newly loaded project instance and its default fallback view dataframe.
        Returns (None, None, None) if the file path is invalid.
    """
    if not os.path.exists(pack_path):
        print(f"❌ File not found: {pack_path}")
        print(f"   Current working directory: {os.getcwd()}")

        search_dirs = ['.', 'data', 'notebooks/data', os.path.dirname(pack_path) or '.']

        # dedupe by absolute path, preserving first-seen order
        seen = set()
        unique_dirs = []
        for d in search_dirs:
            abs_d = os.path.abspath(d)
            if abs_d not in seen:
                seen.add(abs_d)
                unique_dirs.append(d)

        found_any = False
        for d in unique_dirs:
            if os.path.isdir(d):
                candidates = [f for f in os.listdir(d) if f.endswith('.tar.gz')]
                if candidates:
                    found_any = True
                    print(f"   Found .tar.gz files in '{d}':")
                    for c in candidates:
                        print(f"     {os.path.join(d, c)}")

        if not found_any:
            print(f"   No .tar.gz files found in any of: {unique_dirs}")

        print("\nGo back to the User input cell and set a valid pack_path.\n")
        return None, None, None

    else:
        if current_project is not None:
            #>> auto-save previous project
            print("WARNING: Previous project data will be lost!")
            #>> force saving of previous project!
            # Example placeholder: current_project.save("auto_backup/")

        project = Project.load_pack(pack_path)
        view_name = project.list_inventory()['name'].iloc[0]
        view = project.views[view_name]

        return project, view, view_name

def validate_pipeline_state(*required_variable_names, glbs, verbose=True):
    """
    Surgically checks that the notebook state matches the loaded project context 
    and verifies that required downstream variables exist and are fresh.

    glbs = globals() - namespace in calling module

    # in the notebook
    validate_pipeline_state("PARENT_VIEW", "TRUSTED_VIEW", glbs=globals())
    """

    # 1. Base Project Check
    if "project" not in glbs or glbs["project"] is None:
        raise RuntimeError("❌ Pipeline Error: No active 'project' found. Run Act I first.")

    # 2. Base Dataframe Check
    if "FULL_DF" not in glbs or glbs["FULL_DF"] is None:
        raise RuntimeError("❌ Pipeline Error: 'FULL_DF' is not initialized. Run your Data Loading/Clipping cells.")

    if glbs["FULL_DF"].empty:
        raise RuntimeError("❌ Pipeline Error: 'FULL_DF' is empty.")

    # 3. Cross-Contamination Guardrail (Ensures CURRENT_DF belongs to this specific project pack)
    master_ids = set(glbs["project"].master_metadata["alignment_id"])
    current_ids = set(glbs["project"].views[glbs["FULL_VIEW"]]["alignment_id"]) if glbs["FULL_VIEW"] in glbs["project"].views else set(glbs["FULL_DF"]["alignment_id"])
    if verbose: print(len(master_ids), len(current_ids))

    stale_view_ids = set(glbs["project"].views[glbs["FULL_VIEW"]]["alignment_id"])
    full_df_ids = set(glbs["FULL_DF"]["alignment_id"])

    if not set(glbs["FULL_DF"]["alignment_id"]).issubset(master_ids):
        raise RuntimeError(
            f"❌ STATE CORRUPTION: The data in FULL_DF does not belong to the currently loaded project pack ({glbs['project']}).\n"
            "You likely loaded a new pack without resetting the pipeline variables. Please re-run the initialization cells."
        )

    # 4. Dynamic Dependency Check
    for var_name in required_variable_names:
        if var_name not in glbs or glbs[var_name] is None:
            raise RuntimeError(f"❌ Missing Dependency: '{var_name}' is not defined. Execute the preceding cells first.")
        if hasattr(glbs[var_name], 'empty') and glbs[var_name].empty:
            raise RuntimeError(f"❌ Empty Dependency: '{var_name}' contains an empty DataFrame. Check your filter thresholds.")
