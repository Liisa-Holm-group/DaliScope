"""
communities.py
--------------------------------
Pipeline for detecting structural communities in domain fold populations using
structural fingerprint indexing, KNN graph construction, and Infomap partitioning.
"""

import gc
import time
from infomap import Infomap
import numpy as np
import pandas as pd
import daliscope
import matplotlib.pyplot as plt
import seaborn as sns

import daliscope.analysis.structural_fingerprints
import daliscope.analysis.structural_fingerprint_patch
import daliscope.analysis.voxel

def extract_infomap_hierarchy(im: Infomap) -> pd.DataFrame:
    """Extracts node mapping and hierarchical paths from an Infomap instance.

    Parameters
    ----------
    im : Infomap
        An Infomap instance after calling `im.run()`.

    Returns
    -------
    pd.DataFrame
        DataFrame mapping node_id/target_idx to hierarchical module labels.
    """
    records = []
    for node in im.tree:
        if node.is_leaf:
            path_tuple = node.path
            path_str = ".".join(map(str, path_tuple))

            records.append(
                {
                    "target_idx": node.node_id,
                    "tree_path": path_str,
                    "module_level_1": (
                        path_tuple[0] if len(path_tuple) > 0 else np.nan
                    ),
                    "module_level_2": (
                        path_tuple[1] if len(path_tuple) > 1 else np.nan
                    ),
                    "module_level_3": (
                        path_tuple[2] if len(path_tuple) > 2 else np.nan
                    ),
                    "flow": node.flow,
                }
            )

    return pd.DataFrame(records)

def run_community_detection(
    view_name: str,
    domain_range: tuple[int, int] | list[tuple[int, int]],
    project,
    k_neighbors: int = 200,
    cell_size: float = 4.5,
    radial_r0: float = 20.0,
    infomap_args: str = "--two-level --num-trials 10 --seed 42 --no-self-links",
    max_workers: int = 1,
    verbose: bool = True,
) -> tuple[pd.DataFrame, int]:
    """Executes structural fingerprint indexing, KNN graph construction, and 
    Infomap community detection for a given project view.

    Parameters
    ----------
    view_name : str
        The target view name in `project.views` (e.g., 'dom_FOLD').
    domain_range : tuple[int, int] or list[tuple[int, int]]
        Boundary tuple or list of discontiguous domain segment tuples 
        (e.g., (23, 123) or [(23, 50), (60, 123)]).
    project : object
        The DaliScope project instance containing coordinates and views.
    k_neighbors : int, default=200
        Number of nearest neighbors to retrieve per target in `topk_threaded`.
    cell_size : float, default=4.5
        Grid cell size for voxel indexing.
    radial_r0 : float, default=20.0
        Auxiliary domain shell radius size.
    infomap_args : str, default='--two-level --num-trials 10 --seed 42 --no-self-links'
        Command-line flags passed directly to the Infomap instance.
    max_workers : int, default=1
        Maximum threads for structural alignment topk retrieval.
    verbose : bool, default=True
        Whether to print timing, throughput, and progress logs.

    Returns
    -------
    tuple[pd.DataFrame, int]
        - modules_df: DataFrame containing target IDs and assigned module IDs.
        - n_modules: Total count of identified structural modules.
    """
    # 1. Normalize domain_range input
    if isinstance(domain_range, tuple) and len(domain_range) == 2 and isinstance(domain_range[0], int):
        domain_ranges = [domain_range]
    else:
        domain_ranges = list(domain_range)

    global_anchor_q_start = min(seg[0] for seg in domain_ranges)
    global_anchor_q_end = max(seg[1] for seg in domain_ranges)

    # 2. Prepare target slice dataframe
    df = project.views[view_name][
        [
            "target_id",
            "R",
            "t",
            "start_idx",
            "end_idx",
            "alignment_id",
            "target_length",
            "sequ_pileup",
        ]
    ].copy()

    # 3. Build Query Mask & Center of Mass (COM) across all segments
    query_mask = (
        daliscope.analysis.structural_fingerprints.build_query_mask(
            project.query_length, domain_ranges
        )
    )
    query_com = daliscope.analysis.structural_fingerprints.compute_query_com(
        project.query_ca_coords, query_mask
    )

    # 4. Build Voxel Anchor Masks using global boundaries
    anchor_data = daliscope.analysis.voxel.build_anchor_mask(
        project,
        df,
        anchor_q_start=global_anchor_q_start,
        anchor_q_end=global_anchor_q_end,
    )

    target_labels = list(anchor_data.keys())
    all_anchor_masks = [anchor_data[t][0] for t in target_labels]
    all_coords = [anchor_data[t][1] for t in target_labels]

    # 4. Fingerprint Index Construction
    t0 = time.time()
    if verbose:
        print("Running build_structural_fingerprint_index ...")

    idx = (
        daliscope.analysis.structural_fingerprints.build_structural_fingerprint_index(
            target_labels,
            all_coords,
            all_anchor_masks,
            query_com=query_com,
            cell_size=cell_size,
            radial_r0=radial_r0,
            idf_alpha=1.0,
            idf_cap=2.0,
        )
    )

    # 5. Top-K Alignment & Graph Construction
    t1 = time.time()
    if verbose:
        print(f"Running {len(target_labels) * k_neighbors} pairwise STRUCTAL alignments ...")
        print("Runtime depends on the selected population, structure sizes, and hardware.")

    neighbors, scores = (
        daliscope.analysis.structural_fingerprint_patch.topk_threaded(
            idx,
            k=k_neighbors,
            max_workers=max_workers,
        )
    )

    t2 = time.time()
    source, target, score = (
        daliscope.analysis.structural_fingerprint_patch.build_knn_graph_fast(
            neighbors, scores
        )
    )

    links = np.column_stack((source, target, score))
    t3 = time.time()

    if verbose:
        print(
            f"Timings (sec): Index={t1-t0:.2f}s, Alignment={t2-t1:.2f}s, "
            f"Graph={t3-t2:.2f}s | Total={t3-t0:.2f}s"
        )

    # 6. Infomap Partitioning
    im = Infomap(infomap_args)
    im.add_links(links)
    im.run()

    # Extract base modules
    base_modules = pd.DataFrame(
        [(node.node_id, node.module_id) for node in im.nodes],
        columns=["target_idx", "module_raw"],
    )

    # Extract hierarchical paths and merge
    hierarchy_df = extract_infomap_hierarchy(im)
    modules_df = base_modules.merge(hierarchy_df, on="target_idx", how="left")
    modules_df["module"] = modules_df["module_level_1"].astype("Int64")

    # Map target_id strings back using original target_df ordering
    target_df = project.views[view_name]
    modules_df["target_id"] = (
        target_df["target_id"].iloc[modules_df["target_idx"]].values
    )

    # Clean Infomap object from memory
    del im
    gc.collect()

    # 7. Add Subset Views back into project
    n_modules = len(modules_df["module"].value_counts())
    parent_df = project.views[view_name]

    for module_id in range(1, n_modules + 1):
        target_ids = modules_df.loc[
            modules_df["module"] == module_id, "target_id"
        ]
        mask = parent_df["target_id"].isin(target_ids)
        sub_name = f"{view_name}_module_{module_id}"

        # Overwrite subset view if present
        if sub_name in project.views:
            del project.views[sub_name]

        project.add_subset(
            name=sub_name,
            parent=view_name,
            mask=mask,
            function="structural_infomap",
        )

    # Add alignment_id, clan, and pfam by indexing target_df directly via positional indices
    modules_df["description"] = parent_df["description"].iloc[modules_df["target_idx"]].values
    modules_df["target_id"] = parent_df["target_id"].iloc[modules_df["target_idx"]].values
    modules_df["clan"] = parent_df["clan"].iloc[modules_df["target_idx"]].values
    modules_df["pfam"] = parent_df["pfam"].iloc[modules_df["target_idx"]].values
    modules_df["clan_domains"] = parent_df["clan_domains"].iloc[modules_df["target_idx"]].values
    modules_df["pfam_domains"] = parent_df["pfam_domains"].iloc[modules_df["target_idx"]].values

    if verbose:
        print(f"Successfully generated {n_modules} module subsets in project.views.")

    return modules_df, n_modules

#------------------

def plot_module_score_distributions(
    raw_scores: dict[str | int, dict[str, list[float] | np.ndarray]],
    figsize: tuple[float, float] = (10, 6),
    showfliers: bool = False,
    title: str = "Within- and between-module STRUCTAL scores",
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Plots side-by-side boxplots comparing 'within' vs. 'between' module STRUCTAL scores.

    Parameters
    ----------
    raw_scores : dict
        Dictionary structured as {module_id: {"within": [...], "between": [...]}}.
    figsize : tuple[float, float], default=(10, 6)
        Width and height of the generated matplotlib figure (ignored if `ax` is supplied).
    showfliers : bool, default=False
        Whether to show outliers beyond the whiskers in the boxplot.
    title : str, default="Within- and between-module STRUCTAL scores"
        Title for the plot.
    ax : plt.Axes or None, default=None
        An existing Matplotlib Axes instance. If None, a new figure and axis are created.

    Returns
    -------
    tuple[plt.Figure, plt.Axes]
        The Matplotlib Figure and Axes objects containing the generated plot.
    """
    modules = sorted(raw_scores.keys())

    data = []
    positions = []
    positions_within = []

    for i, module in enumerate(modules):
        x = i * 2.0 + 1
        data.extend([
            raw_scores[module]["within"],
            raw_scores[module]["between"],
        ])
        positions.extend([x, x + 0.65])
        positions_within.append(x)

    positions_within = np.array(positions_within)

    # Manage figure creation vs passing existing axis
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    ax.boxplot(
        data,
        positions=positions,
        widths=0.5,
        showfliers=showfliers,
    )

    # Center label ticks halfway between the 'within' and 'between' boxes
    ax.set_xticks(positions_within + 0.325)
    ax.set_xticklabels(modules)
    ax.set_xlabel("Module")
    ax.set_ylabel("STRUCTAL score")
    ax.set_title(title)

    fig.tight_layout()

    return fig, ax

def module_pfam_table(project, modules_df, row_group='clab', min_count=0):
    # display colored counts table of pfam category vs. module
    if row_group == "pfam": 
        desc_col = "name"
    else:
        desc_col = "clan_short"
    df = (
        modules_df.value_counts(["module", "clan", "pfam"])
        .reset_index(name="count")
        .merge(project._pfam_names[["clan", "clan_short"]], on="clan", how="left")
        .merge(project._pfam_names[["pfam", "short", "name"]], on="pfam", how="left")
        .sort_values(by=["module", "count"], ascending=[True, False])
        .drop_duplicates()
    )
    df = df[ df['count'] > min_count]

    # 1. Create a composite label column: e.g., "CL0001 (Vp1_capsid)"
    df["row_label"] = (
        df[row_group].astype(str) + " (" + df[desc_col].fillna("").astype(str) + ")"
    )

    # 2. Pivot using the composite row label
    pivot_counts = df.pivot_table(
        index="row_label",
        columns="module",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    log_counts = np.log1p(pivot_counts)

    # 3. Plot Seaborn Heatmap
    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot_counts) * 0.4)))

    sns.heatmap(
        log_counts,
        annot=pivot_counts,
        fmt="g",
        cmap="YlGnBu",
        linewidths=0.5,
        #cbar_kws={"label": "log(1 + count)"},
        ax=ax,
        cbar=False,
    )

    # 4. Move y-axis labels to the right side
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    plt.yticks(rotation=0)  # Keep text horizontal

    plt.title(
        f"Module Counts Stratified by {row_group.capitalize()}", pad=20
    )
    plt.xlabel("Module")
    plt.ylabel(f"{row_group.capitalize()} ({desc_col})")
    plt.tight_layout()

    plt.show()
