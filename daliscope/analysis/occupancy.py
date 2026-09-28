import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def compute_occupancy_profile(df: pd.DataFrame, col: str = 'dssp_pileup') -> np.ndarray:
    """
    Computes a 1D per-position occupancy profile array from a column of aligned strings.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the pileup string column.
    col : str
        Column name containing string alignments of equal length.

    Returns
    -------
    occupancy : np.ndarray (float64)
        1D array of shape (query_length,) containing fraction of non-gap residues [0.0, 1.0].
    """
    # 1. Filter out empty/null strings
    valid_series = df[col].dropna()
    if valid_series.empty:
        return np.array([], dtype=np.float64)

    # 2. Convert pandas Series directly to a contiguous 1D uint8 NumPy byte array
    # ASCII value of '-' is 45
    raw_bytes = np.frombuffer("".join(valid_series).encode("ascii"), dtype=np.uint8)

    # 3. Reshape into a 2D matrix: (n_targets, query_length)
    n_targets = len(valid_series)
    query_length = len(valid_series.iloc[0])
    matrix = raw_bytes.reshape(n_targets, query_length)

    # 4. Vectorized non-gap boolean mask (dash '-' is ASCII 45)
    non_gap_mask = (matrix != 46)

    # 5. Average non-gap occurrences vertically down columns
    occupancy = non_gap_mask.mean(axis=0)

    return occupancy

def plot_occupancy_profile(occupancy: np.ndarray, title: str = "Residue Occupancy Profile") -> None:
    """
    Plots a column/bar chart of the 1D residue occupancy profile.
    
    Parameters
    ----------
    occupancy : np.ndarray
        1D array of fraction values [0.0, 1.0] representing position-specific occupancy.
    title : str
        Chart title.
    """
    if len(occupancy) == 0:
        print("Empty occupancy array provided.")
        return

    residues = np.arange(1, len(occupancy) + 1)  # 1-based residue numbering

    plt.figure(figsize=(12, 4), dpi=100)
    plt.bar(residues, occupancy, width=1.0, color="#2b5c8f", edgecolor="none", alpha=0.85)

    plt.title(title, fontsize=12, fontweight="bold", pad=10)
    plt.xlabel("Query Residue Position", fontsize=10)
    plt.ylabel("Occupancy Fraction", fontsize=10)

    plt.ylim(0, 1.05)
    plt.xlim(1, len(occupancy))
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()

def filter_by_occupancy(
    df: pd.DataFrame,
    occupancy_profile: np.ndarray,
    cutoff: float = 0.8,
    pileup_col: str = "dssp_pileup",
    gamma: float = 0.4,
    percentile_cap: float = 80.0,
) -> tuple[pd.DataFrame, float, np.ndarray]:
    """Scores targets based on position-specific occupancy profile weights.

    Occupancy profile values are exponentially re-weighted (occupancy^gamma)
    and then capped at a percentile threshold to prevent highly over-represented
    short structural motifs (e.g., helix hairpins) from dominating the profile score.

    Parameters
    ----------
    df : pd.DataFrame
        Target DataFrame containing the pileup strings.
    occupancy_profile : np.ndarray
        1D array of shape (query_length,) from round 1.
    cutoff : float
        Normalized score threshold (0.0 to 1.0) for filtering.
    pileup_col : str
        Column name containing string alignments.
    gamma : float
        Exponent for non-linear weight dampening (default 0.4 ~ square root).
    percentile_cap : float
        Percentile rank (0-100) at which to cap the profile weights.

    Returns
    -------
    filtered_df : pd.DataFrame
        Subset of targets exceeding the score cutoff.
    max_score : float
        Maximum theoretical score (sum of transformed occupancy profile).
    norm_scores : np.ndarray
        Normalized scores [0.0, 1.0] for all rows in df.
    """
    #print("This is filter_by_occupancy", np.sum(occupancy_profile))
    valid_series = df[pileup_col].dropna()
    if valid_series.empty or len(occupancy_profile) == 0:
        return df.iloc[0:0], 0.0, np.array([])

    # 1. Exponential dampening (occupancy^gamma)
    weighted_profile = np.power(occupancy_profile, gamma)

    # 2. Cap transformed profile at the 80th percentile
    non_zero = weighted_profile[weighted_profile > 0]

    if len(non_zero) > 0:
        p_cap = np.percentile(non_zero, percentile_cap)
    else:
        p_cap = np.percentile(weighted_profile, percentile_cap)

    # Fallback if p_cap is still zero (e.g., all zeros)
    if p_cap == 0:
        capped_profile = weighted_profile
    else:
        capped_profile = np.minimum(weighted_profile, p_cap)

    # 3. Theoretical maximum score using the transformed profile
    #print("capped_profile", capped_profile)
    max_score = float(np.sum(capped_profile))
    #print("max_score", max_score)

    if max_score == 0:
        norm_scores = np.zeros(len(df), dtype=np.float64)
        return df.iloc[0:0], 0.0, norm_scores

    # 4. Vectorized 2D matrix conversion (n_targets, query_length)
    raw_bytes = np.frombuffer(
        "".join(valid_series).encode("ascii"), dtype=np.uint8
    )
    n_targets = len(valid_series)
    query_length = len(capped_profile)
    matrix = raw_bytes.reshape(n_targets, query_length)

    # 5. Non-gap mask: ignore '-' (45) and '.' (46)
    non_gap_mask = (matrix != 45) & (matrix != 46)
    #print("non_gap_mask", non_gap_mask)

    # 6. Target-wise score calculation (Matrix-Vector dot product)
    raw_scores = non_gap_mask @ capped_profile

    # 7. Normalize score (scale 0.0 to 1.0)
    norm_scores = raw_scores / max_score

    # 8. Apply cutoff mask
    mask = norm_scores > cutoff
    filtered_df = df.iloc[mask].copy()

    return filtered_df, max_score, norm_scores


def old_plot_score_distribution(scores: np.ndarray, cutoff: float = 0.8) -> None:
    """
    Plots a histogram of target normalized occupancy scores with a cutoff line.
    """
    if len(scores) == 0:
        print("No scores to plot.")
        return

    n_pass = np.sum(scores > cutoff)
    pct_pass = (n_pass / len(scores)) * 100

    plt.figure(figsize=(9, 4.5), dpi=100)

    # Histogram plot
    n, bins, patches = plt.hist(
        scores,
        bins=50,
        range=(0.0, 1.0),
        color="#34495e", 
        edgecolor="white", 
        alpha=0.75
    )

    # Highlight bars above cutoff in a distinct color
    for i in range(len(patches)):
        if bins[i] >= cutoff:
            patches[i].set_facecolor("#27ae60")
            patches[i].set_alpha(0.85)

    # Cutoff line
    plt.axvline(
        x=cutoff,
        color="#e74c3c", 
        linestyle="--", 
        linewidth=2,
        label=f"Cutoff ({cutoff}) — {n_pass:,} Targets ({pct_pass:.1f}%)"
    )

    plt.title("Target Normalized Occupancy Score Distribution", fontsize=12, fontweight="bold", pad=10)
    plt.xlabel("Normalized Occupancy Score [0.0 - 1.0]", fontsize=10)
    plt.ylabel("Target Count", fontsize=10)

    plt.xlim(0.0, 1.0)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.legend(loc="upper left", frameon=True, fontsize=10)

    plt.tight_layout()
    plt.show()


def old_filter_by_occupancy(
    df: pd.DataFrame,
    occupancy_profile: np.ndarray,
    cutoff: float = 0.8,
    pileup_col: str = 'dssp_pileup'
) -> tuple[pd.DataFrame, float, np.ndarray]:
    """
    Scores targets based on position-specific occupancy profile weights.
    
    Parameters
    ----------
    df : pd.DataFrame
        Target DataFrame containing the pileup strings.
    occupancy_profile : np.ndarray
        1D array of shape (query_length,) from round 1.
    cutoff : float
        Normalized score threshold (0.0 to 1.0) for filtering.
    pileup_col : str
        Column name containing string alignments.
        
    Returns
    -------
    filtered_df : pd.DataFrame
        Subset of targets exceeding the score cutoff.
    max_score : float
        Maximum theoretical score (sum of occupancy profile).
    norm_scores : np.ndarray
        Normalized scores [0.0, 1.0] for all rows in df.
    """
    valid_series = df[pileup_col].dropna()
    if valid_series.empty or len(occupancy_profile) == 0:
        return df.iloc[0:0], 0.0, np.array([])

    # 1. Theoretical maximum score
    max_score = float(np.sum(occupancy_profile))

    if max_score == 0:
        norm_scores = np.zeros(len(df), dtype=np.float64)
        return df.iloc[0:0], 0.0, norm_scores

    # 2. Vectorized 2D matrix conversion (n_targets, query_length)
    raw_bytes = np.frombuffer("".join(valid_series).encode("ascii"), dtype=np.uint8)
    n_targets = len(valid_series)
    query_length = len(occupancy_profile)
    matrix = raw_bytes.reshape(n_targets, query_length)

    # 3. Non-gap mask: ignore '-' (45) and '.' (46)
    non_gap_mask = (matrix != 45) & (matrix != 46)

    # 4. Target-wise score calculation (Matrix-Vector dot product)
    raw_scores = non_gap_mask @ occupancy_profile

    # 5. Normalize score (scale 0.0 to 1.0)
    norm_scores = raw_scores / max_score

    # 6. Apply cutoff mask
    mask = norm_scores > cutoff
    filtered_df = df.iloc[mask].copy()

    return filtered_df, max_score, norm_scores

def plot_score_distribution(scores: np.ndarray, cutoff: float = 0.8, title: str = None) -> None:
    """
    Plots a histogram of target normalized occupancy scores with a cutoff line.
    """
    if len(scores) == 0:
        print("No scores to plot.")
        return

    n_pass = np.sum(scores > cutoff)
    pct_pass = (n_pass / len(scores)) * 100

    plt.figure(figsize=(4,2), dpi=100)

    # Histogram plot
    n, bins, patches = plt.hist(
        scores,
        bins=50,
        range=(0.0, max(scores)),
        color="#34495e", 
        edgecolor="white", 
        alpha=0.75
    )

    # Highlight bars above cutoff in a distinct color
    for i in range(len(patches)):
        if bins[i] >= cutoff:
            patches[i].set_facecolor("#27ae60")
            patches[i].set_alpha(0.85)

    # Cutoff line
    plt.axvline(
        x=cutoff,
        color="#e74c3c", 
        linestyle="--", 
        linewidth=2,
        label=f"Cutoff ({cutoff}) — {n_pass:,} Targets ({pct_pass:.1f}%)"
    )

    plt.title(title, fontsize=12, fontweight="bold", pad=10)
    plt.xlabel("Normalized Occupancy Score [0.0 - 1.0]", fontsize=10)
    plt.ylabel("Target Count", fontsize=10)

    #plt.xlim(0.0, 1.0)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.legend(loc="upper left", frameon=True, fontsize=10)
    plt.yscale('log')

    plt.tight_layout()
    plt.show()

def create_fold_subset(project, view_name, show_plot=False, cutoff=0.5):
    "Generate occupancy profile of parent view, score targets, register _FOLD subset"
    parent_df = project.views[view_name]
    fold_view_name = f"{view_name}_FOLD"

    # Run scoring pass
    occupancy = compute_occupancy_profile(parent_df, col='dssp_pileup')
    if show_plot: plot_occupancy_profile(occupancy, f"{view_name} occupancy")
    fold_df, max_possible_score, target_norm_scores = filter_by_occupancy(
        df=parent_df,
        occupancy_profile=occupancy,
        cutoff=cutoff,
        pileup_col='dssp_pileup'
    )

    # Store subset into project views dictionary
    #print("add subset ",view_name, fold_view_name)
    project.add_subset(
        name = fold_view_name,
        mask = parent_df.index.isin(fold_df.index),
        parent = view_name,
        function = "filter_by_occupancy",
        parameters = { "cutoff": cutoff},
    )

    # Reporting - DISABLED
    #print(f"Maximum Possible Score (Sum of Occupancy Profile): {max_possible_score:.3f}")
    #print(f"Total Targets Evaluated: {len(FULL_DF)}")
    #print(f"Targets passing cutoff (> {cutoff}): {len(fold_df)} ({(len(fold_df)/len(FULL_DF))*100:.2f}%)")

    # --- Run Plot ---
    if show_plot: plot_score_distribution(target_norm_scores, cutoff=cutoff, title=fold_view_name)

def get_fold_matrix(project, view_list):
    """Generate count matrix of domain_FOLD co-occurrences"""
    # Create labels and formatted target names
    labels = [f"d_{i}" for i in range(len(view_list))]
    fold_names = [f"{name}_FOLD" for name in view_list]

    # Print Legend
    print("Dataset Key:")
    for label, fold_name in zip(labels, fold_names):
        print(f"  {label} -> {fold_name}", len(project.views[fold_name]))
    print()

    # Build set representation for each _FOLD view
    sets = [
        set(project.views[fold_name]["target_id"].dropna())
        for fold_name in fold_names
    ]

    # Compute pairwise overlap matrix
    shared_matrix = pd.DataFrame(
        [[len(s1 & s2) for s2 in sets] for s1 in sets],
        index=labels,
        columns=labels,
    )

    return shared_matrix

def plot_shared_matrix_heatmap(
    shared_matrix, view_key=None, figsize=(7,6), cmap="Blues"
):
    """Plots the co-occurrence matrix as an annotated, hierarchical clustered heatmap.

    Parameters:
    -----------
    shared_matrix : pd.DataFrame
        Symmetric count matrix from get_fold_matrix.
    view_key : dict, optional
        Mapping dictionary { 'df_0': 'dom_0-182_FOLD', ... } for legend logging.
    figsize : tuple, default=(7, 6)
        Figure size dimensions.
    cmap : str, default='Blues'
        Seaborn color palette ('Blues', 'YlGnBu', 'viridis', etc.).
    """
    # 1. Print legend mapping if key provided
    if view_key:
        print("Heatmap Key:")
        for label, full_name in view_key.items():
            print(f"  {label} -> {full_name}")
        print()

    # 2. Render Clustered Heatmap
    g = sns.clustermap(
        shared_matrix,
        annot=True,  # Display raw hit counts in cells
        fmt="d",  # Integer formatting
        cmap=cmap,
        linewidths=0.8,
        linecolor="white",
        figsize=figsize,
        #cbar_kws={"label": "Shared Target Count"},
        tree_kws={"linewidths": 1.2},
    )

    # 3. Styling adjustments
    #g.ax_heatmap.set_title(
    #    "Domain _FOLD Co-occurrence", fontsize=13, pad=15
    #)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=0, fontweight="bold")
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontweight="bold")

    plt.show()
