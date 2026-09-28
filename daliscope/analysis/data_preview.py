import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from sklearn.svm import SVC
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import seaborn as sns

import daliscope
import daliscope.analysis.signature
import daliscope.analysis.domain_comparison
import daliscope.widgets.factory
from IPython.display import display, HTML
from html import escape

def show_heading(title, level=3):
    display(HTML(f"<h{level}>{escape(str(title))}</h{level}>"))

# Section 1 helper functions
def _pretty_top_hit_table(df, k=3):
    # Create a mapping for shorter column names
    short_cols = {
        'target_id': 'Target',
        'z_score': 'Z-Score',
        'rmsd': 'RMSD',
        'sequence_identity': '% Id',
        'description': 'Description',
        'alignment_length': 'aln_len',
        'seq_id_pct': '% Id',
    }

    # 2. Select columns, scale sequence identity, sort, and format
    df_view = (
        df[
            [
                'target_id',
                'z_score',
                'rmsd',
                'alignment_length',
                'sequence_identity',
                'description',
            ]
        ]
        .assign(seq_id_pct=lambda x: x['sequence_identity'] * 100)
        .drop(columns=['sequence_identity'])
        .sort_values('z_score', ascending=False)
        .head(k)
        .rename(columns=short_cols)
    )

    # 3. Format floats: 2 decimals for RMSD, 0 decimals (integer) for % Id
    display(
        df_view.style.format(
            {
                'RMSD': '{:.1f}', 
                'aln_len': '{:.0f}',
                'Z-Score': '{:.1f}',
                '% Id': '{:.0f}',  # Displays as an integer percentage (e.g., 85%)
            }
        )
    )

def load_and_inspect_project(pack_path):
    """Initializes project and displays top hit summary."""
    show_heading("Loading data")
    project, full_df, full_view = daliscope.core.project.initialize_project(
        pack_path
    )

    show_heading("Top hits")
    _pretty_top_hit_table(full_df)

    return project, full_df, full_view


def build_and_register_subsets(project, full_df, full_view):
    """Fits Dali boundaries, orders filtered subsets by size, registers them

    to the project, and returns the boundary results.
    """
    show_heading("Creating filtered subsets")
    _, res_list = daliscope.analysis.data_preview.fit_robust_dali_boundary(
        full_df
    )

    # Extract hits & masks
    df_0, mask_0 = (
        daliscope.analysis.data_preview.extract_selected_hits_and_mask(
            full_df, res_list[0][0], res_list[0][1]
        )
    )
    df_1, mask_1 = (
        daliscope.analysis.data_preview.extract_selected_hits_and_mask(
            full_df, res_list[1][0], res_list[1][1]
        )
    )

    # Ensure FILTERED_0 is always the larger subset
    is_0_larger = len(df_0) >= len(df_1)
    name_0 = "FILTERED_0" if is_0_larger else "FILTERED_1"
    name_1 = "FILTERED_1" if is_0_larger else "FILTERED_0"

    project.add_subset(
        name_0, mask_0, full_view, "fit_robust_dali_boundary", {"res": res_list[0]}
    )
    project.add_subset(
        name_1, mask_1, full_view, "fit_robust_dali_boundary", {"res": res_list[1]}
    )

    return res_list

def fit_robust_dali_boundary(
    df, z_col="z_score", cov_col="query_coverage", seq_col=None
):
    """Fits robust bisection boundaries based on low-coverage KDE valleys.

    Guarantees two classes for SVM fitting even if data is homogeneous.
    """
    work_df = df.copy()

    # 1. Non-redundant filtering (if sequence identity available)
    if seq_col and seq_col in work_df.columns:
        work_df = work_df[work_df[seq_col] < 0.90]

    work_df["log_z"] = np.log1p(work_df[z_col])

    # 2. Clip upper extreme outliers for stable KDE density estimation
    p95_z = np.percentile(work_df["log_z"], 95)

    # 3. Separate low-coverage background region
    low_cov_mask = work_df[cov_col] < 0.5
    low_cov_z = work_df.loc[low_cov_mask, "log_z"]

    valley_cutoffs = []

    if len(low_cov_z) > 20 and low_cov_z.nunique() > 1:
        kde = gaussian_kde(low_cov_z)
        z_grid = np.linspace(low_cov_z.min(), p95_z, 300)
        densities = kde(z_grid)

        # Invert densities to find local minima (valleys)
        inv_densities = -densities
        valleys_idx, _ = find_peaks(inv_densities)

        # Filter valleys after the main background noise peak
        main_peak_idx = np.argmax(densities)
        valid_valleys = [v for v in valleys_idx if v > main_peak_idx]

        for v_idx in valid_valleys[:2]:
            valley_cutoffs.append(z_grid[v_idx])

    # Fallback cutoffs if KDE fails or yields < 2 valleys
    if len(valley_cutoffs) == 0:
        valley_cutoffs = [
            np.percentile(work_df["log_z"], 25),
            np.percentile(work_df["log_z"], 75),
        ]
    elif len(valley_cutoffs) == 1:
        valley_cutoffs.append(valley_cutoffs[0] + 0.5)

    results = []
    res_list = []

    for idx, cutoff_log_z in enumerate(valley_cutoffs):
        cov_threshold = 0.5 if idx == 0 else 0.75

        # Initial signal labeling
        is_signal = (
            (work_df["log_z"] > cutoff_log_z)
            | (work_df[cov_col] > cov_threshold)
        ).astype(int)

        # --- SAFETY CHECK: Ensure exactly 2 classes (0 and 1) ---
        unique_classes = np.unique(is_signal)

        if len(unique_classes) < 2:
            # Fallback to median split on log_z if labeling collapsed to 1 class
            median_z = work_df["log_z"].median()
            is_signal = (work_df["log_z"] >= median_z).astype(int)

            # If median split still yields 1 class (e.g. all values identical), split by index
            if len(np.unique(is_signal)) < 2:
                is_signal = np.zeros(len(work_df), dtype=int)
                is_signal[: len(work_df) // 2] = 1

        X = work_df[["log_z", cov_col]]
        clf = SVC(kernel="linear", C=1.0)
        clf.fit(X, is_signal)

        w = clf.coef_[0]
        b = clf.intercept_[0]

        results.append(
            {
                "valley_log_z": cutoff_log_z,
                "clf": clf,
                "weights": (w[0], w[1]),
                "intercept": b,
            }
        )

        res_list.append([(w[0], w[1]), b])

    return results, res_list

def extract_selected_hits_and_mask(
    df, weights, intercept, z_col="z_score", cov_col="query_coverage", verbose=False
):
    """Evaluates the SVM decision boundary against the input dataframe and returns

    both the filtered DataFrame and a boolean pd.Series mask aligned with df.
    """
    w1, w2 = weights
    b = intercept

    # 1. Compute decision score directly aligned with input index
    decision_score = w1 * np.log1p(df[z_col]) + w2 * df[cov_col] + b

    # 2. Boolean mask aligned with df.index (True for selected hits)
    mask = decision_score > 0

    # 3. Create filtered DataFrame
    df_classified = df.copy()
    df_classified["boundary_score"] = decision_score
    df_classified["is_selected"] = mask

    selected_hits_df = (
        df_classified[mask]
        .sort_values(by=z_col, ascending=False)
        .reset_index(drop=True)
    )

    # 4. Summary metrics
    total_hits = len(df)
    selected_count = mask.sum()
    if verbose:
      print(
        f"Selected {selected_count:,} out of {total_hits:,} total hits ({selected_count / total_hits:.1%})"
      )

    # Return both the filtered DataFrame and the index-aligned boolean mask
    return selected_hits_df, mask

def plot_dali_hits_with_boundaries(df, res_list, z_col="z_score", cov_col="query_coverage"):
    _df = df[[z_col,cov_col,'sequence_identity']].dropna()
    x = _df[z_col]
    y = _df[cov_col]
    seq = _df["sequence_identity"]

    # 1. SET UP THE GRID DESIGN
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.toolbar_visible = False

    gs = fig.add_gridspec(
        2, 2,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        wspace=0.05,
        hspace=0.05
    )

    # Define the 3 axes (and lock their shared scaling)
    ax = fig.add_subplot(gs[1, 0])                     # Bottom-Left (Main Scatter)
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)    # Top (Marginal X KDE)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)    # Right (Marginal Y KDE)

    # 2. BACKGROUND REGIONS (z-score bands on Main Plot)
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    # We clip the lower limit of x at 2 to match the final scatter plot set_xlim
    plot_x_min = 2

    ax.axvspan(0, 4,   color="lightgray", alpha=0.3, label="likely unrelated")
    ax.axvspan(4, 8,   color="orange",    alpha=0.2, label="marginal")
    ax.axvspan(8, 12,  color="yellow",    alpha=0.2, label="likely fold")
    ax.axvspan(12, 16, color="green",  alpha=0.15, label="likely superfamily")
    ax.axvspan(16, x_max, color="blue", alpha=0.1, label="likely family")

    # 3. CATEGORICAL MARKERS (Sequence Identity on Main Plot)
    homolog = seq >= 0.3
    twilight = (seq > 0.15) & (seq < 0.3)
    dark = seq <= 0.15

    ax.scatter(
        x[homolog], y[homolog],
        c="red", marker="o", s=40,
        label="homolog (seq_id ≥ 0.3)", edgecolor="red"
    )

    ax.scatter(
        x[twilight], y[twilight],
        c="cyan", marker="s", s=40,
        label="twilight (seq_id 0.15–0.3)", edgecolor="lightblue"
    )

    ax.scatter(
        x[dark], y[dark],
        c="darkgrey", marker="^", s=40,
        label="dark (seq_id ≤ 0.15)", alpha=0.7
    )

    # 4. COMPUTE AND PLOT KDE CURVES
    # Smooth Top KDE (z_score)
    x_eval = np.linspace(plot_x_min, x_max, 300)
    kde_x = gaussian_kde(x)
    ax_histx.plot(x_eval, kde_x(x_eval), color="slategrey", lw=1.5)
    ax_histx.fill_between(x_eval, 0, kde_x(x_eval), color="slategrey", alpha=0.3)

    # Smooth Right KDE (query_coverage) - Rotated 90 degrees
    y_eval = np.linspace(y_min, y_max, 300)
    kde_y = gaussian_kde(y)
    ax_histy.plot(kde_y(y_eval), y_eval, color="slategrey", lw=1.5)
    # fill_betweenx fills horizontally from x=0 to the KDE curve values across the y_eval range
    ax_histy.fill_betweenx(y_eval, 0, kde_y(y_eval), color="slategrey", alpha=0.3)

    # 5. CLEAN UP MARGINAL AXES (and control height scaling)
    ax_histx.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_histx.tick_params(axis="y", left=False, labelleft=False)
    ax_histx.set_ylim(bottom=0)

    # --- ADD THIS TO SQUISH THE TOP KDE VERTICALLY ---
    # e.g., Make the limit 1.5x larger than the peak, pushing the peak down
    ax_histx.set_ylim(top=max(kde_x(x_eval)) * 3.0)

    for spine in ["top", "left", "right"]:
        ax_histx.spines[spine].set_visible(False)

    ax_histy.tick_params(axis="y", labelleft=False, left=False)
    ax_histy.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_histy.set_xlim(left=0)

    # --- ADD THIS TO SQUISH THE RIGHT KDE HORIZONTALLY ---
    # e.g., Make the limit 1.5x larger than the peak, pushing the peak to the left
    ax_histy.set_xlim(right=max(kde_y(y_eval)) * 3.0)

    for spine in ["top", "bottom", "right"]:
        ax_histy.spines[spine].set_visible(False)

    # 1. Extract raw linear Z-scores and coverage
    z_vals = df[z_col]
    cov = df[cov_col]

    # 4. Generate curve evaluation points on linear Z domain
    z_min = max(0, z_vals.min())
    x_line = np.linspace(z_min, z_vals.max(), 300)
    y_max = max(1.05, cov.max() * 1.05 if cov.max() > 1.0 else 1.05)

    # Palette for multiple curved boundaries
    line_colors = [
        "#d9534f",
        "#0275d8",
        "#5cb85c",
        "#f0ad4e",
        "#6f42c1",
        "#17a2b8",
    ]

    # 5. Plot each boundary curve in res_list
    for idx, res in enumerate(res_list):
        w1, w2 = res[0][0], res[0][1]
        b = res[1]

        # Calculate curve coordinates in linear Z space
        # Log-space line: w1 * ln(1 + Z) + w2 * Cov + b = 0
        # Solved for Coverage: Cov = (-w1 * ln(1 + Z) - b) / w2
        y_line = (-w1 * np.log1p(x_line) - b) / w2

        color = line_colors[idx % len(line_colors)]
        label = f"Cutoff {idx}: {w1:.2f} ln(1+Z) + {w2:.2f} Cov = {-b:.2f}"

        ax.plot(
            x_line,
            y_line,
            color=color,
            linestyle="--",
            linewidth=2.2,
            zorder=4 + idx,
            label=label,
        )

    # 6. AXES LIMITS + MAIN FORMATTING
    ax.set_xlabel("z_score")
    ax.set_ylabel("query_coverage")

    ax_histx.set_title("Naive hit classification (heuristic thresholds)", fontweight='bold', pad=10)

    ax.set_xlim(plot_x_min, x_max)
    ax.set_ylim(y_min, y_max)

    #ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        fontsize=9,
        framealpha=0.9,
    )

    fig.subplots_adjust(bottom=0.30)

    display(fig)
    plt.close(fig)

    return plt

def get_hc_profile(seed_df, nr_threshold=0.4):
    nr_df = daliscope.analysis.signature.filter_nonredundant_sequences(seed_df, threshold=nr_threshold)
    logo = daliscope.analysis.signature._get_cached_msa(nr_df, "sequ_pileup")
    if logo._info_df is None:
        logo._initialize_everything()

    # 1. Extract fundamental metrics per position (FIXED: peak height is total row sum)
    total_heights = logo._info_df.sum(axis=1)    # Total height of the sequence logo stack
    max_aa = logo._info_df.idxmax(axis=1)     # The most dominant letter at this position

    # Global profile metrics based on stack heights
    global_avg_height = total_heights[total_heights > 0].mean()
    global_max_height = total_heights.max()

    return global_max_height, global_avg_height, len(nr_df)

def generate_summary_table(
    df_map: dict, seq_id_col: str = "sequence_identity"
) -> pd.DataFrame:
    summary_data = {}

    for name, df in df_map.items():
        max_ic, median_ic, nr_size =  get_hc_profile(df)

        # Standardize identity scaling (0-1 vs 0-100%)
        seq_id = df[seq_id_col]
        if seq_id.max() <= 1.0:
            seq_id = seq_id * 100.0

        signal_strength = (
            (max_ic - median_ic)
            if (pd.notnull(max_ic) and pd.notnull(median_ic) and median_ic != 0)
            else np.nan
        )

        # Count Pfam annotations
        pfam_count = 0
        if "pfam" in df.columns:
            pfam_count = (
                df["pfam"] 
                .ne("Unassigned")
                .sum()
            )

        # Format values as strict strings to prevent pandas type coercion
        total = len(df)
        n35 = int((seq_id > 35.0).sum())
        n15 = int(((seq_id >= 15.0) & (seq_id <= 35.0)).sum())
        n0 = int((seq_id < 15.0).sum())
        npfam = int(pfam_count)
        summary_data[name] = {
            "(a) Total number of targets": f"{total:,d}",
            "(a.1) Non-redundant targets": f"{nr_size:,d}",
            "(a.1) percentage": f"{nr_size/total*100:.0f}%",
            "(a.2) Targets in >35% identity regime": f"{n35:,d}",
            "(a.2) percentage": f"{n35/total*100:.0f}%",            
            "(a.3) Targets in 15-35% identity regime": f"{n15:,d}",
            "(a.3) percentage": f"{n15/total*100:.0f}%",            
            "(a.4) Targets in <15% identity regime": f"{n0:,d}",
            "(a.4) percentage": f"{n0/total*100:.0f}%",            
            "(b) Alignments with Pfam annotation": f"{npfam:,d}",
            "(b) percentage": f"{npfam/total*100:.0f}%",            
            "(c) Signature strength (c.1 - c.2)": f"{signal_strength:.2f}",
            "(c.1) Maximum IC peak": f"{max_ic:.2f}",
            "(c.2) Average IC value": f"{median_ic:.2f}",
        }

    return pd.DataFrame(summary_data)

def apply_summary_highlights(df: pd.DataFrame):
    """Applies conditional cell highlighting to summary_df based on string values

    and propagates colors to corresponding target/IC parent rows in the same column.
    """
    GREEN = "background-color: #d4edda; color: #155724; font-weight: bold;"
    YELLOW = "background-color: #fff3cd; color: #856404; font-weight: bold;"

    # Initialize empty style matrix matching df dimensions
    styles = pd.DataFrame("", index=df.index, columns=df.columns)

    def parse_num(val):
        try:
            return float(str(val).replace("%", "").replace(",", "").strip())
        except (ValueError, TypeError):
            return None

    # Maps row index substrings to their parent/associated row index substrings
    parent_map = {
        "(b) percentage": "(b) Alignments with Pfam annotation",
        "(a.3) percentage": "(a.3) Targets in 15-35% identity regime",
        "(a.4) percentage": "(a.4) Targets in <15% identity regime",
    }

    # Helper to find exact matching row index by substring
    def find_row(pattern):
        for idx in df.index:
            if pattern in str(idx):
                return idx
        return None

    for col in df.columns:
        # ---------------------------------------------------------------------
        # Rule 1: Signature strength (c) -> propagates to (c.1) and (c.2)
        # ---------------------------------------------------------------------
        row_c = find_row("(c) Signature strength")
        if row_c is not None:
            val_c = parse_num(df.loc[row_c, col])
            if val_c is not None:
                color_c = None
                if val_c > 3.0:
                    color_c = GREEN
                elif val_c > 2.0:
                    color_c = YELLOW

                if color_c:
                    # Highlight (c), (c.1), and (c.2)
                    styles.loc[row_c, col] = color_c

                    #row_c1 = find_row("(c.1) Maximum IC peak")
                    #if row_c1 is not None:
                    #    styles.loc[row_c1, col] = color_c

                    #row_c2 = find_row("(c.2) Average IC value")
                    #if row_c2 is not None:
                    #    styles.loc[row_c2, col] = color_c

        # ---------------------------------------------------------------------
        # Rules 2-4: Percentage rows -> propagate to the row directly above
        # ---------------------------------------------------------------------
        for pct_pattern, parent_pattern in parent_map.items():
            row_pct = find_row(pct_pattern)
            if row_pct is not None:
                val_pct = parse_num(df.loc[row_pct, col])
                if val_pct is not None:
                    color_pct = None

                    # Rule 2: Pfam percentage
                    if pct_pattern == "(b) percentage":
                        if val_pct < 25.0:
                            color_pct = GREEN
                        elif val_pct < 50.0:
                            color_pct = YELLOW

                    # Rule 3: Twilight regime percentage
                    elif pct_pattern == "(a.3) percentage":
                        if val_pct > 50.0:
                            color_pct = GREEN

                    # Rule 4: Dark regime percentage
                    elif pct_pattern == "(a.4) percentage":
                        if val_pct > 50.0:
                            color_pct = YELLOW

                    if color_pct:
                        # Highlight the percentage cell
                        styles.loc[row_pct, col] = color_pct

                        # Highlight the parent count row directly above it
                        row_parent = find_row(parent_pattern)
                        if row_parent is not None:
                            styles.loc[row_parent, col] = color_pct

    return df.style.apply(lambda _: styles, axis=None)

def composite_domain_coverage(df, query_length):
    "plot side-by-side msa-heatmap and coverage contingency table"
    res = daliscope.analysis.domain_comparison.run_domain_comparison_wizard(
        df, query_length, query_length
    )

    contingency_matrix = np.array(res["contingency"]) * 100

    # 1. Build visualization layout
    fig, (ax_kde, ax_heat) = plt.subplots(
        1, 2, figsize=(12, 4), dpi=80, gridspec_kw={"width_ratios": [6, 3]}
    )

    # 2. Render MSA RGB array directly into ax_kde
    msa = daliscope.widgets.factory._get_cached_msa(df, "dssp_pileup")

    # Trigger lazy initialization of _rgb_array if not already computed
    if msa._rgb_array is None:
        lookup = msa.config.get("rgb", {})
        matrix = msa.df[msa.pileup_col].astype(str).tolist()
        msa._rgb_array = np.array(
            [[lookup.get(char, (255, 255, 255)) for char in row] for row in matrix],
            dtype=np.uint8,
        )

    ax_kde.imshow(msa._rgb_array, aspect="auto")
    ax_kde.set_title("DSSP Pileup Alignment (E red, H blue, L green)")
    ax_kde.set_xlabel("Query Position")
    ax_kde.set_ylabel("Target Index")

    # 3. Contingency Heatmap on the right
    sns.heatmap(
        contingency_matrix,
        rasterized=True,
        annot=True,
        fmt=".1f",
        cmap="Blues",
        cbar_kws={"label": "Percentage of Total Hits (%)"},
        ax=ax_heat,
        annot_kws={"size": 12, "weight": "bold"},
        xticklabels=["t_single", "t_multi"],
        yticklabels=["q_single", "q_multi"],
    )

    ax_heat.set_title("Coverage < 0.5 or > 0.5")
    ax_heat.set_ylabel("Query Coverage")
    ax_heat.set_xlabel("Target Coverage")

    plt.tight_layout()
    plt.show()
