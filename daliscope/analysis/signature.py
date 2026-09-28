from daliscope.mechanics.metrics import filter_nonredundant_sequences
from daliscope.widgets.factory import widget_view
from daliscope.optics.msa import _get_cached_msa
from daliscope.analysis.wrappers import next_subset_name
from daliscope.optics.plotly import plot_sequence_profile_heatmap
import math
import pandas as pd
import numpy as np

def derive_signature_profile(project, seed_view, nr_threshold=0.40, top_ic=12, verbose=True):
    """
    Locks a profile seed view, reduces to non-redundant subset, derives the
    top-k conserved-position profile (hc_profile), registers it as a model
    with provenance, and produces the profile heatmap + report table.

    Returns
    -------
    seed_df, nr_df, hc_profile, logo, model_name, report
    """
    seed_df = project.views[seed_view]
    nr_df = filter_nonredundant_sequences(seed_df, threshold=nr_threshold)

    logo = _get_cached_msa(nr_df, "sequ_pileup")
    if logo._info_df is None:
        logo._initialize_everything()

    hc_profile = get_hc_profile(
        logo._info_df, top_ic=top_ic, verbose=verbose,
    )

    # Register the profile as a named model, with a fallback-numbered default name
    model_name = next_subset_name(project, seed_view, 'model')
    project.add_model(
        name       = model_name,
        artifact   = hc_profile,
        parent     = seed_view,
        function   = "get_hc_profile",
        parameters = {"nr_threshold": nr_threshold, "top_ic": top_ic},
    )

    # Plot + report
    positions, matrix = plot_sequence_profile_heatmap(hc_profile)
    report = hc_report(hc_profile, project.query_sequence).sort_values(
        "max_value", ascending=False
    )

    return seed_df, nr_df, hc_profile, logo, model_name, report

def get_hc_profile(
    info_df,
    top_ic=12,
    verbose=False
):
    """
    Analyzes an Information-Content DataFrame to provide automated biological insights.
    """
    # Standard biological groupings for amino acids
    polar_charged = set(['R', 'K', 'D', 'E', 'N', 'Q', 'S', 'T', 'Y', 'H', 'C'])
    hydrophobic = set(['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'P', 'G'])

    # 1. Extract fundamental metrics per position (FIXED: peak height is total row sum)
    total_heights = info_df.sum(axis=1)    # Total height of the sequence logo stack
    max_aa = info_df.idxmax(axis=1)     # The most dominant letter at this position

    # Global profile metrics based on stack heights
    global_avg_height = total_heights[total_heights > 0].mean()
    global_max_height = total_heights.max()

    if verbose:
        print(f"Analyzed {len(info_df)} positions.")
        print(f"Global average stack height: {global_avg_height:.2f} bits")
        print(f"Global maximum stack height: {global_max_height:.2f} bits\n")

    # Gather the top-k most conserved positions (by information content) for rule C and D
    hc_positions = info_df.loc[total_heights.sort_values(ascending=False).index[:top_ic]]

    # Gather highly conserved positions for rule C and D
    #highly_conserved_mask = total_heights >= low_conservation_cutoff
    #hc_positions = info_df[highly_conserved_mask]

    return hc_positions

def get_sequence_context(seq, i, window=5):
    """
    Returns the amino acid context around 1-based position i in seq:
    (i-window)..(i-1), space, i, space, (i+1)..(i+window).
    Out-of-bounds positions are padded with spaces.
    """
    def residue_at(pos):
        if 1 <= pos <= len(seq):
            return seq[pos - 1]
        return " "

    left = "".join(residue_at(i - offset) for offset in range(window, 0, -1))
    center = residue_at(i)
    right = "".join(residue_at(i + offset) for offset in range(1, window + 1))

    return f"{left} {center} {right}"


def hc_report(hc_positions, query_sequence):
    # 1. Get the column name of the maximum value for each row
    max_cols = hc_positions.idxmax(axis=1)
    # 2. Get the actual maximum value for each row
    max_vals = hc_positions.max(axis=1)
    # 3. Combine them into a clean report DataFrame
    report_df = pd.DataFrame({
        'max_column': max_cols,
        'max_value': max_vals
    })
    # 4. Add sequence context for each position (index of hc_positions)
    report_df['sequence_context'] = [
        get_sequence_context(query_sequence, i+1)
        for i in report_df.index
    ]
    return report_df


##########


def add_nonredundant_subset(project, parent_view: str, threshold: float = 0.40, subset_name: str = None):
    """
    Creates a non-redundant subset of an existing project view and registers it.

    Parameters
    ----------
    project : the project object (must have .views dict and .add_subset method)
    parent_view : name of the existing view to filter (key into project.views)
    threshold : redundancy threshold passed to filter_nonredundant_sequences
    subset_name : optional custom name for the new subset; defaults to f"{parent_view}_nr"

    Returns
    -------
    The result of project.add_subset(...)
    """
    if not parent_view:
        raise ValueError("parent_view must be a non-empty view name.")
    if parent_view not in project.views:
        raise KeyError(f"'{parent_view}' not found in project.views.")

    parent_df = project.views[parent_view]

    name = subset_name or f"{parent_view}_nr"

    threshold = 0.40

    nr_view_name = f"{parent_view}_nr"
    if nr_view_name in project.views:
        existing_params = project.get_view_parameters(nr_view_name)  # however your project exposes this
        existing_threshold = existing_params.get("threshold")
        if math.isclose(existing_threshold, threshold, abs_tol=0.01):
            print(f"'{nr_view_name}' already exists with threshold={threshold} — reusing it.")
        else:
            print(f"'{nr_view_name}' already exists but was built with threshold={existing_threshold}, "
              f"not the requested {threshold}. Choose a custom name, or delete the existing view first: "
              f"del project.views['{nr_view_name}']")
            # don't proceed — force explicit user decision
            return

    # child is a non-redundant subset of parent
    child_df = filter_nonredundant_sequences(parent_df, threshold=threshold)

    # mask must be built against parent_df, not child_df
    mask = parent_df["target_id"].isin(child_df["target_id"])

    return project.add_subset(
        name       = name,
        mask       = mask,
        parent     = parent_view,
        function   = 'add_nonredundant_subset',
        parameters = {"threshold": threshold},
    )

########################

import numpy as np
import pandas as pd
import plotly.express as px

def rescore_targets_safely(hc_profile: pd.DataFrame, targets_df: pd.DataFrame, positions: list[int], score_col_name='profile_score'):
    """
    Safely rescores target sequences against a specified list of profile positions.
    Calculates the score across 'positions' and the half-max-score threshold.
    
    Returns:
    --------
    final_results : pd.DataFrame
        DataFrame with calculated scores, sorted descending.
    """
    if not positions:
        raise ValueError("positions list cannot be empty.")

    # 1. Isolate metadata columns and drop rows missing sequence data
    meta_cols = [col for col in ['target_id', 'z_score', 'sequence_identity', 'description', 'pfam', 'clan'] if col in targets_df.columns]
    clean_targets = targets_df[meta_cols + ['sequ_pileup']].dropna(subset=['sequ_pileup']).copy()

    # 2. Extract profile matrix for valid amino acid columns and selected positions
    aa_cols = [c for c in hc_profile.columns if len(c) == 1 and c.isalpha()]
    profile_matrix = hc_profile.loc[positions, aa_cols]

    # 4. Safe NumPy vectorization using the selected positions
    sequences_np = np.array([list(seq) for seq in clean_targets['sequ_pileup']], dtype='S1')
    pos_indices = np.array(positions)
    target_aas_matrix = sequences_np[:, pos_indices].astype(str)

    # 5. Map target amino acids to column indices and retrieve scores
    aa_to_idx = {aa: i for i, aa in enumerate(aa_cols)}
    mapping_grid = np.vectorize(lambda x: aa_to_idx.get(x, 0))
    target_aa_indices = mapping_grid(target_aas_matrix)

    scoring_lookup = profile_matrix.to_numpy()
    n_targets, n_pos = target_aa_indices.shape
    scoring_matrix = scoring_lookup[np.arange(n_pos), target_aa_indices]

    # 6. Compute total score per target across selected positions
    target_scores = np.sum(scoring_matrix, axis=1)

    # 7. Reconstruct output DataFrame
    clean_targets[score_col_name] = target_scores

    # Sort by score descending
    final_results = clean_targets.sort_values(by=score_col_name, ascending=False)

    # Clean display columns
    display_cols = [col for col in final_results.columns if col != 'sequ_pileup']

    return final_results[display_cols]

def plot_final_results_histogram_interactive(final_results: pd.DataFrame, trusted_df: pd.DataFrame,
                                               initial_threshold: float, n_positions: int,
                                               score_col_name: str = 'profile_score'):
    """
    Interactive version of plot_final_results_histogram: a slider controls
    score_threshold, redrawing the cutoff line live. Returns get_state, so the
    final chosen threshold can be read after interaction via get_state()['score_threshold'].
    """
    trusted_ids = set(trusted_df['target_id'])
    plot_df = final_results.copy()
    is_positive = plot_df[score_col_name] > 0
    is_trusted = plot_df['target_id'].isin(trusted_ids)

    conditions = [~is_positive, is_positive & is_trusted, is_positive & ~is_trusted]
    choices = ['Zero/Background Noise', 'Trusted Positive Hit', 'NEW Positive Hit']
    plot_df['Screening Group'] = np.select(conditions, choices, default='Unknown')

    score_min = float(plot_df[score_col_name].min())
    score_max = float(plot_df[score_col_name].max())

    def render(score_threshold):
        n_above = int((plot_df[score_col_name] >= score_threshold).sum())

        fig = px.histogram(
            plot_df, x=score_col_name, color='Screening Group', nbins=50, log_y=True,
            title=f'Target Screening Profile Distribution ({n_positions} positions) — '
                  f'{n_above} above threshold',
            labels={score_col_name: f'Cumulative Profile Score ({n_positions} positions)', 'count': 'Log(Count)'},
            color_discrete_map={
                'Zero/Background Noise': '#d62728',
                'Trusted Positive Hit': '#1f77b4',
                'NEW Positive Hit': '#2ca02c'
            },
            barmode='overlay'
        )
        fig.add_vline(x=score_threshold, line_dash="dash", line_color="black",
                       annotation_text="Score Cutoff", annotation_position="top right")
        fig.update_layout(template="plotly_white", hovermode="x unified", legend_title_text="Target Status")
        fig.show()
        return None

    root, get_state, controls = widget_view(
        func=render,
        title="Adjust motif-positive score threshold",
        params={
            "score_threshold": {
                "type": "floatslider",
                "value": initial_threshold,
                "min": score_min,
                "max": score_max,
                "step": (score_max - score_min) / 200,
                "label": "Score threshold",
            },
        },
        show_default_buttons=True,
        return_state=True,
    )
    return get_state

def static_plot_final_results_histogram(final_results: pd.DataFrame, trusted_df: pd.DataFrame,
                                  score_threshold: float, n_positions: int,
                                  score_col_name: str = 'profile_score'):
    """
    Plots a log-scaled histogram with three distinct target groups:
    1. Zero/Background Noise (Score <= 0)
    2. Trusted Positive Hits (Score > 0 and in trusted_df)
    3. NEW Positive Hits (Score > 0 and NOT in trusted_df)
    """
    trusted_ids = set(trusted_df['target_id'])

    plot_df = final_results.copy()
    is_positive = plot_df[score_col_name] >0
    is_trusted = plot_df['target_id'].isin(trusted_ids)

    conditions = [~is_positive, is_positive & is_trusted, is_positive & ~is_trusted]
    choices = ['Zero/Background Noise', 'Trusted Positive Hit', 'NEW Positive Hit']
    plot_df['Screening Group'] = np.select(conditions, choices, default='Unknown')

    fig = px.histogram(
        plot_df, x=score_col_name, color='Screening Group', nbins=50, log_y=True,
        title=f'Target Screening Profile Distribution ({n_positions} positions)',
        labels={score_col_name: f'Cumulative Profile Score ({n_positions} positions)', 'count': 'Log(Count)'},
        color_discrete_map={
            'Zero/Background Noise': '#d62728',
            'Trusted Positive Hit': '#1f77b4',
            'NEW Positive Hit': '#2ca02c'
        },
        barmode='overlay'
    )
    fig.add_vline(x=score_threshold, line_dash="dash", line_color="black",
                   annotation_text="Score Cutoff", annotation_position="top right")
    fig.update_layout(template="plotly_white", hovermode="x unified", legend_title_text="Target Status")
    fig.show()

def run_profile_screening_pipeline(
    hc_profile: pd.DataFrame,
    full_df: pd.DataFrame,
    trusted_df: pd.DataFrame,
    positions: list[int],
    initial_threshold: float,
    n_positions: int = None,
):
    """
    Runs the score-independent part of the pipeline:
    1. Rescores targets safely using selected profile positions.
    2. Shows the interactive histogram with a draggable threshold slider.

    Returns (final_results, get_threshold_state) — call
    summarize_screening_results(final_results, trusted_df, get_threshold_state(), ...)
    after the user has settled on a threshold.
    """
    print("Rescoring targets against profile model...")
    score_col_name = 'profile_score'
    final_results = rescore_targets_safely(hc_profile, full_df, positions, score_col_name=score_col_name)

    get_threshold_state = plot_final_results_histogram_interactive(
        final_results, trusted_df, initial_threshold,
        n_positions or len(positions), score_col_name
    )

    return final_results, get_threshold_state


def summarize_screening_results(
    final_results: pd.DataFrame,
    trusted_df: pd.DataFrame,
    score_threshold: float,
    n_tail: int = 10,
    show_model_set: bool = True,
    score_col_name: str = 'profile_score',
):
    """
    Runs the threshold-dependent part of the pipeline, using whatever
    score_threshold the user settled on via the interactive histogram:
    1. Flags targets against the trusted seed model set.
    2. Slices data into an inspection table around the score threshold.
    3. Prints summary metrics.

    Returns the inspection_df actually displayed (matching_targets can be
    derived as final_results[final_results[score_col_name] >= score_threshold]).
    """
    trusted_ids = set(trusted_df['target_id'])
    df_sorted = final_results.sort_values(score_col_name, ascending=False).copy()

    is_trusted = df_sorted['target_id'].isin(trusted_ids).values
    is_above = (df_sorted[score_col_name] >= score_threshold).values

    conditions = [
        is_above & is_trusted,
        is_above & ~is_trusted,
        ~is_above & is_trusted,
        ~is_above & ~is_trusted
    ]
    choices = [
        'Trusted (Model Set)',
        'NEW HIT',
        'Trusted (Below Threshold!)',
        'NEW (Below Threshold)'
    ]
    df_sorted['status'] = np.select(conditions, choices, default='Unknown')

    above_threshold = df_sorted[is_above]
    below_threshold = df_sorted[~is_above].head(n_tail)

    if not show_model_set:
        above_threshold = above_threshold[above_threshold['status'] != 'Trusted (Model Set)']

    inspection_df = pd.concat([above_threshold, below_threshold])

    all_above = df_sorted[is_above]
    new_discoveries = all_above[all_above['status'] == 'NEW HIT']
    model_recovered = all_above[all_above['status'] == 'Trusted (Model Set)']
    model_dropped = df_sorted[df_sorted['status'] == 'Trusted (Below Threshold!)']

    print(f"\n--- Discovery Summary (Threshold: {score_threshold:.1f})")
    print(f"Total Hits Above Threshold   : {len(all_above)}")
    print(f"  └─ Recovered Model Seeds   : {len(model_recovered)}")
    print(f"  └─ NEW Expanded Targets    : {len(new_discoveries)}  ◄ 🌟 Search Expanded!")
    print(f"Model seeds missed below line: {len(model_dropped)}")
    if not show_model_set:
        print("ℹ️ Note: 'Trusted (Model Set)' entries have been suppressed from the display table below.")

    cols_to_show = ['target_id', 'status', score_col_name] + [
        col for col in inspection_df.columns if col not in ['target_id', 'status', score_col_name]
    ]
    with pd.option_context('display.max_rows', None):
        display(inspection_df[cols_to_show])

    return inspection_df

###########################

def compute_max_possible_score(hc_profile: pd.DataFrame, positions: list) -> float:
    """Maximum possible cumulative profile score over the given positions."""
    aa_cols = [c for c in hc_profile.columns if len(c) == 1 and c.isalpha()]
    return hc_profile.loc[positions, aa_cols].max(axis=1).sum()


def plot_final_results_histogram_interactive(
    final_results: pd.DataFrame,
    trusted_df: pd.DataFrame,
    max_possible_score: float,
    initial_threshold_percent: float,
    n_positions: int,
    score_col_name: str = 'profile_score',
    window: int = 10,
):
    """
    Interactive histogram: a slider (1% steps, 0-100%) sets score_threshold
    as a percentage of max_possible_score. Alongside the histogram, shows a
    BOUNDED inspection table (2*window rows max, regardless of threshold
    position) of the targets closest to the current cutoff — the weakest
    passers and the strongest non-passers — for live functional validation
    (checking descriptions) while adjusting the threshold.

    Returns get_threshold_state, a callable returning
    {'score_threshold': <actual score>} after interaction — matching the
    key Cell III(b) reads.
    """
    trusted_ids = set(trusted_df['target_id'])
    plot_df = final_results.copy()
    is_positive = plot_df[score_col_name] > 0
    is_trusted = plot_df['target_id'].isin(trusted_ids)

    conditions = [~is_positive, is_positive & is_trusted, is_positive & ~is_trusted]
    choices = ['Zero/Background Noise', 'Trusted Positive Hit', 'NEW Positive Hit']
    plot_df['Screening Group'] = np.select(conditions, choices, default='Unknown')

    # Sort once, descending by score — reused for the borderline slice on
    # every redraw rather than re-sorting per threshold change.
    df_sorted = plot_df.sort_values(score_col_name, ascending=False).reset_index(drop=True)

    inspection_cols = ['target_id', score_col_name, 'description'] + [
        c for c in df_sorted.columns
        if c not in ('target_id', score_col_name, 'description', 'Screening Group')
    ]

    def pct_to_score(pct):
        return (pct / 100.0) * max_possible_score

    initial_pct = min(max(initial_threshold_percent * 100.0, 0), 100)

    def render(threshold_pct, table_window=window):
        score_threshold = pct_to_score(threshold_pct)
        is_above = df_sorted[score_col_name] >= score_threshold
        n_above = int(is_above.sum())

        # --- Histogram ---
        fig = px.histogram(
            plot_df, x=score_col_name, color='Screening Group', nbins=50, log_y=True,
            title=f'Target Screening Profile Distribution ({n_positions} positions) — '
                  f'{n_above} above threshold ({threshold_pct:.0f}% of max)',
            labels={score_col_name: f'Cumulative Profile Score ({n_positions} positions)', 'count': 'Log(Count)'},
            color_discrete_map={
                'Zero/Background Noise': '#d62728',
                'Trusted Positive Hit': '#1f77b4',
                'NEW Positive Hit': '#2ca02c'
            },
            barmode='overlay'
        )
        fig.add_vline(x=score_threshold, line_dash="dash", line_color="black",
                       annotation_text=f"Score Cutoff ({threshold_pct:.0f}% of max)", annotation_position="top right")
        fig.update_layout(template="plotly_white", hovermode="x unified", legend_title_text="Target Status")
        fig.show()

        # --- Bounded borderline table: closest above + closest below threshold ---
        # is_above/df_sorted share the same (score-descending) row order, so
        # the LAST True rows are the weakest passers, and the FIRST False
        # rows are the strongest non-passers.
        closest_above = df_sorted[is_above].tail(table_window)
        closest_below = df_sorted[~is_above].head(table_window)

        print(f"\n{n_above} targets above threshold. Showing the {len(closest_above)} weakest "
              f"passers and {len(closest_below)} strongest non-passers for inspection:")
        with pd.option_context('display.max_rows', None, 'display.max_colwidth', 60):
            display(pd.concat([
                closest_above.assign(status='ABOVE (borderline)'),
                closest_below.assign(status='BELOW (borderline)'),
            ])[['status'] + inspection_cols])

        return None

    root, get_state, controls = widget_view(
        func=render,
        title="Adjust motif-positive score threshold (% of max possible score)",
        params={
            "threshold_pct": {
                "type": "floatslider",
                "value": initial_pct,
                "min": 0,
                "max": 100,
                "step": 1,
                "label": "Threshold (% of max)",
            },
        },
        show_default_buttons=True,
        return_state=True,
    )

    def get_threshold_state():
        """Returns {'score_threshold': <actual score>} for the slider's current position."""
        return {"score_threshold": pct_to_score(get_state()["threshold_pct"])}

    return get_threshold_state

####

def run_profile_screening_pipeline(
    hc_profile: pd.DataFrame,
    full_df: pd.DataFrame,
    trusted_df: pd.DataFrame,
    positions: list,
    initial_threshold: float,
    window: int = 10,
    score_col_name: str = 'profile_score',
):
    """
    Runs the score-independent part of the pipeline:
    1. Rescores targets safely using selected profile positions.
    2. Shows the interactive histogram + bounded borderline inspection table.

    initial_threshold : an actual score value (not a percent) seeding the
        slider's starting position — converted internally to a percentage
        of max_possible_score.
    window : Power User control — number of rows shown above AND below the
        threshold in the live inspection table (table is always <= 2*window
        rows, regardless of how many targets actually pass).

    Returns (final_results, get_threshold_state) — call
    summarize_screening_results(final_results, trusted_df, get_threshold_state()['score_threshold'], ...)
    after the user has settled on a threshold.
    """
    print("Rescoring targets against profile model...")
    final_results = rescore_targets_safely(hc_profile, full_df, positions, score_col_name=score_col_name)

    max_possible_score = compute_max_possible_score(hc_profile, positions)
    initial_threshold_percent = 0.0 if max_possible_score == 0 else initial_threshold / max_possible_score

    get_threshold_state = plot_final_results_histogram_interactive(
        final_results, trusted_df, max_possible_score,
        initial_threshold_percent=initial_threshold_percent,
        n_positions=len(positions),
        score_col_name=score_col_name,
        window=window,
    )

    return final_results, get_threshold_state

def compute_score_threshold(hc_profile: pd.DataFrame, positions: list, threshold_percent: float) -> float:
    """Cumulative profile score cutoff = threshold_percent * max possible score over the given positions."""
    aa_cols = [c for c in hc_profile.columns if len(c) == 1 and c.isalpha()]
    pos_maxes = hc_profile.loc[positions, aa_cols].max(axis=1)
    max_possible_score = pos_maxes.sum()
    return max_possible_score * threshold_percent
