import matplotlib.pyplot as plt
import seaborn as sns
from daliscope.widgets.factory import launch_interactive_view, with_refresh
from daliscope.mechanics.metrics import get_domain_ranges
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import numpy as np

def select_representatives(
    original_descriptions: list[str],
    tfidf_scores: dict[str, float],
    top_k: int = 20,
    processed_descriptions: list[str] = None,
    decay_factor: float = 0.0,  # 0.0 = zero-out covered n-grams; 0.5 = halve their weight
) -> list[tuple[str, float]]:
    """
    Greedily selects top representative descriptions using submodular coverage decay.
    Prevents repeated selection of the same functional keyphrase.
    """
    if processed_descriptions is None:
        processed_descriptions = original_descriptions

    # Make a mutable copy of n-gram weights so we don't mutate external state
    current_tfidf = dict(tfidf_scores)

    selected_reps = []
    selected_indices = set()

    for _ in range(min(top_k, len(original_descriptions))):
        # 1. Rescore remaining candidate descriptions with current (damped) TF-IDF weights
        candidates = []
        for idx, (orig, proc) in enumerate(zip(original_descriptions, processed_descriptions)):
            if idx in selected_indices:
                continue

            # Score using current damped TF-IDF
            score = _score_single_description(proc, current_tfidf)
            candidates.append((score, orig, proc, idx))

        if not candidates:
            break

        # 2. Sort deterministically by (score, original_text)
        candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
        best_score, best_orig, best_proc, best_idx = candidates[0]

        if best_score <= 0.0:
            break  # Stop if no remaining informative descriptions

        selected_reps.append((best_orig, best_score))
        selected_indices.add(best_idx)

        # 3. DAMPING STEP: Decay weights of all n-grams covered by the chosen description
        formatted_proc = " " + best_proc.strip()
        for ngram in list(current_tfidf.keys()):
            if current_tfidf[ngram] > 0 and ngram in formatted_proc:
                current_tfidf[ngram] *= decay_factor

    return selected_reps


def _score_single_description(processed_desc: str, tfidf_dict: dict[str, float]) -> float:
    """Helper to compute mean max-TF-IDF for a single processed description."""
    formatted_desc = " " + processed_desc.strip()
    words = formatted_desc.split()
    if not words:
        return 0.0

    char_scores = np.zeros(len(formatted_desc))

    for ngram, val in tfidf_dict.items():
        if val <= 0:
            continue
        start = 0
        while True:
            idx = formatted_desc.find(ngram, start)
            if idx == -1:
                break
            char_scores[idx : idx + len(ngram)] = np.maximum(
                char_scores[idx : idx + len(ngram)], val
            )
            start = idx + 1

    word_scores = []
    curr_pos = 0
    for word in words:
        w_idx = formatted_desc.find(word, curr_pos)
        w_len = len(word)
        if w_idx != -1:
            word_scores.append(char_scores[w_idx : w_idx + w_len].max())
            curr_pos = w_idx + w_len

    return float(np.mean(word_scores)) if word_scores else 0.0

def get_query_domain_length(domain_ranges: list) -> int:
    """
    Compute the total query length covered by a list of domain ranges.
    Sum of (end - start) for each non-overlapping range.

    Parameters
    ----------
    domain_ranges : list of (start, end) tuples, 0-based half-open

    Returns
    -------
    int — total number of query positions covered
    """
    return sum(end - start for start, end in domain_ranges)

def run_domain_comparison_wizard(df, query_domain_size, query_length, z_score_low_threshold=2.0):
    """
    Analyzes and compares query vs target coverage using KDE peaks 
    and a 2x2 structural contingency table. No plotting, pure data analysis.
    """
    cols = ['query_coverage', 'target_coverage', 'alignment_length', 'z_score']
    clean_df = df[cols].dropna()
    clean_df = clean_df[clean_df['z_score'] >= z_score_low_threshold]

    if len(clean_df) == 0:
        return {
            "clean_df_size": 0, "q_peak": 0, "t_peak": 0, "median_length": 0,
            "peak_inference": "⚠️ No valid coverage data available for this view.",
            "contingency": pd.DataFrame(), "table_inference": "No data.",
            "eval_range": np.linspace(0, 1.0, 500), "q_densities": np.zeros(500), "t_densities": np.zeros(500)
        }

    q_cov = clean_df['query_coverage'] / query_domain_size * query_length
    t_cov = clean_df['target_coverage']
    lengths = clean_df['alignment_length']

    # Eval KDE over the standard 0 to 1 coverage range
    eval_range = np.linspace(0, 1.0, 500)

    # Catch edge case where variance is 0 (e.g., identical values)
    try:
        kde_q = gaussian_kde(q_cov)
        q_densities = kde_q(eval_range)
    except Exception:
        q_densities = np.zeros_like(eval_range)

    try:
        kde_t = gaussian_kde(t_cov)
        t_densities = kde_t(eval_range)
    except Exception:
        t_densities = np.zeros_like(eval_range)

    # Locate major peak x-coordinates
    q_peak = eval_range[np.argmax(q_densities)] if np.any(q_densities) else 0.0
    t_peak = eval_range[np.argmax(t_densities)] if np.any(t_densities) else 0.0
    peak_diff = t_peak - q_peak

    # Calculate alignment length metrics
    median_length = int(lengths.median())

    if median_length < 50:
        length_context = f"definite fragments (median length: {median_length} aa, which is < 50 aa)"
    elif median_length > 100:
        length_context = f"definitely valid domains (median length: {median_length} aa, which is > 100 aa)"
    else:
        length_context = f"borderline/intermediate domains (median length: {median_length} aa)"

    # Peak Inference Rules
    if abs(peak_diff) < 0.2 and q_peak > 0.5 and t_peak > 0.5:
        peak_inference = (
            "✅ BOTH PEAKS HIGH & CLOSE:\n"
            f"   Query peak ({q_peak:.2f}) and Target peak ({t_peak:.2f}) are aligned and high.\n"
            "   Inference: The query and most targets are likely a single, mutually matched domain."
        )
    elif peak_diff < -0.4:
        peak_inference = (
            "⚠️ TARGET PEAK SIGNIFICANTLY SMALLER THAN QUERY PEAK:\n"
            f"   Target peak ({t_peak:.2f}) is much lower than Query peak ({q_peak:.2f}) (diff: {peak_diff:.2f}).\n"
            "   Inference: Target structures are likely multidomain or contain fragmentary matches.\n"
            f"   (Supporting context: Alignments are characterized as {length_context})."
        )
    else:
        peak_inference = (
            "ℹ️ MIXED PEAK DISTRIBUTION:\n"
            f"   Query peak: {q_peak:.2f}, Target peak: {t_peak:.2f} (diff: {peak_diff:.2f}).\n"
            "   No dominant single-domain or strictly fractional target pattern observed."
        )

    # Contingency Table
    q_high = 'q_single' # 'q_cov > 0.5'
    q_low = 'q_multi' # 'q_cov ≤ 0.5'
    t_high = 't_single' # 't_cov > 0.5'
    t_low = 't_multi' # 't_cov ≤ 0.5'
    q_cat = np.where(q_cov > 0.5, q_high, q_low)
    t_cat = np.where(t_cov > 0.5, t_high, t_low)

    contingency = pd.crosstab(index=q_cat, columns=t_cat, normalize='all')

    for idx in [q_low, q_high]:
        if idx not in contingency.index:
            contingency.loc[idx] = 0.0
    for col in [t_low, t_high]:
        if col not in contingency.columns:
            contingency[col] = 0.0

    contingency = contingency.reindex(
        index=[q_high, q_low],
        columns=[t_low, t_high]
    )

    p_high_q_high_t = contingency.loc[q_high, t_high]
    p_high_q_low_t  = contingency.loc[q_high, t_low]
    p_low_q_high_t  = contingency.loc[q_low, t_high]
    p_low_q_low_t   = contingency.loc[q_low, t_low]

    mass_threshold = 0.25
    table_inference = "📊 CONTINGENCY MASS INFERENCE:\n"

    if p_high_q_high_t > mass_threshold:
        table_inference += (
            f"   • Major mass ({p_high_q_high_t*100:.1f}%) in [High Query, High Target].\n"
            "     Inference: Confirmed single-domain structural match on both sides."
        )
    elif p_high_q_low_t > mass_threshold:
        table_inference += (
            f"   • Major mass ({p_high_q_low_t*100:.1f}%) in [High Query, Low Target].\n"
            "     Inference: Multidomain targets. The query is fully covered, but target proteins have extra domain architectures."
        )
    elif p_low_q_high_t > mass_threshold:
        table_inference += (
            f"   • Major mass ({p_low_q_high_t*100:.1f}%) in [Low Query, High Target].\n"
            "     Inference: Multidomain query. The target is fully covered, but the query protein has extra domain architectures."
        )
    elif p_low_q_low_t > mass_threshold:
        table_inference += (
            f"   • Major mass ({p_low_q_low_t*100:.1f}%) in [Low Query, Low Target].\n"
            "     Inference: Fragmentary matches. Both query and targets share only local structural motifs."
        )
    else:
        table_inference += "   • Mass is distributed across multiple quadrants. Structurally heterogeneous population."

    return {
        "clean_df_size": len(clean_df),
        "q_peak": q_peak,
        "t_peak": t_peak,
        "median_length": median_length,
        "peak_inference": peak_inference,
        "contingency": contingency,
        "table_inference": table_inference,
        "eval_range": eval_range,
        "q_densities": q_densities,
        "t_densities": t_densities
    }


def launch_domain_wizard_widget(project, title="Domain Coverage Profile Wizard"):
    """
    Launches the combined Domain Comparison widget with view-refresh dropdowns,
    z-score adjustment, and bookmarks.
    """
    wizard_params = {
        'view_name': {
            'type':    'dropdown',
            'options': list(project.views.keys()),
            'value':   list(project.views.keys())[0],
            'label':   'View Selection',
        },
        'z_score_low_threshold': {
            # Corrected parameter type to match your widget_view mappings:
            'type':    'floatslider',  # Or 'floattext' if you want direct input
            'value':   2.0,
            'min':     2.0,
            'max':     32.0,
            'step':    0.5,
            'label':   'Z-Score Min',
        }
    }

    launch_interactive_view(
        func=make_domain_comparison_fn(project),
        title=title,
        params=wizard_params,
        project=project,
        control_wrappers={
            'view_name': with_refresh(lambda: list(project.views.keys())),
        }
    )


def plot_domain_comparison(df, query_domain_size, query_length, z_score_low_threshold=2.0):
    """
    Extracts analytics using run_domain_comparison_wizard, generates the 1x2 Matplotlib
    figure, forces it to render first, and then prints the text report underneath.
    """
    # 1. Run the decoupled analysis wizard
    res = run_domain_comparison_wizard(df, query_domain_size, query_length, z_score_low_threshold)

    if res["clean_df_size"] == 0:
        print("No hits found matching current criteria.")
        return None

    # 2. Build the visualization layout
    fig, (ax_kde, ax_heat) = plt.subplots(
        1, 2,
        figsize=(10, 3),
        gridspec_kw={'width_ratios': [7, 3]}
    )

    # A. KDE plot
    ax_kde.plot(res["eval_range"], res["q_densities"], label=f"Query Coverage (Peak: {res['q_peak']:.2f})", color="royalblue", lw=2)
    ax_kde.fill_between(res["eval_range"], 0, res["q_densities"], color="royalblue", alpha=0.15)

    ax_kde.plot(res["eval_range"], res["t_densities"], label=f"Target Coverage (Peak: {res['t_peak']:.2f})", color="orange", lw=2)
    ax_kde.fill_between(res["eval_range"], 0, res["t_densities"], color="orange", alpha=0.15)

    ax_kde.set_xlim(0, 1.0)
    ax_kde.set_ylim(bottom=0)
    ax_kde.set_xlabel("Coverage Value")
    ax_kde.set_ylabel("Density (KDE)")
    ax_kde.set_title("Coverage Distribution Smooth Profiles")
    ax_kde.legend(loc="upper left")
    ax_kde.grid(True, linestyle="--", alpha=0.5)

    # B. Contingency Heatmap
    sns.heatmap(
        res["contingency"] * 100, 
        annot=True,
        fmt=".1f", 
        cmap="Blues", 
        cbar_kws={'label': 'Percentage of Total Hits (%)'},
        ax=ax_heat,
        annot_kws={"size": 12, "weight": "bold"}
    )
    ax_heat.set_title("Structural Contingency Layout")
    ax_heat.set_ylabel("Query Coverage")
    ax_heat.set_xlabel("Target Coverage")

    plt.tight_layout()

    # === THE FLIP MECHANISM ===
    # Force the current figure to draw right now inside the widget container
    plt.show()

    if True: return None

    # 3. Print out the report (will now land strictly below the forced canvas)
    print("="*80)
    print("                    DOMAIN ARCHITECTURE WIZARD REPORT                    ")
    print("="*80)
    print(f"Dataset Size : {res['clean_df_size']} hits analyzed")
    print(f"Median alignment length: {res['median_length']}")
    print("-"*80)
    print(res['peak_inference'])
    print("-"*80)
    print(res['table_inference'])
    print("="*80)

    # Return None so widget_view doesn't try to draw a second duplicate copy
    return None

def make_domain_comparison_fn(project):
    """
    Returns a widget_view-compatible runner. 
    Matches current view selections and resolves physical domain sizes dynamically.
    """
    def runner(view_name, z_score_low_threshold, **kwargs):
        df = project.views[view_name]

        # Trace back to nearest Silver view
        domain_ranges = get_domain_ranges(project, view_name)

        query_length = project.query_length

        if domain_ranges:
            # Reuses the domain_range parser
            query_domain_size = get_query_domain_length(domain_ranges)
        else:
            query_domain_size = query_length

        # Suppress automatic drawing so widget_view manages display(fig) cleanly
        plt.ioff()
        fig = plot_domain_comparison(
            df=df,
            query_domain_size=query_domain_size,
            query_length=query_length,
            z_score_low_threshold=z_score_low_threshold
        )
        plt.ion()
        return fig

    return runner
