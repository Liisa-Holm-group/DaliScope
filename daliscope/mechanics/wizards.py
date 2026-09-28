import numpy as np
import pandas as pd

#############################

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _passes_effect(var: str, median_diff: float, min_effect: dict) -> bool:
    if not min_effect or var not in min_effect:
        return True
    return abs(median_diff) >= min_effect[var]


def group_distribution_test(
    df: pd.DataFrame,
    k: int = 100,
    variables: list = None,
    min_group_size: int = 3,
    alpha: float = 1e-5,
    min_effect: dict = None,
) -> pd.DataFrame:
    """Two-pass distribution test that identifies the primary Home Clan first,

    then purges it from the background of alternative groups to protect
    statistical power.
    """
    if variables is None:
        variables = ["z_score", "sequence_identity", "query_coverage"]

    required = {"pfam", "clan", "pfam_rank", "clan_rank"} | set(variables)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    clean_df = df[["pfam", "clan"]].dropna()
    cl = dict(zip(clean_df["pfam"], clean_df["clan"]))
    cl.update({clan: clan for clan in clean_df["clan"].unique()})

    # ---------------------------------------------------------------- #
    # PASS 1: Locate the Home Clan Anchor via Clan-level Z-Scores      #
    # ---------------------------------------------------------------- #
    home_clan = None
    top_k_clans = df[df["clan_rank"] <= k]["clan"].dropna().unique()

    if len(top_k_clans) > 0 and "z_score" in df.columns:
        clan_z_medians = {}
        for clan in top_k_clans:
            c_vals = df[df["clan"] == clan]["z_score"].dropna()
            if len(c_vals) >= min_group_size:
                clan_z_medians[clan] = c_vals.median()

        if clan_z_medians:
            # The clan with the highest median z_score is our structural home
            home_clan = max(clan_z_medians, key=clan_z_medians.get)

    # ---------------------------------------------------------------- #
    # PASS 2: Main Evaluation Loop with Unbiased Backgrounds          #
    # ---------------------------------------------------------------- #
    group_configs = [("pfam", "pfam_rank"), ("clan", "clan_rank")]
    records = []

    for group_col, rank_col in group_configs:
        top_k_groups = df[df[rank_col] <= k][group_col].dropna().unique()
        if len(top_k_groups) == 0:
            continue

        for var in variables:
            p_values = []
            valid_groups = []

            for group in top_k_groups:
                grp_vals = df[df[group_col] == group][var].dropna()
                if len(grp_vals) < min_group_size:
                    continue

                # --- THE FIX: Isolate the standard background slice ---
                bg_mask = df[group_col] != group

                # If this group isn't the home clan (or a family inside it),
                # drop the dominant home clan entirely out of their background population!
                if home_clan and group != home_clan and cl.get(group) != home_clan:
                    bg_mask = bg_mask & (df["clan"] != home_clan)

                bg_vals = df[bg_mask][var].dropna()
                # -------------------------------------------------------

                stat, p = mannwhitneyu(
                    grp_vals, bg_vals, alternative="greater"
                )
                p_values.append(p)
                valid_groups.append((group, stat, grp_vals, bg_vals))

            if not p_values:
                continue

            _, q_values, _, _ = multipletests(
                p_values, alpha=alpha, method="fdr_bh"
            )

            for idx, (group, stat, grp_vals, bg_vals) in enumerate(
                valid_groups
            ):
                grp_med = grp_vals.median()
                grp_max = grp_vals.max()
                bg_med = bg_vals.median()
                records.append(
                    {
                        "group_type": group_col,
                        "group": group,
                        "clan": cl[group],
                        "variable": var,
                        "n": len(grp_vals),
                        "mean": round(grp_vals.mean(), 3),
                        "median": round(grp_med, 3),
                        "max": round(grp_max, 3),
                        "bg_median": round(bg_med, 3),
                        "median_diff": round(grp_med - bg_med, 3),
                        "fold_change": (
                            round(grp_med / bg_med, 3) if bg_med != 0 else np.nan
                        ),
                        "u_stat": round(stat, 1),
                        "p_value": round(p_values[idx], 4),
                        "q_value": round(q_values[idx], 4),
                        "significant": (
                            bool(q_values[idx] < alpha)
                            and _passes_effect(
                                var, grp_med - bg_med, min_effect
                            )
                        ),
                    }
                )

    result = pd.DataFrame(records)
    if result.empty:
        return result

    return result.sort_values(
        ["group_type", "variable", "significant", "q_value"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)


# ===================================================================== #
#  PART 2: THE INTERPRETATION WIZARD                                    #
# ===================================================================== #


def get_significant_hits(
    results: pd.DataFrame, target_variables: list
) -> pd.DataFrame:
    return results[
        results["significant"] & results["variable"].isin(target_variables)
    ].copy()


def determine_home_hierarchy(sig_hits: pd.DataFrame) -> tuple:
    z_hits = sig_hits[sig_hits["variable"] == "z_score"]
    clan_z = z_hits[z_hits["group_type"] == "clan"].sort_values(
        "median", ascending=False
    )
    home_clan = clan_z.iloc[0]["group"] if not clan_z.empty else None

    pfam_z = z_hits[z_hits["group_type"] == "pfam"].sort_values(
        "median", ascending=False
    )
    if home_clan and not pfam_z.empty:
        clan_families = pfam_z[pfam_z["clan"] == home_clan]
        home_family = (
            clan_families.iloc[0]["group"]
            if not clan_families.empty
            else pfam_z.iloc[0]["group"]
        )
    else:
        home_family = pfam_z.iloc[0]["group"] if not pfam_z.empty else None

    return home_clan, home_family


def extract_group_stats(
    results_df: pd.DataFrame, group_id: str, variables: list
) -> list:
    group_rows = results_df[results_df["group"] == group_id]
    stats = []
    if group_rows.empty:
        return ["No stat metrics available"]

    n_val = int(group_rows["n"].max())
    stats.append(f"n={n_val}")

    for var in variables:
        var_row = group_rows[group_rows["variable"] == var]
        if not var_row.empty:
            med = var_row["median"].values[0]
            diff = var_row["median_diff"].values[0]
            q = var_row["q_value"].values[0]
            name = var.replace("_", " ")

            if q < 0.001:
                star = " ***"
            elif q < 0.01:
                star = " **"
            elif q < 0.05:
                star = " *"
            else:
                star = " ns"

            stats.append(
                f"{name}: med={med:.2f} (Δ={diff:+.2f}, q={q:.1e}{star})"
            )

    return stats


def build_home_section(
    sig_hits: pd.DataFrame, home_clan: str, home_family: str, variables: list
) -> list:
    lines = ["═" * 60, "PRIMARY STRUCTURAL HOME REFERENCE", "─" * 40]
    if home_clan:
        clan_stats = extract_group_stats(sig_hits, home_clan, variables)
        lines.append(f"HOME CLAN:   {home_clan}")
        for stat in clan_stats:
            lines.append(f"      · {stat}")
    else:
        lines.append("HOME CLAN:   None Significant")

    if home_family:
        fam_stats = extract_group_stats(sig_hits, home_family, variables)
        lines.append(f"HOME FAMILY: {home_family}")
        for stat in fam_stats:
            lines.append(f"      · {stat}")
    else:
        lines.append("HOME FAMILY: None Significant")
    return lines


def build_candidates_section(
    sig_hits: pd.DataFrame, home_clan: str, home_family: str, variables: list
) -> list:
    lines = ["", "SIGNIFICANT NOVEL / OUTLIER CLAN CANDIDATES", "─" * 40]
    all_sig_clans = sig_hits[sig_hits["group_type"] == "clan"]["group"].unique()
    candidate_clans = [c for c in all_sig_clans if c != home_clan]

    if not candidate_clans:
        lines.append(
            "  No alternative clans reached significance across the tested criteria."
        )
        return lines

    lines.append(
        f"  Found {len(candidate_clans)} candidates for unification into {home_clan}:"
    )
    for clan_id in sorted(candidate_clans):
        clan_stats = extract_group_stats(sig_hits, clan_id, variables)
        lines.append(f"\n  ✓ Clan: {clan_id}")
        for stat in clan_stats:
            lines.append(f"      · {stat}")
    return lines


def clan_discovery_wizard(results: pd.DataFrame, df: pd.DataFrame) -> str:
    target_variables = ["z_score", "query_coverage", "sequence_identity"]
    sig_hits = get_significant_hits(results, target_variables)
    if sig_hits.empty:
        return "No significant group-variable associations found."

    home_clan, home_family = determine_home_hierarchy(sig_hits)

    report_lines = []
    report_lines.extend(
        build_home_section(sig_hits, home_clan, home_family, target_variables)
    )
    report_lines.extend(
        build_candidates_section(sig_hits, home_clan, home_family, target_variables)
    )
    return "\n".join(report_lines)

#############################

def classify_seq_identity(pid: float) -> str:
    if pid > 0.30:
        return "homolog"
    elif pid >= 0.15:
        return "twilight"
    else:
        return "dark"


def classify_z_score(z: float) -> str:
    if z > 16:
        return "family"
    elif z >= 12:
        return "superfamily"
    elif z >= 8:
        return "fold"
    elif z >= 4:
        return "marginal"
    else:
        return "noise"


def domain_wizard(master_metadata: pd.DataFrame,
                  seq_id_col: str = 'sequence_identity',
                  z_score_col: str = 'z_score',
                  coverage_col: str = 'query_coverage',
                  coverage_threshold: float = 0.5,
                  coverage_fraction_threshold: float = 0.3,
                  ) -> str:
    """
    Generate plain-language orientation commentary for Naive User,
    based on the overall hit distribution in master_metadata.

    Logic
    -----
    1. Look at the BEST hit (highest z_score) — its seq identity
       and z_score bands jointly determine the homology classification
       message and the recommended next step.
    2. Look at overall query_coverage distribution to suggest whether
       the query is likely single-domain or multidomain.
    3. Always mention Signature Motif as a tool that works regardless
       of how remote the homology is.

    Parameters
    ----------
    master_metadata : the full hit table (or any filtered view of it)
    seq_id_col, z_score_col, coverage_col : column names
    coverage_threshold : query_coverage value counted as "high coverage"
    coverage_fraction_threshold : fraction of hits that must clear
                                  coverage_threshold for the query to
                                  be called "probably single-domain"

    Returns
    -------
    Multi-paragraph string — ready to print() or display in a notebook
    """
    df = master_metadata.dropna(subset=[seq_id_col, z_score_col])
    if df.empty:
        return ("No scored hits available — cannot generate orientation "
                "commentary. Check that the pack loaded correctly and "
                "alignments were scored.")

    lines = []

    # ---------------------------------------------------------------- #
    # 1. Best hit — homology classification                            #
    # ---------------------------------------------------------------- #
    best = df.loc[df[z_score_col].idxmax()]
    best_pid   = best[seq_id_col]
    best_z     = best[z_score_col]
    best_id    = best.get('target_id', '?')

    pid_class = classify_seq_identity(best_pid)
    z_class   = classify_z_score(best_z)

    if False:
      lines.append(
        f"Best hit: {best_id}  "
        f"(Z-score={best_z:.1f}, sequence identity={best_pid*100:.0f}%)"
      )
      lines.append("")

    if pid_class == "homolog":
        lines.append(
            "Your best hit has sequence identity above 30% — this is "
            "a clear homolog. The query's classification at the family "
            "level is essentially trivial: a standard sequence search "
            "(BLAST/HMMER) would likely have found this same relationship "
            "without needing structural comparison at all."
        )

    elif pid_class == "twilight":
        if z_class in ("family", "superfamily"):
            lines.append(
                f"Your best hit sits in the 'twilight zone' for sequence "
                f"identity (15-30%), but the structural Z-score "
                f"({best_z:.1f}) is strong enough to suggest a "
                f"{z_class}-level relationship. Sequence comparison alone "
                f"would likely have missed or under-supported this call — "
                f"this is exactly the regime where structural comparison "
                f"earns its place. Recommendation: check the Pfam "
                f"annotation of this and related top hits — if they agree "
                f"on a specific Pfam family or clan, that gives you a "
                f"precise functional assignment to anchor on."
            )
        else:
            lines.append(
                f"Your best hit has twilight-zone sequence identity "
                f"(15-30%) but only a {z_class}-level structural Z-score "
                f"({best_z:.1f}). The structural signal alone is not "
                f"strongly diagnostic here. Treat any family-level claim "
                f"with caution until corroborated by Pfam or other "
                f"independent evidence."
            )

    else:  # dark
        if z_class in ("family", "superfamily"):
            lines.append(
                f"Your best hit has very low sequence identity "
                f"({best_pid*100:.0f}%) — invisible to sequence search — "
                f"but the structural Z-score ({best_z:.1f}) indicates a "
                f"likely {z_class}-level relationship. This is the core "
                f"value case for structural comparison: a real "
                f"evolutionary relationship recovered purely from fold "
                f"geometry. Recommendation: check the Pfam annotation of "
                f"this and other strong hits — if several independently "
                f"point to the same family, that's strong corroborating "
                f"evidence despite the sequence divergence."
            )
        elif z_class == "fold":
            lines.append(
                f"Your best hit has very low sequence identity "
                f"({best_pid*100:.0f}%) and a fold-level Z-score "
                f"({best_z:.1f}) — this suggests structural similarity "
                f"without necessarily implying common ancestry (could be "
                f"convergent fold usage). Recommendation: check whether "
                f"top hits share an unambiguous Pfam CLAN association — "
                f"clan-level agreement across independent hits would "
                f"support a real (if very remote) evolutionary link, "
                f"while disagreement suggests the structural similarity "
                f"may be coincidental fold convergence rather than homology."
            )
        else:
            lines.append(
                f"Your best hit has very low sequence identity "
                f"({best_pid*100:.0f}%) and only a {z_class}-level "
                f"Z-score ({best_z:.1f}). This combination does not "
                f"provide strong evidence of homology on its own. "
                f"Treat results from this dataset cautiously — consider "
                f"whether the query or domain selection needs revisiting."
            )

    lines.append("")
    lines.append(
        "Before trusting the homology classification above too far, run "
        "the Rigid Core analysis — it checks whether the query's anchor "
        "domain is a complete fold or just a partial fragment, and does "
        "the same for your top hits. A partial fold on either side can "
        "make a real relationship look weaker than it is, or a weak one "
        "look stronger."
    )

    lines.append("")
    lines.append(
        "Signature Motif analysis (derived from "
        "the structurally conserved positions across your trusted hit "
        "set) can reveal likely enzymatic or functional residues even "
        "when overall homology is too remote for confident family "
        "assignment — it does not depend on the classification above."
    )


    # ---------------------------------------------------------------- #
    # 2. Single-domain vs multidomain hint, from coverage distribution  #
    # ---------------------------------------------------------------- #
    lines.append("")
    cov = df[coverage_col].dropna()
    if len(cov) > 0:
        frac_high_cov = (cov >= coverage_threshold).mean()

        if frac_high_cov >= coverage_fraction_threshold:
            lines.append(
                f"{frac_high_cov*100:.0f}% of hits cover at least "
                f"{coverage_threshold*100:.0f}% of the query length — "
                f"this pattern is consistent with a SINGLE-DOMAIN query. "
                f"Verify by looking at the query structure in the py3Dmol "
                f"viewer and the DSSP heatmap: a single compact domain "
                f"should appear as one continuous folded unit with no "
                f"obvious break in secondary structure pattern."
            )
        else:
            lines.append(
                f"Only {frac_high_cov*100:.0f}% of hits cover at least "
                f"{coverage_threshold*100:.0f}% of the query length — "
                f"this pattern is consistent with a MULTIDOMAIN query, "
                f"where most hits only match one part of it. "
                f"Verify by looking at the query structure in the py3Dmol "
                f"viewer and the DSSP heatmap: look for visually distinct "
                f"globular regions connected by a linker, or breaks in "
                f"the secondary structure pattern that suggest separate "
                f"domains. If confirmed, consider clipping to individual "
                f"domains (Step 2) before further analysis — combining "
                f"unrelated domains in one view can obscure real signal "
                f"for each domain individually."
            )
    else:
        lines.append(
            "Query coverage data unavailable — cannot assess "
            "single- vs multi-domain status from this view."
        )

    return "\n".join(lines)

##############################

def analyze_sequence_logo_wizard(
    info_df,
    low_conservation_cutoff=2.0,
    diversity_tolerance=0.15,
    enzymatic_height_multiplier=1.5,
    verbose=False
):
    """
    Analyzes an Information-Content DataFrame to provide automated biological insights.
    """
    # Standard biological groupings for amino acids
    polar_charged = set(['R', 'K', 'D', 'E', 'N', 'Q', 'S', 'T', 'Y', 'H', 'C'])
    hydrophobic = set(['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'P', 'G'])

    # 1. Extract fundamental metrics per position (FIXED: peak height is total row sum)
    total_heights = info_df.sum(axis=1)  # Total height of the sequence logo stack
    max_aa = info_df.idxmax(axis=1)      # The most dominant letter at this position

    # Global profile metrics based on stack heights
    global_avg_height = total_heights.mean()
    global_max_height = total_heights.max()

    if verbose:
        print("--- Sequence Logo Wizard Report ---")
        print(f"Analyzed {len(info_df)} positions.")
        print(f"Global average stack height: {global_avg_height:.2f} bits")
        print(f"Global maximum stack height: {global_max_height:.2f} bits\n")

    # 2. Rule Evaluation
    insights = []

    # Rule A: Overall Conservation Check
    if global_max_height < low_conservation_cutoff:
        insights.append(
            f"⚠️ LOW CONSERVATION: The absolute maximum stack height ({global_max_height:.2f} bits) "
            f"is below {low_conservation_cutoff}. The entire alignment shows poor sequence conservation."
        )

    # Rule B: Diversity Check (Homogeneous conservation profile)
    height_variance = total_heights.std()
    if height_variance < (global_avg_height * diversity_tolerance):
        insights.append(
            f"⚖️ LOW DIVERSITY PROFILE: Stack heights are uniformly flat across the profile "
            f"(std dev: {height_variance:.2f}). This suggests a lack of distinct functional hyper-conserved regions."
        )

    # Gather highly conserved positions for rule C and D
    highly_conserved_mask = total_heights >= low_conservation_cutoff
    hc_positions = info_df[highly_conserved_mask]

    if not hc_positions.empty:
        hc_aas = max_aa[highly_conserved_mask]
        hc_heights = total_heights[highly_conserved_mask]

        # Count polar vs hydrophobic in the top conserved positions
        polar_count = sum(1 for aa in hc_aas if aa in polar_charged)
        hydrophobic_count = sum(1 for aa in hc_aas if aa in hydrophobic)

        # Rule C: Enzymatic function likely
        prominent_polar_peaks = [
            (pos, aa, height) for pos, aa, height in zip(hc_positions.index, hc_aas, hc_heights)
            if aa in polar_charged and height > (global_avg_height * enzymatic_height_multiplier)
        ]

        if 1 <= len(prominent_polar_peaks) <= 5:
            pos_summary = ", ".join([f"{aa}{pos} ({h:.2f}b)" for pos, aa, h in prominent_polar_peaks])
            insights.append(
                f"🧬 ENZYMATIC FUNCTION LIKELY: Found a small cluster of highly prominent, "
                f"polar/charged residues significantly above average heights: {pos_summary}. "
                f"This looks like a catalytic active site or binding pocket."
            )

        # Rule D: Enzymatic function unlikely
        if hydrophobic_count > polar_count and (hydrophobic_count / len(hc_aas)) > 0.60:
            insights.append(
                f"🛑 ENZYMATIC FUNCTION UNLIKELY: Conserved positions are heavily dominated by "
                f"non-polar/hydrophobic amino acids ({hydrophobic_count}/{len(hc_aas)}). "
                f"This strongly points to structural core preservation or a transmembrane domain rather than enzymatic activity."
            )
    else:
        if global_max_height >= low_conservation_cutoff:
            insights.append("ℹ️ No single positions cross the strict conservation threshold rules.")

    # 3. Print Results
    if insights:
        for insight in insights:
            return(insight)
    else:
        return("✅ Profile looks standard. No extreme structural or catalytic anomalies detected.")


########################

from scipy import stats
from scipy.signal import find_peaks


def distribution_wizard(df: pd.DataFrame,
                        variables: list = None) -> str:
    """
    Examines distributions of alignment quality variables and
    generates data-driven plain-language insights.

    Parameters
    ----------
    df        : DataFrame with alignment statistics
    variables : columns to examine — defaults to the standard set
    """
    if variables is None:
        variables = ['z_score', 'alignment_length', 'sequence_identity',
                     'query_coverage', 'target_coverage', 'target_length']

    variables = [v for v in variables if v in df.columns]
    n         = len(df)
    lines     = []

    lines.append(f"DISTRIBUTION ANALYSIS  ({n} alignments)")
    lines.append("═" * 60)

    # ---------------------------------------------------------------- #
    # Helper: detect modality                                           #
    # ---------------------------------------------------------------- #
    def modality(series: pd.Series) -> str:
        """Rough modality detection via kernel density peaks."""
        s = series.dropna()
        if len(s) < 20:
            return 'insufficient data'
        kde_x = np.linspace(s.min(), s.max(), 300)
        kde_y = stats.gaussian_kde(s, bw_method=0.2)(kde_x)
        peaks, props = find_peaks(kde_y,
                                  prominence=kde_y.max() * 0.15,
                                  distance=20)
        n_peaks = len(peaks)
        if n_peaks == 0:
            return 'flat'
        elif n_peaks == 1:
            return f'unimodal (peak ≈ {kde_x[peaks[0]]:.2f})'
        elif n_peaks == 2:
            return (f'bimodal (peaks ≈ {kde_x[peaks[0]]:.2f} '
                    f'and {kde_x[peaks[1]]:.2f})')
        else:
            return f'{n_peaks}-modal'

    def pct(mask) -> str:
        return f"{100 * mask.sum() / n:.0f}%"

    # ---------------------------------------------------------------- #
    # Z-score                                                           #
    # ---------------------------------------------------------------- #
    if 'z_score' in df.columns:
        z       = df['z_score'].dropna()
        z_max   = z.max()
        z_med   = z.median()
        z_mod   = modality(z)
        frac_noise      = pct(z < 4)
        frac_marginal   = pct((z >= 4)  & (z < 8))
        frac_fold       = pct((z >= 8)  & (z < 12))
        frac_superfam   = pct((z >= 12) & (z < 16))
        frac_family     = pct(z >= 16)

        lines.append("\nZ-SCORE")
        lines.append(f"  Distribution: {z_mod}")
        lines.append(f"  Median={z_med:.1f}, max={z_max:.1f}")
        lines.append(f"  Noise (<4): {frac_noise}  |  Marginal (4-8): {frac_marginal}  |  "
                     f"Fold (8-12): {frac_fold}  |  Superfamily (12-16): {frac_superfam}  |  "
                     f"Family (>16): {frac_family}")

        if z_max > 30:
            lines.append(f"  ✓ Very high maximum Z-score ({z_max:.1f}) — "
                         f"strong structural homologs present.")
        if 'bimodal' in z_mod:
            lines.append(f"  ⚠ Bimodal Z-score distribution suggests two "
                         f"structurally distinct populations — consider "
                         f"domain clipping if one peak is below 8.")
        if float(frac_noise.rstrip('%')) > 30:
            lines.append(f"  ⚠ {frac_noise} of hits are likely noise (Z<4) — "
                         f"consider filtering before analysis.")

    # ---------------------------------------------------------------- #
    # Alignment length                                                  #
    # ---------------------------------------------------------------- #
    if 'alignment_length' in df.columns:
        al      = df['alignment_length'].dropna()
        al_med  = al.median()
        al_max  = al.max()
        al_mod  = modality(al)
        frac_short = pct(al < 50)

        lines.append("\nALIGNMENT LENGTH")
        lines.append(f"  Distribution: {al_mod}")
        lines.append(f"  Median={al_med:.0f} residues, max={al_max:.0f}")

        if 'bimodal' in al_mod:
            lines.append(f"  ⚠ Bimodal alignment length suggests a mix of "
                         f"full-domain and partial hits. Short alignments at "
                         f"low Z-scores are likely noise.")
        if float(frac_short.rstrip('%')) > 20:
            lines.append(f"  ⚠ {frac_short} of alignments are very short "
                         f"(<50 residues) — may indicate fragmentary matches "
                         f"or noisy hits at the periphery of the fold family.")
        if al_max > 400:
            lines.append(f"  ℹ Maximum alignment length is {al_max:.0f} residues — "
                         f"some hits cover a large portion of the query, "
                         f"consistent with full-domain or multidomain matches.")

    # ---------------------------------------------------------------- #
    # Sequence identity                                                 #
    # ---------------------------------------------------------------- #
    if 'sequence_identity' in df.columns:
        si       = df['sequence_identity'].dropna()
        si_max   = si.max()
        si_med   = si.median()
        si_mod   = modality(si)
        frac_dark      = pct(si < 0.15)
        frac_twilight  = pct((si >= 0.15) & (si < 0.30))
        frac_homolog   = pct(si >= 0.30)

        lines.append("\nSEQUENCE IDENTITY")
        lines.append(f"  Distribution: {si_mod}")
        lines.append(f"  Median={si_med*100:.1f}%, max={si_max*100:.1f}%")
        lines.append(f"  Dark zone (<15%): {frac_dark}  |  "
                     f"Twilight (15-30%): {frac_twilight}  |  "
                     f"Clear homologs (>30%): {frac_homolog}")

        if si_med < 0.15:
            lines.append(f"  ✓ Median sequence identity is {si_med*100:.1f}% — "
                         f"the majority of hits are in the midnight zone. "
                         f"These are remote structural homologs invisible to "
                         f"sequence-based search — exactly where DALI adds value.")
        elif si_med < 0.30:
            lines.append(f"  ℹ Median sequence identity ({si_med*100:.1f}%) is in "
                         f"the twilight zone — structural alignment is essential "
                         f"for reliable homology assignment here.")
        else:
            lines.append(f"  ℹ Median sequence identity ({si_med*100:.1f}%) is above "
                         f"the twilight zone — many hits could also have been "
                         f"found by sequence search alone.")
        if 'bimodal' in si_mod:
            lines.append(f"  ⚠ Bimodal identity distribution may indicate "
                         f"two distinct homolog classes — e.g. close paralogs "
                         f"plus a separate set of remote homologs.")

    # ---------------------------------------------------------------- #
    # Query coverage                                                    #
    # ---------------------------------------------------------------- #
    if 'query_coverage' in df.columns:
        qc      = df['query_coverage'].dropna()
        qc_med  = qc.median()
        qc_mod  = modality(qc)
        frac_low  = pct(qc < 0.3)
        frac_high = pct(qc > 0.7)

        lines.append("\nQUERY COVERAGE")
        lines.append(f"  Distribution: {qc_mod}")
        lines.append(f"  Median={qc_med*100:.1f}%")
        lines.append(f"  Low coverage (<30%): {frac_low}  |  "
                     f"High coverage (>70%): {frac_high}")

        if qc_med > 0.7:
            lines.append(f"  ✓ Most hits cover >70% of the query — consistent "
                         f"with a single-domain query with well-represented "
                         f"homologs.")
        elif qc_med < 0.4:
            lines.append(f"  ⚠ Median query coverage is only {qc_med*100:.1f}% — "
                         f"the query is likely multidomain and hits are "
                         f"matching individual domains. Consider domain "
                         f"clipping (Step 2) before further analysis.")
        if 'bimodal' in qc_mod:
            lines.append(f"  ⚠ Bimodal query coverage strongly suggests "
                         f"a multidomain query — one population of hits "
                         f"covers a single domain, another covers more. "
                         f"Inspect the DSSP heatmap and py3Dmol viewer "
                         f"to identify the domain boundary.")

    # ---------------------------------------------------------------- #
    # Target coverage                                                   #
    # ---------------------------------------------------------------- #
    if 'target_coverage' in df.columns:
        tc      = df['target_coverage'].dropna()
        tc_med  = tc.median()
        tc_mod  = modality(tc)
        frac_module = pct(tc < 0.3)

        lines.append("\nTARGET COVERAGE")
        lines.append(f"  Distribution: {tc_mod}")
        lines.append(f"  Median={tc_med*100:.1f}%")

        if float(frac_module.rstrip('%')) > 20:
            lines.append(f"  ℹ {frac_module} of targets have low coverage (<30%) — "
                         f"these are likely mobile module matches: a conserved "
                         f"functional domain present in much larger, unrelated "
                         f"multidomain proteins. The shared module may be "
                         f"structurally identical but the host proteins have "
                         f"different architectures and biological roles.")
        if tc_med > 0.7:
            lines.append(f"  ✓ Most hits cover >70% of the target — "
                         f"suggesting hits are similar-sized proteins, "
                         f"not just domain insertions.")

    # ---------------------------------------------------------------- #
    # Target length                                                     #
    # ---------------------------------------------------------------- #
    if 'target_length' in df.columns:
        tl      = df['target_length'].dropna()
        tl_med  = tl.median()
        tl_mod  = modality(tl)
        frac_large = pct(tl > 500)

        lines.append("\nTARGET LENGTH")
        lines.append(f"  Distribution: {tl_mod}")
        lines.append(f"  Median={tl_med:.0f} residues")

        if 'bimodal' in tl_mod:
            lines.append(f"  ℹ Bimodal target length suggests the hit set "
                         f"contains both compact single-domain proteins and "
                         f"larger multidomain proteins carrying the same fold.")
        if float(frac_large.rstrip('%')) > 30:
            lines.append(f"  ℹ {frac_large} of targets are large (>500 residues) — "
                         f"many hits are likely multidomain proteins where "
                         f"only one domain matches the query.")

    lines.append("\n" + "═" * 60)
    return "\n".join(lines)

########################
