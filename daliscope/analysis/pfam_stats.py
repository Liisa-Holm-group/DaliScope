"""
Pfam/clan distribution analysis.

One-stop function for the full workflow:
  1. violin plots     (optics)
  2. Mann-Whitney U   (mechanics/wizards)
  3. clan discovery   (mechanics/wizards)
"""

from __future__ import annotations
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import HTML, display

from daliscope.mechanics.wizards import (
    group_distribution_test,
    clan_discovery_wizard,
)
from daliscope.optics.matplot import plot_violins

import plotly.express as px

#==============

def top_category_label(
    df: pd.DataFrame,
    column: str,
    k: int = 10,
    other_label: str = "Other",
    output_name: str | None = None,
) -> pd.Series:
    """
    Return a Series that keeps the top-k categories in `column` and
    collapses all other categories into `other_label`.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the source column.
    column : str
        Name of the categorical column to summarize.
    k : int, default 10
        Number of most frequent categories to retain.
    other_label : str, default "Other"
        Label assigned to categories outside the top-k.
    output_name : str, optional
        Name of the returned Series. Defaults to `column + "_label"`.

    Returns
    -------
    pd.Series
        Series aligned with `df.index`, suitable for direct assignment
        to the dataframe.

    Examples
    --------
    df["pfam_label"] = top_category_label(df, "pfam", k=10)

    df["clan_label"] = top_category_label(
        df, "clan", k=5, output_name="clan_label"
    )
    """
    if output_name is None:
        output_name = f"{column}_label"

    top_k = df[column].value_counts().iloc[:k].index

    return (
        df[column]
        .where(df[column].isin(top_k), other=other_label)
        .rename(output_name)
    )

# Base URLs
CLAN_BASE_URL = "https://www.ebi.ac.uk/interpro/set/pfam/"
PFAM_BASE_URL = "https://www.ebi.ac.uk/interpro/entry/pfam/"


def add_pfam_links_vectorized(df: pd.DataFrame) -> pd.DataFrame:
    """Pure NumPy vectorized HTML link generation for 'clan' and 'pfam' columns."""
    df = df.copy()

    def _vectorized_link(series, base_url, prefix):
        # Convert series to string representation cleanly
        s_str = series.astype(str)

        # Vectorized prefix check using Pandas .str accessor
        mask = s_str.str.startswith(prefix).fillna(False)

        # Vectorized string concatenation
        linked = (
            '<a href="' + base_url + s_str + '/" target="_blank">' + s_str + "</a>"
        )

        # Assign back formatted links where valid, keeping original series values elsewhere
        return np.where(mask, linked, series)

    if "clan" in df.columns:
        df["clan"] = _vectorized_link(df["clan"], CLAN_BASE_URL, "CL")

    if "pfam" in df.columns:
        df["pfam"] = _vectorized_link(df["pfam"], PFAM_BASE_URL, "PF")

    return df

def display_html_df(df: pd.DataFrame):
    """Renders DataFrame containing raw HTML markup directly in Jupyter."""
    display(HTML(df.to_html(escape=False, index=False)))


import sys
from IPython.display import display

def pfam_counts_df(
    df,
    pfam_names,
    view_name=None,
    k=5,
    print_tables=True,
):

    if pfam_names is None or pfam_names.empty:
        print("Warning: No pfam_names table supplied.")
        sys.stdout.flush()
        return None

    # 2. Always work on isolated copies
    df = df.copy()
    pfam_names_local = pfam_names.copy()

    # 3. Group counts
    clan_counts_df = (
        df.groupby("clan", dropna=False)
        .size()
        .reset_index(name="count")
    )

    pfam_counts_df = (
        df.groupby("pfam", dropna=False)
        .size()
        .reset_index(name="count")
    )

    # 4. Merge Pfam details
    pfam_cols = [
        c for c in ["clan", "clan_short", "pfam", "short", "name"]
        if c in pfam_names_local.columns
    ]

    pfam_counts_merged = pfam_counts_df.merge(
        pfam_names_local[pfam_cols].drop_duplicates("pfam"),
        on="pfam",
        how="left",
    )

    # 5. Merge clan details
    clan_cols = [
        c for c in ["clan", "clan_short"]
        if c in pfam_names_local.columns
    ]

    clan_counts_merged = clan_counts_df.merge(
        pfam_names_local[clan_cols].drop_duplicates("clan"),
        on="clan",
        how="left",
    )

    # 6. Render outputs
    if print_tables:
        population_label = (
            f"'{view_name}' population"
            if view_name is not None
            else "population"
        )

        print(f"Largest families in {population_label} (n={len(df)}):")

        top_families = (
            pfam_counts_merged
            .sort_values("count", ascending=False)
            .head(k)
        )

        formatted_families = add_pfam_links_vectorized(
            top_families.copy()
        )

        if formatted_families is not None and not formatted_families.empty:
            display_html_df(formatted_families)
        else:
            print("No family records to display.")

        print(f"\nLargest superfamilies in {population_label}:")

        top_clans = (
            clan_counts_merged
            .sort_values("count", ascending=False)
            .head(k)
        )

        formatted_clans = add_pfam_links_vectorized(
            top_clans.copy()
        )

        if formatted_clans is not None and not formatted_clans.empty:
            display_html_df(formatted_clans)
        else:
            print("No superfamily records to display.")

    sys.stdout.flush()
    return None

def pfam_counts(
    view_name,
    project,
    pfam_names=None,
    k=5,
    sort_by=None,
    ascending=False,
    print_tables=True,
):
    # 1. Fetch current view DataFrame
    df = project.views.get(view_name, None)
    if df is None or df.empty:
        print(f"View '{view_name}' is empty or does not exist.")
        sys.stdout.flush()
        return None

    # 2. Safely resolve pfam_names
    if pfam_names is None:
        pfam_names = getattr(project, "_pfam_names", None)

    if pfam_names is None or (isinstance(pfam_names, pd.DataFrame) and pfam_names.empty):
        print(f"Warning: No pfam_names table found in project for view '{view_name}'.")
        sys.stdout.flush()
        return None

    # 3. Always work on isolated copies to prevent state corruption
    pfam_names_local = pfam_names.copy()

    # 4. Groupby counts
    clan_counts_df = (
        df.groupby(["clan"], dropna=False).size().reset_index(name="count")
    )
    pfam_counts_df = (
        df.groupby(["pfam"], dropna=False).size().reset_index(name="count")
    )

    # 5. Merge pfam details
    pfam_cols = [c for c in ["clan", "clan_short", "pfam", "short", "name"] if c in pfam_names_local.columns]
    pfam_counts_merged = pfam_counts_df.merge(
        pfam_names_local[pfam_cols].drop_duplicates("pfam"),
        on="pfam",
        how="left",
    )

    # 6. Merge clan details
    clan_cols = [c for c in ["clan", "clan_short"] if c in pfam_names_local.columns]
    clan_counts_merged = clan_counts_df.merge(
        pfam_names_local[clan_cols].drop_duplicates("clan"),
        on="clan",
        how="left",
    )

    # 7. Render outputs
    if print_tables:
        print(f"Largest families in '{view_name}' population (n={len(df)}):")
        top_families = pfam_counts_merged.sort_values("count", ascending=False).head(k)

        # Apply link formatting cleanly
        formatted_families = add_pfam_links_vectorized(top_families.copy())
        if formatted_families is not None and not formatted_families.empty:
            display_html_df(formatted_families)
        else:
            print("No family records to display.")

        print(f"\nLargest superfamilies in '{view_name}' population:")
        top_clans = clan_counts_merged.sort_values("count", ascending=False).head(k)

        formatted_clans = add_pfam_links_vectorized(top_clans.copy())
        if formatted_clans is not None and not formatted_clans.empty:
            display_html_df(formatted_clans)
        else:
            print("No superfamily records to display.")

    # Explicitly flush stdout so ipywidgets output capture displays immediately
    sys.stdout.flush()
    return None

def print_pfam_names_table_with_counts(pfam_names_df: pd.DataFrame, plot_df: pd.DataFrame):
    """
    Generates a metadata table where count values match the sunburst plot 1:1,
    properly preserving Pfam subcategories (like 'Other Pfam') within their specific Clan.
    """
    # 1. Aggregate counts by BOTH clan_clean and pfam_clean (matching the sunburst path)
    counts = (
        plot_df.groupby(['clan_clean', 'pfam_clean'], as_index=False)
        .size()
        .rename(columns={'clan_clean': 'clan', 'pfam_clean': 'pfam', 'size': 'count'})
    )

    # 2. Separate standard Pfams from synthetic ones ('Other Pfam', 'Unassigned Pfam', 'Other Clan')
    known_metadata = pfam_names_df.copy()

    # Merge known Pfam metadata on 'pfam' ID
    table_df = counts.merge(known_metadata, on='pfam', how='left', suffixes=('', '_meta'))

    # 3. Fill in metadata for synthetic/collapsed rows dynamically
    # If 'clan_meta' exists from pfam_names_df, prioritize it; otherwise use the collapsed 'clan' from plot_df
    if 'clan_meta' in table_df.columns:
        table_df['clan'] = table_df['clan_meta'].fillna(table_df['clan'])
        table_df.drop(columns=['clan_meta'], inplace=True)

    # Fill metadata fields for 'Other Pfam', 'Unassigned Pfam', etc.
    table_df['clan_short'] = table_df['clan_short'].fillna('N/A')

    def fill_short(row):
        if pd.notna(row['short']):
            return row['short']
        return row['pfam']

    def fill_name(row):
        if pd.notna(row['name']):
            return row['name']
        if row['pfam'] == 'Other Pfam':
            return f"Aggregated Pfams (rank > k) in {row['clan']}"
        if row['pfam'] == 'Unassigned Pfam':
            return "Unannotated Pfam Domain"
        return row['pfam']

    table_df['short'] = table_df.apply(fill_short, axis=1)
    table_df['name'] = table_df.apply(fill_name, axis=1)

    # 4. Arrange columns and sort
    cols = ["pfam", "count", "clan", "clan_short", "short", "name"]
    table_df = table_df[[c for c in cols if c in table_df.columns]]

    return table_df

def plot_pfam_clan_sunburst(df: pd.DataFrame, k: int = 15):
    """
    Groups dataframe into Clan -> Pfam hierarchy based on top-k rank cutoffs.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Must contain ['pfam', 'pfam_rank', 'clan', 'clan_rank'].
    k : int
        Ranks > k are collapsed into 'Other Clan' and 'Other Pfam'.
    """
    plot_df = df.copy()

    # Fill NaN values explicitly
    plot_df['clan'] = plot_df['clan'].fillna('Unassigned Clan').astype(str)
    plot_df['pfam'] = plot_df['pfam'].fillna('Unassigned Pfam').astype(str)
    plot_df['clan_rank'] = plot_df['clan_rank'].fillna(999999)
    plot_df['pfam_rank'] = plot_df['pfam_rank'].fillna(999999)

    # 1. Apply rank cutoffs
    plot_df['clan_clean'] = plot_df.apply(
        lambda r: r['clan'] if r['clan'] == 'Unassigned Clan' 
        else (r['clan'] if r['clan_rank'] <= k else 'Other Clan'), axis=1
    )

    plot_df['pfam_clean'] = plot_df.apply(
        lambda r: r['pfam'] if r['pfam'] == 'Unassigned Pfam' 
        else (r['pfam'] if r['pfam_rank'] <= k else 'Other Pfam'), axis=1
    )

    # 2. Aggregate counts across collapsed categories
    df_grouped = (
        plot_df.groupby(['clan_clean', 'pfam_clean'], as_index=False)
        .size()
        .rename(columns={'size': 'count'})
    )

    # 3. Create Plotly Sunburst
    fig = px.sunburst(
        df_grouped,
        path=['clan_clean', 'pfam_clean'],
        values='count',
        title=f'Pfam Hierarchy (Top-{k} Cutoff)',
        color='clan_clean',
        color_discrete_sequence=px.colors.qualitative.Pastel
    )

    fig.update_traces(
        textinfo="label+value",
        texttemplate="%{label}<br>n: %{value}"
    )

    fig.update_layout(margin=dict(t=40, l=0, r=0, b=0))

    # Pass the transformed plot dataframe so the table can compute exact matching counts
    return fig, plot_df

#===========

def pfam_distribution_analysis(
    df: pd.DataFrame,
    k: int = 5,
    variables: list = None,
    group:str = 'clan',
    alpha: float = 1e-5,
    min_effect: dict = None,
    show_plots: bool = True,
    show_table: bool = False,
    return_results: bool = True,
) -> pd.DataFrame | None:
    """
    Full Pfam/clan distribution analysis in one call.

    Steps
    -----
    1. Violin plots of variable distributions by top-k clan
    2. Mann-Whitney U tests with FDR correction
    3. Clan discovery wizard commentary

    Parameters
    ----------
    df           : view dataframe with pfam, clan, pfam_rank, clan_rank
                   and the variables to test
    k            : top-k groups to highlight
    group        : categories to test - default 'clan'
    variables    : variables to test — default ['z_score',
                   'sequence_identity', 'query_coverage']
    alpha        : significance threshold (default 1e-5 — tighter than
                   conventional 0.05 due to large-N sensitivity)
    min_effect   : optional {variable: min median_diff} for effect-size
                   gating — see group_distribution_test()
    show_plots   : display violin plots inline (default True)
    show_table   : display stats table (default False)
    return_results: return the results DataFrame (default True)

    Returns
    -------
    results DataFrame if return_results=True, else None
    """
    if variables is None:
        variables = ['z_score', 'sequence_identity', 'query_coverage']

    # --- step 1: violin plots ---
    if show_plots:
        # plot k pfam-ranks
        df_violin = df[df['clan_rank']<=k]
        fig = plot_violins( df_violin, group, ['z_score', 'query_coverage', 'sequence_identity'], cols=3)       
        plt.close(fig)

    # --- step 2: statistical tests ---
    results = group_distribution_test(
        df         = df,
        k          = k,
        variables  = variables,
        alpha      = alpha,
        min_effect = min_effect,
    )

    if results.empty:
        print("No groups with sufficient data to test.")
        return results if return_results else None

    if show_table:
        print(f"\n{results['significant'].sum()} significant associations "
          f"(q < {alpha}):\n")
        display(results[results['significant']].reset_index(drop=True))

    # --- step 3: wizard commentary ---
    print()
    print(clan_discovery_wizard(results, df=df))

    return results if return_results else None

##################################

def logo_significance_test(
    df: pd.DataFrame,
    pileup_col: str = 'sequ_pileup',
    group_col: str = 'clan',
    group_value: str = 'CL0209',
    query_length: int = None,
    gap_char: str = '.',
    min_occupancy: float = 0.1,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Per-position Mann-Whitney U test comparing residue conservation
    between a group (e.g. home clan) and the background (all others).

    For each position i, the test variable is binary:
        1 if the target has a non-gap residue at position i
        0 if it has a gap

    This tests whether the group is more consistently COVERED at
    each position — i.e. whether structural alignment places a
    residue there more reliably in the group than in background.

    For the dominant-residue identity (not just presence), a
    second pass compares the frequency of the consensus residue
    between group and background.

    Returns
    -------
    DataFrame with one row per position:
        position, occupancy_group, occupancy_background,
        consensus_residue, freq_group, freq_background,
        u_stat, p_value, q_value, significant, stars
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    import numpy as np

    valid    = df.dropna(subset=[pileup_col, group_col])
    in_group = valid[valid[group_col] == group_value][pileup_col].tolist()
    out_group= valid[valid[group_col] != group_value][pileup_col].tolist()

    if not in_group or not out_group:
        raise ValueError(f"No sequences for group {group_value!r} "
                         f"or background in column {group_col!r}")

    L = query_length or len(in_group[0])

    def to_matrix(seqs):
        return np.array([[0 if c == gap_char else 1
                          for c in s[:L]] for s in seqs],
                        dtype=np.float32)

    M_grp = to_matrix(in_group)   # (n_group, L)
    M_bg  = to_matrix(out_group)  # (n_bg, L)

    records  = []
    p_values = []
    positions = []

    for i in range(L):
        grp_occ = M_grp[:, i]
        bg_occ  = M_bg[:, i]

        # skip positions with very low overall occupancy — noise
        total_occ = (grp_occ.sum() + bg_occ.sum()) / (len(grp_occ) + len(bg_occ))
        if total_occ < min_occupancy:
            continue

        stat, p = mannwhitneyu(grp_occ, bg_occ, alternative='greater')

        # consensus residue in group at this position
        if grp_occ.sum() > 0:
            residues_grp = [s[i] for s in in_group
                           if i < len(s) and s[i] != gap_char]
            residues_bg  = [s[i] for s in out_group
                           if i < len(s) and s[i] != gap_char]
            if residues_grp:
                from collections import Counter
                consensus = Counter(residues_grp).most_common(1)[0][0]
                freq_grp  = residues_grp.count(consensus) / max(len(grp_occ), 1)
                freq_bg   = (residues_bg.count(consensus) /
                             max(len(bg_occ), 1)) if residues_bg else 0.0
            else:
                consensus = gap_char
                freq_grp  = freq_bg = 0.0
        else:
            consensus = gap_char
            freq_grp  = freq_bg = 0.0

        records.append({
            'position':             i,
            'occupancy_group':      float(grp_occ.mean()),
            'occupancy_background': float(bg_occ.mean()),
            'consensus_residue':    consensus,
            'freq_group':           round(freq_grp, 3),
            'freq_background':      round(freq_bg,  3),
            'u_stat':               round(stat, 1),
            'p_value':              p,
        })
        p_values.append(p)
        positions.append(i)

    if not records:
        return pd.DataFrame()

    _, q_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')

    for rec, q in zip(records, q_values):
        rec['q_value']     = round(q, 6)
        rec['significant'] = bool(q < alpha)
        rec['stars']       = (
            '***' if q < 1e-5 else
            '**'  if q < 1e-3 else
            '*'   if q < alpha else ''
        )

    return pd.DataFrame(records)

##################################
