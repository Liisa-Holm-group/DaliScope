import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
import matplotlib.patches as mpatches

import math
import seaborn as sns

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from daliscope.mechanics.metrics import compute_family_presence
from daliscope.analysis.wrappers import apply_rank_pooling

def volcano_plots(results: pd.DataFrame,
                  variables: list = None,
                  alpha: float = 0.05,
                  x_col: str = 'median_diff',
                  label_col: str = 'group') -> go.Figure:
    """
    Interactive Plotly volcano plots — one panel per variable.
    Hover shows all stats; significant points labelled on hover,
    not cluttered on the plot itself.

    Parameters
    ----------
    results   : output of group_distribution_test()
    variables : which variables to plot — defaults to all unique in results
    alpha     : significance threshold line on y axis
    x_col     : 'median_diff' or 'fold_change'
    label_col : column used for hover label
    """
    if variables is None:
        variables = results['variable'].unique().tolist()

    n_vars = len(variables)
    fig = make_subplots(
        rows=1, cols=n_vars,
        subplot_titles=[v.replace('_', ' ') for v in variables],
        horizontal_spacing=0.08,
    )

    colors      = {'pfam': '#185FA5', 'clan': '#1D9E75'}
    sig_thresh  = -np.log10(alpha)
    legend_seen = set()

    for col_idx, var in enumerate(variables, start=1):
        sub = results[results['variable'] == var].copy()
        sub = sub.dropna(subset=[x_col, 'q_value'])
        sub['neg_log10_q'] = -np.log10(sub['q_value'].clip(lower=1e-10))

        for group_type, grp in sub.groupby('group_type'):
            color = colors.get(group_type, '#888888')

            for is_sig, marker_size, marker_opacity, suffix in [
                (False, 8,  0.4, ' (ns)'),
                (True,  12, 0.9, ''),
            ]:
                subset = grp[grp['significant'] == is_sig]
                if subset.empty:
                    continue

                legend_label = f"{group_type}{suffix}"
                show_legend  = legend_label not in legend_seen
                if show_legend:
                    legend_seen.add(legend_label)

                hover = [
                    f"<b>{row[label_col]}</b><br>"
                    f"group_type: {row['group_type']}<br>"
                    f"{x_col}: {row[x_col]:.3f}<br>"
                    f"median: {row['median']:.3f}<br>"
                    f"bg median: {row['bg_median']:.3f}<br>"
                    f"fold_change: {row.get('fold_change', np.nan):.3f}<br>"
                    f"n: {row['n']}<br>"
                    f"q: {row['q_value']:.4f}<br>"
                    f"{'✓ significant' if row['significant'] else 'ns'}"
                    for _, row in subset.iterrows()
                ]

                fig.add_trace(
                    go.Scatter(
                        x=subset[x_col],
                        y=subset['neg_log10_q'],
                        mode='markers+text' if is_sig else 'markers',
                        text=subset[label_col] if is_sig else None,
                        textposition='top center',
                        textfont=dict(size=8, color=color),
                        marker=dict(
                            size=marker_size,
                            color=color,
                            opacity=marker_opacity,
                            line=dict(width=1, color='white')
                            if is_sig else dict(width=0),
                        ),
                        name=legend_label,
                        legendgroup=legend_label,
                        showlegend=show_legend,
                        hovertext=hover,
                        hoverinfo='text',
                    ),
                    row=1, col=col_idx,
                )

        # significance threshold line
        x_range = sub[x_col]
        x_pad   = (x_range.max() - x_range.min()) * 0.1 or 0.1
        fig.add_shape(
            type='line',
            x0=x_range.min() - x_pad, x1=x_range.max() + x_pad,
            y0=sig_thresh, y1=sig_thresh,
            line=dict(color='red', dash='dash', width=1),
            row=1, col=col_idx,
        )
        # zero line
        fig.add_shape(
            type='line',
            x0=0, x1=0,
            y0=0, y1=sub['neg_log10_q'].max() * 1.1,
            line=dict(color='grey', dash='dot', width=1),
            row=1, col=col_idx,
        )

        # axis labels
        fig.update_xaxes(title_text=x_col.replace('_', ' '),
                         row=1, col=col_idx)
        fig.update_yaxes(title_text='-log₁₀(q)' if col_idx == 1 else '',
                         row=1, col=col_idx)

    fig.update_layout(
        title='Volcano plots — group vs background',
        height=500,
        width=420 * n_vars,
        template='simple_white',
        legend=dict(
            title='group type',
            font=dict(size=10),
        ),
        hoverlabel=dict(bgcolor='white', font_size=11),
    )

    return fig

def plot_violins(df, x_category, y_variables, cols=2, palette='viridis'):
    """
    Creates a grid of violin plots for multiple Y variables against one X category.
    """
    n_vars = len(y_variables)
    rows = math.ceil(n_vars / cols)

    # Create the figure and axes
    fig, axes = plt.subplots(rows, cols, figsize=(9, 3), squeeze=False)
    axes_flat = axes.flatten()

    for i, y_var in enumerate(y_variables):
        sns.violinplot(
            data=df,
            x=x_category,
            y=y_var,
            ax=axes_flat[i], # This tells Seaborn WHERE to plot
            hue=x_category,
            palette=palette,
            inner='quartile',
            legend=False
        )

        # Formatting individual subplots
        axes_flat[i].set_title(f'{y_var} by {x_category}', fontsize=12)
        axes_flat[i].tick_params(axis='x', rotation=45)
        axes_flat[i].set_xlabel('') # Hide x-label to save space

    # Hide empty subplots
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis('off')

    plt.tight_layout()
    plt.show()

import ast
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plot_domain_architecture(
    df: pd.DataFrame,
    pfam_col: str = "pfam_domains",
    max_targets: int = 50,
    group_by_architecture: bool = True,
    sort_col: str = None,
    ascending: bool = False,
) -> "plt.Figure":
    """
    Domain architecture diagram — one horizontal track per representative target.

    Architecture definition
    -----------------------
    An architecture is defined as the unique set of domain identifiers:

        PF1, PF2, PF2, PF3

    is equivalent to:

        PF3, PF1, PF2

    because domain order and copy number are ignored.

    Representative selection
    ------------------------
    When group_by_architecture=True, one representative target is selected
    for each unique architecture. The representative is the row with the
    highest z_score.

    The resulting representatives are then sorted by z_score descending.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.

    pfam_col : str
        Column containing domain annotations.

        For pfam_domains, identifiers are taken from:

            domain["pfam_id"]

        For clan_domains, identifiers are taken from:

            domain["clan"]

        if available.

    max_targets : int
        Maximum number of representative targets to plot.

    group_by_architecture : bool
        If True, select one highest-Z-score representative per architecture.

        If False, do not perform representative selection.

    sort_col : str or None
        Sorting column when group_by_architecture=False.

    ascending : bool
        Sorting direction when sort_col is used.
    """

    # -----------------------------------------------------------------
    # Validate input
    # -----------------------------------------------------------------

    if pfam_col not in df.columns:
        raise KeyError(
            "Column not found: "
            + repr(pfam_col)
        )


    plot_df = (
        df
        .dropna(
            subset=[pfam_col]
        )
        .copy()
    )


    if plot_df.empty:

        fig, ax = plt.subplots(
            figsize=(8, 2)
        )

        ax.text(
            0.5,
            0.5,
            "No targets with domain annotations",
            ha="center",
            va="center",
        )

        ax.axis("off")

        return fig


    # -----------------------------------------------------------------
    # Normalize domain annotations
    #
    # Handles:
    #
    #   list of dictionaries
    #
    # and stringified lists such as:
    #
    #   "[{'pfam_id': 'PF00959', ...}]"
    # -----------------------------------------------------------------

    def normalize_domains(value):

        if isinstance(
            value,
            (list, tuple)
        ):

            return list(value)


        if isinstance(
            value,
            str
        ):

            value = value.strip()

            if not value:
                return []

            try:

                parsed = ast.literal_eval(
                    value
                )

            except (
                ValueError,
                SyntaxError,
            ):

                return []


            if isinstance(
                parsed,
                (list, tuple)
            ):

                return list(parsed)


            return []


        return []


    plot_df[pfam_col] = (
        plot_df[pfam_col]
        .apply(
            normalize_domains
        )
    )


    plot_df = plot_df[
        plot_df[pfam_col].apply(
            lambda domains:
                len(domains) > 0
        )
    ].copy()


    if plot_df.empty:

        fig, ax = plt.subplots(
            figsize=(8, 2)
        )

        ax.text(
            0.5,
            0.5,
            "No valid domain annotations",
            ha="center",
            va="center",
        )

        ax.axis("off")

        return fig


    # -----------------------------------------------------------------
    # Domain identifier
    #
    # pfam_domains:
    #
    #     PF00959
    #
    # clan_domains:
    #
    # Prefer "clan", otherwise allow "pfam_id" as fallback so the
    # function remains tolerant of partially normalized data.
    # -----------------------------------------------------------------

    def domain_identifier(domain):

        if not isinstance(
            domain,
            dict
        ):

            return None


        if pfam_col == "clan_domains":

            identifier = (
                domain.get("clan")
                or domain.get("pfam_id")
            )

        else:

            identifier = (
                domain.get("pfam_id")
            )


        if identifier is None:
            return None


        identifier = str(
            identifier
        ).strip()


        if not identifier:
            return None


        return identifier


    # -----------------------------------------------------------------
    # Canonical architecture key
    #
    # IMPORTANT:
    #
    # sorted(set(...))
    #
    # deliberately ignores:
    #
    #   - sequence order
    #   - copy number
    #
    # Thus:
    #
    #   [PF1, PF2, PF2]
    #
    # and:
    #
    #   [PF2, PF1]
    #
    # have the same architecture.
    # -----------------------------------------------------------------

    def architecture_key(domains):

        identifiers = []

        for domain in domains:

            identifier = (
                domain_identifier(
                    domain
                )
            )

            if identifier is not None:

                identifiers.append(
                    identifier
                )


        if not identifiers:
            return None


        return tuple(
            sorted(
                set(
                    identifiers
                )
            )
        )


    plot_df["_architecture_key"] = (
        plot_df[pfam_col]
        .apply(
            architecture_key
        )
    )


    plot_df = (
        plot_df
        .dropna(
            subset=[
                "_architecture_key"
            ]
        )
        .copy()
    )


    # -----------------------------------------------------------------
    # Representative selection
    # -----------------------------------------------------------------

    if group_by_architecture:

        if "z_score" not in plot_df.columns:

            raise KeyError(
                "group_by_architecture=True requires "
                "a 'z_score' column."
            )


        plot_df["_plot_z_score"] = (
            pd.to_numeric(
                plot_df["z_score"],
                errors="coerce",
            )
        )


        # Missing Z-scores are allowed but rank below real Z-scores.
        #
        # Stable mergesort ensures deterministic selection when two rows
        # have identical architecture and identical Z-score.
        plot_df = (
            plot_df
            .sort_values(
                "_plot_z_score",
                ascending=False,
                kind="mergesort",
                na_position="last",
            )
        )


        # Because the table is already sorted highest -> lowest,
        # drop_duplicates keeps the highest-Z representative.
        plot_df = (
            plot_df
            .drop_duplicates(
                subset=[
                    "_architecture_key"
                ],
                keep="first",
            )
        )


        # Explicitly retain descending Z-score order after selection.
        plot_df = (
            plot_df
            .sort_values(
                "_plot_z_score",
                ascending=False,
                kind="mergesort",
                na_position="last",
            )
        )


    elif (
        sort_col is not None
        and sort_col in plot_df.columns
    ):

        plot_df = (
            plot_df
            .sort_values(
                sort_col,
                ascending=ascending,
                kind="mergesort",
            )
        )


    # -----------------------------------------------------------------
    # Apply limit AFTER representative selection.
    #
    # This is important: otherwise the first 50 raw rows could contain
    # many duplicates and suppress architectures whose representative
    # occurs later.
    # -----------------------------------------------------------------

    plot_df = (
        plot_df
        .head(
            max_targets
        )
        .copy()
    )


    # -----------------------------------------------------------------
    # Collect displayed domain identifiers
    # -----------------------------------------------------------------

    all_domain_ids = sorted(
        {
            identifier

            for domains in plot_df[pfam_col]

            for domain in domains

            for identifier in [
                domain_identifier(
                    domain
                )
            ]

            if identifier is not None
        }
    )


    cmap = plt.get_cmap(
        "tab20"
    )


    colors = {
        identifier:
            cmap(
                index % cmap.N
            )

        for index, identifier
        in enumerate(
            all_domain_ids
        )
    }


    # -----------------------------------------------------------------
    # Name lookup
    # -----------------------------------------------------------------

    name_lookup = {}


    for domains in plot_df[pfam_col]:

        for domain in domains:

            if not isinstance(
                domain,
                dict
            ):
                continue


            identifier = (
                domain_identifier(
                    domain
                )
            )


            if identifier is None:
                continue


            name = (
                domain.get("name")
                or ""
            )


            if (
                identifier not in name_lookup
                and name
            ):

                name_lookup[
                    identifier
                ] = str(name)


    # -----------------------------------------------------------------
    # Legend labels
    # -----------------------------------------------------------------

    def make_legend_label(identifier):

        name = (
            name_lookup.get(
                identifier,
                ""
            )
        )


        if name:

            return (
                str(identifier)
                + "  "
                + name
            )


        return str(
            identifier
        )


    # -----------------------------------------------------------------
    # Figure
    # -----------------------------------------------------------------

    fig_height = max(
        2,
        0.45 * len(plot_df) + 1.5,
    )


    fig, ax = plt.subplots(
        figsize=(
            12,
            fig_height,
        )
    )


    # -----------------------------------------------------------------
    # Draw rows
    # -----------------------------------------------------------------

    for row_idx, (_, row) in enumerate(
        plot_df.iterrows()
    ):

        domains = row[
            pfam_col
        ]


        # Full protein length.
        target_length = row.get(
            "target_length",
            None,
        )


        if (
            target_length is None
            or pd.isna(
                target_length
            )
        ):

            ends = []

            for domain in domains:

                if isinstance(
                    domain,
                    dict
                ):

                    end = domain.get(
                        "end"
                    )

                    if end is not None:

                        try:

                            ends.append(
                                float(end)
                            )

                        except (
                            TypeError,
                            ValueError,
                        ):

                            pass


            target_length = (
                max(ends)
                if ends
                else 100
            )


        # Base sequence line.
        ax.plot(
            [0, target_length],
            [row_idx, row_idx],
            color="#dddddd",
            linewidth=1.5,
            zorder=1,
        )


        # Domain rectangles.
        for domain in domains:

            if not isinstance(
                domain,
                dict
            ):
                continue


            identifier = (
                domain_identifier(
                    domain
                )
            )


            if identifier is None:
                continue


            try:

                start = float(
                    domain["start"]
                )

                end = float(
                    domain["end"]
                )

            except (
                KeyError,
                TypeError,
                ValueError,
            ):

                continue


            width = (
                end - start
            )


            if width <= 0:
                continue


            ax.barh(
                row_idx,
                width,
                left=start,
                height=0.6,
                color=colors[
                    identifier
                ],
                edgecolor="white",
                linewidth=0.5,
                zorder=2,
            )


        # Alignment-region highlight.
        t_left = row.get(
            "target_left",
            None,
        )

        t_right = row.get(
            "target_right",
            None,
        )


        if (
            t_left is not None
            and t_right is not None
            and not pd.isna(t_left)
            and not pd.isna(t_right)
        ):

            try:

                t_left = float(
                    t_left
                )

                t_right = float(
                    t_right
                )


                if t_right >= t_left:

                    rect = (
                        mpatches.Rectangle(
                            (
                                t_left,
                                row_idx - 0.4,
                            ),
                            (
                                t_right
                                - t_left
                                + 1
                            ),
                            0.8,
                            fill=False,
                            edgecolor="#e05c5c",
                            linewidth=1.5,
                            linestyle="--",
                            zorder=3,
                        )
                    )


                    ax.add_patch(
                        rect
                    )

            except (
                TypeError,
                ValueError,
            ):

                pass


    # -----------------------------------------------------------------
    # Row labels
    # -----------------------------------------------------------------

    labels = []


    for _, row in plot_df.iterrows():

        target_id = str(
            row.get(
                "target_id",
                row.name,
            )
        )


        if (
            group_by_architecture
            and "_plot_z_score"
            in plot_df.columns
        ):

            z_score = row[
                "_plot_z_score"
            ]


            if pd.notna(
                z_score
            ):

                label = (
                    target_id
                    + "  (Z="
                    + format(
                        float(z_score),
                        ".1f",
                    )
                    + ")"
                )

            else:

                label = target_id


        elif (
            not group_by_architecture
            and sort_col is not None
            and sort_col in plot_df.columns
        ):

            value = row[
                sort_col
            ]


            try:

                label = (
                    target_id
                    + "  ("
                    + str(sort_col)
                    + "="
                    + format(
                        float(value),
                        ".1f",
                    )
                    + ")"
                )

            except (
                TypeError,
                ValueError,
            ):

                label = target_id

        else:

            label = target_id


        labels.append(
            label
        )


    ax.set_yticks(
        range(
            len(plot_df)
        )
    )


    ax.set_yticklabels(
        labels,
        fontsize=7,
    )


    # -----------------------------------------------------------------
    # Critical ordering step
    #
    # Row 0 is the highest Z-score because representatives were sorted
    # descending. Inverting the axis puts row 0 visually at the TOP.
    # -----------------------------------------------------------------

    ax.invert_yaxis()


    ax.set_xlabel(
        "Residue position"
    )


    if group_by_architecture:

        title = (
            "Domain architecture "
            "(one highest-Z representative per composition)"
        )

    elif sort_col is not None:

        title = (
            "Domain architecture "
            "(sorted by "
            + str(sort_col)
            + ")"
        )

    else:

        title = (
            "Domain architecture"
        )


    ax.set_title(
        title
    )


    # -----------------------------------------------------------------
    # Legend
    # -----------------------------------------------------------------

    legend_patches = [
        mpatches.Patch(
            color=colors[identifier],
            label=make_legend_label(
                identifier
            ),
        )

        for identifier
        in all_domain_ids
    ]


    if legend_patches:

        ax.legend(
            handles=legend_patches,
            bbox_to_anchor=(
                1.02,
                1,
            ),
            loc="upper left",
            fontsize=8,
        )


    plt.tight_layout()


    return fig


def old_plot_domain_architecture(df: pd.DataFrame,
                             pfam_col: str = 'pfam_domains',
                             max_targets: int = 50,
                             group_by_architecture: bool = True,
                             sort_col: str = None,
                             ascending: bool = False,
                             ) -> 'plt.Figure':
    """
    Domain architecture diagram — one horizontal track per target.

    Row ordering
    ------------
    group_by_architecture=True (default):
        rows sorted so identical architectures are grouped together
        (sort_col and ascending are ignored in this mode)

    group_by_architecture=False:
        rows kept in df's existing order, UNLESS sort_col is given,
        in which case rows are sorted by that column (e.g. 'z_score')
        — useful when the caller has preselected a specific set of
        targets and wants them ordered by some external criterion
        (e.g. z_score descending) rather than grouped by architecture

    Parameters
    ----------
    sort_col   : column to sort by when group_by_architecture=False.
                 If None, df's existing row order is preserved as-is.
    ascending  : sort direction for sort_col (default False, i.e.
                 highest z_score first)
    """

    plot_df = df.dropna(subset=[pfam_col]).head(max_targets).copy()

    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(0.5, 0.5, "No targets with Pfam domain annotations",
               ha='center', va='center')
        ax.axis('off')
        return fig

    if group_by_architecture:
        plot_df['architecture_key'] = plot_df[pfam_col].apply(
            lambda doms: tuple(d['pfam_id'] for d in doms)
        )
        plot_df = plot_df.sort_values('architecture_key')
    elif sort_col is not None:
        plot_df = plot_df.sort_values(sort_col, ascending=ascending)
    # else: keep df's existing row order untouched

    all_pfam_ids = sorted(set(
        d['pfam_id'] for doms in plot_df[pfam_col] for d in doms
    ))
    cmap   = plt.get_cmap('tab20')
    colors = {pid: cmap(i % cmap.N) for i, pid in enumerate(all_pfam_ids)}

    # name and clan lookups straight from the data — no separate
    # pfam_names argument needed, both already resolved at silver
    # build time by build_pfam_domains()
    name_lookup = {
        d['pfam_id']: d['name']
        for doms in plot_df[pfam_col] for d in doms
    }
    clan_lookup = {
        d['pfam_id']: d.get('clan')
        for doms in plot_df[pfam_col] for d in doms
    }

    def make_legend_label(pid: str) -> str:
        name = name_lookup.get(pid, '')
        if pfam_col == 'pfam_domains':
            clan = clan_lookup.get(pid)
            if clan:
                return f"{pid}/{clan}  {name}"
            return f"{pid}  {name}"
        else:
            # clan_domains: pid IS the clan id, no second id to append
            return f"{pid}  {name}"

    fig, ax = plt.subplots(figsize=(12, 0.3 * len(plot_df) + 1))
    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        doms = row[pfam_col]
        # Use actual full protein length for the sequence line (fallback to last domain end if missing)
        target_length = row.get('target_length', None)
        if target_length is None or pd.isna(target_length):
            target_length = max((d['end'] for d in doms), default=100)

        # 1. Base sequence line drawn to full C-terminal sequence length
        ax.plot([0, target_length], [row_idx, row_idx],
            color='#dddddd', linewidth=1.5, zorder=1)

        for dom in doms:
            ax.barh(row_idx, dom['end'] - dom['start'],
                   left=dom['start'], height=0.6,
                   color=colors[dom['pfam_id']],
                   edgecolor='white', linewidth=0.5, zorder=2)

        # 3. Transparent Highlight Box per Target (using target_left and target_right)
        t_left = row.get("target_left", None)
        t_right = row.get("target_right", None)

        if t_left is not None and t_right is not None and not pd.isna(t_left) and not pd.isna(t_right):
            # Calculate full span width (inclusive coordinates)
            width = (t_right - t_left) + 1

            # Single box spanning the whole aligned region for this row
            rect = mpatches.Rectangle(
                (t_left, row_idx - 0.4),  # (x, y) bottom-left origin
                width,                    # total span from target_left to target_right
                0.8,                      # height framing the 0.6 domain bars
                fill=False,               # No fill color inside the box
                edgecolor='#e05c5c',      # Red/coral dashed outline (or '#d97706' for amber)
                linewidth=1.5,
                linestyle='--',
                zorder=3                  # Layers directly over domain bars
            )
            ax.add_patch(rect)

    # row labels: include sort_col value when relevant, so the
    # ordering criterion is visible alongside each row
    if not group_by_architecture and sort_col is not None and sort_col in plot_df.columns:
        labels = [
            f"{row['target_id']}  ({sort_col}={row[sort_col]:.1f})"
            for _, row in plot_df.iterrows()
        ]
    else:
        labels = plot_df['target_id'].tolist()

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Residue position')

    if group_by_architecture:
        title = 'Domain architecture (grouped by composition)'
    elif sort_col is not None:
        title = f'Domain architecture (sorted by {sort_col})'
    else:
        title = 'Domain architecture'
    ax.set_title(title)

    legend_patches = [
        mpatches.Patch(color=colors[pid], label=make_legend_label(pid))
        for pid in all_pfam_ids
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1),
             loc='upper left', fontsize=8)

    plt.tight_layout()
    return fig


def plot_family_presence(population_df: pd.DataFrame, family_col: str, rank_col: str,
                          target_ids, id_col: str = "target_id",
                          k: int = 10, show_zero_present: bool = True,
                          ax=None, title=None):
    """
    Pure rendering function: diverging bar chart of family presence/absence.
 
    Families are selected by rank (pfam_rank/clan_rank — ranked by each
    family's best z_score hit): rank<=k families keep their own bar,
    everything ranked below k is pooled into 'Other', and unranked rows
    into 'Unassigned'. The k named families are sorted by pct_present
    descending (most enriched first); 'Other' and 'Unassigned' are always
    drawn last, since they're aggregate buckets rather than real families.
 
    show_zero_present : if False, named families with zero present-count
        are hidden from the chart (large clans especially can push k into
        dozens of near-empty bars, which is mostly noise) — 'Other' and
        'Unassigned' are always kept regardless, since they summarize the
        rest of the population rather than representing a single family.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 7))
    else:
        fig = ax.figure

    group_series = apply_rank_pooling(population_df, rank_col=rank_col, label_col=family_col, k=k)
    presence = compute_family_presence(population_df, group_series, target_ids, id_col=id_col)

    if presence.empty:
        ax.text(0.5, 0.5, "No data for this category/view combination",
                ha="center", va="center")
        return fig, ax

    pooled_labels = [lbl for lbl in ("Other", "Unassigned") if lbl in presence.index]
    named = presence.drop(index=pooled_labels, errors="ignore").sort_values("pct_present", ascending=False)

    n_hidden = 0
    if not show_zero_present:
        n_hidden = int((named["n_present"] == 0).sum())
        named = named[named["n_present"] > 0]

    shown = pd.concat([named, presence.loc[pooled_labels]])

    ax.bar(shown.index, shown["pct_present"], color="royalblue", label="Present")
    ax.bar(shown.index, shown["pct_absent"], color="indianred", label="Absent")

    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylim(-115, 115)  # room for n= labels above the bars
    ax.set_ylabel("Percentage (%)")

    default_title = f"Presence of target set across top {len(named)} '{family_col}' families (by z_score rank)"
    if n_hidden:
        default_title += f"\n({n_hidden} zero-present families hidden)"
    ax.set_title(title or default_title)

    for i, fam in enumerate(shown.index):
        p_val = shown["pct_present"].loc[fam]
        a_val = shown["pct_absent"].loc[fam]
        total_n = shown["n_total"].loc[fam]
        n_present = shown["n_present"].loc[fam]

        # present/total counts, shown regardless of bar height — unlike the
        # inline percentage labels below, which only fit when the segment
        # is tall enough. This is the number that matters most for large
        # families with a small percentage present.
        ax.text(i, 105, f"{n_present}/{total_n}", ha="center", va="bottom",
                fontsize=9, color="gray", weight="bold")
        if p_val > 5:
            ax.text(i, p_val / 2, f"{p_val:.1f}%", ha="center", color="white", weight="bold")
        if abs(a_val) > 5:
            ax.text(i, a_val / 2, f"{abs(a_val):.1f}%", ha="center", color="white", weight="bold")

    ax.set_xticks(range(len(shown.index)))
    ax.set_xticklabels(shown.index, rotation=45, ha="right")
    fig.tight_layout()

    # Legend placed BELOW the chart (outside the axes, via bbox_to_anchor),
    # rather than ax.legend()'s default 'best' placement — which doesn't
    # reliably avoid the n= count labels sitting at the top of every bar
    # and can end up overlapping the rightmost one.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2, frameon=False)
    fig.subplots_adjust(bottom=0.32)  # reserve room for the legend + rotated x-tick labels

    return fig, ax

def plot_dali_hits(df):
    _df = df[['z_score','query_coverage','sequence_identity']].dropna()
    x = _df["z_score"]
    y = _df["query_coverage"]
    seq = _df["sequence_identity"]

    # 1. SET UP THE GRID DESIGN
    fig = plt.figure(figsize=(10, 6))
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
    # 6. AXES LIMITS + MAIN FORMATTING
    ax.set_xlabel("z_score")
    ax.set_ylabel("query_coverage")

    ax_histx.set_title("Naive hit classification (heuristic thresholds)", fontweight='bold', pad=10)

    ax.set_xlim(plot_x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

    plt.show()
    return plt
