"""
domain_architecture.py

Domain architecture diagrams — one horizontal track per representative
target — with a shared, plotter-agnostic data-prep layer and two thin
rendering functions (matplotlib, Plotly).
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from ..analysis.data_preview import show_heading

# =======================================================================
# Plotter-agnostic data model
# =======================================================================


@dataclass
class DomainSegment:
    identifier: str
    name: str
    start: float
    end: float


@dataclass
class ArchitectureRow:
    row_idx: int
    target_id: str
    label: str
    target_length: float
    segments: list[DomainSegment]
    highlight: Optional[tuple[float, float]]  # (t_left, t_right) or None


@dataclass
class PreparedArchitectureData:
    rows: list[ArchitectureRow]
    all_domain_ids: list[str]
    colors: dict[str, str]  # identifier -> hex color
    name_lookup: dict[str, str]
    display_labels: dict[str, str]
    title: str
    pfam_col: str
    empty_message: Optional[str] = None


# =======================================================================
# Shared helpers
# =======================================================================


def _normalize_domains(value: Any) -> list:
    """Coerce a cell's domain annotations into a list of dicts.

    Handles a native list/tuple, a stringified list, or anything else -> [].
    """
    if isinstance(value, (list, tuple)):
        return list(value)

    if isinstance(value, str):
        value = value.strip()
        if not value:
            return []

        try:
            parsed = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return []

        if isinstance(parsed, (list, tuple)):
            return list(parsed)

        return []

    return []


def _domain_identifier(domain: Any, pfam_col: str) -> Optional[str]:
    """Extract the identifier for one domain dict.

    pfam_domains -> domain["pfam_id"]
    clan_domains -> domain["clan"], falling back to domain["pfam_id"]
    """
    if not isinstance(domain, dict):
        return None

    if pfam_col == "clan_domains":
        identifier = domain.get("clan") or domain.get("pfam_id")
    else:
        identifier = domain.get("pfam_id")

    if identifier is None:
        return None

    identifier = str(identifier).strip()
    return identifier or None


def _architecture_key(domains: list, pfam_col: str) -> Optional[tuple]:
    """Canonical architecture key: sorted, deduplicated identifier set.

    Domain order and copy number are deliberately ignored, so
    [PF1, PF2, PF2] and [PF2, PF1] are the same architecture.
    """
    identifiers = [
        ident
        for d in domains
        if (ident := _domain_identifier(d, pfam_col)) is not None
    ]

    if not identifiers:
        return None

    return tuple(sorted(set(identifiers)))


def _format_display_label(
    identifier: str,
    name: str,
    clan: Optional[str],
    pfam_col: str,
) -> str:
    """Build the legend/hover text for one domain identifier."""
    if pfam_col == "pfam_domains" and clan:
        head = f"{identifier} / {clan}"
    else:
        head = identifier

    return f"{head}  {name}" if name else head


def _get_color_palette(
    all_domain_ids: list[str],
) -> dict[str, str]:
    """Return shared tab20 hex colors."""
    from matplotlib import colormaps
    from matplotlib.colors import to_hex

    cmap = colormaps["tab20"]

    return {
        identifier: to_hex(cmap(i % cmap.N))
        for i, identifier in enumerate(all_domain_ids)
    }


# =======================================================================
# Shared data preparation
# =======================================================================


def prepare_domain_architecture_data(
    df: pd.DataFrame,
    pfam_col: str = "pfam_domains",
    max_targets: int = 50,
    group_by_architecture: bool = True,
    representative_sort: str = "z_score",
    sort_col: Optional[str] = None,
    ascending: bool = False,
) -> PreparedArchitectureData:
    """Build plotter-agnostic data for the domain architecture diagram.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.

    pfam_col : str, default="pfam_domains"
        Column containing domain annotations.

    max_targets : int, default=50
        Maximum number of representative targets to include.

    group_by_architecture : bool, default=True
        If True, select one representative target per unique architecture.

    representative_sort : {"z_score", "architecture_count"}, default="z_score"
        Ordering of representatives when group_by_architecture=True.

        "z_score":
            Select the highest-z_score target for each architecture and
            sort those representatives by z_score descending.

        "architecture_count":
            Select the highest-z_score target for each architecture and
            sort representatives by the number of targets having that
            architecture, descending. Z-score is used as a tie-breaker.

    sort_col : str or None, default=None
        Sort column when group_by_architecture=False.

    ascending : bool, default=False
        Sort direction when sort_col is used.

    Returns
    -------
    PreparedArchitectureData
        Plotter-agnostic representation consumed by both renderers.

    Raises
    ------
    KeyError
        If pfam_col is missing, or if architecture grouping requires
        z_score but the column is absent.

    ValueError
        If representative_sort is invalid.
    """

    if pfam_col not in df.columns:
        raise KeyError(
            "Column not found: " + repr(pfam_col)
        )

    if representative_sort not in {
        "z_score",
        "architecture_count",
    }:
        raise ValueError(
            "representative_sort must be either "
            "'z_score' or 'architecture_count'"
        )

    # -------------------------------------------------------------------
    # Normalize domain annotations
    # -------------------------------------------------------------------

    plot_df = df.dropna(subset=[pfam_col]).copy()

    if plot_df.empty:
        return PreparedArchitectureData(
            rows=[],
            all_domain_ids=[],
            colors={},
            name_lookup={},
            display_labels={},
            title="",
            pfam_col=pfam_col,
            empty_message="No targets with domain annotations",
        )

    plot_df[pfam_col] = plot_df[pfam_col].apply(
        _normalize_domains
    )

    plot_df = plot_df[
        plot_df[pfam_col].apply(lambda d: len(d) > 0)
    ].copy()

    if plot_df.empty:
        return PreparedArchitectureData(
            rows=[],
            all_domain_ids=[],
            colors={},
            name_lookup={},
            display_labels={},
            title="",
            pfam_col=pfam_col,
            empty_message="No valid domain annotations",
        )

    # -------------------------------------------------------------------
    # Build canonical architecture key
    # -------------------------------------------------------------------

    plot_df["_architecture_key"] = plot_df[pfam_col].apply(
        lambda domains: _architecture_key(domains, pfam_col)
    )

    plot_df = plot_df.dropna(
        subset=["_architecture_key"]
    ).copy()

    # -------------------------------------------------------------------
    # Representative selection and ordering
    # -------------------------------------------------------------------

    if group_by_architecture:

        if "z_score" not in plot_df.columns:
            raise KeyError(
                "group_by_architecture=True requires "
                "a 'z_score' column."
            )

        # Numeric representation of Z-score.
        plot_df["_plot_z_score"] = pd.to_numeric(
            plot_df["z_score"],
            errors="coerce",
        )

        # Count the number of original targets belonging to each
        # architecture BEFORE selecting representatives.
        plot_df["_architecture_count"] = (
            plot_df["_architecture_key"]
            .map(
                plot_df["_architecture_key"].value_counts()
            )
        )

        # ---------------------------------------------------------------
        # First determine the representative for each architecture.
        #
        # The representative is ALWAYS the highest-Z target.
        # ---------------------------------------------------------------

        plot_df = plot_df.sort_values(
            "_plot_z_score",
            ascending=False,
            kind="mergesort",
            na_position="last",
        )

        plot_df = plot_df.drop_duplicates(
            subset=["_architecture_key"],
            keep="first",
        )

        # ---------------------------------------------------------------
        # Then order the representatives according to the requested
        # display mode.
        # ---------------------------------------------------------------

        if representative_sort == "z_score":

            plot_df = plot_df.sort_values(
                "_plot_z_score",
                ascending=False,
                kind="mergesort",
                na_position="last",
            )

        else:  # representative_sort == "architecture_count"

            plot_df = plot_df.sort_values(
                [
                    "_architecture_count",
                    "_plot_z_score",
                ],
                ascending=[
                    False,
                    False,
                ],
                kind="mergesort",
                na_position="last",
            )

    elif sort_col is not None and sort_col in plot_df.columns:

        plot_df = plot_df.sort_values(
            sort_col,
            ascending=ascending,
            kind="mergesort",
        )

    # -------------------------------------------------------------------
    # Apply max_targets AFTER representative selection.
    #
    # Otherwise duplicate architectures among the original rows could
    # consume the limit and prevent other architectures from appearing.
    # -------------------------------------------------------------------

    plot_df = plot_df.head(max_targets).copy()

    # -------------------------------------------------------------------
    # Domain IDs and colors
    # -------------------------------------------------------------------

    all_domain_ids = sorted(
        {
            ident
            for domains in plot_df[pfam_col]
            for domain in domains
            if (
                ident := _domain_identifier(
                    domain,
                    pfam_col,
                )
            ) is not None
        }
    )

    colors = _get_color_palette(
        all_domain_ids
    )

    # -------------------------------------------------------------------
    # Domain names / clan lookup
    # -------------------------------------------------------------------

    name_lookup: dict[str, str] = {}
    clan_lookup: dict[str, str] = {}

    for domains in plot_df[pfam_col]:

        for domain in domains:

            if not isinstance(domain, dict):
                continue

            ident = _domain_identifier(
                domain,
                pfam_col,
            )

            if ident is None:
                continue

            name = domain.get("name") or ""

            if ident not in name_lookup and name:
                name_lookup[ident] = str(name)

            if (
                pfam_col == "pfam_domains"
                and ident not in clan_lookup
            ):
                clan = domain.get("clan")

                if clan:
                    clan_lookup[ident] = str(
                        clan
                    ).strip()

    display_labels = {
        ident: _format_display_label(
            ident,
            name_lookup.get(ident, ""),
            clan_lookup.get(ident),
            pfam_col,
        )
        for ident in all_domain_ids
    }

    # -------------------------------------------------------------------
    # Build ArchitectureRow objects
    # -------------------------------------------------------------------

    rows: list[ArchitectureRow] = []

    for row_idx, (_, row) in enumerate(
        plot_df.iterrows()
    ):

        target_id = str(
            row.get(
                "target_id",
                row.name,
            )
        )

        domains = row[pfam_col]

        # ---------------------------------------------------------------
        # Target label
        # ---------------------------------------------------------------

        if (
            group_by_architecture
            and "_plot_z_score" in plot_df.columns
        ):
            z_score = row["_plot_z_score"]

            if pd.notna(z_score):
                label = (
                    f"{target_id}  "
                    f"(Z={float(z_score):.1f})"
                )
            else:
                label = target_id

        elif (
            not group_by_architecture
            and sort_col is not None
            and sort_col in plot_df.columns
        ):
            try:
                label = (
                    f"{target_id}  "
                    f"({sort_col}="
                    f"{float(row[sort_col]):.1f})"
                    f")"
                )
            except (TypeError, ValueError):
                label = target_id

        else:
            label = target_id

        # ---------------------------------------------------------------
        # Target length
        # ---------------------------------------------------------------

        target_length = row.get(
            "target_length",
            None,
        )

        if (
            target_length is None
            or pd.isna(target_length)
        ):
            ends = [
                float(d["end"])
                for d in domains
                if (
                    isinstance(d, dict)
                    and d.get("end") is not None
                )
            ]

            target_length = (
                max(ends)
                if ends
                else 100
            )

        else:
            target_length = float(
                target_length
            )

        # ---------------------------------------------------------------
        # Domain segments
        # ---------------------------------------------------------------

        segments: list[DomainSegment] = []

        for domain in domains:

            if not isinstance(domain, dict):
                continue

            ident = _domain_identifier(
                domain,
                pfam_col,
            )

            if ident is None:
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

            if end - start <= 0:
                continue

            segments.append(
                DomainSegment(
                    identifier=ident,
                    name=name_lookup.get(
                        ident,
                        "",
                    ),
                    start=start,
                    end=end,
                )
            )

        # ---------------------------------------------------------------
        # Alignment-region highlight
        # ---------------------------------------------------------------

        highlight = None

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
                t_left = float(t_left)
                t_right = float(t_right)

                if t_right >= t_left:
                    highlight = (
                        t_left,
                        t_right,
                    )

            except (
                TypeError,
                ValueError,
            ):
                pass

        rows.append(
            ArchitectureRow(
                row_idx=row_idx,
                target_id=target_id,
                label=label,
                target_length=target_length,
                segments=segments,
                highlight=highlight,
            )
        )

    # -------------------------------------------------------------------
    # Title
    # -------------------------------------------------------------------

    if group_by_architecture:

        if representative_sort == "z_score":
            title = (
                "Domain architecture "
                "(one highest-Z representative "
                "per composition, sorted by Z-score)"
            )

        else:
            title = (
                "Domain architecture "
                "(one highest-Z representative "
                "per composition, sorted by composition count)"
            )

    elif sort_col is not None:

        title = (
            f"Domain architecture "
            f"(sorted by {sort_col})"
        )

    else:

        title = "Domain architecture"

    # -------------------------------------------------------------------
    # Return the prepared object
    # -------------------------------------------------------------------

    return PreparedArchitectureData(
        rows=rows,
        all_domain_ids=all_domain_ids,
        colors=colors,
        name_lookup=name_lookup,
        display_labels=display_labels,
        title=title,
        pfam_col=pfam_col,
    )


# =======================================================================
# Matplotlib renderer
# =======================================================================


def plot_domain_architecture_matplotlib(
    df: pd.DataFrame,
    pfam_col: str = "pfam_domains",
    max_targets: int = 50,
    group_by_architecture: bool = True,
    representative_sort: str = "z_score",
    sort_col: Optional[str] = None,
    ascending: bool = False,
):
    """Domain architecture diagram, matplotlib backend.

    See `prepare_domain_architecture_data` for parameter documentation.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    data = prepare_domain_architecture_data(
        df,
        pfam_col=pfam_col,
        max_targets=max_targets,
        group_by_architecture=group_by_architecture,
        representative_sort=representative_sort,
        sort_col=sort_col,
        ascending=ascending,
    )

    if data.empty_message:

        fig, ax = plt.subplots(
            figsize=(8, 2)
        )

        ax.text(
            0.5,
            0.5,
            data.empty_message,
            ha="center",
            va="center",
        )

        ax.axis("off")

        return fig

    fig_height = max(
        2,
        0.45 * len(data.rows) + 1.5,
    )

    fig, ax = plt.subplots(
        figsize=(12, fig_height)
    )

    for row in data.rows:

        # ---------------------------------------------------------------
        # Base sequence line
        # ---------------------------------------------------------------

        ax.plot(
            [0, row.target_length],
            [row.row_idx, row.row_idx],
            color="#dddddd",
            linewidth=1.5,
            zorder=1,
        )

        # ---------------------------------------------------------------
        # Domain rectangles
        # ---------------------------------------------------------------

        for seg in row.segments:

            ax.barh(
                row.row_idx,
                seg.end - seg.start,
                left=seg.start,
                height=0.6,
                color=data.colors[
                    seg.identifier
                ],
                edgecolor="white",
                linewidth=0.5,
                zorder=2,
            )

        # ---------------------------------------------------------------
        # Alignment-region highlight
        # ---------------------------------------------------------------

        if row.highlight is not None:

            t_left, t_right = row.highlight

            ax.add_patch(
                mpatches.Rectangle(
                    (
                        t_left,
                        row.row_idx - 0.4,
                    ),
                    t_right - t_left + 1,
                    0.8,
                    fill=False,
                    edgecolor="#e05c5c",
                    linewidth=1.5,
                    linestyle="--",
                    zorder=3,
                )
            )

    ax.set_yticks(
        range(len(data.rows))
    )

    ax.set_yticklabels(
        [
            row.label
            for row in data.rows
        ],
        fontsize=7,
    )

    # Data rows are already in the desired display order.
    # Invert so row 0 appears at the top.
    ax.invert_yaxis()

    ax.set_xlabel(
        "Residue position"
    )

    ax.set_title(
        data.title
    )

    legend_patches = [
        mpatches.Patch(
            color=data.colors[ident],
            label=data.display_labels[ident],
        )
        for ident in data.all_domain_ids
    ]

    if legend_patches:

        ax.legend(
            handles=legend_patches,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=8,
        )

    plt.tight_layout()

    return fig


# =======================================================================
# Plotly renderer
# =======================================================================


def plot_domain_architecture_plotly(
    df: pd.DataFrame,
    pfam_col: str = "pfam_domains",
    max_targets: int = 50,
    group_by_architecture: bool = True,
    representative_sort: str = "z_score",
    sort_col: Optional[str] = None,
    ascending: bool = False,
    show_legend: bool = True,
):
    """Domain architecture diagram, Plotly backend.

    Same prepared data model as the matplotlib renderer.

    Parameters
    ----------
    show_legend : bool
        If True, add a static color-key legend via invisible dummy
        traces. Hover text remains available on the domain segments.
    """
    import plotly.graph_objects as go

    data = prepare_domain_architecture_data(
        df,
        pfam_col=pfam_col,
        max_targets=max_targets,
        group_by_architecture=group_by_architecture,
        representative_sort=representative_sort,
        sort_col=sort_col,
        ascending=ascending,
    )

    if data.empty_message:

        fig = go.Figure()

        fig.add_annotation(
            text=data.empty_message,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
        )

        fig.update_xaxes(
            visible=False
        )

        fig.update_yaxes(
            visible=False
        )

        return fig

    fig = go.Figure()

    # -------------------------------------------------------------------
    # Base sequence line
    # -------------------------------------------------------------------

    base_x = []
    base_y = []

    for row in data.rows:

        base_x += [
            0,
            row.target_length,
            None,
        ]

        base_y += [
            row.row_idx,
            row.row_idx,
            None,
        ]

    fig.add_trace(
        go.Scatter(
            x=base_x,
            y=base_y,
            mode="lines",
            line=dict(
                color="#dddddd",
                width=1.5,
            ),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    # -------------------------------------------------------------------
    # Domain segments
    # -------------------------------------------------------------------

    bar_x = []
    bar_base = []
    bar_y = []
    bar_color = []
    bar_hover = []

    for row in data.rows:

        for seg in row.segments:

            # display_labels combines identifier/clan and name.
            head = data.display_labels[
                seg.identifier
            ]

            if (
                seg.name
                and head.endswith(seg.name)
            ):
                head = head[
                    :-len(seg.name)
                ].rstrip()

            hover_lines = [
                f"<b>{head}</b>"
            ]

            if seg.name:
                hover_lines.append(
                    seg.name
                )

            hover_lines.append(
                f"Position: "
                f"{int(seg.start)}-"
                f"{int(seg.end)}"
            )

            hover_lines.append(
                f"Target: {row.target_id}"
            )

            bar_x.append(
                seg.end - seg.start
            )

            bar_base.append(
                seg.start
            )

            bar_y.append(
                row.row_idx
            )

            bar_color.append(
                data.colors[
                    seg.identifier
                ]
            )

            bar_hover.append(
                "<br>".join(
                    hover_lines
                )
            )

    fig.add_trace(
        go.Bar(
            x=bar_x,
            base=bar_base,
            y=bar_y,
            orientation="h",
            width=0.6,
            marker=dict(
                color=bar_color,
                line=dict(
                    color="white",
                    width=0.5,
                ),
            ),
            hovertext=bar_hover,
            hovertemplate=(
                "%{hovertext}<extra></extra>"
            ),
            showlegend=False,
        )
    )

    # -------------------------------------------------------------------
    # Alignment-region highlight
    # -------------------------------------------------------------------

    for row in data.rows:

        if row.highlight is None:
            continue

        t_left, t_right = row.highlight

        fig.add_shape(
            type="rect",
            x0=t_left,
            x1=t_right + 1,
            y0=row.row_idx - 0.4,
            y1=row.row_idx + 0.4,
            line=dict(
                color="#e05c5c",
                width=1.5,
                dash="dash",
            ),
            fillcolor="rgba(0,0,0,0)",
        )

    # -------------------------------------------------------------------
    # Optional static legend
    # -------------------------------------------------------------------

    if show_legend:

        for ident in data.all_domain_ids:

            fig.add_trace(
                go.Bar(
                    x=[None],
                    y=[None],
                    orientation="h",
                    marker=dict(
                        color=data.colors[
                            ident
                        ]
                    ),
                    name=data.display_labels[
                        ident
                    ],
                    showlegend=True,
                    hoverinfo="skip",
                )
            )

    # -------------------------------------------------------------------
    # Layout
    # -------------------------------------------------------------------

    fig_height = max(
        300,
        45 * len(data.rows) + 150,
    )

    fig.update_layout(
        title=data.title,
        xaxis_title="Residue position",
        yaxis=dict(
            tickvals=list(
                range(len(data.rows))
            ),
            ticktext=[
                row.label
                for row in data.rows
            ],
            autorange="reversed",
            tickfont=dict(size=10),
        ),
        height=fig_height,
        width=1000,
        barmode="overlay",
        plot_bgcolor="white",
        showlegend=show_legend,
        margin=dict(
            l=10,
            r=10,
            t=50,
            b=40,
        ),
    )

    fig.update_xaxes(
        showgrid=False,
        zeroline=False,
    )

    fig.update_yaxes(
        showgrid=False,
        zeroline=False,
    )

    display(fig)

def plot_domain_architecture(
    df: pd.DataFrame,
    title: str = None,
    max_targets: int = 10,
    pfam_col: str = "pfam_domains",
    renderer: str = "plotly",
    show_legend: bool = True,
    group_by_architecture: bool = True,
    representative_sort: str = "z_score",
    sort_col: str = None,
    ascending: bool = False,
    **kwargs,
):
    """Pure function to plot domain architectures directly from a DataFrame."""
    if title:
        show_heading(title)
    elif "view_name" in df.attrs:
        show_heading(f"{df.attrs['view_name']} (n={len(df)})")

    if renderer == "plotly":
        plot_domain_architecture_plotly(
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

    elif renderer == "matplotlib":
        plot_domain_architecture_matplotlib(
            df=df,
            pfam_col=pfam_col,
            max_targets=max_targets,
            group_by_architecture=group_by_architecture,
            representative_sort=representative_sort,
            sort_col=sort_col,
            ascending=ascending,
            **kwargs,
        )

    else:
        raise ValueError(f"Unknown renderer: {renderer!r}")
