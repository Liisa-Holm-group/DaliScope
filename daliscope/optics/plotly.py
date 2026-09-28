import plotly.express as px
import plotly.io as pio
pio.renderers.default = "iframe"
from IPython.display import display, HTML

import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import math

import plotly.graph_objects as go

def sunburst_trace(df, pfam_col="pfam", clan_col="clan", unclassified_label="(no clan)"):
    """Build just the sunburst trace (not a standalone figure) so it can be
    dropped into a subplot grid."""
    d = df[[pfam_col, clan_col]].copy()
    d[clan_col] = d[clan_col].fillna(unclassified_label)
    counts = d.groupby([clan_col, pfam_col]).size().reset_index(name="count")
    single_fig = px.sunburst(counts, path=[clan_col, pfam_col], values="count")
    return single_fig.data[0]

def sunburst_side_by_side(dataframes, ncol=3):
    """dataframes: dict of {title: df}"""
    n = len(dataframes)
    nrows = math.ceil(n / ncol)

    # specs must be a 2D grid: nrows lists, each containing ncol dicts
    specs = [[{"type": "domain"} for _ in range(ncol)] for _ in range(nrows)]

    # pad titles to the full grid size, or subplot_titles will misalign
    # once nrows*ncol > n
    titles = list(dataframes.keys()) + [""] * (nrows * ncol - n)

    fig = make_subplots(
        rows=nrows, cols=ncol,
        specs=specs,
        subplot_titles=titles,
    )

    for i, (title, d) in enumerate(dataframes.items()):
        row = i // ncol + 1
        col = i % ncol + 1
        fig.add_trace(sunburst_trace(d), row=row, col=col)

    fig.update_layout(title_text="Pfam/clan counts, side by side")
    fig.show()

def plotly_overview(df, x='z_score', y='alignment_length',
                    color='pfam', marker='clan', size=None,
                    hover_cols=None, height=600):

    hover_cols = hover_cols or ['target_id', 'description']
    hover_cols = [c for c in hover_cols if c in df.columns]

    plot_df    = df.copy()  ## BAD!

    # guard against NaN rmsd values (in copy)
    numeric_cols = plot_df.select_dtypes(include='number').columns
    plot_df[numeric_cols] = plot_df[numeric_cols].fillna(0)

    color_col  = color  if (color  and color  in plot_df.columns) else None
    symbol_col = marker if (marker and marker in plot_df.columns) else None
    size_col   = size   if (size   and size   in plot_df.columns) else None

    # detect whether color is categorical or continuous
    is_categorical = (
        color_col is not None and
        plot_df[color_col].dtype == object or
        str(plot_df[color_col].dtype) == 'category'
    ) if color_col else True

    fig = px.scatter(
        plot_df,
        x          = x,
        y          = y,
        color      = color_col,
        symbol     = symbol_col,
        size       = size_col,
        hover_data = hover_cols,
        template   = "simple_white",
        height     = height,
        title      = f"{x} vs {y}  (n={len(df)})",
        color_discrete_sequence = px.colors.qualitative.Dark24
                                  if is_categorical else None,
        color_continuous_scale  = "Viridis"
                                  if not is_categorical else None,
        render_mode = "webgl",
    )

    if size_col is None:
        fig.update_traces(marker=dict(size=12))
    fig.update_traces(marker=dict(opacity=0.7))

    if is_categorical:
        # categorical: colorbar irrelevant, legend handles groups
        fig.update_coloraxes(showscale=False)
        fig.update_layout(
            legend=dict(
                title     = color_col or "",
                font_size = 9,
                x=1.02, y=1,
                xanchor="left",
                yanchor="top",
            ),
            margin=dict(r=180),
        )
    else:
        # continuous: move colorbar below the plot, out of legend area
        fig.update_coloraxes(
            showscale     = True,
            colorbar=dict(
                orientation = "h",        # horizontal bar
                y           = 0.1,       # below the plot
                x           = 0.8,
                xanchor     = "center",
                len         = 0.3,        # 30% of plot width
                thickness   = 12,
                title       = dict(text=color_col, side="bottom"),
            )
        )
        fig.update_layout(
            legend=dict(
                title     = symbol_col or "",
                font_size = 9,
                x=1.02, y=1,
                xanchor="left",
                yanchor="top",
            ),
            margin=dict(r=180, b=200),   # extra bottom margin for colorbar
        )
        # Increase data point size (default is usually around 6)

    fig.show()



def scatter_plot(df, x, y, color='_color_label', symbol=None,
                 size=None, hover_cols=None, height=500, width=800):
    hover_cols = hover_cols or ['target_id', 'description']
    hover_cols = [c for c in hover_cols if c in df.columns]
    size_max   = 20 if size is not None else None

    fig = px.scatter(
        df, x=x, y=y,
        color      = color,
        symbol     = symbol,
        size       = size,
        size_max   = size_max,
        hover_data = hover_cols,
        color_discrete_sequence = px.colors.qualitative.Dark24,
        template   = 'simple_white',
        height     = height,
        title      = f"{x} vs {y}",
        render_mode = "webgl",
    )
    if size is None:
        fig.update_traces(marker=dict(size=12, opacity=0.7))
    else:
        fig.update_traces(marker=dict(opacity=0.7))

    # --- figure size
    fig.update_layout(height=height, width=width)

    # --- ADD EYE-CANDY TOOLBAR CONFIGURATION ---
    plotly_config = {
        'toImageButtonOptions': {
            'format': 'png',
            'filename': f"daliscope_{x}_vs_{y}",
            'height': 500,             # Sharp download size
            'width': 500,
            'scale': 2                 # High-resolution retina scale
        },
        'displaylogo': False,          # Clean up toolbar by hiding Plotly logo
        'modeBarButtonsToRemove': [    # Strip cluttered/redundant actions
            'lasso2d', 
            'select2d', 
#            'autoScale2d'
        ]
    }
    fig.show()
    return None

    # --- deleted bwlow ===#
    # Pass the config dictionary to the HTML exporter
    html_content = fig.to_html(
        full_html=False,
        include_plotlyjs='cdn', 
        config=plotly_config
    )

    display(HTML(html_content))
    return None

######

def plot_sequence_profile_heatmap(df: pd.DataFrame):
    """
    Plots a sequence profile heatmap where the index contains the positions.
    Positions are on the vertical axis (categorical), amino acids on the horizontal axis.
    """
    # 1. Define the BLAST-style property gradient (Polar/Charged -> Hydrophobic/Aromatic)
    aa_property_order = [
        'D', 'E', 'R', 'K', 'H',  # Charged / Highly Polar
        'N', 'Q', 'S', 'T',       # Polar Uncharged
        'G', 'P', 'A',            # Special / Small Aliphatic
        'C', 'V', 'M', 'I', 'L',  # Hydrophobic Aliphatic
        'Y', 'F', 'W'             # Hydrophobic Aromatic
    ]

    # Only keep amino acids that are actually columns in this subset dataframe
    amino_acids = [aa for aa in aa_property_order if aa in df.columns]

    # 2. Extract positions from the index as string categories
    positions = df.index.astype(str).tolist()

    # 3. Pull out the values using our sorted amino acid columns
    matrix_values = df[amino_acids].values

    # 4. Generate the Heatmap
    fig = go.Figure(data=go.Heatmap(
        z=matrix_values,
        x=amino_acids,
        y=positions,
        colorscale='Viridis',
        hovertemplate=(
            "Position: %{y}<br>"
            "Amino Acid: %{x}<br>"
            "Value: %{z:.3f}<extra></extra>"
        )
    ))

    # 5. Layout optimizations for discrete categorical positions
    fig.update_layout(
        title='Sequence Profile Heatmap',
        xaxis_title='Amino Acid (Polar ➔ Hydrophobic Gradient)',
        yaxis_title='Sequence Position (Index)',
        xaxis={'side': 'top'},        # Keep labels at the top for quick scanning
        height=max(400, len(df) * 20)  # Keeps rows from compressing if you have many rows
    )

    # Forces categorical layout and reads top-to-bottom
    fig.update_yaxes(type='category', autorange='reversed')

    fig.show()

    return(positions, matrix_values)

##########

import ast

import pandas as pd
import plotly.graph_objects as go
from matplotlib import colormaps
from matplotlib.colors import to_hex

def plot_domain_architecture_plotly(
    df: pd.DataFrame,
    pfam_col: str = "pfam_domains",
    max_targets: int = 50,
    group_by_architecture: bool = True,
    sort_col: str = None,
    ascending: bool = False,
    show_legend: bool = False,
    fallback_name_col: str = "pfam_domains",
) -> "go.Figure":
    """
    Domain architecture diagram — one horizontal track per representative
    target. Plotly version of plot_domain_architecture(): same data model
    (architecture definition, representative selection, row ordering),
    hover text on each domain segment instead of a static legend.

    NOTE: the normalization / architecture-key / representative-selection
    logic below is copied verbatim from the matplotlib version rather than
    imported from a shared helper. If this goes into production alongside
    the matplotlib version, factor that shared logic out into one place —
    otherwise a future bugfix to (say) domain_identifier() only lands in
    whichever version someone remembers to edit, and the two renderers
    silently drift apart on what counts as "the same architecture".

    Parameters
    ----------
    Same as plot_domain_architecture(), plus:

    show_legend : bool
        If True, add a static color-key legend (built from invisible
        dummy traces) alongside the hover text. Off by default, since
        the point of this version is that hover replaces the legend.

    fallback_name_col : str
        Only used when pfam_col == "clan_domains". Some domains in the
        clan-grouped data have no clan (their identifier falls back to
        the raw pfam_id) and no "name" field either. For those, the
        display name is looked up by pfam_id from this column instead
        (default "pfam_domains", i.e. the per-Pfam-family annotations).
        Set to None to disable the fallback lookup.
    """

    # -----------------------------------------------------------------
    # Validate input
    # -----------------------------------------------------------------

    if pfam_col not in df.columns:
        raise KeyError("Column not found: " + repr(pfam_col))

    plot_df = df.dropna(subset=[pfam_col]).copy()

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No targets with domain annotations",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        return fig

    # -----------------------------------------------------------------
    # Normalize domain annotations (identical to matplotlib version)
    # -----------------------------------------------------------------

    def normalize_domains(value):
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

    plot_df[pfam_col] = plot_df[pfam_col].apply(normalize_domains)
    plot_df = plot_df[plot_df[pfam_col].apply(lambda d: len(d) > 0)].copy()

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No valid domain annotations",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        return fig

    # -----------------------------------------------------------------
    # Domain identifier (identical to matplotlib version)
    # -----------------------------------------------------------------

    def domain_identifier(domain):
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

    # -----------------------------------------------------------------
    # Canonical architecture key (identical to matplotlib version)
    # -----------------------------------------------------------------

    def architecture_key(domains):
        identifiers = [
            domain_identifier(d) for d in domains
            if domain_identifier(d) is not None
        ]
        if not identifiers:
            return None
        return tuple(sorted(set(identifiers)))

    plot_df["_architecture_key"] = plot_df[pfam_col].apply(architecture_key)
    plot_df = plot_df.dropna(subset=["_architecture_key"]).copy()

    # -----------------------------------------------------------------
    # Representative selection (identical to matplotlib version)
    # -----------------------------------------------------------------

    if group_by_architecture:
        if "z_score" not in plot_df.columns:
            raise KeyError(
                "group_by_architecture=True requires a 'z_score' column."
            )

        plot_df["_plot_z_score"] = pd.to_numeric(
            plot_df["z_score"], errors="coerce"
        )

        plot_df = plot_df.sort_values(
            "_plot_z_score", ascending=False,
            kind="mergesort", na_position="last",
        )
        plot_df = plot_df.drop_duplicates(
            subset=["_architecture_key"], keep="first"
        )
        plot_df = plot_df.sort_values(
            "_plot_z_score", ascending=False,
            kind="mergesort", na_position="last",
        )

    elif sort_col is not None and sort_col in plot_df.columns:
        plot_df = plot_df.sort_values(
            sort_col, ascending=ascending, kind="mergesort"
        )

    plot_df = plot_df.head(max_targets).copy()

    # -----------------------------------------------------------------
    # Colors (matplotlib tab20, converted to hex for Plotly)
    # -----------------------------------------------------------------

    all_domain_ids = sorted({
        identifier
        for domains in plot_df[pfam_col]
        for domain in domains
        for identifier in [domain_identifier(domain)]
        if identifier is not None
    })

    cmap = colormaps["tab20"]
    colors = {
        identifier: to_hex(cmap(i % cmap.N))
        for i, identifier in enumerate(all_domain_ids)
    }

    # -----------------------------------------------------------------
    # Name lookup (identical logic to matplotlib version)
    # -----------------------------------------------------------------

    name_lookup = {}
    for domains in plot_df[pfam_col]:
        for domain in domains:
            if not isinstance(domain, dict):
                continue
            identifier = domain_identifier(domain)
            if identifier is None:
                continue
            name = domain.get("name") or ""
            if identifier not in name_lookup and name:
                name_lookup[identifier] = str(name)

    # -----------------------------------------------------------------
    # Row labels (identical logic to matplotlib version)
    # -----------------------------------------------------------------

    labels = []
    for _, row in plot_df.iterrows():
        target_id = str(row.get("target_id", row.name))

        if group_by_architecture and "_plot_z_score" in plot_df.columns:
            z_score = row["_plot_z_score"]
            label = (
                f"{target_id}  (Z={float(z_score):.1f})"
                if pd.notna(z_score) else target_id
            )
        elif (not group_by_architecture and sort_col is not None
              and sort_col in plot_df.columns):
            try:
                label = f"{target_id}  ({sort_col}={float(row[sort_col]):.1f})"
            except (TypeError, ValueError):
                label = target_id
        else:
            label = target_id

        labels.append(label)

    # -----------------------------------------------------------------
    # Build traces
    #
    # Numeric row positions (0, 1, 2, ...) are used throughout, with
    # tick labels applied afterwards via yaxis.tickvals/ticktext - this
    # mirrors the matplotlib version's set_yticks/set_yticklabels split
    # exactly, and sidesteps any risk of duplicate-label collisions
    # that a categorical (label-as-position) y-axis would have.
    # -----------------------------------------------------------------

    fig = go.Figure()

    # --- base sequence line, one Scatter trace for all rows,
    #     segments separated by None so they don't connect ---
    base_x, base_y = [], []
    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        domains = row[pfam_col]

        target_length = row.get("target_length", None)
        if target_length is None or pd.isna(target_length):
            ends = [
                float(d["end"]) for d in domains
                if isinstance(d, dict) and d.get("end") is not None
            ]
            target_length = max(ends) if ends else 100

        base_x += [0, target_length, None]
        base_y += [row_idx, row_idx, None]

    fig.add_trace(go.Scatter(
        x=base_x, y=base_y,
        mode="lines",
        line=dict(color="#dddddd", width=1.5),
        hoverinfo="skip",
        showlegend=False,
    ))

    # --- domain segments, one Bar trace, array-based ---
    bar_x, bar_base, bar_y, bar_color, bar_hover = [], [], [], [], []

    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        target_id = str(row.get("target_id", row.name))
        domains = row[pfam_col]

        for domain in domains:
            if not isinstance(domain, dict):
                continue
            identifier = domain_identifier(domain)
            if identifier is None:
                continue
            try:
                start = float(domain["start"])
                end = float(domain["end"])
            except (KeyError, TypeError, ValueError):
                continue

            width = end - start
            if width <= 0:
                continue

            name = name_lookup.get(identifier, "")
            hover_lines = [f"<b>{identifier}</b>"]
            if name:
                hover_lines.append(name)
            hover_lines.append(f"Position: {int(start)}-{int(end)}")
            hover_lines.append(f"Target: {target_id}")

            bar_x.append(width)
            bar_base.append(start)
            bar_y.append(row_idx)
            bar_color.append(colors[identifier])
            bar_hover.append("<br>".join(hover_lines))

    fig.add_trace(go.Bar(
        x=bar_x, base=bar_base, y=bar_y,
        orientation="h",
        width=0.6,
        marker=dict(color=bar_color, line=dict(color="white", width=0.5)),
        hovertext=bar_hover,
        hovertemplate="%{hovertext}<extra></extra>",
        showlegend=False,
    ))

    # --- alignment-region highlight, dashed outline shapes ---
    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        t_left = row.get("target_left", None)
        t_right = row.get("target_right", None)

        if t_left is None or t_right is None:
            continue
        if pd.isna(t_left) or pd.isna(t_right):
            continue

        try:
            t_left = float(t_left)
            t_right = float(t_right)
        except (TypeError, ValueError):
            continue

        if t_right < t_left:
            continue

        fig.add_shape(
            type="rect",
            x0=t_left, x1=t_right + 1,
            y0=row_idx - 0.4, y1=row_idx + 0.4,
            line=dict(color="#e05c5c", width=1.5, dash="dash"),
            fillcolor="rgba(0,0,0,0)",
        )

    # --- optional static legend, invisible dummy traces ---
    if show_legend:
        for identifier in all_domain_ids:
            name = name_lookup.get(identifier, "")
            legend_label = f"{identifier}  {name}" if name else identifier
            fig.add_trace(go.Bar(
                x=[None], y=[None],
                orientation="h",
                marker=dict(color=colors[identifier]),
                name=legend_label,
                showlegend=True,
                hoverinfo="skip",
            ))

    # -----------------------------------------------------------------
    # Layout
    # -----------------------------------------------------------------

    if group_by_architecture:
        title = "Domain architecture (one highest-Z representative per composition)"
    elif sort_col is not None:
        title = f"Domain architecture (sorted by {sort_col})"
    else:
        title = "Domain architecture"

    fig_height = max(300, 45 * len(plot_df) + 150)

    fig.update_layout(
        title=title,
        xaxis_title="Residue position",
        yaxis=dict(
            tickvals=list(range(len(plot_df))),
            ticktext=labels,
            autorange="reversed",  # row 0 (highest Z) at the top
            tickfont=dict(size=10),
        ),
        height=fig_height,
        width=1000,
        barmode="overlay",
        plot_bgcolor="white",
        showlegend=show_legend,
        margin=dict(l=10, r=10, t=50, b=40),
    )

    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=False, zeroline=False)

    return fig

def old_plot_domain_architecture_plotly(
    df: pd.DataFrame,
    pfam_col: str = "pfam_domains",
    max_targets: int = 50,
    group_by_architecture: bool = True,
    sort_col: str = None,
    ascending: bool = False,
    show_legend: bool = False,
) -> "go.Figure":
    """
    Domain architecture diagram — one horizontal track per representative
    target. Plotly version of plot_domain_architecture(): same data model
    (architecture definition, representative selection, row ordering),
    hover text on each domain segment instead of a static legend.

    NOTE: the normalization / architecture-key / representative-selection
    logic below is copied verbatim from the matplotlib version rather than
    imported from a shared helper. If this goes into production alongside
    the matplotlib version, factor that shared logic out into one place —
    otherwise a future bugfix to (say) domain_identifier() only lands in
    whichever version someone remembers to edit, and the two renderers
    silently drift apart on what counts as "the same architecture".

    Parameters
    ----------
    Same as plot_domain_architecture(), plus:

    show_legend : bool
        If True, add a static color-key legend (built from invisible
        dummy traces) alongside the hover text. Off by default, since
        the point of this version is that hover replaces the legend.
    """

    # -----------------------------------------------------------------
    # Validate input
    # -----------------------------------------------------------------

    if pfam_col not in df.columns:
        raise KeyError("Column not found: " + repr(pfam_col))

    plot_df = df.dropna(subset=[pfam_col]).copy()

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No targets with domain annotations",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        return fig

    # -----------------------------------------------------------------
    # Normalize domain annotations (identical to matplotlib version)
    # -----------------------------------------------------------------

    def normalize_domains(value):
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

    plot_df[pfam_col] = plot_df[pfam_col].apply(normalize_domains)
    plot_df = plot_df[plot_df[pfam_col].apply(lambda d: len(d) > 0)].copy()

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No valid domain annotations",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)
        return fig

    # -----------------------------------------------------------------
    # Domain identifier (identical to matplotlib version)
    # -----------------------------------------------------------------

    def domain_identifier(domain):
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

    # -----------------------------------------------------------------
    # Canonical architecture key (identical to matplotlib version)
    # -----------------------------------------------------------------

    def architecture_key(domains):
        identifiers = [
            domain_identifier(d) for d in domains
            if domain_identifier(d) is not None
        ]
        if not identifiers:
            return None
        return tuple(sorted(set(identifiers)))

    plot_df["_architecture_key"] = plot_df[pfam_col].apply(architecture_key)
    plot_df = plot_df.dropna(subset=["_architecture_key"]).copy()

    # -----------------------------------------------------------------
    # Representative selection (identical to matplotlib version)
    # -----------------------------------------------------------------

    if group_by_architecture:
        if "z_score" not in plot_df.columns:
            raise KeyError(
                "group_by_architecture=True requires a 'z_score' column."
            )

        plot_df["_plot_z_score"] = pd.to_numeric(
            plot_df["z_score"], errors="coerce"
        )

        plot_df = plot_df.sort_values(
            "_plot_z_score", ascending=False,
            kind="mergesort", na_position="last",
        )
        plot_df = plot_df.drop_duplicates(
            subset=["_architecture_key"], keep="first"
        )
        plot_df = plot_df.sort_values(
            "_plot_z_score", ascending=False,
            kind="mergesort", na_position="last",
        )

    elif sort_col is not None and sort_col in plot_df.columns:
        plot_df = plot_df.sort_values(
            sort_col, ascending=ascending, kind="mergesort"
        )

    plot_df = plot_df.head(max_targets).copy()

    # -----------------------------------------------------------------
    # Colors (matplotlib tab20, converted to hex for Plotly)
    # -----------------------------------------------------------------

    all_domain_ids = sorted({
        identifier
        for domains in plot_df[pfam_col]
        for domain in domains
        for identifier in [domain_identifier(domain)]
        if identifier is not None
    })

    cmap = colormaps["tab20"]
    colors = {
        identifier: to_hex(cmap(i % cmap.N))
        for i, identifier in enumerate(all_domain_ids)
    }

    # -----------------------------------------------------------------
    # Name lookup (identical logic to matplotlib version)
    # -----------------------------------------------------------------

    name_lookup = {}
    for domains in plot_df[pfam_col]:
        for domain in domains:
            if not isinstance(domain, dict):
                continue
            identifier = domain_identifier(domain)
            if identifier is None:
                continue
            name = domain.get("name") or ""
            if identifier not in name_lookup and name:
                name_lookup[identifier] = str(name)

    # -----------------------------------------------------------------
    # Row labels (identical logic to matplotlib version)
    # -----------------------------------------------------------------

    labels = []
    for _, row in plot_df.iterrows():
        target_id = str(row.get("target_id", row.name))

        if group_by_architecture and "_plot_z_score" in plot_df.columns:
            z_score = row["_plot_z_score"]
            label = (
                f"{target_id}  (Z={float(z_score):.1f})"
                if pd.notna(z_score) else target_id
            )
        elif (not group_by_architecture and sort_col is not None
              and sort_col in plot_df.columns):
            try:
                label = f"{target_id}  ({sort_col}={float(row[sort_col]):.1f})"
            except (TypeError, ValueError):
                label = target_id
        else:
            label = target_id

        labels.append(label)

    # -----------------------------------------------------------------
    # Build traces
    #
    # Numeric row positions (0, 1, 2, ...) are used throughout, with
    # tick labels applied afterwards via yaxis.tickvals/ticktext - this
    # mirrors the matplotlib version's set_yticks/set_yticklabels split
    # exactly, and sidesteps any risk of duplicate-label collisions
    # that a categorical (label-as-position) y-axis would have.
    # -----------------------------------------------------------------

    fig = go.Figure()

    # --- base sequence line, one Scatter trace for all rows,
    #     segments separated by None so they don't connect ---
    base_x, base_y = [], []
    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        domains = row[pfam_col]

        target_length = row.get("target_length", None)
        if target_length is None or pd.isna(target_length):
            ends = [
                float(d["end"]) for d in domains
                if isinstance(d, dict) and d.get("end") is not None
            ]
            target_length = max(ends) if ends else 100

        base_x += [0, target_length, None]
        base_y += [row_idx, row_idx, None]

    fig.add_trace(go.Scatter(
        x=base_x, y=base_y,
        mode="lines",
        line=dict(color="#dddddd", width=1.5),
        hoverinfo="skip",
        showlegend=False,
    ))

    # --- domain segments, one Bar trace, array-based ---
    bar_x, bar_base, bar_y, bar_color, bar_hover = [], [], [], [], []

    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        target_id = str(row.get("target_id", row.name))
        domains = row[pfam_col]

        for domain in domains:
            if not isinstance(domain, dict):
                continue
            identifier = domain_identifier(domain)
            if identifier is None:
                continue
            try:
                start = float(domain["start"])
                end = float(domain["end"])
            except (KeyError, TypeError, ValueError):
                continue

            width = end - start
            if width <= 0:
                continue

            name = name_lookup.get(identifier, "")
            hover_lines = [f"<b>{identifier}</b>"]
            if name:
                hover_lines.append(name)
            hover_lines.append(f"Position: {int(start)}-{int(end)}")
            hover_lines.append(f"Target: {target_id}")

            bar_x.append(width)
            bar_base.append(start)
            bar_y.append(row_idx)
            bar_color.append(colors[identifier])
            bar_hover.append("<br>".join(hover_lines))

    fig.add_trace(go.Bar(
        x=bar_x, base=bar_base, y=bar_y,
        orientation="h",
        width=0.6,
        marker=dict(color=bar_color, line=dict(color="white", width=0.5)),
        hovertext=bar_hover,
        hovertemplate="%{hovertext}<extra></extra>",
        showlegend=False,
    ))

    # --- alignment-region highlight, dashed outline shapes ---
    for row_idx, (_, row) in enumerate(plot_df.iterrows()):
        t_left = row.get("target_left", None)
        t_right = row.get("target_right", None)

        if t_left is None or t_right is None:
            continue
        if pd.isna(t_left) or pd.isna(t_right):
            continue

        try:
            t_left = float(t_left)
            t_right = float(t_right)
        except (TypeError, ValueError):
            continue

        if t_right < t_left:
            continue

        fig.add_shape(
            type="rect",
            x0=t_left, x1=t_right + 1,
            y0=row_idx - 0.4, y1=row_idx + 0.4,
            line=dict(color="#e05c5c", width=1.5, dash="dash"),
            fillcolor="rgba(0,0,0,0)",
        )

    # --- optional static legend, invisible dummy traces ---
    if show_legend:
        for identifier in all_domain_ids:
            name = name_lookup.get(identifier, "")
            legend_label = f"{identifier}  {name}" if name else identifier
            fig.add_trace(go.Bar(
                x=[None], y=[None],
                orientation="h",
                marker=dict(color=colors[identifier]),
                name=legend_label,
                showlegend=True,
                hoverinfo="skip",
            ))

    # -----------------------------------------------------------------
    # Layout
    # -----------------------------------------------------------------

    if group_by_architecture:
        title = "Domain architecture (one highest-Z representative per composition)"
    elif sort_col is not None:
        title = f"Domain architecture (sorted by {sort_col})"
    else:
        title = "Domain architecture"

    fig_height = max(300, 45 * len(plot_df) + 150)

    fig.update_layout(
        title=title,
        xaxis_title="Residue position",
        yaxis=dict(
            tickvals=list(range(len(plot_df))),
            ticktext=labels,
            autorange="reversed",  # row 0 (highest Z) at the top
            tickfont=dict(size=10),
        ),
        height=fig_height,
        width=1000,
        barmode="overlay",
        plot_bgcolor="white",
        showlegend=show_legend,
        margin=dict(l=10, r=10, t=50, b=40),
    )

    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=False, zeroline=False)

    return fig
