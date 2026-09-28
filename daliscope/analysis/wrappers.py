import daliscope

import pandas as pd
from IPython.display import display, HTML
import matplotlib.pyplot as plt

import ast

import numpy as np

import matplotlib.pyplot as plt

from IPython.display import display

def make_pfam_counts_table(view_name, project, **kwargs):
    return daliscope.analysis.pfam_stats.pfam_counts_df(
        project.views[view_name],
        project._pfam_names,
        view_name=view_name,
        **kwargs,
    )


def make_top_architecture_cartoons(view_name, pfam_col, k, project, **kwargs):
    """
    Plot up to k unique domain architectures from a project view.

    Representative selection is performed by
    daliscope.optics.domain_architecture.prepare_domain_architecture_data():

    - architecture = domain composition
    - domain order ignored
    - duplicate domain copies ignored
    - highest-Z target selected as representative
    - representatives ordered by descending Z-score
    """
    project.viz.architecture(
        view_name,
        max_targets=k,
        pfam_col=pfam_col,
        renderer="plotly",
        group_by_architecture=True,
        **kwargs,
    )

def old_make_top_architecture_cartoons(
    view_name,
    pfam_col,
    k,
    project,
    **kwargs,
):
    """
    Plot up to k unique domain architectures from a project view.

    Representative selection is performed by
    daliscope.optics.matplot.plot_domain_architecture():

    - architecture = domain composition
    - domain order ignored
    - duplicate domain copies ignored
    - highest-Z target selected as representative
    - representatives ordered by descending Z-score
    """

    df = project.views[
        view_name
    ].copy()


    fig = (
        daliscope.optics.plotly.old_plot_domain_architecture_plotly(
        #daliscope.optics.matplot.plot_domain_architecture(
            df,
            pfam_col=pfam_col,
            max_targets=k,
            group_by_architecture=True,
            **kwargs,
        )
    )


    display(
        fig
    )

    # plotly, not plt
    #plt.close(
    #    fig
    #)


    return None

def make_sunburst_plot(view_name, k, project, sort_by=None, ascending=False):
    # 1. Fetch live view
    df = project.views[view_name]

    # 2. Generate plot and retrieve plot_df (containing 'pfam_clean')
    fig, plot_df = daliscope.analysis.pfam_stats.plot_pfam_clan_sunburst(df, k=k)

    # 3. Build table from transformed plot_df
    table_df = daliscope.analysis.pfam_stats.print_pfam_names_table_with_counts(
        pfam_names_df=project._pfam_names,
        plot_df=plot_df
    )

    # 4. Sort table
    if sort_by:
        table_df = table_df.sort_values(by=sort_by, ascending=ascending)

    # 5. Display table and return plot
    display(table_df)
    return fig

# auto-generate names with auto-incremented index
def next_subset_name(project, parent: str, label: str = 'trusted') -> str:
    """
    Generate the next available name like '{parent}_{label}_0',
    '{parent}_{label}_1', etc. by checking existing view names.
    """
    i = 0
    while f"{parent}_{label}_{i}" in project.views:
        i += 1
    return f"{parent}_{label}_{i}"

def prev_make_top_architecture_cartoons(view_name, pfam_col, k, project, **kwargs):
    """
    Finds the highest z-score representative target per unique domain architecture 
    and forwards only those top-k representative cartoons to the plotter.
    """
    # 1. Pull active dataframe from project view
    df = project.views[view_name].copy()
    col_to_use = pfam_col if pfam_col in df.columns else 'architecture'

    # Helper function: Extract unique Pfams (or Clans), exclude 'Other'/'Unassigned', sort alphabetically, create hash key
    def create_arch_hash_key(domain_val):
        if isinstance(domain_val, (list, tuple, pd.Series)):
            raw_items = domain_val
        elif pd.notna(domain_val) and isinstance(domain_val, str):
            # Split string representations if stored as comma/space separated
            raw_items = [x.strip() for x in domain_val.replace(';', ',').split(',')]
        else:
            return None

        # Filter out invalid, 'Other', and 'Unassigned' entries
        valid_items = [
            str(item).strip() for item in raw_items
            if pd.notna(item)
            and str(item).strip() != '' 
            and 'Other' not in str(item) 
            and 'Unassigned' not in str(item)
        ]

        if not valid_items:
            return None

        # Deduplicate and sort alphabetically to ensure invariant architecture hashing
        unique_sorted = sorted(list(set(valid_items)))
        return " | ".join(unique_sorted)

    # 2. Assign the canonical architecture hash key
    df['arch_hash'] = df[col_to_use].apply(create_arch_hash_key)

    # Drop rows without valid Pfam compositions
    df = df.dropna(subset=['arch_hash']).copy()

    # 3. Identify structure ID column
    struct_col = 'target_id' if 'target_id' in df.columns else ('query_id' if 'query_id' in df.columns else None)

    # 4. Group by target structure to find its max Z-score and architecture hash
    if struct_col and 'z_score' in df.columns:
        struct_summary = (
            df.groupby(struct_col, as_index=False)
              .agg({'z_score': 'max', 'arch_hash': 'first'})
        )

        # 5. Pick the single best target structure (highest Z-score) per unique architecture hash
        best_struct_per_arch = (
            struct_summary.sort_values(by='z_score', ascending=False)
                          .groupby('arch_hash', as_index=False)
                          .first()
        )

        # 6. Rank unique architectures by their representative max Z-score and take top-k
        top_k_architectures = best_struct_per_arch.sort_values(by='z_score', ascending=False).head(k)

        # Filter raw dataframe for hits belonging to these top-k representative targets
        df_top = df[df[struct_col].isin(top_k_architectures[struct_col])].copy()

    else:
        # Fallback if target_id column is absent
        df_top = (
            df.sort_values(by='z_score', ascending=False) if 'z_score' in df.columns else df
        ).groupby('arch_hash', as_index=False).first().head(k)

    # 7. Update target_id labels to include Z-score for the left-side legend/y-axis
    if struct_col in df_top.columns and 'z_score' in df_top.columns:
        # Map existing target_id -> "target_id (Z: XX.X)"
        def format_target_label(row):
            target = str(row[struct_col])
            z_val = row['z_score']
            return f"{target} (Z: {z_val:.1f})" if pd.notna(z_val) else target

        # If multiple domain rows exist per target, create a new target_label column
        df_top['target_label'] = df_top.apply(format_target_label, axis=1)

        # Override the target_id column so downstream optics plotter picks it up for the legend
        df_top[struct_col] = df_top['target_label']

    # Clean up internal temporary key
    df_top.drop(columns=['arch_hash', 'target_label'], errors='ignore', inplace=True)

    # 8. Render domain architecture plot
    fig = daliscope.optics.matplot.plot_domain_architecture(df_top, pfam_col=col_to_use, **kwargs)

    if fig is None:
        fig = plt.gcf()

    # Render explicitly and suppress duplicate inline shadow
    display(fig)
    plt.clf()
    plt.close(fig)

    return None

def old_make_top_architecture_cartoons(view_name, pfam_col, project, k):
    """
    Adapter that finds the highest z-score representative per unique architecture
    and forwards only those representative hits to the plotter.
    """
    # 1. Pull the live dataframe from the project based on the dropdown selection
    df = project.views[view_name]

    # 2. Filter / group to get highest z-score representative per architecture
    # Assumes 'architecture' (or pfam_col) identifies the unique architecture
    # and 'z_score' is used for ranking representatives
    if 'z_score' in df.columns:
        # Sort by z-score descending and take top representative per architecture
        group_col = pfam_col if pfam_col in df.columns else 'architecture'
        df_top = (
            df.sort_values(by='z_score', ascending=False)
              .groupby(group_col, as_index=False)
              .head(1)
        )
    else:
        df_top = df.copy()

    # 3. Restrict to top-k architectures by rank or row count
    if 'pfam_rank' in df_top.columns:
        df_top = df_top[df_top['pfam_rank'] <= k]
    else:
        df_top = df_top.head(k)

    # 4. Generate the Matplotlib architecture plot
    fig = daliscope.optics.matplot.plot_domain_architecture(df_top)

    # If the function returns None but renders to current gcf, capture it explicitly
    if fig is None:
        fig = plt.gcf()

    # 5. Render figure explicitly into the ipywidgets output canvas
    display(fig)
    plt.close(fig)  # Prevents double-rendering and memory leaks in Jupyter

    return fig


def make_violin_plot(view_name, group_by, k, project, metrics):
    """
    Wrapper function that pulls the live dataframe dynamically from the project 
    view state, filters it, and returns the Matplotlib figure.
    """
    # 1. Fetch the live dataframe dynamically using the active dropdown key
    df = project.views[view_name]

    # 2. Determine the correct rank column dynamically
    rank_col = 'clan_rank' if group_by == 'clan' else 'pfam_rank'

    # 3. Filter for top-k ranks safely
    df_violin = df[df[rank_col] <= k]

    # 4. Generate the plot
    fig = daliscope.analysis.pfam_stats.plot_violins(
        df_violin,
        group_by,
        metrics,
        cols=3
    )

    return fig

def make_family_presence_plot(view_name, category_col, k, show_zero_present, project, population_view, id_col="target_id"):
    """
    Wrapper function that treats `view_name` as the Target View, pulls 
    dataframes dynamically, and returns the Matplotlib figure.
    """
    population_df = project.views[population_view]
    target_df = project.views[view_name]
    rank_col = f"{category_col}_rank"

    # Validation checks
    if category_col not in population_df.columns or rank_col not in population_df.columns:
        print(f"Required rank/category columns missing in population view '{population_view}'.")
        return None
    if id_col not in target_df.columns:
        print(f"Column '{id_col}' not found in target view '{view_name}'.")
        return None

    # Generate the plot
    fig, ax = daliscope.optics.matplot.plot_family_presence(
        population_df=population_df,
        family_col=category_col,
        rank_col=rank_col,
        target_ids=target_df[id_col],
        id_col=id_col,
        k=k,
        show_zero_present=show_zero_present,
        title=(f"Presence of '{view_name}' across top {k} "
               f"'{category_col}' families in '{population_view}'")
    )

    return fig


def apply_rank_pooling(df: pd.DataFrame,
                       rank_col: str = 'pfam_rank',
                       label_col: str = 'pfam',
                       k: int = 10) -> pd.Series:
    """
    Returns a Series of display labels: rank<=k keeps its pfam/clan
    value, rank>k becomes 'Other', NaN becomes 'Unassigned'.
    Never stored back into the view — computed fresh at each plot call.
    """
    rank  = df[rank_col]
    label = df[label_col].copy().astype(object)

    label[rank.isna()]      = 'Unassigned'
    label[rank > k]         = 'Other'
    # rank<=k keeps its own pfam/clan value (already there)

    return label

def make_view_plot_fn(plot_fn, project):
    """Resolves view_name -> df before delegating."""
    def wrapper(view_name, **kwargs):
        df = project.views[view_name]
        return plot_fn(df=df, **kwargs)
    return wrapper

def make_pfam_plot_fn(plot_fn,
                      color_kwarg:  str = 'color',
                      symbol_kwarg: str = 'symbol',
                      size_kwarg:   str = 'size',
                      hover_cols:   list = None):
    """
    Intercepts color_by, marker_by, size_by and k.
    Copies only the columns actually needed for the plot.
    """
    _hover_cols = hover_cols or ['target_id', 'description']

    def wrapper(df, color_by='clan', marker_by='none',
                size_by='none', k=10, **kwargs):

        # --- build needed_cols from all column-name arguments ---
        x_col      = kwargs.get('x')
        y_col      = kwargs.get('y')
        rank_cols  = ['pfam_rank', 'clan_rank']

        candidates = [
            x_col,
            y_col,
            color_by  if color_by  != 'none' else None,
            marker_by if marker_by != 'none' else None,
            size_by   if size_by   != 'none' else None,
            *rank_cols,
            *_hover_cols,
        ]
        needed_cols = list(dict.fromkeys(
            c for c in candidates
            if c is not None and c in df.columns
        ))

        plot_df = df[needed_cols].copy()   # ← targeted copy

        # --- apply rank pooling to color and marker ---
        if color_by != 'none' and color_by in plot_df.columns:
            rank_col = 'clan_rank' if color_by == 'clan' else 'pfam_rank'
            plot_df['_color_label'] = apply_rank_pooling(
                plot_df, rank_col=rank_col, label_col=color_by, k=k
            )
            color_arg = '_color_label'
        else:
            color_arg = None

        if marker_by != 'none' and marker_by in plot_df.columns:
            rank_col = 'clan_rank' if marker_by == 'clan' else 'pfam_rank'
            plot_df['_marker_label'] = apply_rank_pooling(
                plot_df, rank_col=rank_col, label_col=marker_by, k=k
            )
            symbol_arg = '_marker_label'
        else:
            symbol_arg = None

        size_arg = size_by if (size_by != 'none' and
                               size_by in df.columns) else None

        return plot_fn(
            df = plot_df,
            **{color_kwarg:  color_arg},
            **{symbol_kwarg: symbol_arg},
            **{size_kwarg:   size_arg},
            **kwargs,
        )

    return wrapper

def make_view_pfam_plot_fn(plot_fn, project, color_kwarg='color'):
    """
    Convenience: compose both wrappers in one call.
    widget_view receives view_name, color_by, k as parameters;
    plot_fn receives df and color='color_label', knows nothing else.
    """
    pfam_wrapped = make_pfam_plot_fn(plot_fn, color_kwarg=color_kwarg)
    return make_view_plot_fn(pfam_wrapped, project)

def make_pfam_plot_fn(plot_fn, color_kwarg='color',
                      symbol_kwarg='symbol', size_kwarg='size'):
    """
    Intercepts color_by, marker_by, size_by and k.
    Applies rank pooling to color_by and marker_by (both categorical).
    Passes size_by as a numeric column name or None.
    """
    def wrapper(df, color_by='clan', marker_by='none',
                size_by='none', k=10, **kwargs):
        plot_df = df.copy()

        # color pooling
        if color_by != 'none' and color_by in df.columns:
            rank_col = 'clan_rank' if color_by == 'clan' else 'pfam_rank'
            plot_df['_color_label'] = apply_rank_pooling(
                df, rank_col=rank_col, label_col=color_by, k=k
            )
            color_arg = '_color_label'
        else:
            color_arg = None

        # marker/symbol pooling — same top-k logic
        if marker_by != 'none' and marker_by in df.columns:
            rank_col = 'clan_rank' if marker_by == 'clan' else 'pfam_rank'
            plot_df['_marker_label'] = apply_rank_pooling(
                df, rank_col=rank_col, label_col=marker_by, k=k
            )
            symbol_arg = '_marker_label'
        else:
            symbol_arg = None

        # size — numeric column or None
        size_arg = size_by if (size_by != 'none' and
                               size_by in df.columns) else None

        result = plot_fn(
            df     = plot_df,
            **{color_kwarg:  color_arg},
            **{symbol_kwarg: symbol_arg},
            **{size_kwarg:   size_arg},
            **kwargs,
        )
        print(f"wrapper: plot_fn returned {result}")
        return result
    return wrapper


def make_pfam_plot_fn(plot_fn, color_kwarg='color',
                      symbol_kwarg='symbol', size_kwarg='size'):
    """
    Intercepts color_by, marker_by, size_by and k.
    Applies rank pooling to color_by and marker_by (both categorical).
    Passes size_by as a numeric column name or None.
    """
    def wrapper(df, color_by='clan', marker_by='none',
                size_by='none', k=10, **kwargs):
        plot_df = df.copy()

        # color pooling
        if color_by != 'none' and color_by in df.columns:
            rank_col = 'clan_rank' if color_by == 'clan' else 'pfam_rank'
            plot_df['_color_label'] = apply_rank_pooling(
                df, rank_col=rank_col, label_col=color_by, k=k
            )
            color_arg = '_color_label'
        else:
            color_arg = None

        # marker/symbol pooling — same top-k logic
        if marker_by != 'none' and marker_by in df.columns:
            rank_col = 'clan_rank' if marker_by == 'clan' else 'pfam_rank'
            plot_df['_marker_label'] = apply_rank_pooling(
                df, rank_col=rank_col, label_col=marker_by, k=k
            )
            symbol_arg = '_marker_label'
        else:
            symbol_arg = None

        # size — numeric column or None
        size_arg = size_by if (size_by != 'none' and
                               size_by in df.columns) else None

        return plot_fn(
            df     = plot_df,
            **{color_kwarg:  color_arg},
            **{symbol_kwarg: symbol_arg},
            **{size_kwarg:   size_arg},
            **kwargs,
        )
    return wrapper

################

import pandas as pd
from IPython.display import display, HTML

def make_table_window_view(view_name, start_row, window, project, sort_col='z_score', ascending=False, columns=None, max_cols=None):
    """
    Wrapper function that pulls the live dataframe dynamically from project view state,
    ensures ascending order on the DataFrame index, filters selected columns, 
    slices the window [start_row : start_row + window], and renders the complete table
    view with multi-line text wrapping and no row or column truncation.
    """
    # 1. Fetch live dataframe dynamically using active dropdown view key
    df = project.views[view_name]

    ## 2. Guarantee ascending order on the index
    #df = df.sort_index(ascending=True)
    sd = df.sort_values(sort_col, ascending=ascending)

    total_rows = len(df)

    # 3. Filter requested columns if specified and valid
    if columns and isinstance(columns, (list, tuple)):
        valid_cols = [c for c in columns if c in df.columns]
        if valid_cols:
            df = df[valid_cols]

    # 4. Bound start_row safely within dataframe limits
    start_idx = max(0, min(start_row, total_rows))
    end_idx = min(start_idx + window, total_rows)

    # 5. Slice dataframe window
    window_df = df.iloc[start_idx:end_idx]

    if max_cols and isinstance(max_cols, int):
        window_df = window_df.iloc[:, :max_cols]

    # 6. Render context banner
    banner_html = f"""
    <div style="padding: 6px 12px; background-color: #f1f3f5; border-left: 4px solid #0d6efd; margin-bottom: 10px; font-family: sans-serif;">
        <b>Showing rows {start_idx:,} to {end_idx:,}</b> of <b>{total_rows:,}</b> total rows in <code>{view_name}</code> (sorted ascending by index)
    </div>
    """
    display(HTML(banner_html))

    # 7. Apply CSS for multi-line text wrapping & auto cell width expansion
    style_html = """
    <style>
        .dataframe td, .dataframe th {
            white-space: normal !important;
            word-wrap: break-word !important;
            max-width: 400px;
            text-align: left;
            vertical-align: top;
        }
    </style>
    """
    display(HTML(style_html))

    # 8. Render full window slice without row '...' middle-truncation or column width limits
    with pd.option_context('display.max_rows', window, 'display.max_colwidth', None):
        display(window_df)

    return None
