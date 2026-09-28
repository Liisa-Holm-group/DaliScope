import ipywidgets as widgets
from IPython.display import display, clear_output
from daliscope.optics.plotly import plotly_overview

def plotly_explorer(project, x='z_score', y='query_coverage', color='clan', size=None, marker='clan'):
    view_names = list(project.views.keys())
    if not view_names:
        print("❌ No views registered yet")
        return

    def get_col_options(df):
        numeric     = df.select_dtypes(include='number').columns.tolist()
        categorical = [None] + ['pfam', 'clan']
        # df.select_dtypes(exclude='number').columns.tolist()

        return numeric, categorical

    # --- top-level view dropdown only ---
    view_dropdown = widgets.Dropdown(
        options     = view_names,
        value       = view_names[0],
        description = "View:",
        style       = {"description_width": "50px"},
        layout      = widgets.Layout(width="300px"),
    )

    controls_box = widgets.HBox([])   # rebuilt on each view change
    plot_out     = widgets.Output()
    info_out     = widgets.Output()

    def make_col_dropdowns(df):
        """Create fresh dropdowns for the current df — no mutation."""
        numeric, categorical = get_col_options(df)

        x_val     = x     if x     in numeric     else numeric[0]
        y_val     = y     if y     in numeric      else (
                        numeric[1] if len(numeric) > 1 else numeric[0])
        color_val = color ##if color in categorical  else None
        size_val = size if size in numeric else None
        marker_val = marker if marker in categorical else None

        xd = widgets.Dropdown(
            options=numeric, value=x_val,
            description="X:",
            style={"description_width": "30px"},
            layout=widgets.Layout(width="200px"),
        )
        yd = widgets.Dropdown(
            options=numeric, value=y_val,
            description="Y:",
            style={"description_width": "30px"},
            layout=widgets.Layout(width="200px"),
        )
        cd = widgets.Dropdown(
            options=categorical+numeric, value=color_val,
            description="Color:",
            style={"description_width": "50px"},
            layout=widgets.Layout(width="220px"),
        )
        md = widgets.Dropdown(
            options=categorical, value=marker_val,
            description="Marker:",
            style={"description_width": "50px"},
            layout=widgets.Layout(width="220px"),
        )
        sd = widgets.Dropdown(
            options=numeric, value=size_val,
            description="Size:",
            style={"description_width": "50px"},
            layout=widgets.Layout(width="220px"),
        )
        return xd, yd, cd, md, sd

    def draw(xd, yd, cd, md, sd):
        df = project.views[view_dropdown.value]
        with info_out:
            clear_output(wait=True)
            meta = project.views[view_dropdown.value]
            print(f"View: '{view_dropdown.value}'")
        with plot_out:
            clear_output(wait=True)
            plotly_overview(df,
                            x=xd.value,
                            y=yd.value,
                            color=cd.value,
                            size=sd.value,
                            marker=md.value)

    def setup_view(_=None):
        """Rebuild col dropdowns and rewire callbacks for current view."""
        df       = project.views[view_dropdown.value]
        xd, yd, cd, md, sd = make_col_dropdowns(df)

        def on_col_change(_):
            draw(xd, yd, cd, md, sd)

        xd.observe(on_col_change, names="value")
        yd.observe(on_col_change, names="value")
        cd.observe(on_col_change, names="value")
        md.observe(on_col_change, names="value")
        sd.observe(on_col_change, names="value")

        controls_box.children = [xd, yd, cd, md, sd]
        draw(xd, yd, cd, md, sd)

    view_dropdown.observe(setup_view, names="value")

    ui = widgets.VBox([
        widgets.HTML("<b>Scatter plot explorer</b>"),
        view_dropdown,
        controls_box,
        info_out,
        plot_out,
    ])
    display(ui)
    setup_view()   # initial render

#############

import ipywidgets as widgets
from IPython.display import display
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def launch_interactive_selectable_scatter(
    project,
    default_view_name='FILTERED_1',
):
    """
    Interactive selectable scatter plot for Jupyter / Colab.
    """

    if not hasattr(project, "views"):
        raise AttributeError(
            "project must have a 'views' attribute containing the available views."
        )

    views = project.views

    if default_view_name not in views:
        raise ValueError(
            f"Unknown default view '{default_view_name}'. "
            f"Available views: {list(views.keys())}"
        )

    def get_view_df(view_name):
        df = views[view_name]
        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)
        return df.copy()

    def get_reordered_numeric_cols(df):
        cols = df.select_dtypes(include="number").columns.tolist()
        priority = ["z_score", "query_coverage"]
        return [c for c in priority if c in cols] + [
            c for c in cols if c not in priority
        ]

    initial_df = get_view_df(default_view_name)
    numeric_columns = get_reordered_numeric_cols(initial_df)

    if len(numeric_columns) < 2:
        raise ValueError(
            "The selected view must contain at least two numeric columns "
            "for an X/Y scatter plot."
        )

    default_x = "z_score" if "z_score" in numeric_columns else numeric_columns[0]
    default_y = (
        "query_coverage"
        if "query_coverage" in numeric_columns
        else (numeric_columns[1] if len(numeric_columns) > 1 else numeric_columns[0])
    )

    trace_selections = {}

    # ------------------------------------------------------------------
    # Widgets UI Setup
    # ------------------------------------------------------------------
    view_dropdown = widgets.Dropdown(
        options=list(views.keys()),
        value=default_view_name,
        description="View:",
        layout=widgets.Layout(width="280px"),
    )

    x_dropdown = widgets.Dropdown(
        options=numeric_columns,
        value=default_x,
        description="X:",
        layout=widgets.Layout(width="260px"),
    )

    y_dropdown = widgets.Dropdown(
        options=numeric_columns,
        value=default_y,
        description="Y:",
        layout=widgets.Layout(width="260px"),
    )

    color_dropdown = widgets.Dropdown(
        options=[("Clan", "clan"), ("Pfam", "pfam"), ("None", None)],
        value="clan",
        description="Color:",
        layout=widgets.Layout(width="250px"),
    )

    marker_dropdown = widgets.Dropdown(
        options=[("Pfam", "pfam"), ("Clan", "clan"), ("None", None)],
        value="pfam",
        description="Marker:",
        layout=widgets.Layout(width="250px"),
    )

    size_dropdown = widgets.Dropdown(
        options=[
            ("None", None),
            ("Sequence identity", "sequence_identity"),
            ("Query coverage", "query_coverage"),
            ("RMSD", "rmsd"),
            ("Target coverage", "target_coverage"),
        ],
        value=None,
        description="Size:",
        layout=widgets.Layout(width="250px"),
    )

    k_slider = widgets.IntSlider(
        value=10,
        min=1,
        max=50,
        step=1,
        description="Top-k Pfam:",
        continuous_update=False,
        layout=widgets.Layout(width="330px"),
    )

    subset_name_input = widgets.Text(
        value="",
        placeholder="Optional subset name",
        description="Subset:",
        layout=widgets.Layout(width="330px"),
    )

    save_button = widgets.Button(
        description="Save Subset",
        button_style="success",
        icon="save",
        disabled=True,
        layout=widgets.Layout(width="130px"),
    )

    reset_button = widgets.Button(
        description="Reset",
        button_style="warning",
        icon="refresh",
        layout=widgets.Layout(width="110px"),
    )

    output_box = widgets.Output(
        layout=widgets.Layout(
            width="800px",
            border="1px solid #ddd",
            padding="6px",
        )
    )

    plot_container = widgets.Output(
        layout=widgets.Layout(
            width="820px",
        )
    )

    def prepare_dataframe():
        df = get_view_df(view_dropdown.value)
        df = df.reset_index(drop=False)

        if "_custom_index" not in df.columns:
            df["_custom_index"] = df.index

        if "target_id" in df.columns:
            df["_hover_target_id"] = df["target_id"].astype(str)
        else:
            df["_hover_target_id"] = df["_custom_index"].astype(str)

        if "description" in df.columns:
            df["_hover_desc"] = df["description"].fillna("").astype(str)
        else:
            df["_hover_desc"] = ""

        return df

    def apply_top_k_filter(df):
        if "pfam" not in df.columns:
            return df

        pfam_series = df["pfam"]
        if len(pfam_series) > 0:
            first_non_null = None
            for value in pfam_series:
                if value is not None and not (
                    isinstance(value, float) and pd.isna(value)
                ):
                    first_non_null = value
                    break

            if isinstance(first_non_null, (list, tuple, set)):
                counts = {}
                for value in pfam_series:
                    if isinstance(value, (list, tuple, set)):
                        for pfam in value:
                            counts[pfam] = counts.get(pfam, 0) + 1

                if counts:
                    top_pfams = {
                        pfam
                        for pfam, _ in sorted(
                            counts.items(),
                            key=lambda item: item[1],
                            reverse=True,
                        )[: k_slider.value]
                    }

                    mask = pfam_series.apply(
                        lambda value: bool(set(value).intersection(top_pfams))
                        if isinstance(value, (list, tuple, set))
                        else False
                    )

                    return df.loc[mask].copy()

        counts = pfam_series.value_counts(dropna=True)
        if len(counts) > k_slider.value:
            top_pfams = set(counts.head(k_slider.value).index)
            return df.loc[pfam_series.isin(top_pfams)].copy()

        return df

    def on_trace_select(trace, points, selector, trace_idx):
        if points.point_inds and trace.customdata is not None:
            selected_ids = []
            for point_index in points.point_inds:
                try:
                    selected_ids.append(trace.customdata[point_index][0])
                except Exception:
                    pass
            trace_selections[trace_idx] = selected_ids
        else:
            trace_selections[trace_idx] = []

        all_selected_ids = [
            pid for selected in trace_selections.values() for pid in selected
        ]
        all_selected_ids = list(dict.fromkeys(all_selected_ids))
        count = len(all_selected_ids)

        if count > 0:
            save_button.disabled = False
            if not subset_name_input.value.strip():
                subset_name_input.value = (
                    f"BOX_{view_dropdown.value}_{count}_HITS"
                )

            with output_box:
                output_box.clear_output()
                print(
                    f"Selected {count:,} targets across categories. "
                    f"Click 'Save Subset' to commit."
                )
        else:
            save_button.disabled = True
            with output_box:
                output_box.clear_output()
                print("Selection cleared.")

    def create_styled_figure():
        df = prepare_dataframe()

        if color_dropdown.value == "pfam" or "pfam" in df.columns:
            df = apply_top_k_filter(df)

        x_col = x_dropdown.value
        y_col = y_dropdown.value
        color_col = color_dropdown.value
        marker_col = marker_dropdown.value
        size_col = size_dropdown.value

        if size_col is not None and size_col not in df.columns:
            size_col = None

        if color_col is not None and color_col not in df.columns:
            color_col = None

        if marker_col is not None and marker_col not in df.columns:
            marker_col = None

        for c in (color_col, marker_col):
            if c in ("pfam", "clan") and c in df.columns:
                df[c] = df[c].apply(
                    lambda val: ", ".join(map(str, val))
                    if isinstance(val, (list, tuple, set))
                    else str(val)
                    if pd.notna(val)
                    else "Unclassified"
                )

        express_fig = px.scatter(
            df,
            x=x_col,
            y=y_col,
            color=color_col,
            symbol=marker_col,
            size=size_col,
            size_max=16,
            hover_data=[
                "_hover_target_id",
                "_hover_desc",
                x_col,
                y_col,
            ],
            color_discrete_sequence=px.colors.qualitative.Dark24,
            template="simple_white",
            height=520,
            width=800,
            title=f"<b>{x_col}</b> vs <b>{y_col}</b> ({len(df):,} targets)",
            render_mode="webgl",
        )

        hover_template = (
            f"<b>Target ID:</b> %{{customdata[0]}}<br>"
            f"<b>Description:</b> %{{customdata[1]}}<br>"
            f"<b>{x_col}:</b> %{{x}}<br>"
            f"<b>{y_col}:</b> %{{y}}"
            "<extra></extra>"
        )

        express_fig.update_traces(
            hovertemplate=hover_template,
            selected=dict(marker=dict(opacity=1.0)),
            unselected=dict(marker=dict(opacity=0.30)),
        )

        express_fig.update_layout(
            font=dict(family="Arial, sans-serif", size=12),
            margin=dict(l=60, r=40, t=50, b=50),
            showlegend=True,
            legend=dict(
                title=dict(text=f"<b>{color_col or 'Legend'}</b>"),
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
            ),
            xaxis=dict(showgrid=True, gridcolor="#EAEAEA"),
            yaxis=dict(showgrid=True, gridcolor="#EAEAEA"),
        )

        fig_widget = go.FigureWidget(express_fig)
        fig_widget.update_layout(dragmode="select")

        trace_selections.clear()

        for idx, trace in enumerate(fig_widget.data):
            trace_selections[idx] = []
            trace.on_selection(
                lambda tr, pts, sel, i=idx: on_trace_select(tr, pts, sel, i)
            )

        return fig_widget

    def redraw_plot(message=None):
        trace_selections.clear()
        save_button.disabled = True
        subset_name_input.value = ""

        with plot_container:
            plot_container.clear_output(wait=True)
            new_fig = create_styled_figure()
            display(new_fig)

        with output_box:
            output_box.clear_output()
            if message is not None:
                print(message)

    def on_reset_click(button):
        redraw_plot("Plot reset. Ready for a new selection.")

    reset_button.on_click(on_reset_click)

    def on_save_click(button):
        all_selected_ids = [
            pid for selected in trace_selections.values() for pid in selected
        ]
        all_selected_ids = list(dict.fromkeys(all_selected_ids))

        if not all_selected_ids:
            with output_box:
                output_box.clear_output()
                print("No points are selected.")
            save_button.disabled = True
            return

        parent_view_name = view_dropdown.value
        parent_df = get_view_df(parent_view_name)

        subset_name = subset_name_input.value.strip()
        if not subset_name:
            subset_name = (
                f"BOX_{parent_view_name}_{len(all_selected_ids)}_HITS"
            )

        # Generate boolean mask against parent view
        if "target_id" in parent_df.columns:
            mask = parent_df["target_id"].astype(str).isin([str(x) for x in all_selected_ids])
        elif "_custom_index" in parent_df.columns:
            mask = parent_df["_custom_index"].isin(all_selected_ids)
        else:
            mask = parent_df.index.isin(all_selected_ids)

        try:
            _ = project.add_subset(
                name=subset_name,
                mask=mask,
                parent=parent_view_name,
                function="interactive_box_select",
                parameters={
                    "x_col": x_dropdown.value,
                    "y_col": y_dropdown.value,
                    "selected_count": len(all_selected_ids),
                },
            )

            redraw_plot(
                f"Successfully registered subset '{subset_name}' with "
                f"{len(all_selected_ids):,} targets using parent '{parent_view_name}'."
            )
        except Exception as exc:
            with output_box:
                output_box.clear_output()
                print(f"Could not save subset '{subset_name}':\n{exc}")

    save_button.on_click(on_save_click)

    def update_columns_for_view():
        df = get_view_df(view_dropdown.value)
        numeric_cols = get_reordered_numeric_cols(df)

        if len(numeric_cols) < 2:
            with output_box:
                output_box.clear_output()
                print(
                    f"View '{view_dropdown.value}' has fewer than two numeric columns."
                )
            x_dropdown.options = numeric_cols
            y_dropdown.options = numeric_cols
            return

        current_x = x_dropdown.value
        current_y = y_dropdown.value

        x_dropdown.unobserve(update_plot, names="value")
        y_dropdown.unobserve(update_plot, names="value")

        x_dropdown.options = numeric_cols
        y_dropdown.options = numeric_cols

        fallback_x = (
            "z_score" if "z_score" in numeric_cols else numeric_cols[0]
        )
        fallback_y = (
            "query_coverage"
            if "query_coverage" in numeric_cols
            else (
                numeric_cols[1] if len(numeric_cols) > 1 else numeric_cols[0]
            )
        )

        x_dropdown.value = (
            current_x if current_x in numeric_cols else fallback_x
        )
        y_dropdown.value = (
            current_y if current_y in numeric_cols else fallback_y
        )

        x_dropdown.observe(update_plot, names="value")
        y_dropdown.observe(update_plot, names="value")

    def update_plot(change=None):
        if change is not None and change["owner"] is view_dropdown:
            update_columns_for_view()
        redraw_plot()

    view_dropdown.observe(update_plot, names="value")
    x_dropdown.observe(update_plot, names="value")
    y_dropdown.observe(update_plot, names="value")
    color_dropdown.observe(update_plot, names="value")
    marker_dropdown.observe(update_plot, names="value")
    size_dropdown.observe(update_plot, names="value")
    k_slider.observe(update_plot, names="value")

    controls_row1 = widgets.HBox([view_dropdown, x_dropdown, y_dropdown])
    controls_row2 = widgets.HBox(
        [color_dropdown, marker_dropdown, size_dropdown]
    )
    controls_row3 = widgets.HBox(
        [k_slider, subset_name_input, save_button, reset_button]
    )

    ui = widgets.VBox(
        [
            controls_row1,
            controls_row2,
            controls_row3,
            plot_container,
            output_box,
        ]
    )

    display(ui)
    redraw_plot("Ready. Select points with box/lasso, then click 'Save Subset'.")
