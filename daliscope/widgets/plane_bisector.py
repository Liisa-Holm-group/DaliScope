"""
Plane bisector selection widget.

Public API:
    bisector = plane_bisector(project, df=PARENT_DF)
    CHILD_DF, params = bisector.get_selected()
"""

import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display


def compute_mask(x, y, x_sect, y_sect, angle_deg):
    theta = np.deg2rad(angle_deg)
    nx = -np.sin(theta)
    ny = np.cos(theta)
    return (x - x_sect) * nx + (y - y_sect) * ny > 0


def make_view_plot_fn(plot_fn, project):
    """Resolves view_name -> df before delegating."""
    def wrapper(view_name, **kwargs):
        df = project.views[view_name]
        return plot_fn(df=df, **kwargs)
    return wrapper


# =====================================================================
# 1. Pure Rendering Function (Scatterplot fixed, bisector line movable)
# =====================================================================

def plot_bisector(
    df,
    left,
    right,
    bottom,
    top,
    x_sect,
    y_sect,
    angle,
    category_col="pfam",
    top_k=6,
    marker_cycle=None,
    other_marker=".",
    other_color="cyan",
    ax=None,
    **kwargs,
):
    """
    Pure rendering function. Axis bounds (left/right/bottom/top) remain fixed
    while (x_sect, y_sect) and angle define the position of the bisector line.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    if marker_cycle is None:
        marker_cycle = ["o", "s", "^", "D", "v", "P", "X", "*"]

    # Filter valid coordinates
    valid = np.isfinite(df["z_score"]) & np.isfinite(df["query_coverage"])
    work_df = df.loc[valid].copy()
    x = work_df["z_score"].to_numpy(dtype=float)
    y = work_df["query_coverage"].to_numpy(dtype=float)

    if len(x) == 0:
        ax.text(0.5, 0.5, "No valid points in view", ha="center", va="center")
        return fig, ax

    # Category styling setup
    category_to_marker = {}
    category_to_color = {}
    category_masks = {}

    if category_col is not None and category_col in work_df.columns:
        cat_series = work_df[category_col].fillna("(none)")
        rank_col = f"{category_col}_rank"

        if rank_col in work_df.columns:
            rank_lookup = work_df.groupby(category_col)[rank_col].min()
            top_categories = rank_lookup.sort_values().head(top_k).index.tolist()
        else:
            top_categories = cat_series.value_counts().nlargest(top_k).index.tolist()

        cmap = plt.colormaps["tab10"]
        category_to_marker = {cat: m for cat, m in zip(top_categories, marker_cycle)}
        category_to_color = {cat: cmap(i) for i, cat in enumerate(top_categories)}

        cat_values = cat_series.to_numpy()
        category_masks = {cat: (cat_values == cat) for cat in top_categories}
        category_masks["other"] = ~np.isin(cat_values, top_categories)
    else:
        category_masks = {"other": np.ones(len(x), dtype=bool)}

    # Fix Axes Boundaries
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)

    # Calculate line vector in data space, taking aspect ratio into account
    theta = np.deg2rad(angle)
    p0 = ax.transData.transform((0, 0))
    px = ax.transData.transform((1, 0))
    py = ax.transData.transform((0, 1))
    px_per_x = px[0] - p0[0]
    px_per_y = py[1] - p0[1]

    dx_data = np.cos(theta) / px_per_x
    dy_data = np.sin(theta) / px_per_y
    norm = np.hypot(dx_data, dy_data)
    dx_data, dy_data = dx_data / norm, dy_data / norm
    effective_angle = np.degrees(np.arctan2(dy_data, dx_data))

    # Selection masking based on line position (x_sect, y_sect)
    mask = compute_mask(x, y, x_sect, y_sect, effective_angle)
    in_bounds = (x >= left) & (x <= right) & (y >= bottom) & (y <= top)
    selected_in_view = mask & in_bounds

    # Scatter points
    for cat, cat_mask in category_masks.items():
        if not cat_mask.any():
            continue
        marker = category_to_marker.get(cat, other_marker)
        color = category_to_color.get(cat, other_color)
        sel = selected_in_view[cat_mask]

        ax.scatter(
            x[cat_mask],
            y[cat_mask],
            s=np.where(sel, 48, 20),
            c=[color],
            alpha=np.where(sel, 0.95, 0.75),
            marker=marker,
            zorder=2 if cat != "other" else 1,
            label=str(cat) if cat != "other" else "other",
            edgecolors="black" if sel.any() else "none",
            linewidths=np.where(sel, 0.6, 0),
        )

    # Draw Bisector Line and Drag Origin Handle at (x_sect, y_sect)
    t_vals = np.array([-1e6, 1e6])
    ax.plot(
        x_sect + t_vals * dx_data,
        y_sect + t_vals * dy_data,
        "r-",
        linewidth=2,
        zorder=3,
    )
    ax.plot(x_sect, y_sect, "ro", markersize=8, zorder=4)

    ax.set_title(f"Selected in view: {selected_in_view.sum()} | Angle: {angle:.0f}\u00b0")
    ax.set_xlabel("z_score")
    ax.set_ylabel("query_coverage")

    if category_to_marker:
        ax.legend(loc="best", fontsize=8, title=category_col, markerscale=1.2)

    return fig, ax


# =====================================================================
# 2. Stateful Selection Widget Manager
# =====================================================================

class BisectorSelectionWidget:
    """
    Interactive bisector selector.

    Controls are ordinary ipywidgets displayed separately from the
    matplotlib widget canvas, avoiding layout issues with ipympl.
    """

    HANDLE_HIT_RADIUS_PX = 15  # click tolerance for grabbing the center handle

    def __init__(self, project, wrapped_plot_fn, controls, fixed_params=None):
        self.project = project
        self.plot_fn = wrapped_plot_fn
        self.controls = controls
        self.fixed = fixed_params or {}
        self.current_selection = None
        self.df = None          # current dataframe for the active view
        self.work_df = None     # finite-coordinate subset
        self.x = None
        self.y = None

        self._dragging = False
        self._drag_target = None

        # Register callbacks. `view_name` gets a dedicated handler that also
        # refreshes the axis-range sliders (left/right/bottom/top/x_sect/
        # y_sect) to match the newly selected view's actual data range —
        # otherwise those sliders keep whatever bounds were computed for
        # the ORIGINAL view at construction time, and switching views just
        # shows a near-empty/near-identical plot (new points falling
        # outside the stale axis limits) rather than visibly updating.
        for key, widget in self.controls.items():
            if key == "view_name":
                widget.observe(self._on_view_change, names="value")
            else:
                widget.observe(self._on_control_change, names="value")

        # Create figure with auto-display suppressed. Without plt.ioff(),
        # plt.subplots() triggers matplotlib's interactive-mode auto-display
        # in addition to the explicit canvas display below,
        # producing a duplicate ("shadow") copy of the canvas in the output.
        with plt.ioff():
            self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self._connect_events()
        self.redraw()

        # Display controls
        display(self.controls["view_name"])
        display(widgets.HBox([
            self.controls["left"],
            self.controls["right"],
        ]))
        display(widgets.HBox([
            self.controls["bottom"],
            self.controls["top"],
        ]))
        display(widgets.HBox([
            self.controls["x_sect"],
            self.controls["y_sect"],
        ]))
        display(self.controls["angle"])

        # Embed the ipympl canvas to avoid the direct-display PNG preview
        # involved in the reported detached-artist drawing error.
        canvas = self.fig.canvas
        if isinstance(canvas, widgets.Widget):
            display(widgets.Box([canvas]))
        else:
            display(canvas)

    def _connect_events(self):
        canvas = self.fig.canvas
        canvas.mpl_connect("button_press_event", self._on_press)
        canvas.mpl_connect("motion_notify_event", self._on_motion)
        canvas.mpl_connect("button_release_event", self._on_release)

    def _on_press(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return

        x_sect = self.controls["x_sect"].value
        y_sect = self.controls["y_sect"].value

        # Hit-test in PIXEL space, not raw data units — z_score and
        # query_coverage live on very different scales, so a data-space
        # hypot() would be dominated by whichever axis has larger numbers.
        px_per_x, px_per_y = self._pixel_scale()
        dx_px = (event.xdata - x_sect) * px_per_x
        dy_px = (event.ydata - y_sect) * px_per_y
        dist_px = np.hypot(dx_px, dy_px)

        self._dragging = True
        self._drag_target = "center" if dist_px < self.HANDLE_HIT_RADIUS_PX else "angle"

    def _on_motion(self, event):
        if (
            not self._dragging
            or event.inaxes != self.ax
            or event.xdata is None
            or event.ydata is None
        ):
            return

        if self._drag_target == "center":
            self.controls["x_sect"].value = float(
                np.clip(
                    event.xdata,
                    self.controls["left"].value,
                    self.controls["right"].value,
                )
            )
            self.controls["y_sect"].value = float(
                np.clip(
                    event.ydata,
                    self.controls["bottom"].value,
                    self.controls["top"].value,
                )
            )

        else:
            x0 = self.controls["x_sect"].value
            y0 = self.controls["y_sect"].value

            # Convert the data-space direction to the mouse into the same
            # pixel-space angle convention plot_bisector() uses when drawing
            # (see its `effective_angle` calculation) — otherwise the line's
            # rendered direction won't match where the mouse is dragged to.
            px_per_x, px_per_y = self._pixel_scale()
            dx_screen = (event.xdata - x0) * px_per_x
            dy_screen = (event.ydata - y0) * px_per_y
            angle = np.degrees(np.arctan2(dy_screen, dx_screen))

            self.controls["angle"].value = float(np.clip(angle, -90.0, 90.0))

    def _on_release(self, event):
        self._dragging = False
        self._drag_target = None

    def _on_view_change(self, change):
        """
        Fires when the view dropdown changes. Recomputes left/right/bottom/
        top/x_sect/y_sect for the NEW view's actual data range (same
        calculation plane_bisector() does at construction time), then
        redraws once. Without this, those sliders would keep the previous
        view's bounds, and the plot would appear not to update — or worse,
        show a blank canvas with the bisector line off-screen, if the new
        view's data range doesn't overlap the stale one at all.
        """
        view_name = change["new"]
        df = self.project.views[view_name]
        valid = np.isfinite(df["z_score"]) & np.isfinite(df["query_coverage"])
        work_df = df.loc[valid]

        slider_keys = ["left", "right", "bottom", "top", "x_sect", "y_sect"]

        if len(work_df) == 0:
            print(f"Warning: view '{view_name}' has no valid (finite z_score, "
                  f"query_coverage) points — axis bounds left unchanged.")
            self.redraw()
            return

        x_min, x_max = float(work_df["z_score"].min()), float(work_df["z_score"].max())
        y_min, y_max = float(work_df["query_coverage"].min()), float(work_df["query_coverage"].max())
        step_x = round((x_max - x_min) / 100, 2) or 0.05
        step_y = round((y_max - y_min) / 100, 2) or 0.01

        new_specs = {
            "left":   (x_min, x_max, x_min, step_x),
            "right":  (x_min, x_max, x_max, step_x),
            "bottom": (y_min, y_max, y_min, step_y),
            "top":    (y_min, y_max, y_max, step_y),
            "x_sect": (x_min, x_max, (x_min + x_max) / 2.0, step_x),
            "y_sect": (y_min, y_max, (y_min + y_max) / 2.0, step_y),
        }

        for key in slider_keys:
            self.controls[key].unobserve(self._on_control_change, names="value")

        try:
            for key, (new_min, new_max, new_value, new_step) in new_specs.items():
                w = self.controls[key]
                # Widen to a huge, unconditionally-safe range before
                # assigning value, so the assignment can NEVER land outside
                # [min, max] regardless of how the old and new bounds
                # relate to each other — avoids having to reason about
                # old-vs-new ordering, which is easy to get subtly wrong.
                w.min = -1e12
                w.max = 1e12
                w.value = float(new_value)
                w.step = float(new_step)
                w.min = float(new_min)
                w.max = float(new_max)
        except Exception as e:
            print(f"\u26a0\ufe0f Failed to update axis-range sliders for view '{view_name}': {e}")
        finally:
            # ALWAYS re-observe, even if the update above raised partway
            # through — otherwise a single failure here would silently
            # disable redraw for every future slider interaction, not just
            # this one view switch.
            for key in slider_keys:
                self.controls[key].observe(self._on_control_change, names="value")

        self.redraw()

    def _on_control_change(self, change):
        self.redraw()

    def redraw(self):
        self.ax.clear()

        kwargs = {
            k: w.value
            for k, w in self.controls.items()
            if k != "view_name"
        }
        kwargs.update(self.fixed)

        view_name = self.controls["view_name"].value

        # Current dataframe for the active view (read directly — do NOT
        # route through plot_fn here, which returns a (fig, ax) tuple, not
        # a dataframe).
        self.df = self.project.views[view_name]

        valid = (
            np.isfinite(self.df["z_score"])
            & np.isfinite(self.df["query_coverage"])
        )

        self.work_df = self.df.loc[valid].copy()
        self.x = self.work_df["z_score"].to_numpy(dtype=float)
        self.y = self.work_df["query_coverage"].to_numpy(dtype=float)

        self.plot_fn(
            view_name=view_name,
            ax=self.ax,
            **kwargs,
        )

        self.fig.canvas.draw_idle()

    def get_selected(self):
        """
        Return
        ------
        selection : pandas.DataFrame
            Rows currently selected by the bisector, restricted to the
            visible crop box (left/right/bottom/top) — matching exactly
            what is highlighted on screen.

        params : dict
            Current widget parameters for provenance.
        """
        state = self.get_state()

        theta = np.deg2rad(state["angle"])

        px_per_x, px_per_y = self._pixel_scale()

        dx_data = np.cos(theta) / px_per_x
        dy_data = np.sin(theta) / px_per_y

        norm = np.hypot(dx_data, dy_data)
        dx_data /= norm
        dy_data /= norm

        effective_angle = np.degrees(np.arctan2(dy_data, dx_data))

        mask = compute_mask(
            self.x,
            self.y,
            state["x_sect"],
            state["y_sect"],
            effective_angle,
        )

        in_bounds = (
            (self.x >= state["left"]) & (self.x <= state["right"])
            & (self.y >= state["bottom"]) & (self.y <= state["top"])
        )

        selection = self.work_df.loc[mask & in_bounds].copy()

        return selection, state

    def _pixel_scale(self):
        p0 = self.ax.transData.transform((0, 0))
        px = self.ax.transData.transform((1, 0))
        py = self.ax.transData.transform((0, 1))

        return (
            px[0] - p0[0],
            py[1] - p0[1],
        )

    def get_state(self):
        state = {k: w.value for k, w in self.controls.items()}
        state.update(self.fixed)
        return state


# =====================================================================
# 3. Launcher Wrapper
# =====================================================================

def plane_bisector(
    project,
    df=None,
    category_col="pfam",
    top_k=6,
    title="Plane Bisector View",
    return_state=False,
    initial_params=None,
    **kwargs,
):
    """
    Custom launcher that encapsulates dragging logic while matching the
    (root, get_state, controls) return contract when return_state=True.

    initial_params : dict, optional
        A previously-captured params dict (e.g. from bisector.get_selected()
        in an earlier session) used to pre-seed the sliders' starting
        values, so the widget opens already positioned where you left off
        instead of at its plain defaults. Any keys not present fall back to
        the normal defaults. Purely a starting point — dragging the sliders
        still overrides it, same as usual.
    """
    initial_params = initial_params or {}

    view_keys = list(project.views.keys())
    if df is None:
        df = project.views[view_keys[0]]

    # Default the view dropdown to whichever registered view actually
    # matches the caller-supplied df, so the initial slider bounds (computed
    # from df, below) always agree with the initially-displayed view.
    default_view_name = next(
        (name for name, view_df in project.views.items() if view_df is df),
        view_keys[0],
    )

    valid = np.isfinite(df["z_score"]) & np.isfinite(df["query_coverage"])
    work_df = df.loc[valid]

    x_min, x_max = float(work_df["z_score"].min()), float(work_df["z_score"].max())
    y_min, y_max = float(work_df["query_coverage"].min()), float(work_df["query_coverage"].max())

    step_x = round((x_max - x_min) / 100, 2) or 0.05
    step_y = round((y_max - y_min) / 100, 2) or 0.01

    controls = {
        "view_name": widgets.Dropdown(
            options=view_keys,
            value=default_view_name,
            description="View:",
        ),
        "left": widgets.FloatSlider(
            min=x_min, max=x_max, value=x_min, step=step_x, description="Left"
        ),
        "right": widgets.FloatSlider(
            min=x_min, max=x_max, value=x_max, step=step_x, description="Right"
        ),
        "bottom": widgets.FloatSlider(
            min=y_min, max=y_max, value=y_min, step=step_y, description="Bottom"
        ),
        "top": widgets.FloatSlider(
            min=y_min, max=y_max, value=y_max, step=step_y, description="Top"
        ),
        "x_sect": widgets.FloatSlider(
            min=x_min, max=x_max, value=(x_min + x_max) / 2.0, step=step_x, description="X Center"
        ),
        "y_sect": widgets.FloatSlider(
            min=y_min, max=y_max, value=(y_min + y_max) / 2.0, step=step_y, description="Y Center"
        ),
        "angle": widgets.FloatSlider(
            min=-90.0, max=90.0, value=-90.0, step=1.0, description="Angle"
        ),
    }

    wrapped_plot_fn = make_view_plot_fn(plot_bisector, project)

    widget_inst = BisectorSelectionWidget(
        project=project,
        wrapped_plot_fn=wrapped_plot_fn,
        controls=controls,
        fixed_params={"category_col": category_col, "top_k": top_k},
    )

    return widget_inst
