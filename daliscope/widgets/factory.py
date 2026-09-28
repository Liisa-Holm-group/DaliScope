import ipywidgets as widgets
from IPython.display import display, clear_output
import matplotlib.pyplot as plt  # <--- Ensure plt is accessible
import matplotlib

# Inside daliscope/widgets/factory.py
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display

import daliscope.analysis.wrappers

def launch_interactive_table_viewer(func, title, project, params=None, fixed=None, **kwargs):
    """
    Specialized launcher for windowed dataframe table viewing.
    Configures default sliders for `start_row` and `window` and hands off to the master launcher.
    """
    table_params = {
        'start_row': {
            'type': 'intslider',
            'value': 0,
            'min': 0,
            'max': 10000,
            'step': 10,
            'label': 'Start Row',
        },
        'window': {
            'type': 'intslider',
            'value': 25,
            'min': 5,
            'max': 200,
            'step': 5,
            'label': 'Window Size',
        }
    }

    # Safely unpack and inject user-customized params from notebook
    if params:
        for param_key, param_config in params.items():
            table_params[param_key] = param_config

    # Hand off to master launcher—it automatically attaches the view_name dropdown
    launch_interactive_view(
        func=func,
        title=title,
        params=table_params,
        fixed=fixed,
        project=project
    )


def launch_interactive_sunburst_plot(func, title, project, params=None, fixed=None, **kwargs):
    """
    Specialized wrapper for Sunburst plots. Configures default parameters
    and delegates view management to the master launcher.
    """
    sunburst_params = {
        'k': {
            'type': 'intslider',
            'value': 15,
            'min': 1,
            'max': 50,
            'label': 'Top-k Cutoff',
        }
    }

    # Safely unpack and inject user-customized params from notebook
    if params:
        for param_key, param_config in params.items():
            sunburst_params[param_key] = param_config

    # Hand off to master launcher—it patches in view_name automatically
    launch_interactive_view(
        func=func,
        title=title,
        params=sunburst_params,
        fixed=fixed,
        project=project
    )


def with_refresh(source_fn):
    """Returns a wrap_fn for control_wrappers — places a refresh
    button immediately to the right of the control it wraps."""
    def wrap_fn(ctrl):
        return widgets.HBox([ctrl, refresh_button(ctrl, source_fn)])
    return wrap_fn

def launch_msa_widget(func, query_length, title=""):
    """Convenience function to call widget_view - no dropdowns

    func: msa function (msa.heatmap or msa.logo
    title: text

    Example:
    
    a = msa(CURRENT_DF,'dssp_pileup')
    daliscope.widgets.factory.launch_msa_widget(a.heatmap, project.query_length, title="Stacked alignment colored by DSSP")

    """
    widget_view(
        func = func, # sequ_logo.heatmap, sequ_logo.logo
        title = title,
        params = {
             "left": {
                "type": "intslider", "value": 0,
                "min": 0, "max": query_length-1,
                "less_than": "right",
            },
            "right": {
                "type": "intslider", "value": query_length,
                "min": 1, "max": query_length,
            },
           },
        show_default_buttons = False
    )


def launch_interactive_view(func, title, params, fixed=None, control_wrappers=None, project=None, renderer_name=None, **kwargs):
    """
    Orchestrates the widget lifecycle. Automatically injects the global 'view_name'
    dropdown and standardizes the 'project' context inside the fixed parameters.
    """
    plt.ioff()

    # Initialize dictionaries safely
    params = params.copy() if params else {}
    control_wrappers = control_wrappers or {}

    # GLOBAL 'fixed' INJECTION: Always include project first, then layer extra fixed variables
    ordered_fixed = {}
    if project is not None:
        ordered_fixed['project'] = project

    if fixed:
        ordered_fixed.update(fixed)
    fixed = ordered_fixed

    # GLOBAL 'view_name' INJECTION: Auto-inject view_name dropdown if we have a project context
    if project is not None and 'view_name' not in params:
        # 1. Create a new dict with 'view_name' as the FIRST element
        ordered_params = {
            'view_name': {
                'type': 'dropdown',
                'options': list(project.views.keys()),
                'value': list(project.views.keys())[0],
                'label': 'View',
            }
        }
        # 2. Merge the rest of the existing plot-specific params behind it
        ordered_params.update(params)
        params = ordered_params

        # Auto-inject the matching dynamic refresh rule
        if 'view_name' not in control_wrappers:
            control_wrappers['view_name'] = with_refresh(lambda: list(project.views.keys()))

    # Invoke the core factory engine
    wv, get_params, controls = widget_view(
        func=func,
        title=title,
        fixed=fixed,
        params=params,
        control_wrappers=control_wrappers,
        show_default_buttons=False,
        display_result=False,
        return_state=True
    )

    plt.close('all')
    plt.ion()

    components = [wv]
    if project is not None:
        components.append(
            bookmark_button(
                get_params,
                project,
                renderer_name=renderer_name or f"{func.__module__}.{func.__name__}"
            )
        )

    display(widgets.VBox(components))

def launch_interactive_scatter_plot(func, title, project, params=None, fixed=None, **kwargs):
    """
    Specialized wrapper for scatter plots. Standardizes the 'group_by' resolution
    while letting the master launcher handle the global 'view_name' logic.
    
    Accepts `params`, a dictionary of dictionaries containing extra widget configurations 
    (e.g., params={'k': k_dict}).
    """
    # 1. Get the first view's DataFrame dynamically
    first_view_name = list(project.views.keys())[0]
    first_view_df = project.views[first_view_name]

    # 2. Extract and sort the numeric columns from this specific DataFrame
    numeric_cols = sorted(first_view_df.select_dtypes(include='number').columns.tolist())
    categorical_options = ['none', 'pfam', 'clan']
    size_options      = ['none', 'sequence_identity', 'query_coverage']

    # 1. Initialize with the baseline scatter plot configuration
    scatter_params = {
        'view_name': {
            'type':    'dropdown',
            'options': list(project.views.keys()),
            'value':   list(project.views.keys())[0],
            'label':   'View',
        },
        'x': {
            'type': 'dropdown', 'options': numeric_cols,
            'label': 'X', 'value': 'z_score',
        },
        'y': {
            'type': 'dropdown', 'options': numeric_cols,
            'label': 'Y', 'value': 'alignment_length',
        },
        'color_by': {
            'type': 'dropdown', 'options': categorical_options,
            'value': 'clan', 'label': 'Color by',
        },
        'marker_by': {
            'type': 'dropdown', 'options': categorical_options,
            'value': 'none', 'label': 'Marker by',
        },
        'size_by': {
            'type': 'dropdown', 'options': size_options,
            'value': 'none', 'label': 'Size by',
        },
    }

    # 2. Safely unpack and inject any additional nested dictionaries passed in
    if params:
        for param_key, param_config in params.items():
            scatter_params[param_key] = param_config

    widget_view(
        func   = daliscope.analysis.wrappers.make_view_pfam_plot_fn(daliscope.optics.plotly.scatter_plot, project),
        title = title,
        params = scatter_params,
        show_default_buttons = False,
        control_wrappers = {
            'view_name': with_refresh(lambda: list(project.views.keys())),
        },
    )

def launch_interactive_violin_plot(func, title, project, params=None, fixed=None, **kwargs):
    """
    Specialized wrapper for violin plots. Standardizes the 'group_by' resolution
    while letting the master launcher handle the global 'view_name' logic.
    
    Accepts `params`, a dictionary of dictionaries containing extra widget configurations 
    (e.g., params={'k': k_dict, 'l': l_dict}).
    """
    # 1. Initialize with the baseline violin plot configuration
    violin_params = {
        'group_by': {
            'type': 'dropdown',
            'options': ['clan', 'pfam'],
            'value': 'clan',
            'label': 'Group Resolution',
        }
    }

    # 2. Safely unpack and inject any additional nested dictionaries passed in
    if params:
        for param_key, param_config in params.items():
            violin_params[param_key] = param_config

    # Hand off to master launcher—it automatically detects `project`
    # and patches in the view_name logic behind the scenes!
    launch_interactive_view(
        func=func,
        title=title,
        params=violin_params,
        fixed=fixed,
        project=project
    )

def launch_interactive_family_presence(func, title, project, params=None, fixed=None, **kwargs):
    """
    Specialized wrapper for the family presence barchart. 
    Delegates to the master launcher to handle the target 'view_name' 
    dropdown and native refresh button placement automatically.
    """
    # 1. Initialize baseline widget configuration
    barchart_params = {
        'category_col': {
            'type': 'dropdown', 'options': ['pfam', 'clan'],
            'value': 'pfam', 'label': 'Category'
        },
        'k': {
            'type': 'intslider', 'value': 10,
            'min': 1, 'max': 40, 'label': 'Rank cutoff (k)'
        },
        'show_zero_present': {
            'type': 'checkbox', 'value': True,
            'label': 'Show zero-present'
        }
    }

    # 2. Safely unpack any additional nested dictionaries
    if params:
        for param_key, param_config in params.items():
            barchart_params[param_key] = param_config

    # 3. Auto-guess the population view if not explicitly provided
    fixed = fixed or {}
    if 'population_view' not in fixed:
        view_keys = list(project.views.keys())
        pop_view = next((k for k in view_keys if 'FULL' in k.upper()), view_keys[0] if view_keys else None)
        fixed['population_view'] = pop_view

    # Hand off to master launcher
    launch_interactive_view(
        func=func,
        title=title,
        params=barchart_params,
        fixed=fixed,
        project=project
    )


def launch_interactive_architecture_plot(func, title, project, params=None, fixed=None, **kwargs):
    """
    Specialized wrapper for domain architecture cartoons. Standardizes the layout,
    bookmark rendering, and view resolution while keeping plot parameters flexible.
    
    Accepts `params`, a dictionary of dictionaries containing extra widget configurations 
    (e.g., params={'pfam_col': pfam_dict, 'k': k_dict}).
    """

    # 1. Initialize with the baseline violin plot configuration
    arch_params = {
        "pfam_col": {
            "type": "dropdown",
            "options": ["clan_domains", "pfam_domains"],
            "value": "clan_domains",
            "label": "Resolution",
        },
        'k': {
            'type': 'intslider',
            'value': 10,
            'min': 1,
            'max': 30,
            'label': 'Top-k Ranks',
        }
    }

    # 2. Safely unpack and inject any additional nested dictionaries passed in
    if params:
        for param_key, param_config in params.items():
            arch_params[param_key] = param_config

    # Hand off to master launcher—it automatically detects `project`
    # and patches in the view_name logic behind the scenes!
    launch_interactive_view(
        func=func,
        title=title,
        params=arch_params,
        fixed=fixed,
        project=project
    )

def widget_view(func, params, title=None, auto_update=True, fixed=None,
                display_result=True, return_state=False,
                control_wrappers=None, show_default_buttons=True):
    """
    Generic widget factory. Returns a VBox containing all controls
    and plot output. If display_result=True (default), also displays it.

    Parameters
    ----------
    func             : callable — called as func(**kwargs)
    params           : dict of parameter specifications
    title            : optional header string
    auto_update      : redraw on every control change (default True)
    fixed            : dict of extra kwargs always passed to func
    display_result   : if True, calls display() on the built widget
    return_state     : if True, returns (root, get_state, controls)
                       instead of just root
    control_wrappers : optional dict {param_name: wrap_fn}, where
                       wrap_fn(control_widget) -> a container widget
                       (e.g. HBox with a refresh button) placed in
                       the layout instead of the bare control
    show_default_buttons : if True (default), shows the built-in
                           "Full range" and "Reset" buttons below
                           the controls. Set False for widgets that
                           don't need them (e.g. a single dropdown
                           selector) to keep the layout minimal.

    Extended param spec fields (all optional):
    --------------------------------------------
    "default"      : reset value (falls back to "value")
    "less_than"    : name of another param this must be less than
    "greater_than" : name of another param this must be greater than
    "validate"     : callable(value) -> True if valid
    "validate_msg" : error message shown when validate fails
    """

    controls     = {}
    defaults     = {}
    output       = widgets.Output()
    validate_out = widgets.Output()

    # ---------- build controls ----------
    for name, spec in params.items():
        typ     = spec["type"]
        default = spec.get("default", spec.get("value"))
        defaults[name] = default

        if typ == "intslider":
            controls[name] = widgets.IntSlider(
                value=spec.get("value", spec["min"]), min=spec["min"],
                max=spec["max"], step=spec.get("step", 1),
                description=spec.get("label", name),
                continuous_update=spec.get("continuous_update", False),
                style={"description_width": "80px"},
                layout=widgets.Layout(width="400px"),
            )
        elif typ == "floatslider":
            controls[name] = widgets.FloatSlider(
                value=spec.get("value", spec["min"]), min=spec["min"],
                max=spec["max"], step=spec.get("step", 0.05),
                description=spec.get("label", name),
                continuous_update=spec.get("continuous_update", False),
                style={"description_width": "80px"},
                layout=widgets.Layout(width="400px"),
            )
        elif typ == "dropdown":
            controls[name] = widgets.Dropdown(
                options=spec["options"],
                value=spec.get("value", spec["options"][0]),
                description=spec.get("label", name),
                style={"description_width": "80px"},
                layout=widgets.Layout(width="400px"),
            )
        elif typ == "checkbox":
            controls[name] = widgets.Checkbox(
                value=spec.get("value", False),
                description=spec.get("label", name),
            )
        elif typ == "text":
            controls[name] = widgets.Text(
                value=spec.get("value", ""),
                description=spec.get("label", name),
                style={"description_width": "80px"},
                layout=widgets.Layout(width="400px"),
            )
        else:
            raise ValueError(f"Unknown widget type: {typ!r}")

    # ---------- constraint checker ----------D
    def _check_constraints():

        errors = []
        for name, spec in params.items():
            val = controls[name].value
            if "less_than" in spec:
                other = spec["less_than"]
                if other in controls and val >= controls[other].value:
                    errors.append(
                        f"⚠️  {name} ({val}) must be less than "
                        f"{other} ({controls[other].value})"
                    )
            if "greater_than" in spec:
                other = spec["greater_than"]
                if other in controls and val <= controls[other].value:
                    errors.append(
                        f"⚠️  {name} ({val}) must be greater than "
                        f"{other} ({controls[other].value})"
                    )
            if "validate" in spec:
                if not spec["validate"](val):
                    errors.append(spec.get(
                        "validate_msg", f"⚠️  {name} = {val!r} is invalid"
                    ))
        return errors

    # ---------- renderer ----------
    def update(_=None):
        validate_out.clear_output(wait=False)
        errors = _check_constraints()
        if errors:
            with validate_out:
                for e in errors:
                    print(e)
            return
        kwargs = {name: ctrl.value for name, ctrl in controls.items()}
        if fixed:
            kwargs.update(fixed)
        with output:
            clear_output(wait=True)
            plt.close('all')
            fig = func(**kwargs)
            if fig is not None:
                display(fig)

    # ---------- wire callbacks ----------
    if auto_update:
        for control in controls.values():
            control.observe(update, names="value")

    # ---------- buttons ----------
    full_range_btn = None
    reset_btn      = None

    if show_default_buttons:
        full_range_btn = widgets.Button(
            description="↺ Full range", layout=widgets.Layout(width="110px")
        )
        def on_full_range(_):
            for ctrl in controls.values():
                ctrl.unobserve_all()
            for name, ctrl in controls.items():
                if isinstance(ctrl, (widgets.IntSlider, widgets.FloatSlider)):
                    ctrl.value = ctrl.min if "less_than" in params.get(name, {}) \
                                else ctrl.max
            if auto_update:
                for ctrl in controls.values():
                    ctrl.observe(update, names="value")
            update()
        full_range_btn.on_click(on_full_range)

        reset_btn = widgets.Button(
            description="Reset", button_style="warning", icon="undo",
            layout=widgets.Layout(width="100px"),
        )
        def on_reset(_):
            for ctrl in controls.values():
                ctrl.unobserve_all()
            for name, ctrl in controls.items():
                if defaults[name] is not None:
                    ctrl.value = defaults[name]
            if auto_update:
                for ctrl in controls.values():
                    ctrl.observe(update, names="value")
            update()
        reset_btn.on_click(on_reset)

    # ---------- layout ----------
    control_wrappers = control_wrappers or {}
    rows = []
    for name, ctrl in controls.items():
        if name in control_wrappers:
            rows.append(control_wrappers[name](ctrl))
        else:
            rows.append(ctrl)

    children = []
    if title:
        children.append(widgets.HTML(f"<b>{title}</b>"))
    children.extend(rows)
    if show_default_buttons:
        children.append(widgets.HBox([full_range_btn, reset_btn]))
    children.append(validate_out)
    children.append(output)

    root = widgets.VBox(children)

    # ---------- get_state — simple dict, nothing else ----------
    def get_state() -> dict:
        kwargs = {name: ctrl.value for name, ctrl in controls.items()}
        if fixed:
            kwargs.update(fixed)
        return kwargs

    # ---------- initial render ----------
    update()
    if display_result:
        display(root)

    # --- Cleaned up return contracts ---
    if return_state:
        return root, get_state, controls

    if not display_result:
        return root

    # If it was already displayed explicitly, return None to stop Jupyter from drawing it twice
    return None

#########################

from daliscope.optics.msa import _get_cached_msa

def make_msa_plot_fn(project):
    """
    Returns a widget_view-compatible function for msa heatmap/logo.
    Uses _get_cached_msa() so the expensive _initialize_everything()
    only runs once per (view, pileup_col) combination.
    """

    # This signature will NEVER crash due to extra arguments
    def msa_plot(view_name, pileup_col, plot_type, left, right, **kwargs):
        df  = project.views[view_name]
        viz = _get_cached_msa(df, pileup_col)   # cached — cheap if seen before
        print(viz)

        if plot_type == 'heatmap':
            # 1. Temporarily disable interactive plotting to catch the figure
            plt.ioff()

            # 2. Let the library do the heavy lifting of building and plotting the heatmap
            viz.heatmap(left=left, right=right)

            # 3. Capture the figure that viz.heatmap() just created
            fig = plt.gcf()

            # 4. Re-enable interactive mode so the rest of your notebook behaves normally
            plt.ion()

            return fig   # Hand the clean figure over to your widget container

        else:   # logo
            viz.logo(left=left, right=right)
            return None

    return msa_plot


def launch_interactive_msa_plot(title, project):
    """
    Launches a fully interactive MSA heatmap/logo widget 
    with bookmarking capabilities.
    """
    # 1. Generate the 5-argument coordinator function
    coordinating_func = make_msa_plot_fn(project)

    # 2. Define parameters matching the 5 arguments of coordinating_func
    msa_params = {
        'pileup_col': {
            'type':    'dropdown',
            'options': ['dssp_pileup', 'sequ_pileup'],
            'value':   'dssp_pileup',
            'label':   'Alignment',
        },
        'plot_type': {
            'type':    'dropdown',
            'options': ['heatmap', 'logo'],
            'value':   'heatmap',
            'label':   'Plot type',
        },
        'left': {
            'type':      'intslider', 'value': 0,
            'min': 0,    'max': project.query_length - 1,
            'label':     'Left',
            'less_than': 'right',
        },
        'right': {
            'type':         'intslider', 'value': project.query_length,
            'min': 1,       'max': project.query_length,
            'label':        'Right',
            'greater_than': 'left',
        },
    }

    # 3. Launch the interactive panel
    launch_interactive_view(
        func=coordinating_func, # <--- This handles all 5 inputs dynamically!
        title=title,
        params=msa_params,
        project=project,
        control_wrappers={
            'view_name': with_refresh(lambda: list(project.views.keys())),
        },
        show_default_buttons = True
    )

######################

def refresh_button(dropdown_control: widgets.Dropdown,
                   source_fn) -> widgets.Button:
    """
    A button that re-pulls options for one specific dropdown from
    source_fn() and preserves the current selection if still valid.
    Knows only about: one dropdown, one callable. Nothing else.
    """
    btn = widgets.Button(
        description="↺ Refresh",
        layout=widgets.Layout(width="100px")
    )

    def on_click(_):
        current = dropdown_control.value
        new_options = source_fn()
        dropdown_control.options = new_options

        # Proper inline ternary assignment with clean indentation
        dropdown_control.value = (
            current if current in new_options
            else (new_options[0] if new_options else None)
        )

    btn.on_click(on_click)
    return btn

def bookmark_button(get_current_params_fn,
                    project: 'Project',
                    renderer_name: str) -> widgets.VBox:
    """
    A save-to-report button with name input and caption textarea.
    Knows about: project (for bookmark_plot), current params
    (from a callable, not a static snapshot).
    """
    name_input = widgets.Text(placeholder='figure name',
                             layout=widgets.Layout(width='250px'))
    caption    = widgets.Textarea(
        placeholder='Write your observations — becomes the figure caption.',
        layout=widgets.Layout(width='500px', height='70px')
    )
    btn        = widgets.Button(description='📌 Save to report',
                               button_style='success',
                               layout=widgets.Layout(width='150px'))
    status     = widgets.Output()

    def on_save(_):
        with status:
            clear_output(wait=True)
            name = name_input.value.strip()
            if not name:
                print("❌ Give this figure a name"); return
            if name in project.plots:
                print(f"⚠️  '{name}' already saved"); return
            project.bookmark_plot(
                name     = name,
                renderer = renderer_name,
                source   = get_current_params_fn().get('view_name', ''),
                params   = get_current_params_fn(),
                caption  = caption.value.strip(),
            )
            print(f"✓ Saved '{name}' to project.plots")

    btn.on_click(on_save)
    return widgets.VBox([
        widgets.HTML("<hr style='margin:6px 0'/>"),
        widgets.HBox([name_input, btn]),
        caption,
        status,
    ])


def register_button(get_current_params_fn,
                    project: 'Project',
                    register_fn) -> widgets.VBox:
    """
    A register-view button for filtering widgets.
    Knows about: project, the registration function to call,
    and whatever params the caller wants passed to it.
    """
    name_input = widgets.Text(placeholder='new view name',
                             layout=widgets.Layout(width='250px'))
    btn        = widgets.Button(description='Register view',
                               button_style='warning',
                               layout=widgets.Layout(width='140px'))
    status     = widgets.Output()

    def on_register(_):
        with status:
            clear_output(wait=True)
            name = name_input.value.strip()
            if not name:
                print("❌ Give this view a name"); return
            if name in project.views:
                print(f"⚠️  '{name}' already exists"); return
            register_fn(name=name, **get_current_params_fn())
            print(f"✓ Registered view '{name}'")

    btn.on_click(on_register)
    return widgets.VBox([
        widgets.HTML("<hr style='margin:6px 0'/>"),
        widgets.HBox([name_input, btn]),
        status,
    ])
