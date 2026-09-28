import ipywidgets as widgets
from IPython.display import display, clear_output
from daliscope.mechanics.metrics import (
    clip_segments,
    parse_domain_ranges,
)
from daliscope.optics.py3Dmol_viewer import _make_domain_preview_str

def domain_preview_widget(project) -> dict:
    preview_out = widgets.Output()
    status_out  = widgets.Output()

    domain_input = widgets.Text(
        value       = f"0-{project.query_length}",
        description = "Domains:",
        style       = {"description_width": "80px"},
        layout      = widgets.Layout(width="400px"),
        placeholder = project._query_meta.get("domain_string", "e.g. 30-335, 400-520"),
    )
    preview_btn = widgets.Button(
        description  = "Preview domain",
        button_style = "info",
        icon         = "eye",
        layout       = widgets.Layout(width="160px"),
    )

    def on_preview(_):
        with status_out:
            clear_output(wait=True)
        with preview_out:
            clear_output(wait=True)

        raw = domain_input.value.strip()
        if not raw:
            with status_out:
                print("❌ Enter a domain range first (e.g. 30-335)")
            return

        try:
            ranges = parse_domain_ranges(raw)
        except ValueError as e:
            with status_out:
                print(f"❌ Domain range error: {e}")
            return

        try:
            clipped_df = clip_segments(project.segments, domain_ranges=ranges)
        except Exception as e:
            with status_out:
                print(f"❌ Clip error: {e}")
            return

        with status_out:
            print(f"✓ {clipped_df['alignment_id'].nunique()} alignments in {ranges}")

        with preview_out:
            try:
                _make_domain_preview_str(project.query_pdb_str, domain_ranges=ranges)
            except Exception as e:
                print(f"⚠️  py3Dmol preview failed: {e}")

    preview_btn.on_click(on_preview)

    box = widgets.VBox([
        widgets.HTML(
            "<span style='color:grey;font-size:0.85em'>"
            "Format: start-end, start-end &nbsp;(e.g. 30-335, 400-520). "
            "Leave empty for full query length."
            "</span>"
        ),
        widgets.HBox([domain_input, preview_btn]),
        status_out,
        preview_out,
    ], layout=widgets.Layout(border="1px solid #ccc", padding="12px", border_radius="6px"))

    display(box)

    controls = {'domain_range': domain_input}
    return controls


# ================================================================== #
# DomainSelectorWidget — Step 2, repeatable                          #
# Collects domain range from user, previews in py3Dmol,              #
# registers named view on Accept. Can be shown multiple times        #
# to register multiple domain views.                                 #
# ================================================================== #

class DomainSelectorWidget:
    """
    Step 2 — define and register a domain view.
    Repeatable: Accept can be clicked multiple times with different
    ranges to register multiple named views.

    Requires an active Project instance.

    Usage
    -----
        project = Project.load_pack("path/to/pack.tar.gz")

        selector = DomainSelectorWidget(project)
        selector.show()
        # ... user previews, accepts multiple domains ...

        df_view = project.views["dom_138-335"]
    """

    def __init__(self, project):
        """
        Parameters
        ----------
        project : Project
            An instance of the DaliScope Project class.
        """
        self.project = project

        # ---------------------------------------------------------- #
        # Public state                                               #
        # ---------------------------------------------------------- #
        self.views = {}   # name → metadata tracking dictionary

        # preview state — reset after each accept or range edit
        self._preview_clipped_df = None
        self._preview_ranges     = None

        # ---------------------------------------------------------- #
        # Widgets                                                    #
        # ---------------------------------------------------------- #
        self.domain_input = widgets.Text(
            value="",
            description="Domains:",
            layout=widgets.Layout(width="400px"),
            style={"description_width": "80px"},
            placeholder="e.g. 10-150, 220-340",
        )
        self.view_name_input = widgets.Text(
            value="",
            description="View name:",
            layout=widgets.Layout(width="400px"),
            style={"description_width": "80px"},
            placeholder="unique name for this domain selection",
        )
        self.preview_btn = widgets.Button(
            description="Preview domain",
            button_style="info",
            icon="eye",
            layout=widgets.Layout(width="160px"),
        )
        self.accept_btn = widgets.Button(
            description="Accept & Register",
            button_style="success",
            icon="check",
            layout=widgets.Layout(width="160px"),
            disabled=True,   # locked until preview succeeds
        )
        self.clip_status   = widgets.Output()
        self.accept_status = widgets.Output()
        self.preview_out   = widgets.Output()

        self.domain_input.observe(self._on_domain_changed, names="value")
        self.preview_btn.on_click(self._on_preview)
        self.accept_btn.on_click(self._on_accept)

        # ---------------------------------------------------------- #
        # Layout                                                     #
        # ---------------------------------------------------------- #
        self._box = widgets.VBox([
            widgets.HTML(
                "<h3>Step 2 — Define domain ranges</h3>"
                "<span style='color:grey;font-size:0.85em'>"
                "Format: start-end, start-end &nbsp;(e.g. 10-150, 220-340). "
                "Leave empty to use full query length."
                "</span>"
            ),
            self.domain_input,
            self.preview_btn,
            self.clip_status,
            self.preview_out,
            widgets.HTML("<hr style='margin:8px 0'/>"),
            self.view_name_input,
            self.accept_btn,
            self.accept_status,
        ], layout=widgets.Layout(
            border="1px solid #ccc",
            padding="12px",
            border_radius="6px",
        ))

    # -------------------------------------------------------------- #
    # Callbacks                                                      #
    # -------------------------------------------------------------- #

    def _on_domain_changed(self, change):
        """Auto-update view name; invalidate stale preview."""
        try:
            ranges = parse_domain_ranges(change["new"])
            self.view_name_input.value = \
                "dom_" + "_".join(f"{s}-{e}" for s, e in ranges)
        except ValueError:
            pass
        # any edit invalidates the current preview
        self.accept_btn.disabled  = True
        self._preview_clipped_df  = None
        self._preview_ranges      = None

    def _on_preview(self, _):
        """Clip segments and show py3Dmol preview — no registration."""
        with self.clip_status:
            clear_output(wait=True)
        with self.preview_out:
            clear_output(wait=True)

        with self.clip_status:
            if self.project is None:
                print("❌ Load a project first (Step 1)")
                return

            # default to full length if input is empty
            raw = self.domain_input.value.strip()
            if not raw:
                ranges = [(0, self.project.query_length)]
                self.domain_input.value    = f"0-{self.project.query_length}"
                self.view_name_input.value = "full"
            else:
                try:
                    ranges = parse_domain_ranges(raw)
                except ValueError as e:
                    print(f"❌ Domain range error: {e}")
                    return

            print(f"Previewing {ranges} ...")
            try:
                clipped_df = clip_segments(
                    self.project.segments,
                    domain_ranges=ranges,
                )
            except Exception as e:
                print(f"❌ Clip error: {e}")
                return

            self._preview_clipped_df = clipped_df
            self._preview_ranges      = ranges
            self.accept_btn.disabled = False
            print(f"✓ {clipped_df['alignment_id'].nunique()} alignments in selected domain.")
            print(f"  Happy with the selection? Click Accept & Register.")

        with self.preview_out:
            try:
                view = _make_domain_preview_str(
                    self.project.query_pdb_str,
                    domain_ranges=ranges,
                )
                display(view)
            except Exception as e:
                print(f"⚠️  py3Dmol preview failed: {e}")
                print(f"   Domain clipping succeeded — you can still Accept.")

    def _on_accept(self, _):
        """Register the previewed domain directly into the project's views."""
        with self.accept_status:
            clear_output(wait=True)

            if self._preview_clipped_df is None:
                print("❌ Run Preview first")
                return

            view_name = self.view_name_input.value.strip()
            if not view_name:
                print("❌ View name is empty")
                return

            if view_name in self.project.views:
                print(f"⚠️  '{view_name}' already registered — edit the view name to create a new view")
                return

            try:
                # Direct registration into project view manager via API
                self.project.add_view(
                    name=view_name,
                    df=self._preview_clipped_df,
                    parent="full",  # or whichever parent is standard for clipping
                    function="clip_segments",
                    parameters={"domain_ranges": self._preview_ranges}
                )
            except Exception as e:
                print(f"❌ Registration error: {e}")
                return

            self.views[view_name] = {
                "domain_ranges": self._preview_ranges,
            }

            # reset preview state — force re-preview for next domain
            self._preview_clipped_df = None
            self._preview_ranges      = None
            self.accept_btn.disabled = True

            print(f"✓ Registered view '{view_name}' successfully.")
            print(f"  Active project views: {list(self.project.views.keys())}")
            print(f"\nDefine another domain range above, or proceed to analysis.")

    # -------------------------------------------------------------- #
    # Display                                                        #
    # -------------------------------------------------------------- #

    def show(self):
        """Display the domain selector widget."""
        # set sensible default if project already loaded
        if self.project is not None and not self.domain_input.value:
            self.domain_input.value    = f"0-{self.project.query_length}"
            self.view_name_input.value = "full"
        display(self._box)
