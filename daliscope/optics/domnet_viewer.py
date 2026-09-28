import uuid
import anywidget
import ipywidgets as widgets
import traitlets
from IPython.display import display

from .domnet_viewer_core import DEFAULT_COLORS

import anywidget
import traitlets


class DomNetCanvas(anywidget.AnyWidget):
    """Pure ESM 3D Canvas widget with vertical domain bar, range labels, hover info,

    re-centering, and C-alpha trace sticks with variable radii for segments.
    """

    _esm = """
    export default {
      async render({ model, el }) {
        const width = model.get("width");
        const height = model.get("height");

        el.style.width = `${width}px`;
        el.style.height = `${height}px`;
        el.style.display = "block";
        el.style.position = "relative";

        // Main Wrapper (Horizontal layout: 3D canvas left, vertical sidebar right)
        const mainWrapper = document.createElement("div");
        mainWrapper.style.width = "100%";
        mainWrapper.style.height = "100%";
        mainWrapper.style.display = "flex";
        mainWrapper.style.flexDirection = "row";
        mainWrapper.style.position = "relative";
        el.appendChild(mainWrapper);

        // 3D Canvas Container
        const container = document.createElement("div");
        container.style.flex = "1";
        container.style.height = "100%";
        container.style.position = "relative";
        container.style.overflow = "hidden";
        mainWrapper.appendChild(container);

        // Vertical Domain Centering Sidebar
        const domainSidebar = document.createElement("div");
        domainSidebar.style.width = "200px";
        domainSidebar.style.height = "100%";
        domainSidebar.style.display = "flex";
        domainSidebar.style.flexDirection = "column";
        domainSidebar.style.gap = "6px";
        domainSidebar.style.padding = "10px";
        domainSidebar.style.background = "#f8f9fa";
        domainSidebar.style.borderLeft = "1px solid #ddd";
        domainSidebar.style.overflowY = "auto";
        domainSidebar.style.boxSizing = "border-box";
        mainWrapper.appendChild(domainSidebar);

        // Tooltip element for hovering
        const tooltip = document.createElement("div");
        tooltip.style.position = "absolute";
        tooltip.style.padding = "4px 8px";
        tooltip.style.background = "rgba(0, 0, 0, 0.8)";
        tooltip.style.color = "#ffffff";
        tooltip.style.borderRadius = "4px";
        tooltip.style.fontSize = "12px";
        tooltip.style.fontFamily = "sans-serif";
        tooltip.style.pointerEvents = "none";
        tooltip.style.display = "none";
        tooltip.style.zIndex = "9999";
        container.appendChild(tooltip);

        // Load 3Dmol.js dynamically
        if (!window.$3Dmol) {
          await new Promise((resolve, reject) => {
            const script = document.createElement("script");
            script.src = "https://3Dmol.org/build/3Dmol-min.js";
            script.onload = resolve;
            script.onerror = reject;
            document.head.appendChild(script);
          });
        }

        const viewer = $3Dmol.createViewer(container, {
          backgroundColor: "white"
        });

        // Helper: Parse domain string into range structures
        function parseDomainString(domainStr) {
          if (!domainStr) return [];
          const colors = model.get("colors") || [
            "#3366cc", "#dc3912", "#ff9900", "#109618", "#990099",
            "#0099c6", "#dd4477", "#66aa00", "#b82e2e", "#316395"
          ];

          const domainGroups = domainStr.split(",").map(s => s.trim()).filter(Boolean);
          const domainData = [];

          domainGroups.forEach((group, domainIdx) => {
            const color = colors[domainIdx % colors.length];
            const segments = group.split("_").map(s => s.trim()).filter(Boolean);
            const ranges = [];

            segments.forEach(seg => {
              const parts = seg.split("-").map(Number);
              if (parts.length === 2 && !isNaN(parts[0]) && !isNaN(parts[1])) {
                ranges.push({ start: parts[0], end: parts[1] });
              }
            });

            if (ranges.length > 0) {
              domainData.push({
                id: domainIdx + 1,
                rangeStr: group,
                color: color,
                ranges: ranges
              });
            }
          });

          return domainData;
        }

        // Render Vertical Domain Buttons
        function renderDomainBadges(domainData) {
          domainSidebar.innerHTML = "";
          if (domainData.length === 0) return;

          const label = document.createElement("div");
          label.style.fontSize = "12px";
          label.style.fontWeight = "bold";
          label.style.color = "#444";
          label.style.marginBottom = "4px";
          label.innerText = "Center Domain:";
          domainSidebar.appendChild(label);

          domainData.forEach((dom) => {
            const btn = document.createElement("button");
            btn.innerText = `D${dom.id}: ${dom.rangeStr}`;
            btn.title = `Center on Domain ${dom.id} (${dom.rangeStr})`;
            btn.style.background = dom.color;
            btn.style.color = "#fff";
            btn.style.border = "none";
            btn.style.borderRadius = "4px";
            btn.style.padding = "6px 10px";
            btn.style.fontSize = "11px";
            btn.style.fontWeight = "bold";
            btn.style.textAlign = "left";
            btn.style.cursor = "pointer";
            btn.style.wordBreak = "break-all";

            btn.onclick = () => {
              const sele = dom.ranges.map(r => ({ resi: `${r.start}-${r.end}` }));
              viewer.zoomTo({ or: sele }, 500);
            };

            domainSidebar.appendChild(btn);
          });
        }

        function updateStructure() {
          const pdbText = model.get("pdb_text");
          const domainStr = model.get("domain_string");
          const normalRadius = model.get("normal_radius") ?? 0.20;
          const highlightRadius = model.get("highlight_radius") ?? 0.55;
          const segmentMask = model.get("segment_mask") || [];

          if (!pdbText) return;

          viewer.clear();
          viewer.addModel(pdbText, "pdb");

          // Default fallback base style (thin sticks)
          viewer.setStyle({}, { stick: { color: "#cccccc", radius: normalRadius } });

          const domainData = parseDomainString(domainStr);

          // Apply CA stick coloring and variable radii per residue
          domainData.forEach(dom => {
            dom.ranges.forEach(r => {
              for (let resi = r.start; resi <= r.end; resi++) {
                const resIdx = resi - 1; // Convert 1-based PDB residue ID to 0-based index
                const isSegment = Boolean(segmentMask[resIdx]);
                const radius = isSegment ? highlightRadius : normalRadius;

                viewer.setStyle(
                  { resi: `${resi}` },
                  { stick: { color: dom.color, radius: radius } }
                );
              }
            });
          });

          renderDomainBadges(domainData);

          // Active hover text label handler
          viewer.setHoverable(
            {},
            true,
            (atom, v, event, c) => {
              if (atom) {
                let domainName = "Unassigned";
                domainData.forEach(dom => {
                  dom.ranges.forEach(r => {
                    if (atom.resi >= r.start && atom.resi <= r.end) {
                      domainName = `Domain ${dom.id}`;
                    }
                  });
                });

                const rect = container.getBoundingClientRect();
                const x = event.clientX - rect.left + 12;
                const y = event.clientY - rect.top + 12;

                tooltip.innerText = `${atom.resn}${atom.resi} (${domainName})`;
                tooltip.style.left = `${x}px`;
                tooltip.style.top = `${y}px`;
                tooltip.style.display = "block";
              }
            },
            () => {
              tooltip.style.display = "none";
            }
          );

          // Click handler to store selected residue
          viewer.setClickable({}, true, (atom) => {
            if (atom) {
              model.set("selected_residue", atom.resi);
              model.save_changes();
            }
          });

          viewer.zoomTo();
          viewer.render();
        }

        // Initial render
        updateStructure();

        // Model Listeners
        model.on("change:pdb_text", updateStructure);
        model.on("change:domain_string", updateStructure);
        model.on("change:colors", updateStructure);
        model.on("change:normal_radius", updateStructure);
        model.on("change:highlight_radius", updateStructure);
        model.on("change:segment_mask", updateStructure);

        // Center / Reset view listener
        model.on("change:reset_trigger", () => {
          viewer.zoomTo({}, 500);
        });
      }
    };
    """

    pdb_text = traitlets.Unicode("").tag(sync=True)
    domain_string = traitlets.Unicode("").tag(sync=True)
    width = traitlets.Int(900).tag(sync=True)
    height = traitlets.Int(650).tag(sync=True)
    normal_radius = traitlets.Float(0.20).tag(sync=True)
    highlight_radius = traitlets.Float(0.55).tag(sync=True)
    segment_mask = traitlets.List().tag(sync=True)
    colors = traitlets.List(traitlets.Unicode()).tag(sync=True)
    selected_residue = traitlets.Int(None, allow_none=True).tag(sync=True)
    reset_trigger = traitlets.Int(0).tag(sync=True)


class DomNetViewer:
    """Persistent Jupyter DomNet viewer powered by pure anywidget."""

    def __init__(
        self,
        pdb_text="",
        domain_string="",
        width=900,
        height=450,
        colors=None,
        normal_radius=0.20,
        highlight_radius=0.55,
        segment_mask=None,
        **kwargs,
    ):
        self._default_domain_string = str(domain_string)
        self.colors = colors if colors is not None else DEFAULT_COLORS

        if segment_mask is None:
            segment_mask = []

        # --- Controls UI ---
        self._domain_widget = widgets.Text(
            value=str(domain_string),
            description="Domains:",
            layout=widgets.Layout(width="100%"),
            style={"description_width": "70px"},
        )

        self._apply_button = widgets.Button(
            description="Apply", layout=widgets.Layout(width="80px")
        )
        self._reset_button = widgets.Button(
            description="Reset", layout=widgets.Layout(width="80px")
        )

        self._status = widgets.HTML(
            value="<span style='color:#666'>Ready.</span>"
        )

        buttons = widgets.HBox(
            [self._apply_button, self._reset_button],
            layout=widgets.Layout(margin="4px 0 8px 70px"),
        )

        self._controls = widgets.VBox(
            [self._domain_widget, buttons, self._status]
        )

        # --- Canvas Widget ---
        self.canvas = DomNetCanvas(
            pdb_text=str(pdb_text),
            domain_string=str(domain_string),
            width=int(width),
            height=int(height),
            colors=list(self.colors),
            normal_radius=normal_radius,
            highlight_radius=highlight_radius,
            segment_mask=segment_mask,
        )

        # Wire Callbacks
        self._apply_button.on_click(self._on_apply)
        self._reset_button.on_click(self._on_reset)

        self.canvas.observe(self._on_js_selection, names=["selected_residue"])

    @property
    def domain_string(self):
        return self._domain_widget.value

    @domain_string.setter
    def domain_string(self, value):
        val = str(value)
        self._domain_widget.value = val
        self.canvas.domain_string = val

    @property
    def pdb_text(self):
        return self.canvas.pdb_text

    @pdb_text.setter
    def pdb_text(self, value):
        self.canvas.pdb_text = str(value)

    def display(self):
        """Display controls and canvas."""
        display(self._controls)
        display(self.canvas)

    def load(self, pdb_text, domain_string=None, reset_default=True):
        self.pdb_text = str(pdb_text)
        if domain_string is not None:
            self.domain_string = str(domain_string)
            if reset_default:
                self._default_domain_string = str(domain_string)
        self._set_status("Structure replaced.")

    def _on_apply(self, button):
        self.canvas.domain_string = self._domain_widget.value
        self._set_status("Domain ranges applied.")

    def _on_reset(self, button):
        val = self._default_domain_string
        self._domain_widget.value = val
        self.canvas.domain_string = val
        self.canvas.reset_trigger += 1
        self._set_status("Reset view and domain ranges.")

    def _on_js_selection(self, change):
        new_res = change["new"]
        if new_res is not None:
            self._set_status(f"Clicked Residue in JS: {new_res}")

    def _set_status(self, message):
        self._status.value = f"<span style='color:#666'>{message}</span>"
