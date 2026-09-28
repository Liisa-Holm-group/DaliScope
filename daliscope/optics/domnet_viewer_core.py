"""Shared DomNet HTML/Jupyter viewer core.

HTML/JavaScript templates use token replacement rather than Python formatting.
The JavaScript keeps one 3Dmol viewer/model during domain editing; only an
explicit structure replacement removes and adds a model.
"""

import html as html_module
import json
from pathlib import Path

DEFAULT_COLORS = [
    "#3366CC", "#DC3912", "#FF9900", "#109618",
    "#990099", "#0099C6", "#DD4477", "#66AA00",
]


def parse_domain_string(domain_string):
    """Parse ``1-402_725-761,403-724_762-823`` into domain dictionaries."""
    domains = []
    if domain_string is None:
        return domains
    text = str(domain_string).strip()
    if not text:
        return domains
    for domain_index, raw_domain in enumerate(text.split(",")):
        spans = []
        raw_domain = raw_domain.strip()
        if not raw_domain:
            continue
        for raw_span in raw_domain.split("_"):
            raw_span = raw_span.strip()
            if "-" not in raw_span:
                continue
            parts = raw_span.split("-", 1)
            try:
                start = int(parts[0].strip())
                end = int(parts[1].strip())
            except ValueError:
                continue
            if start < 1 or end < 1:
                continue
            if start > end:
                start, end = end, start
            spans.append((start, end))
        if spans:
            domains.append({"id": domain_index + 1, "spans": spans})
    return domains


def validate_domain_string(domain_string):
    domains = parse_domain_string(domain_string)
    if not domains:
        raise ValueError("No valid domains found in domain string.")
    return domains


def domains_to_json(domain_string, colors=None):
    if colors is None:
        colors = DEFAULT_COLORS
    result = []
    for domain in parse_domain_string(domain_string):
        result.append({
            "id": domain["id"],
            "color": colors[(domain["id"] - 1) % len(colors)],
            "spans": domain["spans"],
        })
    return json.dumps(result)


def build_static_legend_html(domain_string, colors=None):
    if colors is None:
        colors = DEFAULT_COLORS
    items = []
    for domain in parse_domain_string(domain_string):
        color = colors[(domain["id"] - 1) % len(colors)]
        ranges = "_".join(
            str(start) + "-" + str(end)
            for start, end in domain["spans"]
        )
        items.append(
            '<div class="legend-item">'
            '<span class="color-box" style="background:'
            + html_module.escape(color, quote=True)
            + ';">'
            '</span><span>Domain '
            + str(domain["id"])
            + ": "
            + html_module.escape(ranges, quote=True)
            + "</span></div>"
        )
    return "\n".join(items)


build_legend_html = build_static_legend_html


def build_viewer_javascript(viewer_id, pdb_text, domain_string,
                           title=None, colors=None, width=900, height=650):
    """Generate the shared browser-side implementation."""
    del title, width, height
    if colors is None:
        colors = DEFAULT_COLORS

    template = r'''(function() {
    "use strict";

    const VIEWER_ID = __VIEWER_ID__;
    const INITIAL_PDB = __PDB_TEXT__;
    const INITIAL_DOMAIN_STRING = __DOMAIN_STRING__;
    const COLORS = __COLORS__;

    const root = document.getElementById(VIEWER_ID);
    const canvas = document.getElementById(VIEWER_ID + "_canvas");
    const status = document.getElementById(VIEWER_ID + "_status");
    const legend = document.getElementById(VIEWER_ID + "_legend");
    const input = document.getElementById(VIEWER_ID + "_domain_input");
    const errorBox = document.getElementById(VIEWER_ID + "_domain_error");
    const applyButton = document.getElementById(VIEWER_ID + "_apply");
    const resetButton = document.getElementById(VIEWER_ID + "_reset");

    let viewer = null;
    let currentDomains = INITIAL_DOMAIN_STRING;
    let hoverInstalled = false;
    let tooltip = null;

    function load3Dmol() {
        if (window.$3Dmol) {
            return Promise.resolve();
        }
        return new Promise(function(resolve, reject) {
            const existing = document.getElementById("domnet-3dmol-script");
            if (existing) {
                existing.addEventListener("load", function() { resolve(); });
                existing.addEventListener("error", function() {
                    reject(new Error("Unable to load 3Dmol.js"));
                });
                return;
            }
            const script = document.createElement("script");
            script.id = "domnet-3dmol-script";
            script.src = "https://3dmol.org/build/3Dmol-min.js";
            script.onload = function() { resolve(); };
            script.onerror = function() {
                reject(new Error("Unable to load 3Dmol.js"));
            };
            document.head.appendChild(script);
        });
    }

    function parseDomains(text) {
        const domains = [];
        if (!text || !text.trim()) {
            return domains;
        }
        text.split(",").forEach(function(rawDomain, domainIndex) {
            rawDomain = rawDomain.trim();
            if (!rawDomain) {
                return;
            }
            const spans = [];
            rawDomain.split("_").forEach(function(rawSpan) {
                rawSpan = rawSpan.trim();
                if (!rawSpan || rawSpan.indexOf("-") < 0) {
                    return;
                }
                const parts = rawSpan.split("-", 2);
                const start = parseInt(parts[0], 10);
                const end = parseInt(parts[1], 10);
                if (!Number.isFinite(start) || !Number.isFinite(end)) {
                    return;
                }
                spans.push([Math.min(start, end), Math.max(start, end)]);
            });
            if (spans.length) {
                domains.push({id: domainIndex + 1, spans: spans});
            }
        });
        return domains;
    }

    function domainColor(domainId) {
        return COLORS[(domainId - 1) % COLORS.length];
    }

    function getDomainId(residueNumber, chain) {
        if (chain && chain !== "A") {
            return null;
        }
        const domains = parseDomains(currentDomains);
        for (let i = 0; i < domains.length; i++) {
            for (let j = 0; j < domains[i].spans.length; j++) {
                const span = domains[i].spans[j];
                if (residueNumber >= span[0] && residueNumber <= span[1]) {
                    return domains[i].id;
                }
            }
        }
        return null;
    }

    function installHover() {
        if (!viewer || hoverInstalled) {
            return;
        }
        hoverInstalled = true;
        tooltip = document.createElement("div");
        tooltip.id = VIEWER_ID + "_tooltip";
        tooltip.style.position = "absolute";
        tooltip.style.display = "none";
        tooltip.style.pointerEvents = "none";
        tooltip.style.zIndex = "1000";
        tooltip.style.background = "rgba(255,255,255,0.95)";
        tooltip.style.color = "#111";
        tooltip.style.border = "1px solid #888";
        tooltip.style.borderRadius = "4px";
        tooltip.style.padding = "4px 7px";
        tooltip.style.fontSize = "12px";
        tooltip.style.fontFamily = "Arial,sans-serif";
        canvas.style.position = "relative";
        canvas.appendChild(tooltip);

        viewer.setHoverable({}, true, function(atom, viewerInstance, event) {
            if (!atom) {
                tooltip.style.display = "none";
                return;
            }
            const residueNumber = parseInt(atom.resi, 10);
            if (!Number.isFinite(residueNumber)) {
                tooltip.style.display = "none";
                return;
            }
            const domainId = getDomainId(residueNumber, atom.chain || "");
            let text = "Residue " + residueNumber;
            if (atom.resn) {
                text += " (" + atom.resn + ")";
            }
            if (domainId !== null) {
                text += "<br>Domain " + domainId;
            }
            tooltip.innerHTML = text;
            if (event) {
                const rect = canvas.getBoundingClientRect();
                tooltip.style.left = (event.clientX - rect.left + 12) + "px";
                tooltip.style.top = (event.clientY - rect.top + 12) + "px";
            }
            tooltip.style.display = "block";
        });

        canvas.addEventListener("mouseleave", function() {
            tooltip.style.display = "none";
        });
    }

    function buildLegend(domains) {
        if (!legend) {
            return;
        }
        legend.innerHTML = "";
        const heading = document.createElement("div");
        heading.innerHTML = "<b>Domains</b>";
        heading.style.marginBottom = "5px";
        legend.appendChild(heading);
        if (!domains.length) {
            const empty = document.createElement("div");
            empty.textContent = "No domains defined.";
            legend.appendChild(empty);
            return;
        }
        domains.forEach(function(domain) {
            const row = document.createElement("div");
            row.style.display = "flex";
            row.style.alignItems = "center";
            row.style.gap = "7px";
            row.style.padding = "4px";
            row.style.cursor = "pointer";
            row.addEventListener("mouseenter", function() {
                row.style.background = "#f0f0f0";
            });
            row.addEventListener("mouseleave", function() {
                row.style.background = "transparent";
            });
            const swatch = document.createElement("span");
            swatch.style.display = "inline-block";
            swatch.style.width = "16px";
            swatch.style.height = "16px";
            swatch.style.borderRadius = "3px";
            swatch.style.background = domainColor(domain.id);
            const label = document.createElement("span");
            label.textContent = "Domain " + domain.id + ": " + domain.spans.map(
                function(span) { return span[0] + "-" + span[1]; }
            ).join("_");
            row.appendChild(swatch);
            row.appendChild(label);
            row.addEventListener("click", function() {
                centerOnDomain(domain);
            });
            legend.appendChild(row);
        });
    }

    function centerOnDomain(domain) {
        if (!viewer) {
            return;
        }
        const selection = {chain: "A", or: []};
        domain.spans.forEach(function(span) {
            selection.or.push({resi: span[0] + "-" + span[1]});
        });
        viewer.zoomTo(selection);
        viewer.render();
    }

    function applyDomainColors() {
        if (!viewer) {
            return;
        }
        const domains = parseDomains(currentDomains);
        viewer.setStyle(
            {chain: "A"},
            {cartoon: {color: "lightgray"}}
        );
        domains.forEach(function(domain) {
            const color = domainColor(domain.id);
            domain.spans.forEach(function(span) {
                viewer.setStyle(
                    {chain: "A", resi: span[0] + "-" + span[1]},
                    {cartoon: {color: color}}
                );
            });
        });
        buildLegend(domains);
        viewer.render();
    }

    function fitGlobalStructure() {
        if (!viewer) {
            return;
        }
        viewer.zoomTo();
        viewer.render();
    }

    function updateDomains(domainString) {
        const domains = parseDomains(domainString);
        if (!domains.length) {
            if (errorBox) {
                errorBox.textContent = "Invalid domain ranges: no valid ranges found.";
            }
            return false;
        }
        currentDomains = domainString;
        if (input) {
            input.value = domainString;
        }
        if (errorBox) {
            errorBox.textContent = "";
        }
        if (!viewer) {
            return false;
        }
        if (tooltip) {
            tooltip.style.display = "none";
        }
        applyDomainColors();
        return true;
    }

    function resetDomains() {
        currentDomains = INITIAL_DOMAIN_STRING;
        if (input) {
            input.value = INITIAL_DOMAIN_STRING;
        }
        if (errorBox) {
            errorBox.textContent = "";
        }
        if (tooltip) {
            tooltip.style.display = "none";
        }
        if (!viewer) {
            return;
        }
        applyDomainColors();
        fitGlobalStructure();
    }

    function replaceStructure(pdbText, domainString) {
        if (!viewer) {
            return false;
        }
        currentDomains = domainString || "";
        if (input) {
            input.value = currentDomains;
        }
        if (tooltip) {
            tooltip.style.display = "none";
        }
        viewer.removeAllModels();
        viewer.removeAllLabels();
        viewer.removeAllShapes();
        if (!pdbText || !pdbText.trim()) {
            if (status) {
                status.textContent = "No structure supplied.";
            }
            if (legend) {
                legend.innerHTML = "";
            }
            return false;
        }
        viewer.addModel(pdbText, "pdb");
        viewer.setStyle({}, {cartoon: {color: "lightgray"}});
        installHover();
        fitGlobalStructure();
        applyDomainColors();
        if (status) {
            status.textContent = "Structure loaded.";
        }
        return true;
    }

    window[VIEWER_ID + "_update_domains"] = function(domainString) {
        return updateDomains(domainString);
    };

    window[VIEWER_ID + "_replace_structure"] = function(pdbText, domainString) {
        return replaceStructure(pdbText, domainString);
    };

    if (applyButton) {
        applyButton.addEventListener("click", function() {
            if (input) {
                updateDomains(input.value.trim());
            }
        });
    }

    if (resetButton) {
        resetButton.addEventListener("click", function() {
            resetDomains();
        });
    }

    if (input) {
        input.addEventListener("keydown", function(event) {
            if (event.key === "Enter") {
                event.preventDefault();
                updateDomains(input.value.trim());
            }
        });
    }

    load3Dmol().then(function() {
        if (!canvas) {
            throw new Error("Viewer canvas element not found.");
        }
        viewer = $3Dmol.createViewer(canvas, {backgroundColor: "white"});
        replaceStructure(INITIAL_PDB, INITIAL_DOMAIN_STRING);
    }).catch(function(error) {
        console.error(error);
        if (status) {
            status.textContent = "Failed to load 3Dmol.js: " + String(error);
        }
    });
})();
'''

    return (
        template
        .replace("__VIEWER_ID__", json.dumps(str(viewer_id)))
        .replace("__PDB_TEXT__", json.dumps(str(pdb_text)))
        .replace("__DOMAIN_STRING__", json.dumps(str(domain_string)))
        .replace("__COLORS__", json.dumps(list(colors)))
    )


def build_viewer_controls_html(viewer_id, domain_string, include_header=True):
    """Build the common editable controls and viewer containers."""
    safe_id = html_module.escape(str(viewer_id), quote=True)
    safe_domains = html_module.escape(str(domain_string), quote=True)
    header = ""
    if include_header:
        header = (
            "    <h2>DomNet structure viewer</h2>\n"
            "    <div class=\"info-block\">"
            "Interactive DomNet domain decomposition</div>\n"
        )

    template = r'''<div id="__VIEWER_ID__" class="domnet-container">
__HEADER__
    <div class="domain-controls">
        <label for="__VIEWER_ID___domain_input">Domain ranges:</label>
        <input id="__VIEWER_ID___domain_input"
               type="text"
               value="__DOMAIN_STRING__"
               spellcheck="false">
        <button id="__VIEWER_ID___apply" type="button">Apply</button>
        <button id="__VIEWER_ID___reset" type="button">Reset</button>
    </div>
    <div id="__VIEWER_ID___domain_error" class="domain-error"></div>
    <div id="__VIEWER_ID___status" class="viewer-status">Loading...</div>
    <div id="__VIEWER_ID___canvas" class="viewer-canvas"></div>
    <div id="__VIEWER_ID___legend" class="viewer-legend"></div>
</div>
'''
    return (
        template
        .replace("__VIEWER_ID__", safe_id)
        .replace("__DOMAIN_STRING__", safe_domains)
        .replace("__HEADER__", header)
    )


def build_viewer_css(width=900, height=650):
    template = r'''.domnet-container {
    width: 100%;
    font-family: Arial, Helvetica, sans-serif;
}
.domnet-container h2 { margin: 0 0 12px 0; }
.info-block { margin: 5px 0 10px 0; font-size: 14px; }
.domain-controls {
    display: flex; align-items: center; gap: 8px; flex-wrap: wrap;
    margin: 8px 0 4px 0;
}
.domain-controls label { font-weight: bold; }
.domain-controls input {
    flex: 1 1 400px; min-width: 260px; padding: 7px 9px;
    border: 1px solid #aaa; border-radius: 4px;
    font-family: monospace; font-size: 14px;
}
.domain-controls button {
    padding: 7px 14px; border: 1px solid #999; border-radius: 4px;
    background: #fff; cursor: pointer; font-size: 14px;
}
.domain-controls button:hover { background: #eee; }
.domain-error { min-height: 18px; margin: 4px 0; color: #b00020;
    font-family: monospace; font-size: 13px; }
.viewer-status { min-height: 18px; margin: 4px 0; color: #666; font-size: 12px; }
.viewer-canvas {
    width: __WIDTH__px; height: __HEIGHT__px; position: relative;
    border: 1px solid #ccc; border-radius: 4px; overflow: hidden;
    background: white;
}
.viewer-legend {
    margin-top: 8px; padding: 8px; border: 1px solid #ddd;
    border-radius: 4px; font-size: 13px; background: white;
    max-width: __WIDTH__px;
}
@media (max-width: 700px) {
    .viewer-canvas { width: 100%; }
    .viewer-legend { max-width: 100%; }
}
'''
    return (
        template
        .replace("__WIDTH__", str(int(width)))
        .replace("__HEIGHT__", str(int(height)))
    )


def build_standalone_html(pdb_text, domain_string, title="DomNet Viewer",
                          width=900, height=650, colors=None):
    """Build a standalone HTML document with editing controls."""
    if colors is None:
        colors = DEFAULT_COLORS
    viewer_id = "domnet_standalone_viewer"
    controls = build_viewer_controls_html(viewer_id, domain_string, True)
    css = build_viewer_css(width, height)
    javascript = build_viewer_javascript(
        viewer_id, pdb_text, domain_string,
        title=title, colors=colors, width=width, height=height,
    )
    template = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>__TITLE__</title>
<style>
__CSS__
</style>
</head>
<body>
__CONTROLS__
<script>
__JAVASCRIPT__
</script>
</body>
</html>
'''
    return (
        template
        .replace("__TITLE__", html_module.escape(str(title), quote=True))
        .replace("__CSS__", css)
        .replace("__CONTROLS__", controls)
        .replace("__JAVASCRIPT__", javascript)
    )


def build_jupyter_html(viewer_id, pdb_text, domain_string,
                       width=900, height=650, colors=None):
    """Build the no-anywidget HTML/JavaScript fragment for Jupyter."""
    if colors is None:
        colors = DEFAULT_COLORS
    controls = build_viewer_controls_html(viewer_id, domain_string, False)
    css = build_viewer_css(width, height)
    javascript = build_viewer_javascript(
        viewer_id, pdb_text, domain_string,
        colors=colors, width=width, height=height,
    )
    template = r'''<style>
__CSS__
</style>
__CONTROLS__
<script>
__JAVASCRIPT__
</script>
'''
    return (
        template
        .replace("__CSS__", css)
        .replace("__CONTROLS__", controls)
        .replace("__JAVASCRIPT__", javascript)
    )


def js_update_domains(viewer_id, domain_string):
    """Generate a browser command that updates domains without model reload."""
    template = r'''(function() {
    const fn = window[__VIEWER_ID__ + "_update_domains"];
    if (typeof fn === "function") {
        fn(__DOMAIN_STRING__);
    }
})();
'''
    return (
        template
        .replace("__VIEWER_ID__", json.dumps(str(viewer_id)))
        .replace("__DOMAIN_STRING__", json.dumps(str(domain_string)))
    )


def js_replace_structure(viewer_id, pdb_text, domain_string):
    """Generate a browser command that replaces the current PDB model."""
    template = r'''(function() {
    const fn = window[__VIEWER_ID__ + "_replace_structure"];
    if (typeof fn === "function") {
        fn(__PDB_TEXT__, __DOMAIN_STRING__);
    }
})();
'''
    return (
        template
        .replace("__VIEWER_ID__", json.dumps(str(viewer_id)))
        .replace("__PDB_TEXT__", json.dumps(str(pdb_text)))
        .replace("__DOMAIN_STRING__", json.dumps(str(domain_string)))
    )


def save_standalone_html(pdb_text, domain_string, output_filename,
                         title="DomNet Viewer", width=900, height=650,
                         colors=None):
    """Generate and save a standalone viewer."""
    text = build_standalone_html(
        pdb_text=pdb_text,
        domain_string=domain_string,
        title=title,
        width=width,
        height=height,
        colors=colors,
    )
    output_path = Path(output_filename)
    output_path.write_text(text, encoding="utf-8")
    return output_path


def build_html(pdb_text, domain_string, title="DomNet Viewer",
               width=900, height=650, colors=None):
    """Compatibility alias for build_standalone_html."""
    return build_standalone_html(
        pdb_text=pdb_text,
        domain_string=domain_string,
        title=title,
        width=width,
        height=height,
        colors=colors,
    )


__all__ = [
    "DEFAULT_COLORS",
    "parse_domain_string",
    "validate_domain_string",
    "domains_to_json",
    "build_static_legend_html",
    "build_legend_html",
    "build_viewer_javascript",
    "build_viewer_controls_html",
    "build_viewer_css",
    "build_standalone_html",
    "build_jupyter_html",
    "js_update_domains",
    "js_replace_structure",
    "save_standalone_html",
    "build_html",
]
