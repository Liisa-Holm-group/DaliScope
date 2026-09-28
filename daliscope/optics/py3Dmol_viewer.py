import py3Dmol
from IPython.display import display, clear_output, HTML
import ipywidgets as widgets
import math
import numpy as np
from typing import List, Tuple, Dict

############

from collections import defaultdict


def parse_pdb(pdb_str):
    """
    Parse a PDB string into structures suitable for residue-level rendering.

    Returns
    -------
    residues : list of dict
        One dictionary per residue (taken from the CA atom):
            resi      residue number (int)
            resn      residue name (str)
            chain     chain ID
            x,y,z     CA coordinates
            b         B-factor
            ca_line   original ATOM line

    residue_atoms : dict
        (chain,resi) -> list of original ATOM/HETATM lines

    residue_lookup : dict
        (chain,resi) -> residue dictionary
    """

    residues = []
    residue_atoms = defaultdict(list)
    residue_lookup = {}

    for line in pdb_str.splitlines():

        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        atom = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21]
        resi = int(line[22:26])

        residue_atoms[(chain, resi)].append(line)

        if atom != "CA":
            continue

        try:
            b = float(line[60:66])
        except ValueError:
            b = 0.0

        r = {
            "key": (chain, resi),
            "resi": resi,
            "resn": resn,
            "chain": chain,
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "b": b,
            "ca_line": line,
        }

        residues.append(r)
        residue_lookup[(chain, resi)] = r

    return residues, residue_atoms, residue_lookup


def make_highlight_pdb(residues,
                     residue_atoms,
                     important_residues):
    """
    Build a minimal PDB containing

      • every atom of highlighted residues

    Returns
    -------
    pdb_text
    residue_lookup
    """

    # ---------------------------------------------------------
    # Normalize highlighted residues
    # ---------------------------------------------------------

    important = set()

    for r in important_residues:

        if isinstance(r, tuple):
            important.add(r)

        else:
            resi = int(r)

            for res in residues:
                if res["resi"] == resi:
                    important.add(res["key"])

    residue_lookup = {r["key"]: r for r in residues}

    out = []
    conect = []

    current_chain = None
    previous_ca_serial = None
    serial = 1

    # residue key -> new CA atom serial number
    ca_serial = {}

    # ---------------------------------------------------------
    # Write ATOM records
    # ---------------------------------------------------------

    for r in residues:

        chain, resi = r["key"]

        if current_chain is not None and chain != current_chain:
            out.append("TER")
            previous_ca_serial = None

        current_chain = chain

        if r["key"] in important:
            atom_lines = residue_atoms[r["key"]]
        else:
            continue

        for line in atom_lines:

            newline = f"{line[:6]}{serial:5d}{line[11:]}"
            out.append(newline)

            serial += 1

    out.append("TER")

    out.append("END")

    return "\n".join(out) + "\n", residue_lookup

# ----------------------------------------------------------------------

def property_to_color(property_dict,
                       nbins=10,
                       vmin=None,
                       vmax=None,
                       missing="#bdbdbd"):
    """
    Construct a residue -> colour mapping.
    Parameters
    ----------
    property_dict : dict
        Either
            {137:0.81, 138:0.42, ...}
        or
            {('A',137):0.81, ('B',52):0.42, ...}
    nbins : int
        Number of colour bins (default 10).
    vmin,vmax : float or None
        Colour scale limits.
        If None, computed from the supplied values.
    missing : str
        Colour for missing values.
    Returns
    -------
    colour_lookup : dict
        Same keys as property_dict, values are "#RRGGBB".
    colour_function : callable
        colour_function(key) returns colour for any residue.
    """
    values = [v for v in property_dict.values() if v is not None and not math.isnan(v)]
    if len(values) == 0:
        def colour_function(key):
            return missing
        return {}, colour_function

    if vmin is None:
        vmin = min(values)
    if vmax is None:
        vmax = max(values)
    if vmax <= vmin:
        vmax = vmin + 1.0

    # 10-colour blue-white-red palette
    palette = [
        "#084594",  # dark blue
        "#2171b5",
        "#4292c6",
        "#41ab5d",  # green
        "#78c679",
        "#c2e699",  # yellow-green
        "#ffff33",  # bright yellow
        "#fdae61",  # orange
        "#f16913",  # orange-red
        "#cb181d",  # deep red
    ]
    if nbins != 10:
        raise NotImplementedError("Current implementation supplies a fixed 10-bin palette.")

    colour_lookup = {}
    scale = (nbins - 1) / (vmax - vmin)
    for key, value in property_dict.items():
        if value is None or math.isnan(value):
            colour_lookup[key] = missing
            continue
        k = int((value - vmin) * scale)
        if k < 0:
            k = 0
        elif k >= nbins:
            k = nbins - 1
        colour_lookup[key] = palette[k]

    def colour_function(key):
        if key in colour_lookup:
            return colour_lookup[key]
        # if key is (chain,resi), also try resi only
        if isinstance(key, tuple) and key[1] in colour_lookup:
            return colour_lookup[key[1]]
        return missing

    return colour_lookup, colour_function

def show_structure(pdb_str,
                    important_residues,
                    property_dict=(),
                    width=800,
                    height=600,
                    backbone_radius=0.18,
                    ca_radius=0.28,
                    highlight_scale=1.8):
    """
    Hybrid residue viewer.

    Backbone
        neutral cylinders built directly from CA coordinates (no bond/stick
        geometry needed for the trace, which keeps the embedded structure small)
        coloured CA spheres

    Highlighted residues
        full side chains
        coloured by residue property
        labelled

    property_dict may use keys:
        137
    or
        ('A', 137)
    """

    # ---- Parse structure ----
    residues, residue_atoms, residue_lookup = parse_pdb(pdb_str)
    highlight_pdb, residue_lookup = make_highlight_pdb(residues, residue_atoms, important_residues)

    # ---- Colour mapper ----
    _, colour = property_to_color(property_dict)

    # ---- Viewer ----
    view = py3Dmol.view(width=width, height=height)

    # ---- Backbone: cylinders between consecutive CA atoms, colored by conservation ----
    draw_ca_backbone(view, residues, colour)

    # ---- CA spheres, enlarged for highlighted/important residues ----
    important = set()
    for r in important_residues:
        if isinstance(r, tuple):
            important.add(r)
        else:
            resi = int(r)
            for x in residues:
                if x["resi"] == resi:
                    important.add(x["key"])

    for r in residues:
        radius = ca_radius * highlight_scale if r["key"] in important else ca_radius
        view.addSphere({
            "center": {"x": r["x"], "y": r["y"], "z": r["z"]},
            "radius": radius,
            "color": colour(r["key"]),
        })

    # ---- Minimal PDB model containing only CA atoms + highlighted residues' side chains ----
    view.addModel(highlight_pdb, "pdb")

    # Color the highlighted side-chain sticks (CA spheres above already carry residue color)
    for key in important:
        chain, resi = key
        view.setStyle(
            {"chain": chain, "resi": resi},
            {"stick": {"color": colour(key), "radius": backbone_radius * 3.0}},
        )

    # ---- Labels for highlighted residues ----
    for key in important:
        r = residue_lookup[key]
        view.addLabel(
            f"{r['resi']} {r['resn']}",
            {
                "position": {"x": r["x"] + 0.8, "y": r["y"] + 0.8, "z": r["z"] + 0.8},
                "fontSize": 12,
                "fontColor": "black",
                "backgroundColor": "white",
                "backgroundOpacity": 0.85,
                "borderColor": "#888888",
                "borderThickness": 1,
            },
        )

    view.zoomTo()
    view.setBackgroundColor("#f8f9fa")
    view.show()



############

def extract_ca_records(pdb_str: str) -> str:
    """
    Utility function for cartoon views: 
    Extracts only CA (alpha-carbon), C and N ATOM records from a full PDB string.
    Preserves non-ATOM/HETATM structural records (HEADER, TER, END, etc.)
    unchanged; drops all non-[CA,C,N] atom lines.

    Parameters
    ----------
    pdb_str : full PDB-formatted structure, newline-separated.

    Returns
    -------
    A new PDB-formatted string containing only CA atoms plus original
    non-atom records.
    """
    kept_lines = []
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16]
            if atom_name in [" CA "]:
                kept_lines.append(line)
            # else: drop this atom line (not a CA)
        else:
            # Keep structural records (HEADER, TER, END, etc.) untouched
            kept_lines.append(line)

    return "\n".join(kept_lines)

def draw_ca_backbone(view, residues, colour, backbone_radius=0.28):
    """
    Draws the CA trace as cylinders directly from CA coordinates
    (no addModel() — keeps the embedded structure minimal).

    residues : list of residue dicts (from parse_pdb), each with
               x, y, z, chain, key.
    colour   : callable, colour(key) -> "#RRGGBB" hex string.
    """
    for a, b in zip(residues[:-1], residues[1:]):
        if a["chain"] != b["chain"]:
            continue

        mx = 0.5 * (a["x"] + b["x"])
        my = 0.5 * (a["y"] + b["y"])
        mz = 0.5 * (a["z"] + b["z"])

        view.addCylinder({
            "start": {"x": a["x"], "y": a["y"], "z": a["z"]},
            "end":   {"x": mx,     "y": my,     "z": mz},
            "radius": backbone_radius,
            "color": colour(a["key"]),
        })
        view.addCylinder({
            "start": {"x": mx,     "y": my,     "z": mz},
            "end":   {"x": b["x"], "y": b["y"], "z": b["z"]},
            "radius": backbone_radius,
            "color": colour(b["key"]),
        })

def _make_domain_preview_str(pdb_str, domain_ranges=None):
    NEUTRAL_COLOR = "#cccccc"
    HIGHLIGHT_COLOR = "#e05c5c"

    ca_pdb_str = extract_ca_records(pdb_str)
    residues, residue_atoms, residue_lookup = parse_pdb(ca_pdb_str)

    highlighted_resi = set()
    if domain_ranges is not None:
        for start, end in domain_ranges:
            highlighted_resi.update(range(start + 1, end + 1))

    def colour(key):
        chain, resi = key
        return HIGHLIGHT_COLOR if resi in highlighted_resi else NEUTRAL_COLOR

    view = py3Dmol.view(width=450, height=300)
    draw_ca_backbone(view, residues, colour)
    view.zoomTo()
    view.setBackgroundColor("#f8f9fa")
    view.show()

def make_domain_panel(pdb_str, domain_ranges, view_name, width=250, height=250):
    NEUTRAL_COLOR = "#cccccc"
    HIGHLIGHT_COLOR = "#e05c5c"

    # 1. Parse C-alpha residues & coordinates
    ca_pdb_str = extract_ca_records(pdb_str)
    residues, residue_atoms, residue_lookup = parse_pdb(ca_pdb_str)

    # 2. Build set of 1-indexed residue numbers to highlight
    highlighted_resi = set()
    if domain_ranges is not None:
        for start, end in domain_ranges:
            highlighted_resi.update(range(start + 1, end + 1))

    # 3. Define color assignment callback
    def colour(key):
        chain, resi = key
        return HIGHLIGHT_COLOR if resi in highlighted_resi else NEUTRAL_COLOR

    out = widgets.Output(
        layout=widgets.Layout(
            width=f"{width}px",
            height=f"{height}px"
        )
    )

    # 4. Draw backbone using light coordinate primitives (spheres & cylinders)
    with out:
        view = py3Dmol.view(width=width, height=height)
        draw_ca_backbone(view, residues, colour)
        view.zoomTo()
        view.setBackgroundColor("#f8f9fa")
        view.show()

    # 5. Caption Widget
    caption = widgets.HTML(
        f"""
        <div style="
            text-align:center;
            font-size:12px;
            font-family:monospace;
            margin-top:4px;
        ">
            {view_name}
        </div>
        """
    )

    return widgets.VBox([out, caption])


def show_domain_gallery(
    items,
    pdb_str,
    ncols=4,
    width=250,
    height=250
):
    # strip to CA
    pdb_str = extract_ca_records(pdb_str)

    panels = [
        make_domain_panel(
            pdb_str,
            domain_ranges,
            view_name,
            width,
            height
        )
        for domain_ranges, view_name in items if domain_ranges is not None
    ]

    rows = [
        widgets.HBox(
            panels[i:i+ncols],
            layout=widgets.Layout(gap="10px")
        )
        for i in range(0, len(panels), ncols)
    ]

    gallery = widgets.VBox(
        rows,
        layout=widgets.Layout(gap="10px")
    )

    display(gallery)

    return

###########

def make_residue_selector_widget(df, query_pdb_str, property_dict, top_ic=5):
    """Creates an interactive checkbox widget with monospace column headers and a 2-column layout.

    Fixes 0-based vs 1-based indexing offsets between DaliScope and PDB
    structures.
    """
    checkboxes = {}
    checkbox_widgets = []

    # CSS to force Courier/Monospace font across all headers and checkbox labels
    font_style = widgets.HTML("""
    <style>
        .mono-header {
            font-family: 'Courier New', Courier, monospace !important;
            font-size: 12px !important;
            font-weight: bold !important;
            color: #475569 !important;
            white-space: pre !important;
            margin-bottom: 4px;
        }
        .mono-checkbox label {
            font-family: 'Courier New', Courier, monospace !important;
            font-size: 12px !important;
            white-space: pre !important;
        }
    </style>
    """)

    # 1. Dynamic Monospace Column Header String
    index_name = df.index.name if df.index.name else "pdb_pos"
    col_names = df.columns.tolist()

    header_str = f"{index_name:<5s} | {col_names[0]:<7s} ({col_names[1]:<7s}) | {col_names[2]}"

    header_col1 = widgets.HTML(
        f"<div class='mono-header'>&nbsp;&nbsp;&nbsp;&nbsp;{header_str}</div>"
    )
    header_col2 = widgets.HTML(
        f"<div class='mono-header'>&nbsp;&nbsp;&nbsp;&nbsp;{header_str}</div>"
    )

    # 2. Build Monospace Checkboxes with 1-based PDB Position
    i = 0
    for pos, row in df.iterrows():
        pdb_pos = pos + 0  # Map 0-based DaliScope index to 1-based PDB position

        label = f"{pdb_pos:<5d} | {str(row[col_names[0]]):<7s} ({row[col_names[1]]:7.2f}) | {row[col_names[2]]}"

        cb = widgets.Checkbox(
            value= i < top_ic,
            description=label,
            indent=False,
            layout=widgets.Layout(width="100%"),
        )
        cb.add_class("mono-checkbox")
        i += 1

        checkboxes[pos] = cb
        checkbox_widgets.append(cb)

    # 3. Arrange in 2 Columns with Headers
    half = (len(checkbox_widgets) + 1) // 2

    col1 = widgets.VBox(
        [header_col1] + checkbox_widgets[:half],
        layout=widgets.Layout(width="50%"),
    )
    col2 = widgets.VBox(
        [header_col2] + checkbox_widgets[half:],
        layout=widgets.Layout(width="50%"),
    )
    checkbox_grid = widgets.HBox(
        [col1, col2], layout=widgets.Layout(width="100%", margin="5px 0px")
    )

    # 4. Control Buttons & Outputs
    select_all_btn = widgets.Button(
        description="Select All", layout=widgets.Layout(width="110px")
    )
    deselect_all_btn = widgets.Button(
        description="Deselect All", layout=widgets.Layout(width="110px")
    )
    preview_btn = widgets.Button(
        description="Preview Structure",
        button_style="primary",
        layout=widgets.Layout(width="150px"),
    )

    out_3d = widgets.Output()
    out_status = widgets.Output()

    # 5. Callbacks
    def get_selected_positions():
        """Returns list of checked index positions (0-based DaliScope int)."""
        return [pos for pos, cb in checkboxes.items() if cb.value]

    def on_select_all(b):
        for cb in checkboxes.values():
            cb.value = True
        update_status()

    def on_deselect_all(b):
        for cb in checkboxes.values():
            cb.value = False
        update_status()

    def update_status():
        with out_status:
            clear_output(wait=True)
            sel_0based = get_selected_positions()
            sel_1based = [p + 1 for p in sel_0based]
            print(
                f"Selected 0-based (DaliScope): {sel_0based}\nSelected 1-based (PDB): {sel_1based}"
            )

    def on_preview_click(b):
        # Convert 0-based DaliScope positions to 1-based PDB residue positions
        plot_residues_pdb = [p + 1 for p in get_selected_positions()]
        update_status()

        with out_3d:
            clear_output(wait=True)
            if not plot_residues_pdb:
                print(
                    "No residues selected. Check at least one position to preview."
                )
                return

            show_structure(
                query_pdb_str,
                plot_residues_pdb,  # Pass 1-based list to PDB visualizer
                property_dict=property_dict,
                backbone_radius=0.28,
                ca_radius=0.28,
                highlight_scale=1.0,
            )

    # Attach Events
    select_all_btn.on_click(on_select_all)
    deselect_all_btn.on_click(on_deselect_all)
    preview_btn.on_click(on_preview_click)

    for cb in checkboxes.values():
        cb.observe(lambda change: update_status(), names="value")

    # 6. Assemble Layout
    controls_box = widgets.HBox(
        [select_all_btn, deselect_all_btn, preview_btn]
    )

    widget_container = widgets.VBox([
        font_style,
        widgets.HTML("<b>Select Residue Positions:</b>"),
        checkbox_grid,
        controls_box,
        out_status,
        out_3d,
    ])

    update_status()
    return widget_container, get_selected_positions


def old_make_residue_selector_widget(df, query_pdb_str, property_dict):
    """
    Creates an interactive checkbox widget with monospace column headers and a 2-column layout.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with index as residue positions (pos) and columns 
        ['max_column', 'max_value', 'sequence_context'].
    query_pdb_str : str
        PDB string passed to show_structure().
    property_dict : dict
        Property values dictionary passed to show_structure().
        
    Returns:
    --------
    widget_container : ipywidgets.VBox
        The displayed UI widget.
    get_selected_positions : function
        Callable function to retrieve the list of currently checked positions.
    """
    checkboxes = {}
    checkbox_widgets = []

    # CSS to force Courier/Monospace font across all headers and checkbox labels
    font_style = widgets.HTML("""
    <style>
        .mono-header {
            font-family: 'Courier New', Courier, monospace !important;
            font-size: 12px !important;
            font-weight: bold !important;
            color: #475569 !important;
            white-space: pre !important;
            margin-bottom: 4px;
        }
        .mono-checkbox label {
            font-family: 'Courier New', Courier, monospace !important;
            font-size: 12px !important;
            white-space: pre !important;
        }
    </style>
    """)

    # 1. Dynamic Monospace Column Header String
    index_name = df.index.name if df.index.name else "pos"
    col_names = df.columns.tolist()

    # Matches character spacing: "pos   | max_col (max_val) | sequence_context"
    header_str = f"{index_name:<5s} | {col_names[0]:<7s} ({col_names[1]:<7s}) | {col_names[2]}"

    header_col1 = widgets.HTML(f"<div class='mono-header'>&nbsp;&nbsp;&nbsp;&nbsp;{header_str}</div>")
    header_col2 = widgets.HTML(f"<div class='mono-header'>&nbsp;&nbsp;&nbsp;&nbsp;{header_str}</div>")

    # 2. Build Monospace Checkboxes
    for pos, row in df.iterrows():
        # Monospace formatted row matching header character width
        label = f"{pos:<5d} | {str(row[col_names[0]]):<7s} ({row[col_names[1]]:7.2f}) | {row[col_names[2]]}"

        cb = widgets.Checkbox(
            value=True,  # Default to checked
            description=label,
            indent=False,
            layout=widgets.Layout(width='100%')
        )
        cb.add_class('mono-checkbox')

        checkboxes[pos] = cb
        checkbox_widgets.append(cb)

    # 3. Arrange in 2 Columns with Headers
    half = (len(checkbox_widgets) + 1) // 2

    col1 = widgets.VBox([header_col1] + checkbox_widgets[:half], layout=widgets.Layout(width='50%'))
    col2 = widgets.VBox([header_col2] + checkbox_widgets[half:], layout=widgets.Layout(width='50%'))
    checkbox_grid = widgets.HBox([col1, col2], layout=widgets.Layout(width='100%', margin='5px 0px'))

    # 4. Control Buttons & Outputs
    select_all_btn = widgets.Button(description="Select All", layout=widgets.Layout(width='110px'))
    deselect_all_btn = widgets.Button(description="Deselect All", layout=widgets.Layout(width='110px'))
    preview_btn = widgets.Button(description="Preview Structure", button_style='primary', layout=widgets.Layout(width='150px'))

    out_3d = widgets.Output()
    out_status = widgets.Output()

    # 5. Callbacks
    def get_selected_positions():
        """Returns list of checked index positions (int)."""
        return [pos for pos, cb in checkboxes.items() if cb.value]

    def on_select_all(b):
        for cb in checkboxes.values():
            cb.value = True
        update_status()

    def on_deselect_all(b):
        for cb in checkboxes.values():
            cb.value = False
        update_status()

    def update_status():
        with out_status:
            clear_output(wait=True)
            sel = get_selected_positions()
            print(f"Selected Positions ({len(sel)}): {sel}")

    def on_preview_click(b):
        plot_residues = get_selected_positions()
        update_status()

        with out_3d:
            clear_output(wait=True)
            if not plot_residues:
                print("No residues selected. Check at least one position to preview.")
                return

            show_structure(
                query_pdb_str,
                plot_residues,
                property_dict=property_dict,
                backbone_radius=0.28,
                ca_radius=0.28,
                highlight_scale=1.0
            )

    # Attach Events
    select_all_btn.on_click(on_select_all)
    deselect_all_btn.on_click(on_deselect_all)
    preview_btn.on_click(on_preview_click)

    for cb in checkboxes.values():
        cb.observe(lambda change: update_status(), names='value')

    # 6. Assemble Layout
    controls_box = widgets.HBox([select_all_btn, deselect_all_btn, preview_btn])

    widget_container = widgets.VBox([
        font_style,
        widgets.HTML("<b>Select Residue Positions:</b>"),
        checkbox_grid,
        controls_box,
        out_status,
        out_3d
    ])

    update_status()
    return widget_container, get_selected_positions

#########

# DomNet VISUALIZATION

from typing import List, Union
import numpy as np

# Standard mapping from 1-letter amino acid codes to 3-letter PDB codes
AA_1_TO_3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    "U": "SEC", "O": "PYL", "X": "UNK",
}


def xyz_to_ca_pdb(
    xyz: np.ndarray,
    seq: Union[str, np.ndarray, List[str], None] = None,
    res_names: Union[List[str], np.ndarray, None] = None,
) -> str:
    """Converts an (N, 3) CA-only numpy array into a PDB-formatted string including

    CONECT records for backbone connectivity.

    Args:
        xyz: (N, 3) array of C-alpha coordinates.
        seq: Optional sequence as a concatenated string (e.g., "ACDEF..."),
          a numpy array of 1-letter characters, or a list of 1-letter strings.
        res_names: Optional explicit list/array of 3-letter residue names. Takes
          precedence over `seq` if provided.

    Returns:
        PDB-formatted string with CA ATOM records and CONECT bonds.
    """
    n_residues = len(xyz)

    # 1. Resolve 3-letter residue names
    if res_names is not None:
        formatted_res_names = [str(r).upper() for r in res_names]
    elif seq is not None:
        # Convert string, numpy array, or list to 1-letter characters
        if isinstance(seq, (str, np.str_)):
            seq_chars = list(str(seq))
        else:
            seq_chars = [str(aa) for aa in seq]

        # Convert 1-letter codes to 3-letter PDB names
        formatted_res_names = [
            AA_1_TO_3.get(aa.upper(), "UNK") for aa in seq_chars
        ]
    else:
        formatted_res_names = ["UNK"] * n_residues

    lines = []

    # 2. ATOM records
    for i in range(n_residues):
        atom_id = i + 1
        res_id = i + 1
        x, y, z = xyz[i]
        res_name = (
            formatted_res_names[i] if i < len(formatted_res_names) else "UNK"
        )

        # Standard PDB ATOM record format
        line = (
            f"ATOM  {atom_id:5d}  CA  {res_name:>3s} A{res_id:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
        lines.append(line)

    # 3. CONECT records (bonds adjacent CA atoms so sticks form a continuous trace)
    for i in range(1, n_residues):
        lines.append(f"CONECT{i:5d}{(i + 1):5d}")

    lines.append("END")
    return "\n".join(lines)

def visualize_communities_ca_sticks(
    xyz: np.ndarray,
    assignments: np.ndarray,
    segments: List[Tuple[str, int, int]],
    normal_radius: float = 0.20,
    highlight_radius: float = 0.55,
    width: int = 800,
    height: int = 600
):
    """
    Renders CA-only coordinates as continuous sticks in py3Dmol colored by community ID,
    highlighting defined helix/beta segments with a larger stick radius.

    Args:
        xyz: (N, 3) array of C-alpha coordinates.
        assignments: (N,) array of community assignments (-1 for unassigned).
        segments: List of (ch, start_res, end_res) tuples for helices and strands.
        normal_radius: Stick radius for non-segment/loop residues.
        highlight_radius: Larger stick radius for helix/beta segment residues.
    """
    n_residues = len(xyz)

    # 1. Build a fast lookup mask for segment membership
    is_segment_res = np.zeros(n_residues, dtype=bool)
    for ch, start_res, end_res in segments:
        is_segment_res[start_res : end_res + 1] = True

    # 2. Generate PDB string from CA coordinates
    pdb_string = xyz_to_ca_pdb(xyz)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, "pdb")

    palette = [
        "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", 
        "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#008080"
    ]

    # Base fallback style
    view.setStyle({'color': 'lightgray', 'stick': {'radius': normal_radius}})

    # 3. Apply community colors and variable stick radius per residue
    for res_idx in range(n_residues):
        comm_id = assignments[res_idx]

        # Determine radius based on segment membership
        radius = highlight_radius if is_segment_res[res_idx] else normal_radius

        # Determine color based on community assignment
        color = palette[comm_id % len(palette)] if comm_id != -1 else "lightgray"

        # Apply style (PDB residue index is 1-based)
        view.setStyle(
            {'resi': str(res_idx + 1)},
            {'stick': {'color': color, 'radius': radius}}
        )

    view.zoomTo()
    return view.show()
