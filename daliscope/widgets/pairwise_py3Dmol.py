import numpy as np
import pandas as pd
import py3Dmol
import ipywidgets as widgets
from IPython.display import display, clear_output, HTML


import numpy as np
import pandas as pd
import py3Dmol


# ==============================================================================
# 1. CORE GEOMETRY & TRANSFORMATION HELPERS
# ==============================================================================

def apply_transform(coords: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Applies rigid body transformation (X @ R.T) + t to row-vector coordinates."""
    R = np.asarray(R)
    t = np.asarray(t).reshape(1, 3)
    return (coords @ R.T) + t

def draw_variable_width_backbone(
    view,
    coords,
    is_aligned_flags,
    color,
    aligned_radius=0.5,
    unaligned_radius=0.12,
):
    # First: thin complete backbone
    for i in range(len(coords) - 1):
        p1 = {
            "x": float(coords[i][0]),
            "y": float(coords[i][1]),
            "z": float(coords[i][2]),
        }
        p2 = {
            "x": float(coords[i + 1][0]),
            "y": float(coords[i + 1][1]),
            "z": float(coords[i + 1][2]),
        }

        view.addCylinder({
            "start": p1,
            "end": p2,
            "radius": unaligned_radius,
            "color": color,
            "fromCap": 1,
            "toCap": 1,
        })

    # Second: thick overlay for aligned segments
    for i in range(len(coords) - 1):
        if not (is_aligned_flags[i] and is_aligned_flags[i + 1]):
            continue

        p1 = {
            "x": float(coords[i][0]),
            "y": float(coords[i][1]),
            "z": float(coords[i][2]),
        }
        p2 = {
            "x": float(coords[i + 1][0]),
            "y": float(coords[i + 1][1]),
            "z": float(coords[i + 1][2]),
        }

        view.addCylinder({
            "start": p1,
            "end": p2,
            "radius": aligned_radius,
            "color": color,
            "fromCap": 1,
            "toCap": 1,
        })

def draw_alignment_connectors_cylinders(
    view: py3Dmol.view,
    q_coords: np.ndarray,
    t_coords: np.ndarray,
    q2t_map: list,
    color: str = "#616161",
    radius: float = 0.48
):
    """Draws dashed 3D cylinder connectors between matched C-alpha pairs."""
    for query_idx, target_idx in enumerate(q2t_map):
        if not isinstance(target_idx, int) or target_idx < 0:
            continue

        if query_idx < len(q_coords) and target_idx < len(t_coords):
            q_pt = q_coords[query_idx]
            t_pt = t_coords[target_idx]

            view.addCylinder({
                "start": {"x": float(q_pt[0]), "y": float(q_pt[1]), "z": float(q_pt[2])},
                "end":   {"x": float(t_pt[0]), "y": float(t_pt[1]), "z": float(t_pt[2])},
                "radius": radius,
                "color": color,
                "dashed": True,
                "fromCap": 1,
                "toCap": 1
            })


# ==============================================================================
# 2. ALIGNMENT MAP BUILDER
# ==============================================================================

def build_q2t_map_from_segments(
    target_ix: int,
    df: pd.DataFrame,
    project,
    query_len: int,
) -> list:
    """
    Build q2t_map using:
      q_start = 0-based index in q_ca_coords
      s_start = 0-based index in the FULL target protein

    Returned map initially contains full-target indices.
    """

    q2t_map = [-1] * query_len

    target_align_id = df.iloc[target_ix]["alignment_id"]

    if not isinstance(project.segments, pd.DataFrame):
        raise TypeError("Expected project.segments to be a pandas DataFrame")

    matched_segments = project.segments[
        project.segments["alignment_id"] == target_align_id
    ]

    for _, seg in matched_segments.iterrows():
        q_start = int(seg["q_start"])
        s_start = int(seg["s_start"])
        length = int(seg["length"])

        for offset in range(length):
            q_idx = q_start + offset
            t_idx = s_start + offset

            if 0 <= q_idx < query_len:
                q2t_map[q_idx] = t_idx

    return q2t_map

# ==============================================================================
# 3. VISUALIZATION BUILDER (THICK/THIN RADIUS DISTINCTION)
# ==============================================================================

def visualize_pair_alignment_variable_width(
    q_ca_coords: np.ndarray,
    target_ca_coords: np.ndarray,
    q2t_map: list,
    R: np.ndarray,
    t: np.ndarray,
    start_idx: int,
    end_idx: int,
    title: str = "",
    aligned_radius: float = 0.5,
    unaligned_radius: float = 0.10,
    connector_radius: float = 0.38,
    q_color: str = "#1f77b4",  # Query Blue
    t_color: str = "#d62728"   # Target Red
) -> py3Dmol.view:
    """Slices, transforms target coordinates, and renders alignment using single color

    per chain with variable thickness for aligned vs. unaligned regions.
    """
    
    # 1. Slice Target sequence span and apply rigid transformation
    t_ca_slice = target_ca_coords[start_idx:end_idx]
    t_ca_transformed = apply_transform(t_ca_slice, R, t)

    # 2. Derive boolean alignment flags for Query and Target residues
    q_aligned = np.zeros(len(q_ca_coords), dtype=bool)
    t_aligned = np.zeros(len(t_ca_transformed), dtype=bool)

    for q_idx, t_idx in enumerate(q2t_map):
        if isinstance(t_idx, int) and t_idx >= 0:
            if q_idx < len(q_aligned):
                q_aligned[q_idx] = True
            if t_idx < len(t_aligned):
                t_aligned[t_idx] = True

    # 3. Initialize viewer
    view = py3Dmol.view(width=800, height=500)

    # 4. Draw Query Backbone (Blue: Thick = Aligned, Thin = Unaligned)
    draw_variable_width_backbone(
        view, q_ca_coords, q_aligned, 
        color=q_color,
        aligned_radius=aligned_radius,
        unaligned_radius=unaligned_radius
    )

    # 5. Draw Target Backbone (Red: Thick = Aligned, Thin = Unaligned)
    draw_variable_width_backbone(
        view, t_ca_transformed, t_aligned, 
        color=t_color,
        aligned_radius=aligned_radius,
        unaligned_radius=unaligned_radius
    )

    # 6. Draw dashed 3D cylinder connectors for aligned pairs
    if False:
      draw_alignment_connectors_cylinders(
        view, 
        q_ca_coords, 
        t_ca_transformed, 
        q2t_map, 
        color="#616161", 
        radius=connector_radius
      )

    if title:
        print(f"Alignment Title: {title}")

    view.zoomTo()
    return view


# ==============================================================================
# 4. MAIN DATAFRAME WRAPPER
# ==============================================================================

def plot_alignment_by_target_ix_variable_width(
    target_ix: int,
    df: pd.DataFrame,
    project,
    q_ca_coords: np.ndarray,
    target_ca_coords: np.ndarray
):
    """Main execution wrapper for single color per chain with thickness variation."""
    target_row = df.iloc[target_ix]

    start_idx = int(target_row["start_idx"])
    end_idx = int(target_row["end_idx"])
    R = np.array(target_row["R"])
    t = np.array(target_row["t"])
    target_id = str(target_row["target_id"])
    rmsd = str(target_row["rmsd"])
    aln_len = str(target_row["alignment_length"])
    pide = str(target_row["sequence_identity"])
    
    # Build q2t_map adjusting for local slice offset
    q2t_map = build_q2t_map_from_segments(
        target_ix=target_ix,
        df=df,
        project=project,
        query_len=len(q_ca_coords),
    )

    
    # Render view
    view = visualize_pair_alignment_variable_width(
        q_ca_coords=q_ca_coords,
        target_ca_coords=target_ca_coords,
        q2t_map=q2t_map,
        R=R,
        t=t,
        start_idx=start_idx,
        end_idx=end_idx,
        title=f"Target: {target_id} | Align ID: {target_row['alignment_id']} | rmsd: {rmsd} | aln_len: {aln_len} | seq-id: {pide}",
        aligned_radius=0.3,     # Thick cylinders for aligned spans
        unaligned_radius=0.1,   # Thin cylinders for unaligned spans
        #connector_radius=0.1,

        q_color="#1f77b4",       # Blue for Query
        t_color="#d62728"        # Red for Target
    )

    return view.show()

# this is looking good with thick/thin lines - size?
if False:
  view_name = 'dom_610-835_FOLD'
  plot_alignment_by_target_ix_variable_width(
    target_ix = 0,
    df = project.views[view_name],
    project = project,
    q_ca_coords = project.query_ca_coords,
    target_ca_coords = project.coords
  )


# ==============================================================================
# 1. CORE GEOMETRY & MAP HELPERS
# ==============================================================================

def apply_transform(coords: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Applies rigid transformation (X @ R.T) + t to row-vector coordinates."""
    R = np.asarray(R, dtype=float)
    t = np.asarray(t, dtype=float).reshape(1, 3)
    return (coords @ R.T) + t


def build_q2t_map_from_segments(
    target_ix: int,
    df: pd.DataFrame,
    project,
    query_len: int,
    slice_start_offset: int = 0
) -> list:
    """Constructs 0-based q2t_map from project.segments for a target index."""
    q2t_map = [-1] * query_len
    target_align_id = df.iloc[target_ix]["alignment_id"]

    if isinstance(project.segments, pd.DataFrame):
        matched_segments = project.segments[project.segments["alignment_id"] == target_align_id]
        for _, seg in matched_segments.iterrows():
            q_start = int(seg["q_start"])
            s_start = int(seg["s_start"]) - slice_start_offset
            length = int(seg["length"])

            for offset in range(length):
                q_idx = q_start + offset
                s_idx = s_start + offset
                if 0 <= q_idx < query_len:
                    q2t_map[q_idx] = s_idx
    else:
        matched_segments = project.segments[target_align_id]
        for q_start, s_start, length in zip(
            matched_segments["q_start"],
            matched_segments["s_start"],
            matched_segments["length"]
        ):
            q_start = int(q_start)
            s_start = int(s_start) - slice_start_offset
            length = int(length)

            for offset in range(length):
                q_idx = q_start + offset
                s_idx = s_start + offset
                if 0 <= q_idx < query_len:
                    q2t_map[q_idx] = s_idx

    return q2t_map


# ==============================================================================
# 2. CYLINDER RENDERERS WITH TAPERING
# ==============================================================================

def GPT_draw_custom_backbone(
    view: py3Dmol.view,
    coords: np.ndarray,
    is_aligned_flags: np.ndarray,
    aligned_color: str,
    unaligned_color: str,
    aligned_radius: float = 0.5,
    unaligned_radius: float = 0.12,
    enable_taper: bool = False,
    dash_length: float = 1.5,
    gap_length: float = 1.0,
):
    """
    Draw a protein C-alpha backbone with:

      - aligned regions: continuous solid backbone
      - unaligned regions: continuous-phase dashed backbone

    The dashed pattern is carried across residue boundaries, so it does
    not restart at every C-alpha/C-alpha segment.

    Parameters
    ----------
    view
        py3Dmol view.

    coords
        (N, 3) Cartesian coordinates.

    is_aligned_flags
        Boolean array of length N. True means the residue is aligned.

    aligned_color, unaligned_color
        Colors for aligned and unaligned backbone.

    aligned_radius, unaligned_radius
        Cylinder radii.

    enable_taper
        If True, use endpoint radii. Usually False gives a cleaner
        backbone representation.

    dash_length
        Length of each visible dash in Angstroms.

    gap_length
        Length of each invisible gap in Angstroms.
    """

    coords = np.asarray(coords, dtype=float)
    flags = np.asarray(is_aligned_flags, dtype=bool)

    if len(coords) < 2:
        return

    if len(flags) != len(coords):
        raise ValueError(
            f"is_aligned_flags has length {len(flags)}, "
            f"but coords has length {len(coords)}"
        )

    # ------------------------------------------------------------
    # Helper: convert a coordinate to the py3Dmol dictionary format
    # ------------------------------------------------------------
    def point(p):
        return {
            "x": float(p[0]),
            "y": float(p[1]),
            "z": float(p[2]),
        }

    # ------------------------------------------------------------
    # Draw one solid backbone segment
    #
    # fromCap/toCap are disabled so adjacent cylinders visually merge
    # rather than producing obvious spherical/capped joints.
    # ------------------------------------------------------------
    def add_solid_segment(p1, p2, radius, color, radius2=None):
        spec = {
            "start": point(p1),
            "end": point(p2),
            "color": color,
            "radius": float(radius),
            "fromCap": 0,
            "toCap": 0,
        }

        if enable_taper and radius2 is not None:
            spec["radius2"] = float(radius2)

        view.addCylinder(spec)

    # ------------------------------------------------------------
    # Draw a dashed polyline.
    #
    # Crucially, dash phase is preserved across coordinate segments.
    # Thus:
    #
    #   Cα----Cα----Cα----Cα
    #
    # becomes something like:
    #
    #   ====  ====  ====  ==== 
    #
    # rather than:
    #
    #   == == | == == | == ==
    #
    # where every residue bond starts a new dash pattern.
    # ------------------------------------------------------------
    def add_dashed_polyline(points, radius, color):
        points = np.asarray(points, dtype=float)

        if len(points) < 2:
            return

        period = dash_length + gap_length

        # Distance already consumed within the current dash/gap period.
        phase = 0.0

        for i in range(len(points) - 1):
            a = points[i]
            b = points[i + 1]

            vec = b - a
            seg_len = np.linalg.norm(vec)

            if seg_len <= 1e-8:
                continue

            direction = vec / seg_len
            pos = 0.0

            while pos < seg_len:
                phase_in_period = phase % period

                # Are we currently in a visible dash or an invisible gap?
                if phase_in_period < dash_length:
                    remaining = dash_length - phase_in_period
                    visible = True
                else:
                    remaining = period - phase_in_period
                    visible = False

                step = min(remaining, seg_len - pos)

                if visible and step > 1e-8:
                    p1 = a + direction * pos
                    p2 = a + direction * (pos + step)

                    view.addCylinder({
                        "start": point(p1),
                        "end": point(p2),
                        "color": color,
                        "radius": float(radius),
                        "fromCap": 1,
                        "toCap": 1,
                    })

                pos += step
                phase += step

    # ------------------------------------------------------------
    # Split the backbone into aligned/unaligned runs.
    #
    # A segment is considered aligned only when BOTH endpoint residues
    # are aligned.
    # ------------------------------------------------------------
    segment_flags = flags[:-1] & flags[1:]

    i = 0
    n_segments = len(segment_flags)

    while i < n_segments:

        aligned = bool(segment_flags[i])

        j = i + 1
        while j < n_segments and bool(segment_flags[j]) == aligned:
            j += 1

        # Segment range is [i, j).
        #
        # The corresponding coordinate range is [i, j+1].
        run_coords = coords[i:j + 1]

        if True: ##aligned:
            # ----------------------------------------------------
            # Solid aligned backbone.
            #
            # We still need one cylinder per geometric Cα segment
            # because the protein backbone bends.
            #
            # But unlike the original implementation:
            #   - no dashed mode
            #   - no caps
            #   - no unnecessary overlay
            # ----------------------------------------------------
            for k in range(len(run_coords) - 1):

                r1 = aligned_radius
                r2 = aligned_radius

                if enable_taper:
                    # Both endpoints belong to an aligned run.
                    # Keep this here so tapering can be extended later.
                    r1 = aligned_radius
                    r2 = aligned_radius

                add_solid_segment(
                    run_coords[k],
                    run_coords[k + 1],
                    r1,
                    aligned_color,
                    r2,
                )

        else:
            # ----------------------------------------------------
            # Dashed unaligned backbone.
            #
            # One polyline is passed to the dash generator so that
            # the dash pattern continues across residue boundaries.
            # ----------------------------------------------------
            add_dashed_polyline(
                run_coords,
                unaligned_radius,
                unaligned_color,
            )

        i = j
            
def draw_custom_backbone(
    view: py3Dmol.view,
    coords: np.ndarray,
    is_aligned_flags: np.ndarray,
    aligned_color: str,
    unaligned_color: str,
    aligned_radius: float,
    unaligned_radius: float,
    enable_taper: bool = True
):
    """Draws backbone cylinders supporting both thickness variation and color switches."""
    num_pts = len(coords)
    if num_pts < 2:
        return

    node_radii = np.where(is_aligned_flags, aligned_radius, unaligned_radius)

    for i in range(num_pts - 1):
        p1 = {"x": float(coords[i][0]), "y": float(coords[i][1]), "z": float(coords[i][2])}
        p2 = {"x": float(coords[i + 1][0]), "y": float(coords[i + 1][1]), "z": float(coords[i + 1][2])}

        segment_aligned = is_aligned_flags[i] and is_aligned_flags[i + 1]
        color = aligned_color if segment_aligned else unaligned_color
         
        r1 = float(node_radii[i])
        r2 = float(node_radii[i + 1])

        cyl_spec = {
            "start": p1,
            "end": p2,
            "color": color,
            "dashed": False,
            "fromCap": 1,
            "toCap": 1
        }

        if enable_taper and r1 != r2:
            cyl_spec["radius"] = r1
            cyl_spec["radius2"] = r2
        else:
            cyl_spec["radius"] = r1 if segment_aligned else float(unaligned_radius)

        view.addCylinder(cyl_spec)


# ==============================================================================
# 3. VIEW GENERATOR WITH INTEGRATED METRICS TITLE
# ==============================================================================

def generate_superimposition_view(
    target_ix: int,
    df: pd.DataFrame,
    project,
    q_ca_coords: np.ndarray,
    target_ca_coords: np.ndarray,
    view_mode: str = "Thickness (Query Blue, Target Red)"
) -> py3Dmol.view:
    """Slices, transforms target coordinates, displays title metrics, and builds view."""
    if target_ix < 0 or target_ix >= len(df):
        raise IndexError(f"Target index {target_ix} is out of bounds for DataFrame with length {len(df)}.")

    target_row = df.iloc[target_ix]

    start_idx = int(target_row["start_idx"])
    end_idx = int(target_row["end_idx"])
    R = np.array(target_row["R"])
    t = np.array(target_row["t"])
    target_id = target_row.get("target_id", target_ix)
    align_id = target_row.get("alignment_id", "N/A")

    # Format header string from available metadata
    z_score = f"{target_row['z_score']:.2f}" if 'z_score' in target_row else 'N/A'
    rmsd = f"{target_row['rmsd']:.2f} Å" if 'rmsd' in target_row else 'N/A'
    align_len = target_row.get('alignment_length', 'N/A')
    
    seq_id_val = target_row.get('sequence_identity', None)
    if isinstance(seq_id_val, (int, float)):
        seq_id = f"{seq_id_val * 100:.1f}%" if seq_id_val <= 1.0 else f"{seq_id_val:.1f}%"
    else:
        seq_id = 'N/A'
        
    desc = target_row.get('description', '')

    # Print explicit title locked to this 3D canvas instance
    title_html = f"""
    <div style="font-family: monospace; font-size: 13px; background-color: #2b2b2b; color: #ffffff; padding: 10px; border-radius: 4px; margin-bottom: 8px;">
        <span style="color: #61afef; font-weight: bold;">[Rendered Structure Target IX: {target_ix}]</span> 
        <strong>ID:</strong> {target_id} | <strong>Align ID:</strong> {align_id}<br/>
        <strong>Z-score:</strong> {z_score} | <strong>RMSD:</strong> {rmsd} | <strong>Len:</strong> {align_len} | <strong>SeqID:</strong> {seq_id}<br/>
        <em>{desc}</em>
    </div>
    """
    display(HTML(title_html))

    # Slice & Transform Target
    t_ca_slice = target_ca_coords[start_idx:end_idx]
    t_ca_transformed = apply_transform(t_ca_slice, R, t)

    # Alignment map and boolean flags
    q2t_map = build_q2t_map_from_segments(
        target_ix, df, project, len(q_ca_coords) ###, start_idx
    )

    q_aligned = np.zeros(len(q_ca_coords), dtype=bool)
    t_aligned = np.zeros(len(t_ca_transformed), dtype=bool)

    for q_idx, t_idx in enumerate(q2t_map):
        if isinstance(t_idx, int) and t_idx >= 0:
            if q_idx < len(q_aligned):
                q_aligned[q_idx] = True
            if t_idx < len(t_aligned):
                t_aligned[t_idx] = True

    view = py3Dmol.view(width=800, height=500)

    if view_mode == "Thickness (Query Blue, Target Red)":
        draw_custom_backbone(view, q_ca_coords, q_aligned, "#1f77b4", "#1f77b4", 0.8, 0.10, enable_taper=True)
        draw_custom_backbone(view, t_ca_transformed, t_aligned, "#d62728", "#d62728", 0.5, 0.10, enable_taper=True)
    
    elif view_mode == "Dual Color (Blue/Gray & Red/Orange)":
        draw_custom_backbone(view, q_ca_coords, q_aligned, "#1f77b4", "#7f7f7f", 0.25, 0.25, enable_taper=False)
        draw_custom_backbone(view, t_ca_transformed, t_aligned, "#d62728", "#ff7f0e", 0.25, 0.25, enable_taper=False)
        
    elif view_mode == "Uniform (Query Blue, Target Red)":
        draw_custom_backbone(view, q_ca_coords, q_aligned, "#1f77b4", "#1f77b4", 0.25, 0.25, enable_taper=False)
        draw_custom_backbone(view, t_ca_transformed, t_aligned, "#d62728", "#d62728", 0.25, 0.25, enable_taper=False)

    draw_alignment_connectors_cylinders(
        view, q_ca_coords, t_ca_transformed, q2t_map, color="#616161", radius=0.08
    )

    view.zoomTo()
    return view


# ==============================================================================
# 4. WIDGET DASHBOARD
# ==============================================================================

def create_superimposition_widget(
    view_name,
    project,
    q_ca_coords: np.ndarray,
    target_ca_coords: np.ndarray
):
    """Widget interface where slider controls live metadata card, and 3D plot

    embeds the rendered target metadata inside its own viewport title.
    """
    df = project.views[view_name]
    max_idx = max(0, len(df) - 1)

    target_slider = widgets.IntSlider(
        value=0,
        min=0,
        max=max_idx,
        step=1,
        description="Target Index:",
        continuous_update=True,
        layout=widgets.Layout(width="450px")
    )

    view_dropdown = widgets.Dropdown(
        options=[
            "Thickness (Query Blue, Target Red)",
            "Dual Color (Blue/Gray & Red/Orange)",
            "Uniform (Query Blue, Target Red)"
        ],
        value="Thickness (Query Blue, Target Red)",
        description="Style:",
        layout=widgets.Layout(width="380px")
    )

    render_button = widgets.Button(
        description="Update 3-D view",
        button_style="primary",
        icon="play",
        layout=widgets.Layout(width="200px")
    )

    metadata_output = widgets.Output()
    canvas_output = widgets.Output()

    def update_metadata_display(target_ix: int):
        with metadata_output:
            clear_output(wait=True)
            cols = ['z_score', 'rmsd', 'alignment_length', 'sequence_identity', 'description']
            available_cols = [c for c in cols if c in df.columns]
            row_data = df.iloc[target_ix][available_cols].to_frame().T
            
            formatters = {}
            if 'z_score' in row_data: formatters['z_score'] = lambda x: f"{x:.2f}" if isinstance(x, (int, float)) else str(x)
            if 'rmsd' in row_data: formatters['rmsd'] = lambda x: f"{x:.2f} Å" if isinstance(x, (int, float)) else str(x)
            if 'sequence_identity' in row_data: formatters['sequence_identity'] = lambda x: f"{x*100:.1f}%" if isinstance(x, (int, float)) and x <= 1.0 else f"{x:.1f}"
            
            styled_html = row_data.to_html(
                index=False, 
                formatters=formatters, 
                classes="table table-striped table-hover",
                escape=False
            )
            
            display(HTML(f"""
            <div style="margin-top: 8px; margin-bottom: 8px; padding: 6px 12px; border-left: 4px solid #1f77b4; background-color: #f8f9fa;">
                <span style="font-size: 12px; color: #555;">SLIDER PREVIEW (Row {target_ix}):</span>
                {styled_html}
            </div>
            """))

    def on_slider_change(change):
        update_metadata_display(change['new'])

    target_slider.observe(on_slider_change, names='value')

    def on_render_clicked(b):
        with canvas_output:
            clear_output(wait=True)
            target_ix = target_slider.value
            mode = view_dropdown.value
            
            v = generate_superimposition_view(
                target_ix=target_ix,
                df=df,
                project=project,
                q_ca_coords=q_ca_coords,
                target_ca_coords=target_ca_coords,
###                view_mode=mode
            )
            v.show()

    render_button.on_click(on_render_clicked)

    controls = widgets.VBox([
        widgets.HBox([target_slider, view_dropdown]),
        render_button
    ])

    dashboard = widgets.VBox([
        controls,
        metadata_output,
        canvas_output
    ])
    
    display(dashboard)
    update_metadata_display(0)

#---------

from IPython.display import HTML, display
import numpy as np
import pandas as pd
import py3Dmol


def generate_transitive_superimposition_view(
    target_ixes: list[int],
    df: pd.DataFrame,
    project,
    target_ca_coords: np.ndarray,
    q_ca_coords: np.ndarray | None = None,
    colors: list[str] | None = None,
    radius: float = 0.25,
    show_query: bool = True,
    query_color: str = "#888888",
    query_radius: float = 0.45,
) -> py3Dmol.view:
    """Renders superimposition view for a list of target indices over an optional
    subtle Query structural scaffold.

    Parameters
    ----------
    target_ixes : list[int]
        List of row indices in `df` to slice and display.
    df : pd.DataFrame
        DataFrame containing target metadata, slice indices, R, and t.
    project : object
        The DaliScope project instance (providing project.query_length).
    target_ca_coords : np.ndarray
        Array containing C-alpha coordinates for targets.
    q_ca_coords : np.ndarray or None, default=None
        Query C-alpha coordinates. If None, falls back to project.query_ca_coords.
    colors : list[str] or None, default=None
        Custom color palette for rendering targets.
    radius : float, default=0.25
        Cylinder thickness radius for target backbone traces.
    show_query : bool, default=True
        Whether to render the Query backbone as a neutral background scaffold.
    query_color : str, default='#888888'
        Color for the Query backbone scaffold.
    query_radius : float, default=0.15
        Cylinder radius for the Query backbone (thinner than targets for visual hierarchy).

    Returns
    -------
    py3Dmol.view
        Configured py3Dmol view object.
    """
    if not isinstance(target_ixes, (list, tuple, np.ndarray)):
        target_ixes = [target_ixes]

    # Validate indices
    for target_ix in target_ixes:
        if target_ix < 0 or target_ix >= len(df):
            raise IndexError(
                f"Target index {target_ix} is out of bounds for DataFrame with length {len(df)}."
            )

    # Resolve Query coordinates
    if q_ca_coords is None and hasattr(project, "query_ca_coords"):
        q_ca_coords = project.query_ca_coords

    # Default color palette for target overlay
    default_colors = [
        "#d62728",
        "#ff7f0e",
        "#2ca02c",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#bcbd22",
        "#17becf",
    ]
    palette = colors if colors is not None else default_colors

    # Build Header metadata HTML block
    header_rows = []

    if show_query and q_ca_coords is not None:
        header_rows.append(
            f'<div style="margin-bottom: 6px; border-bottom: 1px solid #444; padding-bottom: 4px;">'
            f'<span style="color: {query_color}; font-weight: bold;">[Query Scaffold]</span> '
            f"<strong>Length:</strong> {project.query_length} aa"
            f"</div>"
        )

    for idx, target_ix in enumerate(target_ixes):
        target_row = df.iloc[target_ix]
        target_id = target_row.get("target_id", target_ix)
        align_id = target_row.get("alignment_id", "N/A")

        z_score = (
            f"{target_row['z_score']:.2f}"
            if "z_score" in target_row
            else "N/A"
        )
        rmsd = f"{target_row['rmsd']:.2f} Å" if "rmsd" in target_row else "N/A"
        align_len = target_row.get("alignment_length", "N/A")

        seq_id_val = target_row.get("sequence_identity", None)
        if isinstance(seq_id_val, (int, float)):
            seq_id = (
                f"{seq_id_val * 100:.1f}%"
                if seq_id_val <= 1.0
                else f"{seq_id_val:.1f}%"
            )
        else:
            seq_id = "N/A"

        desc = target_row.get("description", "")
        color = palette[idx % len(palette)]

        header_rows.append(
            f'<div style="margin-bottom: 4px;">'
            f'<span style="color: {color}; font-weight: bold;">[Target IX: {target_ix}]</span> '
            f"<strong>ID:</strong> {target_id} | <strong>Align ID:</strong> {align_id} | "
            f"<strong>Z-score:</strong> {z_score} | <strong>RMSD:</strong> {rmsd} | "
            f"<strong>Len:</strong> {align_len}/{project.query_length} | <strong>SeqID:</strong> {seq_id}<br/>"
            f"<em>{desc}</em>"
            f"</div>"
        )

    title_html = f"""
    <div style="font-family: monospace; font-size: 13px; background-color: #2b2b2b; color: #ffffff; padding: 10px; border-radius: 4px; margin-bottom: 8px;">
        {''.join(header_rows)}
    </div>
    """
    display(HTML(title_html))

    # Initialize py3Dmol view
    view = py3Dmol.view(width=800, height=500)

    # 1. Render Query Scaffold FIRST (places it in background layer)
    if show_query and q_ca_coords is not None:
        q_aligned = np.ones(len(q_ca_coords), dtype=bool)
        draw_custom_backbone(
            view,
            q_ca_coords,
            q_aligned,
            query_color,
            query_color,
            query_radius,
            query_radius,
            enable_taper=False,
        )

    # 2. Render Transformed Target Overlay
    for idx, target_ix in enumerate(target_ixes):
        target_row = df.iloc[target_ix]

        start_idx = int(target_row["start_idx"])
        end_idx = int(target_row["end_idx"])
        R = np.array(target_row["R"])
        t = np.array(target_row["t"])

        t_ca_slice = target_ca_coords[start_idx:end_idx]
        t_ca_transformed = apply_transform(t_ca_slice, R, t)
        t_aligned = np.ones(len(t_ca_transformed), dtype=bool)

        color = palette[idx % len(palette)]

        draw_custom_backbone(
            view,
            t_ca_transformed,
            t_aligned,
            color,
            color,
            radius,
            radius,
            enable_taper=False,
        )

    view.zoomTo()
    return view
