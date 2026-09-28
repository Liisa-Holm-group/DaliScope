"""
plot_alignment_traces.py — Compare structural vs sequence alignment traces
==========================================================================

Visualizes two pairwise alignments of the same protein pair as dot-plot
traces, colored by agreement:
  - structural alignment only  (reference)
  - sequence alignment only
  - both alignments agree

Secondary structure is shown on the margins as colored bars.

Input format: gapped FASTA with 4 records per alignment:
    >seq1_id
    GAPPED-SEQUENCE-1
    >seq2_id
    GAPPED-SEQUENCE-2
    >ss1_id   (secondary structure for seq1, same gapping)
    HHHH--EEE
    >ss2_id
    EEE--HHHH

DSSP three-state: H=helix, E=strand, C/other=coil
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib import gridspec
from typing import Optional


# ------------------------------------------------------------------ #
# DSSP three-state reduction                                          #
# ------------------------------------------------------------------ #

SS_HELIX  = {'H', 'G', 'I'}
SS_STRAND = {'E', 'B'}
SS_COIL   = set('TSCL -.')    # everything else including gaps

SS_COLOR  = {
    'H': '#e05c5c',   # helix  — warm red
    'E': '#5b8dd9',   # strand — steel blue
    'C': '#c8c8c8',   # coil   — light grey
}

def _reduce_ss(char: str) -> str:
    """Map full DSSP alphabet to H / E / C."""
    if char in SS_HELIX:  return 'H'
    if char in SS_STRAND: return 'E'
    return 'C'


# ------------------------------------------------------------------ #
# Gapped FASTA parser                                                 #
# ------------------------------------------------------------------ #

def _parse_gapped_fasta(fasta_str: str) -> dict:
    """
    Parse a gapped FASTA string.
    Returns dict {header: sequence_string} preserving order.
    Accepts string or path.
    """
    import os
    if os.path.exists(fasta_str):
        with open(fasta_str) as f:
            fasta_str = f.read()

    records = {}
    header  = None
    buf     = []
    for line in fasta_str.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                records[header] = ''.join(buf)
            header = line[1:].split()[0]
            buf    = []
        else:
            buf.append(line)
    if header is not None:
        records[header] = ''.join(buf)

    return records

def _parse_offset_from_header(header: str) -> int:
    """
    Extract start offset from FASTA header if encoded as:
        >seqid/start-end   (HMMER style, 1-based)
        >seqid:start-end
        >seqid [start=N]
    Returns 0 if no offset found (assume full sequence).
    """
    import re
    # HMMER: >id/23-187
    m = re.search(r'/(\d+)-\d+', header)
    if m:
        return int(m.group(1)) - 1   # convert to 0-based

    # colon style: >id:23-187
    m = re.search(r':(\d+)-\d+', header)
    if m:
        return int(m.group(1)) - 1

    # explicit tag: >id [start=23]
    m = re.search(r'start=(\d+)', header)
    if m:
        return int(m.group(1)) - 1

    return 0

# ------------------------------------------------------------------ #
# Extract alignment trace from gapped sequences                       #
# ------------------------------------------------------------------ #

def _extract_trace(seq1_gapped: str, seq2_gapped: str,
                   offset1: int = 0, offset2: int = 0) -> np.ndarray:
    """
    Convert a gapped alignment to (pos_in_seq1, pos_in_seq2) pairs.
    Positions are 0-based in the UNGAPPED sequences.
    
    offset1, offset2: starting residue index in the full sequence,
    needed when the alignment covers only a subsequence (e.g. sequence
    alignment may start at residue 10 of a protein that struct alignment
    covers from residue 0).
    """
    if len(seq1_gapped) != len(seq2_gapped):
        raise ValueError(
            f"Gapped sequences must have equal column count, "
            f"got {len(seq1_gapped)} vs {len(seq2_gapped)}"
        )
    pairs = []
    i, j  = offset1, offset2
    for c1, c2 in zip(seq1_gapped, seq2_gapped):
        is_gap1 = c1 == '-'
        is_gap2 = c2 == '-'
        if not is_gap1 and not is_gap2:
            pairs.append((i, j))
        if not is_gap1:
            i += 1
        if not is_gap2:
            j += 1
    return np.array(pairs, dtype=np.int32)

def _ungapped_ss(ss_gapped: str) -> str:
    """Strip gaps from a gapped secondary structure string."""
    return ''.join(c for c in ss_gapped if c != '-')


# ------------------------------------------------------------------ #
# Secondary structure margin bar                                       #
# ------------------------------------------------------------------ #

def _draw_ss_bar(ax, ss_string: str, orientation: str,
                 length: int, lw: float = 6.0) -> None:
    """
    Draw a secondary structure bar along an axis margin.

    Parameters
    ----------
    ax          : the margin axes
    ss_string   : ungapped SS string (H/E/C after reduction)
    orientation : 'horizontal' (x-axis) or 'vertical' (y-axis)
    length      : total sequence length (for scaling)
    lw          : line width of the bar segments
    """
    ax.set_xlim(0, length) if orientation == 'horizontal' \
        else ax.set_ylim(0, length)

    reduced = [_reduce_ss(c) for c in ss_string]

    # draw runs of the same state as single line segments
    i = 0
    while i < len(reduced):
        state = reduced[i]
        j = i
        while j < len(reduced) and reduced[j] == state:
            j += 1
        color = SS_COLOR[state]
        if orientation == 'horizontal':
            ax.plot([i, j], [0.5, 0.5], color=color,
                    linewidth=lw, solid_capstyle='butt')
        else:
            ax.plot([0.5, 0.5], [i, j], color=color,
                    linewidth=lw, solid_capstyle='butt')
        i = j

    ax.set_xlim(0, length) if orientation == 'horizontal' \
        else ax.set_ylim(0, length)
    ax.axis('off')


# ------------------------------------------------------------------ #
# Main visualization function                                          #
# ------------------------------------------------------------------ #

def plot_alignment_traces(
    struct_fasta:   str,
    seq_fasta:      str,
    title:          str   = "Structural vs Sequence Alignment",
    dot_size:       float = 4.0,
    alpha:          float = 0.8,
    figsize:        tuple = (10, 10),
    ss_bar_width:   float = 0.04,   # fraction of figure size
    seq1_label:     Optional[str] = None,
    seq2_label:     Optional[str] = None,
    ax:             Optional[plt.Axes] = None,
) -> plt.Figure:
    """
    Compare two pairwise alignments as dot-plot traces.

    Parameters
    ----------
    struct_fasta : gapped FASTA string or path — structural alignment (reference).
                   Must contain exactly 4 records: seq1, seq2, ss1, ss2.
    seq_fasta    : gapped FASTA string or path — sequence alignment.
                   Must contain seq1 and seq2 records (ss records optional).
    title        : figure title
    dot_size     : marker size for alignment dots
    alpha        : dot transparency
    figsize      : figure size in inches
    ss_bar_width : width of secondary structure margin bars (fraction)
    seq1_label   : x-axis label (default: first record header)
    seq2_label   : y-axis label (default: second record header)
    ax           : if provided, draw into this axes (no SS bars, no title)

    Returns
    -------
    matplotlib Figure

    Color scheme
    ------------
    Blue  : structural alignment only
    Red   : sequence alignment only
    Green : both alignments agree (identical pair)

    Secondary structure bars
    ------------------------
    Red bar segment  : helix (H/G/I)
    Blue bar segment : strand (E/B)
    Grey bar segment : coil/other
    """
    # --- parse inputs ---
    struct_records = _parse_gapped_fasta(struct_fasta)
    seq_records    = _parse_gapped_fasta(seq_fasta)

    keys_s = list(struct_records.keys())
    keys_q = list(seq_records.keys())

    if len(keys_s) < 2:
        raise ValueError("struct_fasta must have at least 2 records")
    if len(keys_q) < 2:
        raise ValueError("seq_fasta must have at least 2 records")

    # sequences
    s1_gapped = struct_records[keys_s[0]]
    s2_gapped = struct_records[keys_s[1]]
    q1_gapped = seq_records[keys_q[0]]
    q2_gapped = seq_records[keys_q[1]]

    # secondary structures (optional in seq_fasta)
    ss1_gapped = struct_records.get(keys_s[2]) if len(keys_s) > 2 else None
    ss2_gapped = struct_records.get(keys_s[3]) if len(keys_s) > 3 else None

    # ungapped sequences for length
    s1 = s1_gapped.replace('-', '')
    s2 = s2_gapped.replace('-', '')
    L1 = len(s1)
    L2 = len(s2)

    ss1 = _ungapped_ss(ss1_gapped) if ss1_gapped else None
    ss2 = _ungapped_ss(ss2_gapped) if ss2_gapped else None



    # --- extract traces ---
    struct_trace = _extract_trace(s1_gapped, s2_gapped)   # (N, 2)
    seq_trace    = _extract_trace(q1_gapped, q2_gapped)   # (M, 2)

    # convert to sets of tuples for O(1) membership test
    struct_set = set(map(tuple, struct_trace.tolist()))
    seq_set    = set(map(tuple, seq_trace.tolist()))

    both_set   = struct_set & seq_set
    only_struct= struct_set - seq_set
    only_seq   = seq_set    - struct_set

    struct_records = _parse_gapped_fasta(struct_fasta)
    seq_records    = _parse_gapped_fasta(seq_fasta)

    # --- labels ---
    label1 = seq1_label or keys_s[0]
    label2 = seq2_label or keys_s[1]

    # ============================================================= #
    # Layout: SS bar (top) | main dot plot | SS bar (right)         #
    #         SS bar lives in thin axes above/right of main         #
    # ============================================================= #

    if ax is not None:
        # simple mode: just draw into provided axes, no margins
        _draw_dotplot(ax, only_struct, only_seq, both_set,
                      L1, L2, dot_size, alpha, label1, label2)
        return ax.get_figure()

    has_ss = ss1 is not None or ss2 is not None

    if has_ss:
        # gridspec: [ss_top | main] x [main | ss_right]
        bar_frac = ss_bar_width
        fig = plt.figure(figsize=figsize)
        gs  = gridspec.GridSpec(
            2, 2,
            width_ratios  = [1 - bar_frac, bar_frac],
            height_ratios = [bar_frac, 1 - bar_frac],
            hspace=0.02, wspace=0.02
        )
        ax_top   = fig.add_subplot(gs[0, 0])   # SS bar for seq1 (x-axis)
        ax_main  = fig.add_subplot(gs[1, 0])   # dot plot
        ax_right = fig.add_subplot(gs[1, 1])   # SS bar for seq2 (y-axis)
        fig.add_subplot(gs[0, 1]).axis('off')  # corner — blank

        # draw SS bars
        if ss1:
            _draw_ss_bar(ax_top,   ss1, 'horizontal', L1)
        else:
            ax_top.axis('off')

        if ss2:
            _draw_ss_bar(ax_right, ss2, 'vertical',   L2)
        else:
            ax_right.axis('off')

    else:
        fig, ax_main = plt.subplots(figsize=figsize)

    # --- main dot plot ---
    _draw_dotplot(ax_main, only_struct, only_seq, both_set,
                  L1, L2, dot_size, alpha, label1, label2)

    # --- legend ---
    legend_elements = [
        mpatches.Patch(color='#3a7bd5', alpha=alpha,
                       label=f'Structural only  (n={len(only_struct)})'),
        mpatches.Patch(color='#e05c5c', alpha=alpha,
                       label=f'Sequence only    (n={len(only_seq)})'),
        mpatches.Patch(color='#2ecc71', alpha=alpha,
                       label=f'Both agree       (n={len(both_set)})'),
    ]
    ax_main.legend(handles=legend_elements, loc='upper left',
                   framealpha=0.85, fontsize=9)

    # --- SS legend (if present) ---
    if has_ss:
        ss_patches = [
            mpatches.Patch(color=SS_COLOR['H'], label='Helix'),
            mpatches.Patch(color=SS_COLOR['E'], label='Strand'),
            mpatches.Patch(color=SS_COLOR['C'], label='Coil'),
        ]
        ax_main.legend(
            handles=legend_elements + ss_patches,
            loc='upper left', framealpha=0.85, fontsize=9
        )

    fig.suptitle(title, fontsize=13, fontweight='bold', y=1.01)
    fig.tight_layout()
    return fig, struct_records, seq_records


def _draw_dotplot(ax, only_struct, only_seq, both_set,
                  L1, L2, dot_size, alpha, label1, label2):
    """Draw the three-color dot plot into ax."""

    def _scatter(point_set, color, zorder):
        if not point_set:
            return
        pts = np.array(list(point_set))
        ax.scatter(pts[:, 0], pts[:, 1],
                   s=dot_size, c=color, alpha=alpha,
                   linewidths=0, zorder=zorder)

    # draw in order: struct (bottom), seq, both (top)
    _scatter(only_struct, '#3a7bd5', zorder=2)   # blue
    _scatter(only_seq,    '#e05c5c', zorder=3)   # red
    _scatter(both_set,    '#2ecc71', zorder=4)   # green

    ax.set_xlim(-1, L1)
    ax.set_ylim(-1, L2)
    ax.set_xlabel(f"{label1}  (residue index)", fontsize=10)
    ax.set_ylabel(f"{label2}  (residue index)", fontsize=10)
    ax.set_aspect('equal')

    # light grid to help read positions
    ax.grid(True, linewidth=0.3, color='#dddddd', zorder=1)
    ax.set_axisbelow(True)

    # diagonal reference line (identity diagonal)
    diag_end = min(L1, L2)
    ax.plot([0, diag_end], [0, diag_end],
            color='#aaaaaa', linewidth=0.6,
            linestyle='--', zorder=1, label='_diagonal')

    # stats text
    n_tot = len(only_struct) + len(only_seq) + len(both_set)
    pct   = 100.0 * len(both_set) / n_tot if n_tot else 0.0
    ax.text(0.98, 0.02,
            f"Agreement: {pct:.1f}%\n"
            f"Struct pairs: {len(only_struct) + len(both_set)}\n"
            f"Seq pairs:    {len(only_seq)    + len(both_set)}",
            transform=ax.transAxes,
            ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='white', alpha=0.8))


# ------------------------------------------------------------------ #
# Convenience: summary statistics only                                #
# ------------------------------------------------------------------ #

def alignment_agreement_stats_from_records(
        struct_records: dict, seq_records: dict) -> dict:
    """
    Compute agreement statistics from already-parsed FASTA dicts.
    Accepts the direct output of _parse_gapped_fasta().
    """
    keys_s = list(struct_records.keys())
    keys_q = list(seq_records.keys())

    #_validate_fasta_records(struct_records, "struct_records", min_records=2)
    #_validate_fasta_records(seq_records,    "seq_records",    min_records=2)

    struct_trace = _extract_trace(struct_records[keys_s[0]],
                                  struct_records[keys_s[1]])
    seq_trace    = _extract_trace(seq_records[keys_q[0]],
                                  seq_records[keys_q[1]])

    struct_set = set(map(tuple, struct_trace.tolist()))
    seq_set    = set(map(tuple, seq_trace.tolist()))
    both_set   = struct_set & seq_set
    union_set  = struct_set | seq_set

    return {
        "n_struct_only":        len(struct_set - seq_set),
        "n_seq_only":           len(seq_set - struct_set),
        "n_both":               len(both_set),
        "n_struct_total":       len(struct_set),
        "n_seq_total":          len(seq_set),
        "pct_agreement":        100.0 * len(both_set) / len(union_set)
                                if union_set else 0.0,
        "pct_struct_recovered": 100.0 * len(both_set) / len(struct_set)
                                if struct_set else 0.0,
    }


def alignment_agreement_stats(struct_fasta: str, seq_fasta: str) -> dict:
    """Thin wrapper — parses files then delegates to the core function."""
    return alignment_agreement_stats_from_records(
        _parse_gapped_fasta(struct_fasta),
        _parse_gapped_fasta(seq_fasta)
    )
