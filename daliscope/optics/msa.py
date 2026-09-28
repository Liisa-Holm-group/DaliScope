import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import sqlite3
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import math
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib import colors
import logomaker as lm
from pathlib import Path
from IPython.display import display, clear_output
from collections import OrderedDict

import tempfile, os
from IPython.display import display
from pymsaviz import MsaViz

def show_domain_alignments(
    df,
    left=None,
    right=None,
    seq_col="sequ_pileup",
    id_col=None,
    max_rows=20,
    **msaviz_kwargs,
):
    """
    Render a pyMSAviz plot for the sequence alignment, optionally sliced
    by 1-based inclusive residue position bounds [left, right].
    """
    sub_df = df.head(max_rows)
    if sub_df.empty:
        print("No sequences to display.")
        return None

    ids = sub_df[id_col].astype(str) if id_col else sub_df.index.astype(str)

    # 1. Fall back to sequence boundaries if left/right are not specified
    if left is None:
        left = 1

    if right is None:
        # Determine total length from the first valid sequence string in pileup
        first_seq = sub_df[seq_col].iloc[0]
        right = len(first_seq)

    # 2. Write FASTA to temp file
    fasta_text = "\n".join(
        f">{sid}\n{seq}" for sid, seq in zip(ids, sub_df[seq_col])
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        f.write(fasta_text)
        fasta_path = f.name

    # 3. Render via MsaViz with explicit start/end boundaries
    try:
        mv = MsaViz(
            fasta_path,
            start=left,
            end=right,
            wrap_length=80,
            show_count=False,
            **msaviz_kwargs,
        )
        fig = mv.plotfig()
        fig.suptitle(f"Segment {left}-{right} (n={len(sub_df)} rows)", y=1.02)
    finally:
        os.unlink(fasta_path)

    return fig

# Helper to convert name dict to RGB dict
def name_to_rgb(name_dict):
    return {
        k: tuple(int(c * 255) for c in colors.to_rgb(v))
        for k, v in name_dict.items()
    }

# Define the raw name mappings
AA_NAMES = {
    "A": "limegreen", "V": "limegreen", "L": "limegreen", "I": "limegreen", "M": "limegreen",
    "P": "cyan", "F": "magenta", "W": "magenta", "Y": "magenta",
    "N": "skyblue", "Q": "skyblue", "S": "skyblue", "T": "skyblue",
    "D": "red", "E": "red", "K": "blue", "R": "blue", "H": "orange",
    "C": "gold", "G": "gray", "X": "black"
}

DSSP_NAMES = {
    "H": "blue",     # Alpha Helix (3-state)
    "E": "red",     # Beta Strand (3-state)
    "C": "limegreen",   # Loop (3-state)
    "L": "limegreen",   # Loop (3-state)
    "T": "blue",    # Turn
    "S": "green",   # Bend
    "G": "orange",  # 3-10 Helix
    "I": "magenta", # Pi Helix
    "B": "cyan",    # Beta Bridge
    "-": "lightgray" # Coil/None
}

# The Master Config Object
COLOR_DATA = {
    'sequ_pileup': {
        'names': AA_NAMES,
        'rgb': name_to_rgb(AA_NAMES)
    },
    'dssp_pileup': {
        'names': DSSP_NAMES,
        'rgb': name_to_rgb(DSSP_NAMES)
    }
}

def rgb_lookup(mode):
    try:
        return COLOR_DATA[mode]["rgb"]
    except KeyError:
        raise ValueError(
            f"Unknown mode '{mode}'. "
            f"Expected one of {list(COLOR_DATA)}"
        )

class msa:
    def __init__(self, df, pileup_col, color_mapper=None, pseudocount=0.1, sort_by_seriation=True):
        # 1. Clean data and run seriation internally if requested
        if sort_by_seriation:
            from daliscope.mechanics.metrics import seriate

            # Drop rows missing the target pileup column or the sister column needed for seriation
            required_cols = list(set([pileup_col, 'dssp_pileup', 'sequ_pileup']))
            clean_df = df.dropna(subset=required_cols)[required_cols].copy()

            # Seriate and sort immediately
            clean_df = seriate(clean_df)
            self.df = clean_df.sort_values(by='dssp_order')
        else:
            self.df = df.dropna(subset=[pileup_col])

        self.pileup_col = pileup_col
        self.color_mapper = color_mapper
        self.pseudocount = pseudocount

        # PERSISTENT OBJECTS
        self._info_df = None
        self._logo_obj = None
        self._fig = None

        self.config = COLOR_DATA.get(pileup_col, {"rgb": {}})
        self._rgb_array = None

        # 2. Length consistency checks now run against the safely cleaned/sorted internal dataframe
        if self.df[pileup_col].isna().any():
            raise ValueError(f"Missing values in {pileup_col}")

        lengths = self.df[pileup_col].str.len()
        if lengths.nunique() != 1:
            raise ValueError(
                f"Inconsistent sequence lengths: {sorted(lengths.unique())}"
            )

    def _initialize_everything(self):
        """Run the heavy math AND the heavy rendering only once."""
        print("Initializing heavy assets...")        

        # 1. Math
        # counts include gaps => downweighting gappy positions
        raw_counts = lm.alignment_to_matrix(self.df[self.pileup_col], to_type='counts', characters_to_ignore='')
        counts_df = raw_counts + self.pseudocount
        p_df = counts_df.div(counts_df.sum(axis=1), axis=0)

        # Background math
        bg = pd.Series({'A': 0.0869, 'Q': 0.0392, 'L': 0.0976, 'S': 0.0720,
                        'R': 0.0581, 'E': 0.0627, 'K': 0.0505, 'T': 0.0559,
                        'N': 0.0486, 'G': 0.0707, 'M': 0.0232, 'W': 0.0130,
                        'D': 0.0544, 'H': 0.0230, 'F': 0.0384, 'Y': 0.0234,
                        'C': 0.0141, 'I': 0.0532, 'P': 0.0515, 'V': 0.0674})

        standard_20 = bg.index
        p_df = p_df[p_df.columns.intersection(standard_20)]
        self._info_df = (p_df * np.log2(p_df.div(bg[p_df.columns], axis=1))).fillna(0).clip(lower=0)

        # 2. Rendering (The 4-second part)
        # Create the figure and the logo object once
        self._fig, ax = plt.subplots(figsize=(6, 2))
        self._fig.canvas.toolbar_visible = False

        colormap = COLOR_DATA[self.pileup_col]["names"]

        self._logo_obj = lm.Logo(self._info_df, ax=ax, color_scheme=colormap, vpad=.1, width=.8)
        self._logo_obj.ax.set_ylabel('bits')
        #self._logo_obj.ax.set_ylim([0, 7.0])

        # Close the figure immediately so it doesn't display yet
        plt.close(self._fig)

    def logo(self, left=-1, right=None, positions=[], threshold=None):
        if self._logo_obj is None:
            self._initialize_everything()
        if right is None:
            right = len(self._info_df)
        ax = self._logo_obj.ax
        ax.set_xlim(left, right)

        # remove any previous threshold line before redrawing
        for line in list(ax.lines):
            if getattr(line, "_is_ic_threshold", False):
                line.remove()

        if threshold is not None:
            tl = ax.axhline(threshold, color="crimson", linewidth=2, linestyle="--")
            tl._is_ic_threshold = True
            position_ic = self._info_df.sum(axis=1)
            n_peaks = int((position_ic >= threshold).sum())
            print(f"{n_peaks} position(s) at or above threshold {threshold:.2f}")

        display(self._fig)
        return None

    def heatmap(self, left=None, right=None, bottom=None, top=None, figsize=(6,4)):
        lookup = rgb_lookup(self.pileup_col)
        if self._rgb_array is None:
            # Access self.config["rgb"] as required
            lookup = self.config.get("rgb", {})

            print("Generating RGB Heatmap array...")
            matrix = self.df[self.pileup_col].astype(str).tolist()

            self._rgb_array = np.array([
                [lookup.get(char, (255, 255, 255)) for char in row]
                for row in matrix
            ], dtype=np.uint8)

        plt.figure(figsize=figsize)
        plt.imshow(self._rgb_array, aspect="auto")
        if left or right: plt.xlim(left=left, right=right)
        if bottom or top: plt.ylim(bottom=bottom, top=top)
        plt.show()

from pandas.util import hash_pandas_object

def _get_cached_msa(df, mode):
    if not hasattr(_get_cached_msa, "cache"):
        _get_cached_msa.cache = OrderedDict()

    data_hash = hash_pandas_object(
        df[[mode]],
        index=True,
    ).sum()

    subset_key = (int(data_hash), mode)

    if subset_key in _get_cached_msa.cache:
        _get_cached_msa.cache.move_to_end(subset_key)
    else:
        if len(_get_cached_msa.cache) >= 20:
            _get_cached_msa.cache.popitem(last=False)

        _get_cached_msa.cache[subset_key] = msa(df, mode)

    return _get_cached_msa.cache[subset_key]

def old__get_cached_msa(df, mode):
    """Hidden logic to handle caching and heavy lifting."""
    if not hasattr(_get_cached_msa, "cache"):
        _get_cached_msa.cache = OrderedDict()

    # Key based on data content AND mode (dssp vs sequ)
    subset_key = (hash(tuple(df.index.values)), mode)

    if subset_key in _get_cached_msa.cache:
        _get_cached_msa.cache.move_to_end(subset_key)
    else:
        if len(_get_cached_msa.cache) >= 20: # Higher limit since we store objects individually
            _get_cached_msa.cache.popitem(last=False)

        # Call the heavy msa constructor
        _get_cached_msa.cache[subset_key] = msa(df, mode)

    return _get_cached_msa.cache[subset_key]

def plot_msa(df, left, right, plot_type="logo", data_col="sequ_pileup"):
    """Render an MSA logo-style plot from a dataframe carrying the
    pileup column directly -- no Project/view_name dependency, fully
    testable on any df in isolation."""
    logo = _get_cached_msa(df, data_col)
    if logo._info_df is None:
        logo._initialize_everything()
    return getattr(logo, plot_type)(left=left, right=right)
