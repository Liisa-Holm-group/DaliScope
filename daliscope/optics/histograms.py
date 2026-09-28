import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import sqlite3
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import math
import seaborn as sns
from matplotlib_venn import venn2, venn3
import matplotlib.gridspec as gridspec
from matplotlib import colors
import plotly.express as pxnn
import plotly.graph_objects as go
import logomaker as lm
from pathlib import Path
from IPython.display import display, clear_output
from collections import OrderedDict

##################################
# Data Input & File Handling
##################################

def read_dali_tsv(tsvfile):
    "reads header-less TSV file assumed to be output of ISS_rundali.py"
    dali_tsv_colnames = [
        "cd1",
        "sbjct",
        "z_score",
        "rmsd",
        "ali-length",
        "sbjct-length-cropped",
        "seq-identity",
        "qstarts",
        "sstarts",
        "lengths",
        "rotation",
        "translation",
        "null1",
        "null2",
    ]
    df = pd.read_csv(tsvfile, sep="\t", header=None, names=dali_tsv_colnames)
    print(f"Read {df.shape[0]} rows with {df.shape[1]} columns from {tsvfile}")
    return df

def write_FASTA(df, outfile, id_col='sbjct', seq_col='sbjct-sequence'): 
    fasta_data = ">" + df[id_col] + "\n" + df[seq_col]
    try:
        with open(outfile, "w") as f:
            f.write("\n".join(fasta_data) + "\n")
        print(f"💾 File saved successfully: {outfile}")
    except Exception as e:
        print(f"❌ Error saving file: {e}") 

def write_TSV(df, filename=None):
    """
    Saves the DataFrame as a DALI-compatible TSV file.
    """
    if filename is None:
        # Default name if they are too lazy to provide one
        filename = "dali_export.tsv"
    
    # Ensure it ends with .tsv
    if not filename.endswith('.tsv'):
        filename += '.tsv'
        
    try:
        df.to_csv(filename, sep='\t', index=False)
        print(f"💾 File saved successfully: {filename}")
        
        # If in Google Colab, trigger the browser download automatically
        try:
            from google.colab import files
            files.download(filename)
        except ImportError:
            pass # Not in Colab, just saved locally
            
    except Exception as e:
        print(f"❌ Error saving file: {e}")

##################################
# Pfam & SQLite Database Operations
##################################

def get_project_root(sentinel_name='pfamdata'):
    """
    Traverses upwards from the current file until it finds a 
    specific folder or file that marks the project root.
    """
    # Start from the current directory
    current_path = Path.cwd()
    
    # Climb up the tree (parents) until the sentinel is found
    for path in [current_path] + list(current_path.parents):
        if (path / sentinel_name).exists():
            return path
            
    # Fallback to current directory if not found
    return current_path

def initialize_db(db_rel_path='pfamdata/pfam_data.db', verbose=False):
    root = get_project_root()
    db_path = (root / db_rel_path).resolve()
    if verbose:
        if db_path.exists():
            print(f"✅ Root found at: {root}")
            print(f"✅ Database mapped: {db_path}")
        else:
            print(f"❌ Error: Could not find database at {db_path}")
    return str(db_path)

# In your notebook:
# db = initialize_db()

def get_pfam_data(db_path, AC_list=None):
    """Replacement for pf.pfam_data_df filtering"""
    conn = sqlite3.connect(db_path)
    # Replaces pf.pfam_data_df[pf.pfam_data_df['AC'].isin([...])]
    placeholders = ','.join(['?'] * len(AC_list))
    query = f"SELECT DISTINCT AC, CL, ML, DE FROM pfam_entries WHERE AC IN ({placeholders}) or CL in ({placeholders})"
    df = pd.read_sql_query(query, conn, params=AC_list + AC_list)
    conn.close()
    return df
    
def get_pfam_composition(db_path, novel_df):
    conn = sqlite3.connect(db_path)
    #print('novel_df',novel_df.head())
    
    # 1. Create and populate the Temporary Table
    conn.execute("CREATE TEMPORARY TABLE temp_sbjcts (sbjct TEXT)")
    unique_sbjcts_input = novel_df[['sbjct']].drop_duplicates()
    unique_sbjcts_input.to_sql('temp_sbjcts', conn, if_exists='append', index=False)
    conn.execute("CREATE INDEX idx_temp_sbjct ON temp_sbjcts(sbjct)")

    # 2. Query all hits for these subjects
    # One subject will return multiple rows here if it has multiple domains
    query = """
    SELECT ts.sbjct, t1.AC as pfam, t1.CL as clan
    FROM temp_sbjcts ts
    JOIN hmmer_hits t2 ON ts.sbjct = t2.target_name
    JOIN pfam_entries t1 ON t2.query_accession = t1.AC
    """
    mapping_df = pd.read_sql_query(query, conn)
    conn.close()
    #print('mapping_df',mapping_df.head())

    # 3. Aggregation Logic (Crucial for Multi-Domain Subjects)
    # We use a set to ensure ['PF00001', 'PF00001'] becomes ['PF00001']
    agg_df = (
        mapping_df.groupby('sbjct')
        .agg({
            'pfam': lambda x: list(set(x.dropna())) if not x.isnull().all() else [],
            'clan': lambda x: list(set(x.dropna())) if not x.isnull().all() else []
        })
        .rename(columns={'pfam': 'pfamlist', 'clan': 'clanlist'})
        .reset_index()
    )
    #print('agg_df',agg_df.head())

    # 4. Final Merge
    # We left-join the aggregated lists back to the original novel_df
    final_df = pd.merge(novel_df, agg_df, on='sbjct', how='left')
    #print('final_df', final_df.head())
    
    # 5. Clean up NaNs in the lists for subjects with NO hits
    # Ensures 'stupid user' gets an empty list [] instead of NaN
    final_df['pfamlist'] = final_df['pfamlist'].apply(lambda x: x if isinstance(x, list) else [])
    final_df['clanlist'] = final_df['clanlist'].apply(lambda x: x if isinstance(x, list) else [])
    
    return final_df
    
def set_true_label(db_path, df, clan_id, database='AFDB2'):
    """
    Identifies subjects in the input DataFrame that belong to a specific Pfam Clan.
    Optimized via indexed temporary 'bridge' tables.
    """
    conn = sqlite3.connect(db_path)
    try:
        # 1. Create the table structure manually (Instant)
        conn.execute("DROP TABLE IF EXISTS temp_input")
        conn.execute("CREATE TEMPORARY TABLE temp_input (sbjct TEXT)")

        # 2. Convert to list of tuples for sqlite3
        unique_sbjcts = df[['sbjct']].drop_duplicates()
        data = unique_sbjcts[['sbjct']].values.tolist()

        # 3. Use a transaction for bulk upload
        conn.execute("BEGIN TRANSACTION")
        conn.executemany("INSERT INTO temp_input (sbjct) VALUES (?)", data)
        conn.execute("COMMIT")

        # 4. Create index AFTER upload (much faster than indexing while uploading)
        conn.execute("CREATE INDEX idx_temp_input_sbjct ON temp_input(sbjct)")

        # --- STEP 2: Create the AC Bridge Table ---
        # This narrows 5M+ hits down to just the ones the user cares about
        conn.execute("DROP TABLE IF EXISTS temp_ac")
        conn.execute(f"""
            CREATE TEMPORARY TABLE temp_ac AS 
            SELECT DISTINCT ti.sbjct, h.query_accession AS AC
            FROM temp_input ti
            INNER JOIN hmmer_hits h ON ti.sbjct = h.target_name
            WHERE h.database = '{database}'
        """)

        # --- STEP 3: Index the Bridge ---
        conn.execute("CREATE INDEX idx_temp_ac_join ON temp_ac(AC)")
        conn.execute("ANALYZE temp_ac")

        # --- STEP 4: Final Filter on Clan ---
        # We join our small bridge to pfam_entries to find the Clan members
        query = """
            SELECT DISTINCT t1.sbjct
            FROM temp_ac t1
            INNER JOIN pfam_entries t2 ON t1.AC = t2.AC
            WHERE t2.CL = ?
        """
        
        # Get the list of subjects that passed the Clan filter
        matching_sbjcts = pd.read_sql_query(query, conn, params=(clan_id,))['sbjct'].tolist()

        # If subject is in the list, assign "TRUE", else assign "Unassigned"
        df['TRUE'] = np.where(df['sbjct'].isin(matching_sbjcts), "TRUE", "Unassigned")
        
        return df

    except Exception as e:
        print(f"❌ Error during label processing: {e}")
        return df

    finally:
        # CLEANUP: Remove temp tables to keep the DB light
        conn.execute("DROP TABLE IF EXISTS temp_input")
        conn.execute("DROP TABLE IF EXISTS temp_ac")

    conn.close()
    
##################################
# Data Wrangling & Filtering
##################################

def calculate_penetrance(df, group_col='family', patterns=[]):
    """
    Calculates motif penetrance (%), absolute motif count, 
    and total group size for each family/clan.
    """
    # 1. Flag rows that match any of our target motifs
    required_cols = ['has_motif', 'motif']
    df_temp = df[required_cols].copy()
    df_temp['has_motif'] = df_temp['motif'].isin(patterns)

    # 2. Aggregate: sum (count of True), size (total), and mean (percentage)
    penetrance = (
        df_temp.groupby(group_col)['has_motif']
        .agg(
            motif_count='sum',
            group_size='size',
            penetrance_pct=lambda x: x.mean() * 100
        )
        .reset_index()
    )
    
    # 3. Cast motif_count to int (sum of booleans results in float sometimes)
    penetrance['motif_count'] = penetrance['motif_count'].astype(int)
    
    return penetrance.sort_values(['penetrance_pct', 'motif_count'], ascending=False)

def overwrite_hovertext(df, columns_to_join):
    """
    Concatenates specified columns into a single 'hovertext' column.
    Only processes columns that actually exist in the dataframe.
    """
    # 1. Filter out columns that don't exist to prevent KeyErrors
    valid_cols = [c for c in columns_to_join if c in df.columns]
    
    missing = set(columns_to_join) - set(valid_cols)
    if missing:
        print(f"⚠️ Warning: Columns not found and skipped: {missing}")

    if not valid_cols:
        print("❌ Error: No valid columns found to create hovertext.")
        df['hovertext'] = ''
        return df

    # 2. Functional concatenation: 
    # Convert to string, replace 'nan' with empty string, and join
    def join_row(row):
        # Only join elements that aren't empty/null
        parts = [str(val) for val in row if val is not None and str(val).lower() != 'nan']
        return " | ".join(parts)

    if len(valid_cols) == 1:
        print("valid_cols",valid_cols)
        df['hovertext'] = df[valid_cols[0]]
    else:
        df['hovertext'] = df[valid_cols].apply(join_row, axis=1)
        
    return df

def frequency(df, col):
    return df[col].value_counts()
    
def get_intersection(priority_list, available_columns):
    """
    Keeps the order of priority_list, but only includes 
    items found in available_columns.
    """
    available_set = set(available_columns) # for O(1) lookup speed
    return [item for item in priority_list if item in available_set]

def print_full_table(df, cols, sorts, max_rows=None):
    # check keys are present
    actual_cols = get_intersection(cols, df.columns)
    actual_sorts = get_intersection(sorts, df.columns)
    # This forces absolute visibility of every single cell
    with pd.option_context('display.max_rows', max_rows, 
                           'display.min_rows', None, 
                           'display.max_colwidth', None):
        display(df[actual_cols].sort_values(actual_sorts))
        
def summarize_results(df):
    rows, cols = df.shape
    print(f"Your data set consists of {rows:,} structural matches, each with {cols} descriptive attributes (columns).")
    
def filter_plane(df, x='z_score', y='ali-length', x_cutoff=None, y_cutoff=None, p1=None, p2=None):
    """
    Filters points based on three possible scenarios.
    p1 and p2 should be tuples like (x1, y1).
    """
    # Create a mask of True values (initially keep everything)
    mask = pd.Series([True] * len(df), index=df.index)

    # Scenario (i): X Cutoff
    if x_cutoff is not None:
        print(f"Applying X-cutoff: {x} > {x_cutoff}")
        mask &= (df[x] > x_cutoff)

    # Scenario (ii): Y Cutoff
    if y_cutoff is not None:
        print(f"Applying Y-cutoff: {y} > {y_cutoff}")
        mask &= (df[y] > y_cutoff)

    # Scenario (iii): Line Boundary
    if p1 is not None and p2 is not None:
        x1, y1 = p1
        x2, y2 = p2
        print(f"Applying Line-filter: Above line through {p1} and {p2}")
        
        # Avoid division by zero for vertical lines
        if x2 - x1 == 0:
            mask &= (df[x] > x1)
        else:
            # Calculate slope (m)
            m = (y2 - y1) / (x2 - x1)
            # Calculate expected Y on the line for each X in the dataframe: y = y1 + m(x - x1)
            expected_y = y1 + m * (df[x] - x1)
            mask &= (df[y] > expected_y)

    return df[mask].copy()

def group_top_labels(df, k=10, label_in="clan", label_out="hue", verbose=False):
    """
    Sorts df by z_score, keeps the top k labels from label_in, 
    and groups all others into 'other'.
    """
    # 1. Validation & Initialization
    if label_in not in df.columns:
        print(f"Warning: {label_in} not in columns. Setting to 'Unassigned'.")
        df[label_in] = "Unassigned"

    # 2. Identify top k unique categories based on z-score
    # We sort by z-score, take the label column, drop duplicates, and grab top k
    top_categories = (
        df.sort_values("z_score", ascending=False)[label_in]
        .dropna()
        .drop_duplicates()
        .head(k)
        .values
    )

    # 3. Vectorized assignment
    # Start everyone as 'other'
    df[label_out] = "other"
    
    # Update only those that are in the top_categories
    mask = df[label_in].isin(top_categories)
    df.loc[mask, label_out] = df[label_in]

    # Optional: Print for debugging
    if verbose:
        for cat in top_categories:
            print(f"Assigning top category: {cat}")

    return df

def sort_by_frequency(df, col):
    if df.empty:
        return df
    # Calculate counts
    counts = df[col].value_counts()
    # Map counts and sort: High counts (Unassigned) first -> Low counts last (on top)
    return df.assign(_tmp_counts=df[col].map(counts)) \
             .sort_values("_tmp_counts", ascending=False) \
             .drop(columns=["_tmp_counts"])

import pandas as pd

def set_motif(df, positions, k=4, verbose=True, base=1):
    """Creates a motif Series by concatenating characters from 'sequ_pileup' at specified indices.
    Lumps infrequent motifs (outside top-k) into an 'other' category.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a 'sequ_pileup' column.
    positions : list[int]
        0-indexed sequence positions to extract for the motif.
    k : int, default=4
        Number of top most frequent motifs to retain before pooling remaining into 'other'.
    verbose : bool, default=True
        Whether to print raw and pooled motif counts.

    Returns
    -------
    raw_counts : pd.Series
        Value counts of all unpooled raw extracted motifs.
    pooled_counts : pd.Series
        Value counts of top-k motifs + 'other'.
    motif_series : pd.Series
        The resulting pooled motif Series (aligned with df.index).
    """
    # 1. Generate raw motifs safely and quickly
    def _extract_motif(s):
        if not isinstance(s, str):
            return "-" * len(positions)
        n = len(s)
        return "".join(s[i-base] if 0 <= i < n else "-" for i in positions)

    raw_motifs = df["sequ_pileup"].apply(_extract_motif)
    raw_counts = raw_motifs.value_counts()

    if verbose:
        print("\n--- Frequent Motifs (Top 10 Raw) ---")
        print(raw_counts.head(10))

    # 2. Identify top-k motif groups
    topk = raw_counts.nlargest(k).index

    # 3. Consolidate infrequent motifs into 'other'
    motif_series = raw_motifs.where(raw_motifs.isin(topk), other="other")
    motif_series.name = "motif"
    
    pooled_counts = motif_series.value_counts()

    if verbose:
        print(f"\n--- Selected Groups (Top-{k} + 'other') ---")
        print(pooled_counts)

    return raw_counts, pooled_counts, motif_series
    

def old_set_motif(df, positions, k=4, verbose=True):
    """
    Creates a 'motif' column by concatenating characters from 'sequ_pileup' at specified indices.
    Lumps infrequent motifs (outside top-k) into an 'other' category.
    """
    # 1. Generate the motif strings
    # We use a list comprehension inside the lambda for speed
    df["raw_motif"] = df["sequ_pileup"].apply(
        lambda s: "".join([s[i] if len(s) > i else "-" for i in positions])
    )

    if verbose:
        print("\n--- Frequent Motifs (Top 10 Raw) ---")
        print(df["raw_motif"].value_counts().head(10))

    # 2. Identify top-k motif groups
    topk = df["raw_motif"].value_counts().nlargest(k).index

    # 3. Consolidate infrequent motifs into 'other'
    # .where(condition, other_value) keeps values where condition is True
    df["motif"] = df["raw_motif"].where(df["raw_motif"].isin(topk), other="other")

    if verbose:
        print(f"\n--- Selected Groups (Top-{k} + 'other') ---")
        print(df["motif"].value_counts())
        
    return df

def print_motif_table(db_path, df, group_col='family', motif_col='motif'):
    # 1. Get the list of family IDs
    family_list = sorted(df[group_col].unique())

    # 2. Pivot the motif data
    # Create the raw counts first
    motif_counts = pd.crosstab(df[group_col], df[motif_col])
    
    # NEW: Calculate total counts per group (family)
    group_totals = motif_counts.sum(axis=1).to_frame(name='Total_Hits')
    
    # Create percentages
    motif_pcts = motif_counts.div(motif_counts.sum(axis=1), axis=0) * 100
    motif_pcts = motif_pcts.round(2)

    # 3. Rename columns and combine with totals
    motif_pcts.columns = [f'%_{c}' for c in motif_pcts.columns]
    
    # Join the counts and percentages together
    # This puts 'Total_Hits' right next to the family ID
    motif_data = group_totals.join(motif_pcts).reset_index()

    # 4. Fetch the Pfam metadata
    tmp = get_pfam_data(db_path, family_list)

    # 5. Merge the Metadata with the combined data
    final_table = tmp.merge(
        motif_data, 
        left_on='AC', 
        right_on='family', 
        how='right'
    ).fillna(0)

    # 6. Define the display columns
    # Added 'Total_Hits' to the display list
    motif_cols = [col for col in final_table.columns if col.startswith('%_')]
    display_cols = ['family', 'CL', 'ML', 'DE', 'Total_Hits'] + motif_cols

    # 7. Print the full table
    print_full_table(final_table, display_cols, ['CL','DE'])
    
##################################
# Statistical Distributions & Comparison
##################################

def plot_distribution_ax(ax, data_frame, column_name, bins=25, color='seagreen'):
    """Helper: Plots histogram + KDE onto a specific axis."""
    data = data_frame[column_name].dropna()
    
    # Histogram
    ax.hist(data, bins=bins, density=True, alpha=0.5, 
            color=color, edgecolor='black', label='Hist')
    
    # KDE
    kde = gaussian_kde(data)
    x_range = np.linspace(data.min(), data.max(), 500)
    ax.plot(x_range, kde(x_range), color='darkred', lw=1.5, label='KDE')
    
    ax.set_title(f'Dist: {column_name}', fontsize=10)
    ax.grid(axis='y', linestyle='--', alpha=0.5)

def plot_distributions(df, columns, cols=2, bins=25, color='seagreen'):
    """
    Creates a grid of distribution plots for a list of columns.
    """
    n_vars = len(columns)
    rows = math.ceil(n_vars / cols)
    
    # Create the grid
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), squeeze=False)
    axes = axes.flatten() # Flatten to 1D array for easy iteration
        
    for i, col_name in enumerate(columns):
        plot_distribution_ax(axes[i], df, col_name, bins=bins, color=color)
    
    # Hide any unused subplots (e.g., if you have 7 plots in an 8-slot grid)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    plt.show()
    
def plot_violins(df, x_category, y_variables, cols=2, palette='viridis'):
    """
    Creates a grid of violin plots for multiple Y variables against one X category.
    """
    n_vars = len(y_variables)
    rows = math.ceil(n_vars / cols)
    
    # Create the figure and axes
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 6, rows * 5), squeeze=False)
    axes_flat = axes.flatten()
    
    for i, y_var in enumerate(y_variables):
        sns.violinplot(
            data=df,
            x=x_category,
            y=y_var,
            ax=axes_flat[i], # This tells Seaborn WHERE to plot
            hue=x_category,
            palette=palette,
            inner='quartile',
            legend=False
        )
        
        # Formatting individual subplots
        axes_flat[i].set_title(f'{y_var} by {x_category}', fontsize=12)
        axes_flat[i].tick_params(axis='x', rotation=45)
        axes_flat[i].set_xlabel('') # Hide x-label to save space
        
    # Hide empty subplots
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis('off')
        
    plt.tight_layout()
    plt.show()

def plot_multi_comparative_distributions(df_list, labels, columns, bars=True, cols=3, bins=25, palette='viridis'):
    """
    Plots a grid of overlaid histograms + KDEs for any number of dataframes.
    df_list: List of DataFrames [df1, df2, df3...]
    labels: List of strings ['Set A', 'Set B', 'Set C'...]
    """
    num_plots = len(columns)
    rows = math.ceil(num_plots / cols)
    
    # Get a color map to handle multiple series automatically
    cmap = plt.get_cmap(palette)
    colors = [cmap(i) for i in np.linspace(0, 0.9, len(df_list))]
    
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4))
    axes = axes.flatten()

    for i, col_name in enumerate(columns):
        ax = axes[i]
        
        for df, label, color in zip(df_list, labels, colors):
            if col_name in df.columns:
                data = df[col_name].dropna()
                if not data.empty:
                    # Histogram
                    if bars:
                        ax.hist(data, bins=bins, density=True, alpha=0.25, 
                            color=color, edgecolor='none', label=label)
                    
                    # KDE
                    if len(data) > 1:
                        kde = gaussian_kde(data)
                        x_range = np.linspace(data.min(), data.max(), 200)
                        ax.plot(x_range, kde(x_range), color=color, lw=2, label=label)
        ax.set_title(f'Distribution: {col_name}', fontweight='bold')
        ax.legend(prop={'size': 8}, frameon=False)
        ax.grid(axis='y', linestyle='--', alpha=0.2)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    plt.show() # Handles the display here so you don't need to return fig

def draw_family_pie_grid(df, family_col='family', motif_col='motif', cols=5, target=None):
    # 1. Get unique families and determine grid size
    families = sorted(df[family_col].unique())
    n_fams = len(families)
    rows = math.ceil(n_fams / cols)
    
    # 2. Setup colors for all unique motifs
    all_motifs = sorted(df[motif_col].unique())
    # Using a qualitative color map
    prop_cycle = plt.rcParams['axes.prop_cycle']
    color_list = prop_cycle.by_key()['color']
    motif_colors = {m: color_list[i % len(color_list)] for i, m in enumerate(all_motifs)}
    
    # Force 'CH' to be RoyalBlue if it exists
    if target in motif_colors:
        motif_colors[target] = colors.to_rgb('royalblue')

    # 3. Create the Figure
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
    axes = axes.flatten() if n_fams > 1 else [axes]
    
    # 4. Loop through each family and plot
    for i, fam in enumerate(families):
        ax = axes[i]
        fam_data = df[df[family_col] == fam]
        counts = fam_data[motif_col].value_counts()
        total_n = len(fam_data)
        
        # Prepare data for this specific pie
        labels = counts.index
        sizes = counts.values
        colors_for_pie = [motif_colors[m] for m in labels]
        
        # Draw the pie
        wedges, texts, autotexts = ax.pie(
            sizes, 
            autopct='%1.1f%%', 
            startangle=140, 
            colors=colors_for_pie,
            pctdistance=0.75,
            wedgeprops={'edgecolor': 'white', 'linewidth': 1}
        )
        
        # Style labels inside the pie
        plt.setp(autotexts, size=9, weight="bold", color="white")
        ax.set_title(f"{fam}\n(n={total_n})", fontsize=11, fontweight='bold')

    # 5. Cleanup: remove empty subplots if grid is larger than family count
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
        
    # 6. Add a Global Legend at the bottom
    legend_handles = [plt.Rectangle((0,0),1,1, color=motif_colors[m]) for m in all_motifs]
    fig.legend(legend_handles, all_motifs, loc='lower center', ncol=4, 
               title="Motif Categories", bbox_to_anchor=(0.5, 0.02))
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.show()

def draw_diverging_percentage_plot(
    df,
    hits,
    grouping_col='module',
    target_id_col='target_id',
    display='percentage',
    yscale='linear',
    max_groups=None,
    show_zeros=False
):
    """
    Plot MOTIF presence/absence as a diverging bar plot.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe.

    hits : list
        List of target_id values that are hits for the MOTIF.

    grouping_col : str, default='module'
        Column used to group the data.

    target_id_col : str, default='target_id'
        Column containing target IDs.

    display : {'percentage', 'count'}, default='percentage'
        Whether to display percentages of each group or raw counts.

    yscale : {'linear', 'log'}, default='linear'
        Y-axis scaling. 'log' uses a symmetric logarithmic scale
        (symlog) so that the zero-line remains meaningful.

    max_groups : int or None, default=None
        Maximum number of groups to display. Groups are selected
        according to the displayed positive quantity.

    show_zeros : bool, default=True
        Whether to show groups with zero MOTIF occurrences.
    """

    # Validate arguments
    if display not in {'percentage', 'count'}:
        raise ValueError(
            "display must be either 'percentage' or 'count'"
        )

    if yscale not in {'linear', 'log'}:
        raise ValueError(
            "yscale must be either 'linear' or 'log'"
        )

    if max_groups is not None and max_groups < 1:
        raise ValueError(
            "max_groups must be >= 1 or None"
        )

    # Work on a copy
    df = df.copy()

    # 1. Identify MOTIF presence/absence
    df['is_motif'] = df[target_id_col].isin(hits)

    # 2. Count presence/absence within each group
    counts = (
        df.groupby([grouping_col, 'is_motif'])
          .size()
          .unstack(fill_value=0)
    )

    # Ensure both columns exist
    if True not in counts.columns:
        counts[True] = 0

    if False not in counts.columns:
        counts[False] = 0

    # 3. Raw counts
    present = counts[True]
    absent = counts[False]
    totals = present + absent

    # 4. Optionally remove groups with zero MOTIF occurrences
    if not show_zeros:
        keep = present > 0

        present = present[keep]
        absent = absent[keep]
        totals = totals[keep]

    # 5. Convert to requested display quantity
    if display == 'percentage':
        positive = present / totals * 100
        negative = -(absent / totals * 100)
    else:
        positive = present
        negative = -absent

    # 6. Sort by displayed positive quantity
    sorted_idx = positive.sort_values(
        ascending=False
    ).index

    # 7. Limit number of groups
    if max_groups is not None:
        sorted_idx = sorted_idx[:max_groups]

    # 8. Plot
    fig, ax = plt.subplots(figsize=(12, 4))

    x = range(len(sorted_idx))

    ax.bar(
        x,
        positive.loc[sorted_idx],
        color='royalblue',
        label='MOTIF present'
    )

    ax.bar(
        x,
        negative.loc[sorted_idx],
        color='indianred',
        label='MOTIF absent'
    )

    # 9. Zero line
    ax.axhline(
        0,
        color='black',
        linewidth=1
    )

    # 10. Y-axis scaling
    if yscale == 'log':
        ax.set_yscale('symlog', linthresh=1)

    # 11. Axis labels and title
    if display == 'percentage':
        ax.set_ylabel(f'Percentage of {grouping_col} (%)')
        ax.set_title(
            f'MOTIF penetrance by {grouping_col}',
            pad=35
        )
    else:
        ax.set_ylabel('Number of rows')
        ax.set_title(
            f'MOTIF occurrences by {groupg_col}',
            pad=35
        )

    ax.set_xlabel(grouping_col)

    ax.set_xticks(list(x))
    ax.set_xticklabels(
        sorted_idx,
        rotation=45,
        ha='right'
    )

    # 12. Group counts: n_positive / n_total
    for i, group in enumerate(sorted_idx):

        n_positive = present.loc[group]
        n_total = totals.loc[group]

        ax.text(
            i,
            1.02,
            #f'{n_positive} / {n_total}',
            f'{n_positive}',
            transform=ax.get_xaxis_transform(),
            ha='center',
            va='bottom',
            fontsize=9,
            color='gray',
            weight='bold'
        )

        # Positive label
        pos = positive.loc[group]

        if pos != 0:
            label = (
                f'{abs(pos):.1f}%'
                if display == 'percentage'
                else f'{int(abs(pos))}'
            )

            ax.text(
                i,
                pos / 2,
                label,
                ha='center',
                va='center',
                color='white',
                weight='bold'
            )

        # Negative label
        neg = negative.loc[group]

        if neg != 0:
            label = (
                f'{abs(neg):.1f}%'
                if display == 'percentage'
                else f'{int(abs(neg))}'
            )

            ax.text(
                i,
                neg / 2,
                label,
                ha='center',
                va='center',
                color='white',
                weight='bold'
            )

    # Leave room for group-count annotations
    plt.subplots_adjust(top=0.88)

    ax.legend()

    plt.show()

def find_pattern_targets(
    df,
    patterns,
    *,
    target_col="target_id",
    sequence_col="sequ_pileup",
    position_base=1,
):
    """
    Return a nonredundant list of target_ids matching any pattern.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain `target_col` and `sequence_col`.

    patterns : list
        Each pattern is a list of (position, allowed_residues) pairs.

        Example:
            [
                [(34, ["H"]), (36, ["D"])],
                [(35, ["H"]), (37, ["H", "D"])],
            ]

        A pattern matches when ALL of its position/residue constraints
        are satisfied. Any matching pattern is sufficient.

    position_base : int, default=1
        Whether sequence positions are 1-based or 0-based.
        Use 1 for biological sequence numbering.

    Returns
    -------
    list
        Nonredundant target_ids, in order of first occurrence in df.
    """

    if position_base not in (0, 1):
        raise ValueError("position_base must be 0 or 1")

    required_cols = {target_col, sequence_col}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Normalize patterns once.
    normalized_patterns = []
    for pattern in patterns:
        constraints = []

        for position, residues in pattern:
            # Allow "D" as well as ["D"].
            if isinstance(residues, str):
                residues = [residues]

            residues = {str(r).upper() for r in residues}

            if not residues:
                raise ValueError(
                    f"Empty residue set at position {position}"
                )

            index = position - position_base

            if index < 0:
                raise ValueError(
                    f"Invalid position {position} for position_base={position_base}"
                )

            constraints.append((index, residues))

        normalized_patterns.append(constraints)

    matched = []

    for row in df.itertuples(index=False):
        target_id = getattr(row, target_col)
        sequence = getattr(row, sequence_col)

        if not isinstance(sequence, str):
            continue

        sequence = sequence.upper()

        # OR over patterns
        for pattern in normalized_patterns:

            # AND over constraints within a pattern
            if all(
                index < len(sequence) and sequence[index] in residues
                for index, residues in pattern
            ):
                matched.append(target_id)
                break

    # Remove duplicates while preserving dataframe order.
    return list(dict.fromkeys(matched))

def old_draw_diverging_percentage_plot(df, family_col='family', motif_col='motif', targetlist=None):
    # 1. Calculate Presence/Absence counts
    df['is_target'] = df[motif_col].isin(targetlist)
    counts = df.groupby([family_col, 'is_target']).size().unstack(fill_value=0)
    
    # Ensure both True (Present) and False (Absent) columns exist
    if True not in counts.columns: counts[True] = 0
    if False not in counts.columns: counts[False] = 0
    
    # 2. Convert to Percentages
    totals = counts.sum(axis=1)
    pct_present = (counts[True] / totals) * 100
    pct_absent = -(counts[False] / totals) * 100  # Negative for downward bar
    
    # 3. Sort by percentage of presence
    sorted_idx = pct_present.sort_values(ascending=False).index
    
    # 4. Plot
    fig, ax = plt.subplots(figsize=(12, 7))
    
    ax.bar(sorted_idx, pct_present.loc[sorted_idx], color='royalblue', label=f'Present')
    ax.bar(sorted_idx, pct_absent.loc[sorted_idx], color='indianred', label='Absent')
    
    # 5. Styling & Labels
    ax.axhline(0, color='black', linewidth=1)
    ax.set_ylim(-115, 115) # Room for text
    ax.set_ylabel('Percentage (%)')
    ax.set_title(f"Normalized Presence of '{targetlist}' Motif by Family")
    
    # Add n=XX labels and percentage text
    for i, fam in enumerate(sorted_idx):
        p_val = pct_present.loc[fam]
        a_val = pct_absent.loc[fam]
        total_n = totals.loc[fam]
        
        # Total sample size label
        ax.text(i, 105, f'n={total_n}', ha='center', va='bottom', fontsize=9, color='gray', weight='bold')
        # Percentage labels
        if p_val > 5: ax.text(i, p_val/2, f'{p_val:.1f}%', ha='center', color='white', weight='bold')
        if abs(a_val) > 5: ax.text(i, a_val/2, f'{abs(a_val):.1f}%', ha='center', color='white', weight='bold')

    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
    plt.show()

##################################
# Core Visualization Engines
##################################

# Global Plotting Configuration
_SETTINGS = {
    'engine': 'plotly',             # Options: 'plotly', 'matplotlib', 'seaborn'
    'theme': 'plotly_white',        # Default theme/style
    'template': 'ggplot2',          # Aesthetic template
    'default_figsize': (10, 6),     # For Matplotlib/Seaborn engines
    'color_continuous_scale': 'Viridis',
    'show_grid': True,
    'precision': 2                  # Decimal precision for hovertext/labels
}

def resolve_setting(key, override_val=None, default='plotly'):
    """Helper to prioritize: Local Override > Global Setting > Hard Default."""
    if override_val is not None:
        return override_val
    return _SETTINGS.get(key, default)

def plot_scatter(df, x='z_score', y='ali-length', color=None, marker=None, size=None, hover_col=None, engine=None, **kwargs):
    # 1. Resolve which engine to use (Consulting our global _SETTINGS)
    selected_engine = resolve_setting('engine', engine)

    # 1. UNIVERSAL SORTING: Smallest categories on top
    plot_df = df.copy() # BAD!
    sort_col = color if color else marker
    if sort_col and sort_col in plot_df.columns:
        if not pd.api.types.is_numeric_dtype(plot_df[sort_col]):
            counts = plot_df[sort_col].value_counts()
            plot_df = plot_df.assign(_tmp=plot_df[sort_col].map(counts)) \
                             .sort_values("_tmp", ascending=False) \
                             .drop(columns="_tmp")

    # 2. Route to the correct renderer
    if selected_engine == 'plotly':
        return _render_plotly(plot_df, x, y, color, marker, size, hover_col=hover_col, **kwargs)
    else:
        # Pass hovertext to MPL as well, in case we want to use it for labels later
        return _render_mpl(plot_df, x, y, color, marker, size, hover_name=hover_col, **kwargs)

def _render_plotly(df, x, y, color, marker, size, hover_col=None, **kwargs):
    import plotly.express as px
    
    # 1. Categorical check for color
    is_cat = False
    if color and color in df.columns:
        if df[color].dtype == 'object' or df[color].dtype == 'bool':
            is_cat = True

    # 2. Build parameter dictionary
    params = {
        "data_frame": df, "x": x, "y": y, "color": color, "symbol": marker,
        "size": size, "title": kwargs.get('title'), "opacity": kwargs.get('opacity', 0.8),
        "template": "plotly_white", "size_max": 12
    }

    # 3. Handle Hovertext Logic
    # Use the user-defined column if it exists, otherwise check for a default
    hover_col = hover_col if hover_col in df.columns else None
    
    if hover_col:
        params["hover_name"] = hover_col
        # hover_data ensures the technical values are still visible in the box
        params["hover_data"] = {col: True for col in [x, y, color, marker] if col}
    
    # 4. Color Scale Logic
    if is_cat:
        params["color_discrete_sequence"] = px.colors.qualitative.D3
    else:
        params["color_continuous_scale"] = "viridis"

    fig = px.scatter(render_mode="webgl", **params)

    # 5. UI Polishing
    fig.update_traces(hoverlabel=dict(namelength=0)) 
    if not is_cat and marker:
        fig.update_layout(legend=dict(yanchor="bottom", y=0.02, xanchor="right", x=0.98, bgcolor="rgba(255,255,255,0.5)"))
        
    return fig
    
def _render_mpl(df, x, y, color, marker, size, **kwargs):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 1. Determine if color is continuous
    is_cont = False
    if color and pd.api.types.is_numeric_dtype(df[color]):
        is_cont = True

    # 2. Plotting
    if is_cont:
        # Continuous Branch
        pts = ax.scatter(df[x], df[y], c=df[color], cmap='viridis', 
                         s=df[size] if size else 20, alpha=0.8, edgecolors='w', lw=0.5)
        fig.colorbar(pts, ax=ax, label=color)
        
        if marker:
            # Create manual legend for markers only
            unique_m = df[marker].unique()
            m_styles = ['o', 's', 'D', '^', 'v', 'p', 'P', '*']
            m_map = dict(zip(unique_m, m_styles))
            h = [Line2D([0], [0], marker=m_map[m], color='w', label=str(m), 
                        markerfacecolor='gray', markersize=8) for m in unique_m]
            ax.legend(handles=h, title=marker, loc='lower right')
    else:
        # Categorical Branch (Seaborn)
        sns.scatterplot(data=df, x=x, y=y, hue=color, style=marker, size=size, palette="tab10", ax=ax)
        
        # --- THE FIX: FILTER THE LEGEND ---
        handles, labels = ax.get_legend_handles_labels()
        
        # We keep handles only if their label is NOT the 'size' column name 
        # and not a numeric value associated with size scaling.
        new_handles, new_labels = [], []
        
        # Stop-words: titles we want to ignore in the legend
        forbidden_titles = [size] if size else []
        
        # Iterate through current legend and skip the size section
        skip_mode = False
        for h, l in zip(handles, labels):
            if l == size:
                skip_mode = True # Found the 'size' header, start skipping
                continue
            if skip_mode:
                # In Seaborn, size values are usually numeric strings in the legend
                # If we hit a new non-numeric label, we've likely left the size section
                try:
                    float(l.replace('<', '').replace('>', ''))
                    continue # It's a size value, skip it
                except ValueError:
                    skip_mode = False # Not a size value, stop skipping
            
            new_handles.append(h)
            new_labels.append(l)

        ax.legend(new_handles, new_labels, bbox_to_anchor=(1.05, 1), loc='upper left')

    ax.set_title(kwargs.get('title'))
    fig.tight_layout()
    plt.close(fig) 
    return fig

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2


def plot_safe_comparison_with_marginal(
    df1,
    df2,
    on_category,
    x_var,
    y_var,
    labels=None,
    max_points=15000,
    sort_by_frequency=None,
):
    """Plots two DataFrames with marginal distributions.

    All Venn diagram counts and scatter points represent total ROW INSTANCES.
    """
    if callable(sort_by_frequency):
        df1_sorted = sort_by_frequency(df1, on_category)
        df2_sorted = sort_by_frequency(df2, on_category)
    else:
        df1_sorted, df2_sorted = df1.copy(), df2.copy()

    if labels is None:
        labels = ["df1", "df2"]

    palette = {
        "Common": "#0072B2",  # Deep Blue
        labels[0]: "#D55E00",  # Vermillion
        labels[1]: "#CC79A7",  # Purple
    }

    # 1. Category Sets for Overlap Matching
    cat_set1 = set(df1_sorted[on_category].dropna().unique())
    cat_set2 = set(df2_sorted[on_category].dropna().unique())
    common_cats = cat_set1.intersection(cat_set2)

    # 2. Partition Rows by Category Membership
    # Rows in df1 whose category is shared with df2
    df1_common_rows = df1_sorted[df1_sorted[on_category].isin(common_cats)]
    # Rows in df1 whose category exists ONLY in df1
    df1_only_rows = df1_sorted[~df1_sorted[on_category].isin(common_cats)]

    # Rows in df2 whose category is shared with df1
    df2_common_rows = df2_sorted[df2_sorted[on_category].isin(common_cats)]
    # Rows in df2 whose category exists ONLY in df2
    df2_only_rows = df2_sorted[~df2_sorted[on_category].isin(common_cats)]

    # 3. Exact Row Instance Counts for Venn Regions
    n_only1 = len(df1_only_rows)
    n_only2 = len(df2_only_rows)
    n_common = len(df1_common_rows) + len(df2_common_rows)

    # 4. Build Plotting DataFrame for Scatter/Marginals
    d1_sub = df1_sorted[[on_category, x_var, y_var]].copy()
    d1_sub["Series"] = d1_sub[on_category].apply(
        lambda x: "Common" if x in common_cats else labels[0]
    )

    d2_only = df2_only_rows[[on_category, x_var, y_var]].copy()
    d2_only["Series"] = labels[1]

    plot_df = pd.concat([d1_sub, d2_only], ignore_index=True)

    if len(plot_df) > max_points:
        plot_df = plot_df.sample(n=max_points, random_state=42)

    # 5. Layout & Venn Plotting
    fig = plt.figure(figsize=(18, 8))
    gs = gridspec.GridSpec(4, 4, hspace=0.3, wspace=0.3)

    ax_v = fig.add_subplot(gs[1:3, 0])

    # Pass explicit disjoint row counts as a 3-tuple: (Ab, aB, AB)
    v = venn2(
        subsets=(n_only1, n_only2, n_common),
        set_labels=(
            f"{labels[0]}\n(rows: {len(df1):,})",
            f"{labels[1]}\n(rows: {len(df2):,})",
        ),
        ax=ax_v,
    )

    if v:
        if v.get_patch_by_id("11"):
            v.get_patch_by_id("11").set_color(palette["Common"])
        if v.get_patch_by_id("10"):
            v.get_patch_by_id("10").set_color(palette[labels[0]])
        if v.get_patch_by_id("01"):
            v.get_patch_by_id("01").set_color(palette[labels[1]])

    # 6. Scatter & Marginals
    ax_main = fig.add_subplot(gs[1:4, 1:3])
    ax_x = fig.add_subplot(gs[0, 1:3], sharex=ax_main)
    ax_y = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    sns.scatterplot(
        data=plot_df,
        x=x_var,
        y=y_var,
        hue="Series",
        style="Series",
        palette=palette,
        alpha=0.6,
        s=60,
        ax=ax_main,
        legend=True,
    )

    sns.kdeplot(
        data=plot_df,
        x=x_var,
        hue="Series",
        palette=palette,
        ax=ax_x,
        legend=False,
        fill=True,
        alpha=0.3,
    )
    sns.kdeplot(
        data=plot_df,
        y=y_var,
        hue="Series",
        palette=palette,
        ax=ax_y,
        legend=False,
        fill=True,
        alpha=0.3,
    )

    ax_x.axis("off")
    ax_y.axis("off")
    ax_main.grid(True, linestyle=":", alpha=0.6)

    plt.show()

def old_plot_safe_comparison_with_marginal(df1, df2, on_category, x_var, y_var, labels=None, max_points=15000):
    """
    Plots two dataframes with marginal distributions, 
    ensuring rare categories are drawn on top.
    """
        # Apply sorting to both DataFrames
    df1_sorted = sort_by_frequency(df1, on_category)
    df2_sorted = sort_by_frequency(df2, on_category)

    # Series labels
    if labels is None: labels=['df', 'df_ref']
    
    # 1. High-Contrast Palette (Accessible)
    palette = {
        'Common': '#0072B2',        # Deep Blue
        labels[0]: '#D55E00',       # Orange/Vermillion
        labels[1]: '#F0E442'        # Bright Yellow
    }

    # --- [Set Logic: Identical to previous high-speed version] ---
    set1_ids = set(df1_sorted[on_category].dropna().unique())
    set2_ids = set(df2_sorted[on_category].dropna().unique())
    common_ids = set1_ids.intersection(set2_ids)
    only1_ids = set1_ids - set2_ids
    only2_ids = set2_ids - set1_ids

    d1_sub = df1_sorted[df1_sorted[on_category].isin(set1_ids)][[on_category, x_var, y_var]].copy()
    d2_sub = df2_sorted[df2_sorted[on_category].isin(only2_ids)][[on_category, x_var, y_var]].copy()
    d1_sub['Series'] = d1_sub[on_category].apply(lambda x: 'Common' if x in common_ids else labels[0])
    d2_sub['Series'] = labels[1]
    plot_df = pd.concat([d1_sub, d2_sub], ignore_index=True)

    if len(plot_df) > max_points:
        plot_df = plot_df.sample(n=max_points, random_state=42)
    # -----------------------------------------------------------

    # 2. Create complex layout using GridSpec
    fig = plt.figure(figsize=(18, 8))
    gs = gridspec.GridSpec(4, 4, hspace=0.3, wspace=0.3)

    # Venn Diagram (Left side, centered vertically)
    ax_v = fig.add_subplot(gs[1:3, 0])
    v = venn2(subsets=(len(only1_ids), len(only2_ids), len(common_ids)), 
              set_labels=(labels[0], labels[1]), ax=ax_v)
    if v:
        if v.get_patch_by_id('11'): v.get_patch_by_id('11').set_color(palette['Common'])
        if v.get_patch_by_id('10'): v.get_patch_by_id('10').set_color(palette[labels[0]])
        if v.get_patch_by_id('01'): v.get_patch_by_id('01').set_color(palette[labels[1]])

    # Main Scatter Plot
    ax_main = fig.add_subplot(gs[1:4, 1:3])
    # Marginal X (Top)
    ax_x = fig.add_subplot(gs[0, 1:3], sharex=ax_main)
    # Marginal Y (Right)
    ax_y = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    # 3. Plotting the Scatter
    sns.scatterplot(data=plot_df, x=x_var, y=y_var, hue='Series', style='Series',
                    palette=palette, alpha=0.6, s=60, ax=ax_main, legend=True)

    # 4. Plotting the Marginals
    # Common practice: use KDE or Hist with 'fill=True'
    sns.kdeplot(data=plot_df, x=x_var, hue='Series', palette=palette, 
                ax=ax_x, legend=False, fill=True, alpha=0.3)
    sns.kdeplot(data=plot_df, y=y_var, hue='Series', palette=palette, 
                ax=ax_y, legend=False, fill=True, alpha=0.3)

    # Cleanup Marginal Axes
    ax_x.axis('off')
    ax_y.axis('off')
    ax_main.grid(True, linestyle=':', alpha=0.6)
    
    plt.show()

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3


def plot_3way_bright_comparison(
    df1,
    df2,
    df3,
    on_cat,
    x_var,
    y_var,
    labels=("set1", "set2", "set3"),
    plot_KDE=False,
):
    """Plots a 3-way comparison with high-saturation scatter points using distinct marker

    shapes for 3-way, 2-way, and 1-way intersection tiers.
    """
    intersection_map = {
        "100": f"Unique to {labels[0]}",
        "010": f"Unique to {labels[1]}",
        "001": f"Unique to {labels[2]}",
        "110": f"{labels[0]} & {labels[1]}",
        "101": f"{labels[0]} & {labels[2]}",
        "011": f"{labels[1]} & {labels[2]}",
        "111": "Common to All",
    }

    # 1. Marker Map: Distinct shapes by intersection complexity
    marker_map = {
        f"Unique to {labels[0]}": "o",  # Circle
        f"Unique to {labels[1]}": "o",  # Circle
        f"Unique to {labels[2]}": "o",  # Circle
        f"{labels[0]} & {labels[1]}": "X",  # Filled Cross (capital X)
        f"{labels[0]} & {labels[2]}": "^",  # Triangle Up
        f"{labels[1]} & {labels[2]}": "v",  # Triangle Down
        "Common to All": "s",  # Square
    }
    
    # 2. High-Saturation Palette
    saturated_colors = plt.cm.Set1.colors
    color_palette = {
        name: saturated_colors[i]
        for i, name in enumerate(intersection_map.values())
    }

    # Extract categories per dataframe for overlap evaluation
    cats1 = set(df1[on_cat].dropna())
    cats2 = set(df2[on_cat].dropna())
    cats3 = set(df3[on_cat].dropna())

    # 3. Build Row-Instance DataFrames
    dfs = []
    for df in [df1, df2, df3]:
        sub = df[[on_cat, x_var, y_var]].copy()
        sub["code"] = sub[on_cat].apply(
            lambda v: f"{'1' if v in cats1 else '0'}{'1' if v in cats2 else '0'}{'1' if v in cats3 else '0'}"
        )
        sub["Intersection"] = sub["code"].map(intersection_map)
        dfs.append(sub)

    plot_df = pd.concat(dfs, ignore_index=True)

    # 4. Calculate Exact Disjoint Row Instance Counts for Venn Regions
    venn_counts = (
        len(plot_df[plot_df["code"] == "100"]),  # 100
        len(plot_df[plot_df["code"] == "010"]),  # 010
        len(plot_df[plot_df["code"] == "110"]),  # 110
        len(plot_df[plot_df["code"] == "001"]),  # 001
        len(plot_df[plot_df["code"] == "101"]),  # 101
        len(plot_df[plot_df["code"] == "011"]),  # 011
        len(plot_df[plot_df["code"] == "111"]),  # 111
    )

    fig = plt.figure(figsize=(22, 11))
    gs = gridspec.GridSpec(
        4, 5, hspace=0.4, wspace=0.3, width_ratios=[2, 1, 1, 1, 0.5]
    )

    # 5. Venn Diagram with Saturated Patches
    ax_v = fig.add_subplot(gs[0:4, 0])
    set_labels = [
        f"{name}\n(rows: {len(df):,})"
        for name, df in zip(labels, [df1, df2, df3])
    ]

    v = venn3(subsets=venn_counts, set_labels=set_labels, ax=ax_v)

    if v:
        for label in v.subset_labels:
            if label:
                label.set_fontsize(14)
        for label in v.set_labels:
            if label:
                label.set_fontsize(16)

        for code, name in intersection_map.items():
            patch = v.get_patch_by_id(code)
            if patch:
                patch.set_color(color_palette[name])
                patch.set_alpha(0.65)

    # 6. Scatter Plot with Dynamic Markers & High Saturation
    ax_main = fig.add_subplot(gs[1:4, 1:4])
    ax_x = fig.add_subplot(gs[0, 1:4], sharex=ax_main)
    ax_y = fig.add_subplot(gs[1:4, 4], sharey=ax_main)

    sns.scatterplot(
        data=plot_df,
        x=x_var,
        y=y_var,
        hue="Intersection",
        style="Intersection",  # Applies marker_map shapes
        markers=marker_map,
        palette=color_palette,
        s=85,  # Slightly larger for clear marker visibility
        alpha=0.95,
        edgecolor="black",
        linewidth=0.5,
        ax=ax_main,
        hue_order=list(intersection_map.values()),
    )

    if plot_KDE:
        sns.kdeplot(
            data=plot_df,
            x=x_var,
            hue="Intersection",
            palette=color_palette,
            ax=ax_x,
            fill=True,
            alpha=0.35,
            legend=False,
        )
        sns.kdeplot(
            data=plot_df,
            y=y_var,
            hue="Intersection",
            palette=color_palette,
            ax=ax_y,
            fill=True,
            alpha=0.35,
            legend=False,
        )

    # Formatting
    ax_x.axis("off")
    ax_y.axis("off")
    ax_main.grid(True, linestyle=":", alpha=0.6)
    ax_main.legend(
        loc="lower right", fontsize="small", frameon=True, facecolor="white"
    )

    plt.show()
    
def old_plot_3way_bright_comparison(df1, df2, df3, on_cat, x_var, y_var, labels=('set1', 'set2', 'set3'), plot_KDE=False):
    # 1. Define Intersection Map
    intersection_map = {
        '100': f'Unique to {labels[0]}', '010': f'Unique to {labels[1]}', '001': f'Unique to {labels[2]}',
        '110': f'{labels[0]} & {labels[1]}', '101': f'{labels[0]} & {labels[2]}', '011': f'{labels[1]} & {labels[2]}',
        '111': 'Common to All'
    }

    # 2. Get Bright Standard Colors (Set2 is excellent for this)
    # We grab 7 colors from the 'Set2' colormap
    bright_colors = plt.cm.Set2.colors 
    color_palette = {name: bright_colors[i] for i, name in enumerate(intersection_map.values())}

    # 3. Extract Sets and Build Plotting DataFrame
    s1, s2, s3 = [set(df[on_cat].dropna().unique()) for df in [df1, df2, df3]]
    all_ids = s1 | s2 | s3
    
    coords = pd.concat([df[[on_cat, x_var, y_var]] for df in [df1, df2, df3]]).drop_duplicates(on_cat)
    plot_df = coords[coords[on_cat].isin(all_ids)].copy()
    plot_df['Intersection'] = plot_df[on_cat].apply(
        lambda v: intersection_map.get(f"{'1' if v in s1 else '0'}{'1' if v in s2 else '0'}{'1' if v in s3 else '0'}", 'NA')
    )

    # 4. Layout
    fig = plt.figure(figsize=(22, 11))
    gs = gridspec.GridSpec(4, 5, hspace=0.4, wspace=0.3, width_ratios=[2, 1, 1, 1, 0.5])

    # 5. Styled Venn Diagram
    ax_v = fig.add_subplot(gs[0:4, 0])
    v = venn3([s1, s2, s3], set_labels=labels, ax=ax_v)
    
    if v:
        # Scale up the font sizes
        for label in v.subset_labels:
            if label: label.set_fontsize(16)  # Numbers inside circles
        for label in v.set_labels:
            if label: label.set_fontsize(18)  # A, B, C labels

        # Apply the bright palette to patches
        for code, name in intersection_map.items():
            patch = v.get_patch_by_id(code)
            if patch:
                patch.set_color(color_palette[name])
                patch.set_alpha(0.7)

    # 6. Scatter and Marginals
    ax_main = fig.add_subplot(gs[1:4, 1:4])
    ax_x = fig.add_subplot(gs[0, 1:4], sharex=ax_main)
    ax_y = fig.add_subplot(gs[1:4, 4], sharey=ax_main)

    sns.scatterplot(
        data=plot_df, x=x_var, y=y_var, hue='Intersection', 
        palette=color_palette, s=70, alpha=0.8, ax=ax_main,
        hue_order=list(intersection_map.values())
    )

    # Add marginal KDEs with matching bright colors
    # if one set is much larger than the others, distribution plot is uninformative
    if plot_KDE:
        sns.kdeplot(data=plot_df, x=x_var, hue='Intersection', palette=color_palette, ax=ax_x, fill=True, legend=False)
        sns.kdeplot(data=plot_df, y=y_var, hue='Intersection', palette=color_palette, ax=ax_y, fill=True, legend=False)

    # Formatting
    ax_x.axis('off')
    ax_y.axis('off')
    ax_main.grid(True, linestyle=':', alpha=0.6)
    ax_main.legend(loc='lower right', fontsize='small', frameon=True, facecolor='white')

    plt.show()

def compare_sets(*dfs, labels=None, on_column='sbjct'):
    """
    Compare 2 or 3 DataFrames using Venn diagrams and scatter plots.
    
    Parameters:
    -----------
    *dfs : DataFrames
        Pass 2 or 3 pandas DataFrames as positional arguments.
    labels : list of str, optional
        Names for the sets (e.g., ['Set A', 'Set B']). 
        Defaults to ['DF1', 'DF2', 'DF3'].
    on_column : str
        The column used to identify common members (default: 'sbjct').
    """
    num_dfs = len(dfs)
    for x in dfs: 
        if x.empty:
            print(f"❌ Error: You provided an empty DataFrame.")
            return
    
    # 1. Enforce the "2 or 3" rule
    if num_dfs < 2 or num_dfs > 3:
        print(f"❌ Error: compare_sets() needs 2 or 3 DataFrames. You provided {num_dfs}.")
        return

    # 2. Handle default labels
    if labels is None:
        labels = [f"Set {i+1}" for i in range(num_dfs)]
    elif len(labels) != num_dfs:
        print(f"⚠️ Warning: Provided {len(labels)} labels for {num_dfs} sets. Using defaults.")
        labels = [f"Set {i+1}" for i in range(num_dfs)]

    # 3. Route to the appropriate engine
    if num_dfs == 2:
        print(f"📊 Generating 2-way comparison: {labels[0]} vs {labels[1]}")
        return plot_safe_comparison_with_marginal(dfs[0], dfs[1], on_category=on_column, labels=labels, x_var='z_score', y_var='ali-length')
    
    elif num_dfs == 3:
        print(f"📊 Generating 3-way comparison: {', '.join(labels)}")
        return plot_3way_bright_comparison(dfs[0], dfs[1], dfs[2], on_cat=on_column, labels=labels, x_var='z_score', y_var='ali-length')

##################################
# Sequence Analysis & Heatmaps Class
##################################

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

class msa:
    def __init__(self, df, pileup_col, color_mapper=None, pseudocount=0.1):
        self.df = df
        self.pileup_col = pileup_col
        self.color_mapper = color_mapper
        self.pseudocount = pseudocount
        
        # PERSISTENT OBJECTS
        self._info_df = None 
        self._logo_obj = None
        self._fig = None

        # 1. Grab the specific config for this column
        # Falls back to an empty dict with an empty 'rgb' map if not found
        self.config = COLOR_DATA.get(pileup_col, {"rgb": {}})
        self._rgb_array = None

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
        self._fig, ax = plt.subplots(figsize=(12, 3))
        colormap = COLOR_DATA[self.pileup_col]["names"]
        
        self._logo_obj = lm.Logo(self._info_df, ax=ax, color_scheme=colormap, vpad=.1, width=.8)
        self._logo_obj.ax.set_ylabel('bits')
        #self._logo_obj.ax.set_ylim([0, 7.0])
        
        # Close the figure immediately so it doesn't display yet
        plt.close(self._fig)

    def logo(self, left=-1, right=None, positions=[]):
        """This will now be instantaneous."""
        if self._logo_obj is None:
            self._initialize_everything()

        # Update the viewport of the EXISTING axes
        if right is None: 
            right = len(self._info_df)
        self._logo_obj.ax.set_xlim(left, right)

        # Use IPython display to show the pre-rendered figure
        display(self._fig)
        return None
           
    def heatmap(self, left=None, right=None, bottom=None, top=None):
        if self._rgb_array is None:
            # Access self.config["rgb"] as required
            rgb_lookup = self.config.get("rgb", {})
            print("Generating RGB Heatmap array...")
            matrix = self.df[self.pileup_col].astype(str).tolist()
            
            self._rgb_array = np.array([
                [rgb_lookup.get(char, (255, 255, 255)) for char in row]
                for row in matrix
            ], dtype=np.uint8)
 
        plt.figure(figsize=(5, 3), dpi=300)
        plt.imshow(self._rgb_array, aspect="auto")
        if left or right: plt.xlim(left=left, right=right)
        if bottom or top: plt.ylim(bottom=bottom, top=top)
        plt.show()

def _get_cached_msa(df, mode):
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

# --- User Facing Functions ---

def plot_dssp(df):
    """Start here to analyze secondary structure and recurrent domains."""
    return _get_cached_msa(df, 'dssp_pileup')

def plot_sequ(df):
    """Start here to analyze sequence conservation and motifs."""
    return _get_cached_msa(df, 'sequ_pileup')

##################################
# Workflows
##################################

def global_plots(df):
    print("""
    =========== BASIC PLOTS ON UNFILTERED DATASET
    """)
    # plot selected variables
    print("""
    ----------- GLOBAL DISTRIBUTION OF SELECTED VARIABLES

    DALI ranks results by Z-score. Alignment length can be used to separate good 
    hits from noise at low Z-scores. Pay attention to the maximum value on the x-axis.
    
    Hits with < 20% sequence identity fall into the "Midnight Zone". At this level, 
    standard sequence searches fail to show a clear relationship, yet these proteins 
    may still be remote homologs with identical folds. Pay attention to the maximum
    value on the x-axis.

    Low target_coverage indicates a match to a mobile module - a conserved functional 
    domain that evolution "copy-pasted" into different multidomain structures. While 
    the shared module is highly similar, the overall proteins likely have different 
    architectures and distinct biological roles.
    """)
    plot_distributions(df, ['z_score', 'rmsd', 'sequence_identity', 'target_length', 'query_coverage', 'target_coverage'], cols=3, bins=25, color='seagreen')

    # scatter plot
    print("""
    ----------- INTERACTIVE SCATTERPLOT
    
    DALI results are annotated with Pfam clans / families. The closest neighbors of the Query 
    are at the top right. 

    Use the interactive legend to toggle categories on and off.
    """)

    print("""
    ----------- VIOLIN PLOTS
    """)
    print("""
    ----------- SEQUENCE LOGO
    
    Pay attention to scale of Information (bits). The theoretical maximum is > 6 bits 
    for the rarest amino acids (W, C) and ~4 bits for common amino acids. Note that 
    the logo only displays amino acids with frequencies exceeding the database 
    background, highlighting true conservation over random noise.

    In highly diverse datasets, global conservation often looks 'blurry'. Downstream
    filtering for coherent structural similarity can reveal sharper, clade-specific 
    signals that are otherwise diluted in the global average. This can even uncover 
    residues essential for a specific sub-family's unique function.
    
    """)
    print("""
    ----------- STRUCTURAL COHERENCE
     
    Stacked structural alignment (helix / blue, strand / red, loop / green). Insertions 
    relative to the Query are omitted. Rows are reordered to maximize similarity between 
    adjacent rows.

    In large, multidomain proteins, "mobile modules" appear as rectangular blocks in the 
    stacked alignment, marking where different proteins share the same 3D fold despite 
    having different overall architectures.
    """)

def plot_overlaid_distributions(series_list, bins=25, alpha=0.5, figsize=(6,3), title="Overlaid Distributions"):
    """
    Plots overlaid histograms for a list of pandas Series in a single plot.
    
    Parameters:
    -----------
    series_list : list of pd.Series
        A list of pandas Series to plot (e.g., [df1['rmsd'], df2['rmsd']]).
    bins : int, default 25
        Number of histogram bins.
    alpha : float, default 0.5
        Transparency level so overlaps are visible.
    figsize : tuple, default (10, 6)
        Dimensions of the output plot.
    title : str, default "Overlaid Distributions"
        The title of the plot.
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Iterate directly over the list of Series
    for i, series in enumerate(series_list):
        # Drop NaN values to avoid breaking the histogram binning pipeline
        data = series.dropna()
        
        # Pull the series name for the legend; fallback to index if name is missing
        label = series.name if series.name is not None else f"Series {i + 1}"
        
        ax.hist(
            data, 
            bins=bins, 
            alpha=alpha, 
            label=label,
            histtype='bar'
        )
    
    # Style the single visual canvas
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlabel("Value", fontsize=8)
    ax.set_ylabel("Frequency", fontsize=8)
    
    # Show the legend to identify each Series
    #ax.legend(loc='upper right', frameon=True)
    # Place the upper-left corner of the legend slightly to the right (1.05) of the top-right corner (1.0)
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), frameon=True)
    
    plt.tight_layout()
    plt.show()