import pandas as pd
import numpy as np
import re
import sys
import ast
import os
import difflib
from io import StringIO

# --- Configuration & Aliases ---

NAMESPACE_MAP = {
    "PDB": 1,
    "AFDB2": 2, "AFDB": 2, "ALPHA_FOLD": 2,
    "BFVD": 3,
    "VIRO3D": 4, "VIRO": 4, "VAD": 4, "TOXO": 4
}

def get_namespace_id(name):
    """Maps namespace strings (and aliases) to fixed integers."""
    if not isinstance(name, str) or name.strip() == "":
        return 0
    key = name.strip().upper()
    # Explicitly catch "QUERY" to avoid it slipping into the ledger
    if key == "QUERY":
        return 0
    return NAMESPACE_MAP.get(key, 0)

def safe_int(val):
    """Handles the '.' quirk and float-strings."""
    if val == '.' or val is None or (isinstance(val, float) and pd.isna(val)):
        return None
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None

# --- Parsing Engines ---

# ------------------------
# Column handling
# ------------------------

COLUMN_ALIASES = {
    "z_score": ["z_score", "zscore", "z-score", "z", 'z_score', 'Zscore'],
    "target_id": ["target_id", 'Target', 'Subject', "target", "sbjct", "Shjct", "subject"],
    "namespace": ["namespace", "database", "Database", 'division'],
    "qstarts": ["qstarts", "q_start", "q_starts"],
    "sstarts": ["sstarts", "s_start", "s_starts"],
    "lengths": ["lengths", "length", 'segment_lengths'],
}


def normalize(name):
    return name.lower().replace("-", "_").replace(" ", "_")


def resolve_col(df, canonical):
    aliases = COLUMN_ALIASES.get(canonical, [canonical])
    aliases = [normalize(a) for a in aliases]

    cols = list(df.columns)

    # exact match
    for col in cols:
        if col in aliases:
            return col

    # fuzzy fallback
    match = difflib.get_close_matches(canonical, cols, n=1, cutoff=0.8)
    return match[0] if match else None


# ------------------------
# Fast parser
# ------------------------

def parse_modern_tsv(file_path):
    with open(file_path, "r") as f:
        cleaned_text = "".join(
            line for line in f
            if not line.lstrip().startswith("#")
        )

    df = pd.read_csv(StringIO(cleaned_text), sep="\t")

    # --- normalize columns ---
    df.columns = [normalize(c) for c in df.columns]

    # --- resolve columns ---
    c_t = resolve_col(df, "target_id")
    c_n = resolve_col(df, "namespace")
    c_qs = resolve_col(df, "qstarts")
    c_ss = resolve_col(df, "sstarts")
    c_ln = resolve_col(df, "lengths")

    # --- filter valid targets (vectorized) ---
    target = df[c_t].astype(str).str.strip()
    mask = (target != ".") & (target != "") & (target.str.lower() != "query")
    df = df.loc[mask].copy()

    # --- assign alignment_id ---
    df["alignment_id"] = np.arange(1, len(df) + 1)

    # --- summary table ---
    summary = pd.DataFrame({
        "alignment_id": df["alignment_id"],
        "target_id": target[mask].values,
        "z_score": pd.to_numeric(df["z_score"], errors="coerce"),
        "namespace_raw": df[c_n].astype(str) if c_n else None,
    })

    # --- segments (vectorized-ish) ---
    if all([c_qs, c_ss, c_ln]):

        # parse lists once
        qs = df[c_qs].map(_safe_literal)
        ss = df[c_ss].map(_safe_literal)
        ln = df[c_ln].map(_safe_literal)

        seg_df = pd.DataFrame({
            "alignment_id": df["alignment_id"],
            "q": qs,
            "s": ss,
            "l": ln,
        })

        # explode all three together
        seg_df = seg_df.explode(["q", "s", "l"], ignore_index=True)

        # convert to numeric
        seg_df["q"] = pd.to_numeric(seg_df["q"], errors="coerce")
        seg_df["s"] = pd.to_numeric(seg_df["s"], errors="coerce")
        seg_df["l"] = pd.to_numeric(seg_df["l"], errors="coerce").fillna(0)

        # drop invalid
        seg_df = seg_df.dropna(subset=["q", "s"])

        segments = pd.DataFrame({
            "alignment_id": seg_df["alignment_id"].astype(np.int32),
            "q_start": seg_df["q"].astype(np.int32),
            "s_start": seg_df["s"].astype(np.int32),
            "length": seg_df["l"].astype(np.int32),
        })

    else:
        segments = pd.DataFrame(columns=["alignment_id", "q_start", "s_start", "length"])

    return summary, segments


# ------------------------
# helpers
# ------------------------

def _safe_literal(x):
    try:
        s = str(x).strip()
        if not s or s in ('[]', 'nan', 'None'):
            return []
        # try ast first for properly formatted lists
        if ',' in s:
            return ast.literal_eval(s)
        # fallback: extract all integers from the string
        nums = re.findall(r'-?\d+', s)
        return [int(n) for n in nums] if nums else []
    except Exception:
        return []

# [Logic for parse_legacy_txt remains same as previous unified script]
def parse_legacy_txt(file_path):
    """Parses the DALI text file into Summary and Segments lists."""
    summary_data = []
    segments_data = []
    current_section = None
    fields = {"idx": (0, 5), "qs": (21, 26), "qe": (27, 32), "ss": (36, 41), "se": (42, 47)}

    with open(file_path, 'r') as f:
        for line in f:
            clean_line = line.strip()
            if not clean_line:
                continue

            if line.startswith("# No:"):
                current_section = "SUMMARY"
                continue
            elif line.startswith("# Structural equivalences"):
                current_section = "EQUIVALENCES"
                continue
            elif line.startswith("-") or (current_section and line.startswith("#")):
                current_section = None
                continue

            if current_section == "SUMMARY":
                parts = line.split()
                if len(parts) < 3: continue
                try:
                    summary_data.append({
                        'alignment_id': int(re.sub(r'\D', '', parts[0])),
                        'target_id': re.sub(r'[^A-Za-z0-9]', '', parts[1]),
                        'z_score': float(parts[2]),
                        'namespace_raw': None
                    })
                except ValueError:
                    continue

            elif current_section == "EQUIVALENCES":
                if True: #try:
                    d = {k: line[a:b].strip().rstrip(':') for k, (a, b) in fields.items()}
                    segments_data.append({
                        'alignment_id': int(d['idx']),
                        'q_start': int(d['qs']), 's_start': int(d['ss']),
                        'length': int(d['qe']) - int(d['qs']) + 1
                    })
                #except: continue

    return pd.DataFrame(summary_data), pd.DataFrame(segments_data)

def get_namespace_id(name):
    """Maps namespace strings or numeric strings to fixed integers."""
    if name is None or str(name).strip() == "" or str(name).upper() == "QUERY":
        return 0

    val = str(name).strip().upper()

    # 1. Check if it's already a valid numeric string ("1", "2", etc.)
    if val in ["1", "2", "3", "4"]:
        return int(val)

    # 2. Check the map for strings ("PDB", "AFDB2", etc.)
    return NAMESPACE_MAP.get(val, 0)

def process_dali(input_file, namespace_override=None):
    ext = os.path.splitext(input_file)[1].lower()

    # Parse the raw data
    df_sum_raw, df_seg = parse_modern_tsv(input_file) if ext == '.tsv' else parse_legacy_txt(input_file)

    if df_sum_raw.empty:
        return print(f"No valid data found in {input_file}")

    # --- 1. Build the ID Mapping Table ---
    unique_targets = df_sum_raw[['target_id', 'namespace_raw']].drop_duplicates('target_id').copy()

    if namespace_override:
        # override everything (fast path)
        ns_series = pd.Series(namespace_override, index=unique_targets.index)
    else:
        ns_series = unique_targets['namespace_raw'].astype(str).str.strip().str.upper()

    # vectorized namespace mapping
    df_mapping = unique_targets.copy()

    # handle numeric namespaces directly
    is_numeric = ns_series.isin(["1", "2", "3", "4"])
    df_mapping["namespace_id"] = 0

    df_mapping.loc[is_numeric, "namespace_id"] = ns_series[is_numeric].astype(int)

    # map string namespaces
    df_mapping.loc[~is_numeric, "namespace_id"] = (
        ns_series[~is_numeric].map(NAMESPACE_MAP).fillna(0).astype(int)
    )

    # assign protein_id
    df_mapping.insert(0, 'protein_id', np.arange(1, len(df_mapping) + 1))

    # --- 2. Build and Export Summary ---
    # Merge to bring BOTH protein_id AND namespace_id into the summary
    df_sum = df_sum_raw.merge(
        df_mapping[['target_id', 'protein_id', 'namespace_id']], 
        on='target_id', 
        how='left'
    )

    # CRITICAL: Added 'namespace_id' to the export list
    summary_cols = ['alignment_id', 'protein_id', 'namespace_id', 'z_score']
    df_sum[summary_cols].to_csv("summary.tsv", sep="\t", index=False)

    # --- 3. Export Mappings and Segments ---
    df_mapping.to_csv("id_mapping.tsv", sep="\t", index=False)
    df_seg.to_csv("segments.tsv", sep="\t", index=False)

    print(f"Processed {len(df_mapping)} targets.")
    print(f"Summary saved with columns: {summary_cols}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"USAGE: {sys.argv[0]} <dali_output_txt_or_tsv> <PDB|AFDB2|BFVD|VIRO3D>")
        sys.exit(1)
    process_dali(sys.argv[1], namespace_override=sys.argv[2])
