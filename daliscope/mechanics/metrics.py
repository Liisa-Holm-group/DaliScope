# ------------------------------------------------------------------ #
# Pure function — outside the class                                  #
# ------------------------------------------------------------------ #

from typing import Optional
from numba import njit
import pandas as pd
import numpy as np
import re

from numba import njit
import daliscope.analysis.occupancy

def compute_family_presence(population_df: pd.DataFrame, group_series: pd.Series,
                             target_ids, id_col: str = "target_id") -> pd.DataFrame:
    """
    For each distinct value in group_series (aligned to population_df by
    index — typically the output of apply_rank_pooling), computes how many
    rows have id_col in target_ids ("present") vs not ("absent"), as counts
    and percentages.
 
    Returns a DataFrame indexed by group label, with columns:
        n_present, n_absent, n_total, pct_present, pct_absent
    pct_absent is stored as a NEGATIVE percentage, for direct use in a
    diverging bar chart.
    """
    target_id_set = set(target_ids)
    is_target = population_df[id_col].isin(target_id_set)

    counts = (
        population_df.groupby(group_series)[id_col].count().rename("n_total").to_frame()
    )
    present_counts = population_df.loc[is_target].groupby(group_series[is_target])[id_col].count()
    counts["n_present"] = present_counts.reindex(counts.index).fillna(0).astype(int)
    counts["n_absent"] = counts["n_total"] - counts["n_present"]

    counts["pct_present"] = (counts["n_present"] / counts["n_total"]) * 100
    counts["pct_absent"] = -(counts["n_absent"] / counts["n_total"]) * 100

    return counts


def get_domain_ranges(project, view_name: str) -> list:
    """
    Return domain_ranges for a view by tracing back through
    parent_objects until a view with domain_ranges is found.
    Returns None if no ancestor has domain_ranges.
    """
    visited = set()
    current = view_name

    while current is not None and current not in visited:
        visited.add(current)
        prov = project.provenance.get(current)
        if prov is None:
            break

        # found it
        if 'domain_ranges' in prov.parameters:
            return prov.parameters['domain_ranges']

        # climb to parent
        current = prov.parent_objects[0] if prov.parent_objects else None

    return None   # no ancestor has domain_ranges — e.g. full view


def add_pfam_clan_ranks(df: pd.DataFrame,
                        pfam_col: str = 'pfam',
                        clan_col: str = 'clan') -> pd.DataFrame:
    """
    Add pfam_rank and clan_rank columns to a view dataframe.

    Rank = order of first appearance when targets are sorted by
    z_score descending. The Pfam/clan seen in the highest-z_score
    target gets rank 1, the next new one gets rank 2, etc.
    Targets sharing the same Pfam/clan share the same rank.
    NaN stays NaN (displayed as 'Unassigned' at plot time).

    Called once in register_domain_view() so ranks are stable
    and consistent across all downstream plot calls on this view.
    """
    df = df.sort_values('z_score', ascending=False).copy()

    for col, rank_col in [(pfam_col, 'pfam_rank'), (clan_col, 'clan_rank')]:
        rank    = 1
        seen    = {}
        ranks   = []
        for val in df[col]:
            if pd.isna(val) or val == 'Unassigned':
                ranks.append(np.nan)
                continue
            if val not in seen:
                seen[val] = rank
                rank += 1
            ranks.append(seen[val])
        df[rank_col] = ranks

    return df

##############################################
# --- PERFORMANCE-OPTIMIZED SERIATION (DSSP) ---
##############################################
@njit
def encode_dssp_optimized(pileups, n_rows, n_pos):
    # Create a LUT for ASCII: E=69, H=72, L=76, others=0
    # Values represent the two bits for the bitmat
    lut = np.zeros(256, dtype=np.uint8)
    lut[69] = 1  # E -> 10 (binary)
    lut[72] = 2  # H -> 01 (binary)
    lut[76] = 3  # L -> 11 (binary)

    bitmat = np.zeros((n_rows, 2 * n_pos), dtype=np.uint8)

    for i in range(n_rows):
        for j in range(n_pos):
            val = lut[pileups[i, j]]
            if val == 1:      # E
                bitmat[i, 2*j] = 1
            elif val == 2:    # H
                bitmat[i, 2*j+1] = 1
            elif val == 3:    # L
                bitmat[i, 2*j] = 1
                bitmat[i, 2*j+1] = 1
    return bitmat

@njit
def nearest_neighbor_reordering_fast(bitmat):
    n_rows, n_bits = bitmat.shape

    # Calculate similarity matrix via dot product
    # (equivalent to popcount of bitwise AND for uint8 0/1 matrices)
    # We cast to float32 for the dot product and back to int32 for the scores
    sim = (bitmat.astype(np.float32) @ bitmat.T.astype(np.float32)).astype(np.int32)

    visited = np.zeros(n_rows, dtype=np.uint8)
    order = np.empty(n_rows, dtype=np.int32)

    order[0] = 0
    visited[0] = 1

    for step in range(1, n_rows):
        last = order[step - 1]
        best_score = -1
        best_idx = -1

        # We only need the row for 'last'
        current_similarities = sim[last]

        for j in range(n_rows):
            if not visited[j]:
                score = current_similarities[j]
                if score > best_score:
                    best_score = score
                    best_idx = j

        order[step] = best_idx
        visited[best_idx] = 1

    return order

def seriate(df, force=False):
    if not 'dssp_pileup' in df.columns:
        print('[-] dssp_pileup missing, seriation skipped')
        df['dssp-order'] = np.arange(len(df))
    if not force and 'dssp_order' in df.columns:
        print('[-] dssp-order column present, not overwritten')
        return df
    pileups = df["dssp_pileup"].to_list()
    n_rows = len(pileups)
    n_pos = len(pileups[0])

    # convert to fixed-width uint8 array for Numba
    pileup_arr = np.zeros((n_rows, n_pos), dtype=np.uint8)
    for i, row in enumerate(pileups):
        pileup_arr[i, :] = np.frombuffer(row.encode("ascii"), dtype=np.uint8)

    bitmat = encode_dssp_optimized(pileup_arr, n_rows, n_pos)
    order = nearest_neighbor_reordering_fast(bitmat)

    df = df.copy()
    df["dssp_order"] = order
    return df

##############################################
# Non-redundant subset generation
##############################################

def nonredundant_subset(pileup_series: pd.Series,
                        threshold: float = 0.9,
                        gap_char: str = '.') -> pd.Index:
    """
    Greedy non-redundant subset of gapped pileup strings.
    Accepts a sequence if its identity to ALL already-accepted
    sequences is below threshold.

    Parameters
    ----------
    pileup_series : pd.Series of gapped strings, all same length
    threshold     : maximum allowed pairwise identity (default 0.9)
    gap_char      : gap character (default '.')

    Returns
    -------
    pd.Index of accepted row indices (subset of pileup_series.index)

    Performance
    -----------
    Uses vectorized numpy — compares each candidate against all
    accepted sequences at once. O(N * n_accepted * L / 64) with
    bitwise operations. Fast for thousands of rows.
    """
    # drop NaN rows — can't compare them
    valid = pileup_series.dropna()
    if valid.empty:
        return pd.Index([])

    L = len(valid.iloc[0])

    # encode as uint8 ASCII array — (N, L)
    X = np.frombuffer(
        b''.join(s.encode('ascii') for s in valid),
        dtype=np.uint8
    ).reshape(len(valid), L)

    gap_code = ord(gap_char)

    # sort by decreasing occupancy (non-gap count) — best first
    occupancy = (X != gap_code).sum(axis=1)
    sort_order = np.argsort(occupancy)[::-1]

    X_sorted      = X[sort_order]
    original_idx  = valid.index[sort_order]

    # greedy acceptance
    accepted_mask = np.zeros(len(valid), dtype=bool)
    accepted_rows = []   # list of accepted X rows for vectorized comparison

    for i in range(len(X_sorted)):
        candidate = X_sorted[i]   # (L,)

        if not accepted_rows:
            # first candidate always accepted
            accepted_mask[i] = True
            accepted_rows.append(candidate)
            continue

        # stack accepted sequences: (n_accepted, L)
        A = np.array(accepted_rows, dtype=np.uint8)

        # positions where BOTH candidate and accepted have a residue
        cand_present = candidate != gap_code          # (L,)
        acc_present  = A != gap_code                  # (n_accepted, L)
        both_present = cand_present[None, :] & acc_present  # (n_accepted, L)

        matches = ((A == candidate[None, :]) & both_present).sum(axis=1)  # (n_accepted,)

        # len_shorter = min(n_residues_in_candidate, n_residues_in_accepted)
        n_candidate = cand_present.sum()                        # scalar
        n_accepted_lens = acc_present.sum(axis=1)               # (n_accepted,)
        len_shorter = np.minimum(n_candidate, n_accepted_lens)  # (n_accepted,)

        with np.errstate(invalid='ignore', divide='ignore'):
            identity = np.where(len_shorter > 0,
                                matches / len_shorter,
                                0.0)

        if identity.max() < threshold:
            accepted_mask[i] = True
            accepted_rows.append(candidate)

    return original_idx[accepted_mask]

def filter_nonredundant_sequences(data_frame, threshold=0.40, verbose=False):
    """
    Wrapper: Filters the dataframe for non-redundant sequences based on a threshold.
    """
    # 1. Get the indices for the non-redundant subset
    nr_indices = nonredundant_subset(
        data_frame['sequ_pileup'], 
        threshold=threshold
    )

    # 2. Slice the dataframe
    filtered_df = data_frame.loc[nr_indices]

    # 3. Log the results
    if verbose: print(f"Reduced {len(data_frame):,} → {len(filtered_df):,} sequences at {threshold} identity")

    return filtered_df

##############################################

def architecture_frequency_table(df: pd.DataFrame,
                                 pfam_col: str = 'pfam_domains'
                                 ) -> pd.DataFrame:
    """
    One row per distinct domain architecture (ordered tuple of
    domain IDs), with count and representative target_id, sorted
    by frequency descending.

    Architecture labels are concise ID-only tuples:
        'pfam_domains' → "PF00001/CL0001, PF00002/CL0002"
        'clan_domains' → "CL0001, CL0002"
    """
    plot_df = df.dropna(subset=[pfam_col]).copy()

    def architecture_label(doms) -> str:
        if not doms:
            return 'no domains annotated'

        parts = []
        for d in doms:
            pid = d['pfam_id']
            if pfam_col == 'pfam_domains':
                clan = d.get('clan')
                parts.append(f"{pid}/{clan}" if clan else pid)
            else:
                parts.append(pid)
        return ', '.join(parts)

    plot_df['architecture'] = plot_df[pfam_col].apply(architecture_label)

    summary = (
        plot_df.groupby('architecture')
        .agg(n_targets=('target_id', 'count'),
             example=('target_id', 'first'),
             mean_z_score=('z_score', 'mean'))
        .sort_values('n_targets', ascending=False)
        .reset_index()
    )
    return summary

def build_pfam_domains(pfam_hits: pd.DataFrame,
                       pfam_names: pd.DataFrame = None,
                       level: str = 'pfam',
                       e_value_threshold: float = 1e-5
                       ) -> pd.DataFrame:
    """
    Convert flat Pfam hit table into ordered domain-architecture
    lists per target — at either Pfam-family or clan resolution.

    Parameters
    ----------
    pfam_hits : raw project._pfam dataframe with columns
                target_name, query_accession, e_value,
                env_from, env_to, database, clan
    pfam_names : project._pfam_names dataframe with columns
                 pfam, clan, clan_short, short, name
    level      : 'pfam' (default) — group by individual Pfam family.
                 Each domain dict carries 'pfam_id' (the accession)
                 and 'name' (from pfam_names.name).
                 'clan' — group by clan instead, merging adjacent/
                          overlapping hits from the same clan into
                          one architecture entry (since a clan-level
                          view shouldn't show two falsely-distinct
                          blocks for two Pfam IDs that are really
                          the same evolutionary unit). Each domain
                          dict carries 'pfam_id' (the clan id) and
                          'name' (from pfam_names.clan_short).
    e_value_threshold : drop weak/spurious hits before building
                        architecture

    Returns
    -------
    pd.DataFrame with columns: target_id, <output_col>
        where output_col is 'pfam_domains' (level='pfam')
        or 'clan_domains' (level='clan')
    """
    if level not in ('pfam', 'clan'):
        raise ValueError(f"level must be 'pfam' or 'clan', got {level!r}")

    output_col = 'pfam_domains' if level == 'pfam' else 'clan_domains'

    # --- lookups from pfam_names ---
    name_map = dict(zip(pfam_names['pfam'], pfam_names['name'])) \
               if pfam_names is not None else {}
    clan_map = dict(zip(pfam_names['pfam'], pfam_names['clan'])) \
               if pfam_names is not None else {}
    clan_short_map = (
        pfam_names.drop_duplicates('clan')
        .set_index('clan')['clan_short'].to_dict()
        if pfam_names is not None and 'clan' in pfam_names.columns
        else {}
    )

    filtered = pfam_hits[pfam_hits['e_value'] <= e_value_threshold].copy()

    if level == 'clan':
        # attach clan id to each hit (prefer the hit's own 'clan'
        # column if present and non-null, else look up via pfam_names)
        filtered['group_id'] = filtered.apply(
            lambda r: r['clan'] if pd.notna(r.get('clan'))
                      else clan_map.get(r['query_accession'],
                                       r['query_accession']),
            axis=1
        )
    else:
        filtered['group_id'] = filtered['query_accession']

    filtered = filtered.sort_values(['target_name', 'env_from'])

    records = []
    for target_id, group in filtered.groupby('target_name'):

        if level == 'pfam':
            # pfam/name pairs
            domains = [
                {
                    'pfam_id': row['query_accession'],
                    'name':    name_map.get(row['query_accession'],
                                            row['query_accession']),
                    'start':   int(row['env_from']),
                    'end':     int(row['env_to']),
                    'clan':    row.get('clan'),
                    'e_value': float(row['e_value']),
                }
                for _, row in group.iterrows()
            ]
        else:
            # clan/clan_short pairs — merge adjacent/overlapping hits
            # sharing the same group_id into a single block, taking
            # the envelope min(start)/max(end) and best (lowest) e-value
            domains = []
            for clan_id, clan_group in group.groupby('group_id'):
                domains.append({
                    'pfam_id': clan_id,
                    'name':    clan_short_map.get(clan_id, clan_id),
                    'start':   int(clan_group['env_from'].min()),
                    'end':     int(clan_group['env_to'].max()),
                    'clan':    clan_id,
                    'e_value': float(clan_group['e_value'].min()),
                })
            domains.sort(key=lambda d: d['start'])

        records.append({'target_id': target_id, output_col: domains})

    return pd.DataFrame(records)


def compute_mask(x, y, x_sect, y_sect, angle_deg):
    theta = np.deg2rad(angle_deg)
    nx = -np.sin(theta)
    ny =  np.cos(theta)
    return (x - x_sect) * nx + (y - y_sect) * ny > 0

def interpret_hits(df, query_name="query"):
    """
    Inspects top DALI hits and writes interpretations to a new column.
    Rules:
      (a) same family:   sequence_identity > 0.4
      (b) same family:   sequence_identity > 0.2 AND  z_score > 12
      (c) same fold:     sequence_identity <= 0.4 AND z_score > 8
      (d) marginal hit:  z_score between 2 and 8
      (e) no hit:        z_score <= 2

    Twilight zone:
      >30% sequence identity: sequence already tells much of the story.
      15–30%: genuine twilight zone where structural evidence adds substantial value.
      <15%: rely primarily on structural alignment quality and conserved functional features rather than sequence identity alone.
    """
    def classify(row):
        z     = row['z_score']
        seqid = row['sequence_identity']
        lali  = row['alignment_length']
        qcov  = row['query_coverage']
        tid   = row['target_id']

        if seqid > 0.3 and (lali > 100 or qcov > 0.5):
            return f"{query_name} likely in same family (high seq_id)"
        elif seqid > 0.15 and z > 12:
            return f"{query_name} likely in same superfamily (twilight seq_id, high Z)"
        elif z > 8:
            return f"{query_name} likely in same fold (strong Z)"
        elif seqid > 0.15 and Z > 4:
            return f"{query_name} remote similarity (twilight seq_id, weak Z)"
        elif z > 2:
            return f"{query_name} marginal similarity (midnight seq_id, low Z)"
        else:
            return f"no significant similarity"

    cols = ['target_id', 'z_score', 'alignment_length', 'rmsd', 'query_coverage',
            'sequence_identity', 'description']

    mini_df = df[cols].head().copy()
    mini_df['rough call'] = mini_df.apply(classify, axis=1)

    return mini_df

def compute_clipped_pileups(
    clipped_df: pd.DataFrame,
    master_metadata: pd.DataFrame,
    target_seqs: np.ndarray,
    target_dssp: Optional[np.ndarray] = None,
    query_length: int = 0,
    domain_ranges: Optional[list] = None,
) -> pd.DataFrame:
    """
    Re-builds sequence and DSSP pileups aligned to the query, restricted
    to the active clipped domain ranges.
    """
    from daliscope.core.project import build_pileups  # or wherever build_pileups resides
    # Determine effective query length for the clipped domain
    if domain_ranges:
        effective_length = sum(end - start + 1 for start, end in domain_ranges)
    else:
        effective_length = query_length

    # Fall back directly to build_pileups using the clipped segments
    pileup_df = build_pileups(
        segments        = clipped_df,
        master_metadata = master_metadata,
        target_seqs     = target_seqs,
        target_dssp     = target_dssp,
        query_length    = effective_length,
    )
    return pileup_df

import pandas as pd
import numpy as np

def clip_segments_to_domain_ranges(segments_df: pd.DataFrame, domain_ranges: list) -> pd.DataFrame:
    """
    Clips alignment segments so that only positions falling inside 
    the Query domain_ranges are kept. Preserves target residue alignment.
    """
    clipped_rows = []

    for _, seg in segments_df.iterrows():
        q_start = seg['q_start']
        q_end = q_start + seg['length']  # half-open
        s_start = seg['s_start']

        for d_start, d_end in domain_ranges:
            # Find overlap between segment [q_start, q_end) and domain [d_start, d_end)
            overlap_start = max(q_start, d_start)
            overlap_end = min(q_end, d_end)

            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                q_offset = overlap_start - q_start

                new_seg = seg.copy()
                new_seg['q_start'] = overlap_start
                new_seg['s_start'] = s_start + q_offset
                new_seg['length'] = overlap_len
                clipped_rows.append(new_seg)

    if not clipped_rows:
        return pd.DataFrame(columns=segments_df.columns)

    return pd.DataFrame(clipped_rows)


def register_domain_view(
    project,
    domain_ranges: list,
    view_name: str
) -> pd.DataFrame:
    """
    Clips alignment segments to query domain ranges, re-evaluates Pfam/Clan 
    overlap strictly for the clipped region, and registers the view.
    """
    if view_name in project.views:
        print(f"View '{view_name}' already exists — choose a different name or delete it first.")
        return

    # 1. Clip segments strictly against Query coordinates
    clipped_df = clip_segments_to_domain_ranges(project.segments, domain_ranges)
    if clipped_df.empty:
        print("Empty dataframe after clipping — view not registered.")
        return

    # Effective query length across domain ranges
    domain_query_length = sum(end - start for start, end in domain_ranges)

    # 2. Compute structural and sequence metrics on clipped segments
    master_stats = compute_metrics(
        clipped_segs    = clipped_df,
        master_metadata = project.master_metadata,
        query_ascii     = project.query_ascii,
        target_seqs     = project.target_seqs,
        query_length    = project.query_length,  # Keep full query coordinate space
        query_ca_coords = project.query_ca_coords,
        target_ca_coords= project.coords,
    )

    # 3. Merge remaining silver columns from master_metadata
    already_present = set(master_stats.columns)
    extra_cols = [
        c for c in project.master_metadata.columns
        if c not in already_present or c == "alignment_id"
    ]
    master_stats = master_stats.merge(
        project.master_metadata[extra_cols],
        on="alignment_id",
        how="inner",
    )

    # 4. Re-evaluate Pfam and Clan for the CLIPPED domain segments
    if getattr(project, "_pfam", None) is not None:
        from ..mechanics.metrics import best_pfam_per_alignment

        best_pfam_df = best_pfam_per_alignment(
            pfam_df         = project._pfam,
            segments_df     = clipped_df,
            master_metadata = project.master_metadata,
        )

        # Drop old parent pfam/clan assignments
        master_stats = master_stats.drop(columns=[c for c in ["pfam", "clan"] if c in master_stats.columns])
        master_stats = master_stats.merge(
            best_pfam_df[["alignment_id", "pfam", "clan"]],
            on="alignment_id",
            how="left"
        )
        master_stats["pfam"] = master_stats["pfam"].fillna("Unassigned")
        master_stats["clan"] = master_stats["clan"].fillna("Unassigned")
    else:
        master_stats["pfam"] = "Unassigned"
        master_stats["clan"] = "Unassigned"

    # 5. Re-generate pileups anchored to the full query length
    pileup_df = build_pileups(
        segments        = clipped_df,
        master_metadata = project.master_metadata,
        target_seqs     = project.target_seqs,
        target_dssp     = getattr(project, 'target_dssp', None),
        query_length    = project.query_length,  # Preserves correct Query column alignment
    )

    for col in ["sequ_pileup", "dssp_pileup"]:
        if col in master_stats.columns:
            master_stats = master_stats.drop(columns=[col])

    master_stats = master_stats.merge(pileup_df, on="alignment_id", how="left")
    master_stats = add_pfam_clan_ranks(master_stats)

    master_stats = master_stats.assign(
        rmsd=master_stats['rmsd_new'],
        R=master_stats['rotation'],
        t=master_stats['translation']
    ).drop(columns=['rmsd_new', 'rotation', 'translation'])

    # 6. Register view
    project.add_view(
        view_name,
        master_stats,
        function="register_domain_view",
        parameters={"domain_ranges": domain_ranges}
    )

    print(f"✓ Registered '{view_name}': {len(master_stats)} alignments")
    return master_stats


def old_register_domain_view(
    project,
    domain_ranges: list,
    view_name: str
) -> pd.DataFrame:
    """
    Clips segments, computes metrics, merges ALL silver columns
    from master_metadata, registers view.
    """
    # Check if view exists already
    if view_name in project.views:
        print(
            f"View '{view_name}' already exists - "
            f"choose a different name or delete it first: "
            f"del project.views['{view_name}']"
        )
        return

    # 1. clip segments using non-contiguous domain ranges
    clipped_df = clip_segments(
        project.segments,
        domain_ranges=domain_ranges,
    )
    if clipped_df.empty:
        print("Empty dataframe after clipping - not registered")
        return

    # Compute effective total query length across all ranges in the domain
    domain_query_length = sum(end - start + 1 for start, end in domain_ranges)

    # 2. compute metrics — passing domain length and coordinate arrays
    master_stats = compute_metrics(
        clipped_segs    = clipped_df,
        master_metadata = project.master_metadata,
        query_ascii     = project.query_ascii,
        target_seqs     = project.target_seqs,
        query_length    = domain_query_length,
        query_ca_coords = project.query_ca_coords,
        target_ca_coords= project.coords, # or project.target_ca_coords
    )

    # 3. merge ALL remaining columns from master_metadata
    already_present = set(master_stats.columns)

    extra_cols = [
        c for c in project.master_metadata.columns
        if c not in already_present or c == "alignment_id"
    ]

    master_stats = master_stats.merge(
        project.master_metadata[extra_cols],
        on="alignment_id",
        how="inner",
    )

    # 4. recompute pileups for the clipped domain
    pileup_df = compute_clipped_pileups(
        clipped_df      = clipped_df,
        master_metadata = project.master_metadata,
        target_seqs     = project.target_seqs,
        target_dssp     = getattr(project, 'target_dssp', None),
        query_length    = project.query_length,
        domain_ranges   = domain_ranges,
    )

    for col in ["sequ_pileup", "dssp_pileup"]:
        if col in master_stats.columns:
            master_stats = master_stats.drop(columns=[col])

    master_stats = master_stats.merge(
        pileup_df, on="alignment_id", how="left"
    )

    master_stats = add_pfam_clan_ranks(master_stats)

    master_stats = master_stats.assign(
        rmsd=master_stats['rmsd_new'],
        R=master_stats['rotation'],
        t=master_stats['translation']
    ).drop(columns=['rmsd_new', 'rotation', 'translation'])

    # 5. register
    project.add_view(
        view_name,
        master_stats,
        function="register_domain_view",
        parameters={"domain_ranges": domain_ranges}
    )

    print(f"✓ Registered '{view_name}': {len(master_stats)} alignments")

def register_multiple_domains(domain_string, project, custom_name: str = ""):
    """Parses domain_string from viewer, creates domain views, and registers them to project."""

    dom_strings = [
        d.strip() for d in domain_string.split(",") if d.strip()
    ]

    for i, domstr in enumerate(dom_strings):

        domain_ranges = [
            tuple(map(int, part.split("-"))) for part in domstr.split("_")
        ]
        default_name = "dom_" + "_".join(f"{s}-{e}" for s, e in domain_ranges)

        # Logical resolution for custom naming:
        # If single domain: use custom_name directly (if provided)
        # If multiple domains: append index if custom_name provided to avoid collision
        if custom_name:
            dom_name = (
                f"{custom_name}_{i+1}"
                if len(dom_strings) > 1
                else custom_name
            )
        else:
            dom_name = default_name

        register_domain_view(
            project=project,
            domain_ranges=domain_ranges,
            view_name=dom_name,
        )

def build_pileups(segments: pd.DataFrame,
                  master_metadata: pd.DataFrame,
                  target_seqs: np.ndarray,
                  target_dssp: np.ndarray,
                  query_length: int) -> pd.DataFrame:
    """
    Build per-alignment sequence and DSSP pileup strings.
    Each string has length query_length:
        '.' at query positions not covered by this alignment
        AA char / DSSP char where the target has an aligned residue

    Parameters
    ----------
    segments        : normalised (0-based) segments dataframe
    master_metadata : provides start_idx per alignment_id
    target_seqs     : flat (N_total,) int32 array of AA codes
    target_dssp     : flat (N_total,) int32 array of DSSP codes, or None
    query_length    : L

    Returns
    -------
    pd.DataFrame with columns: alignment_id, sequ_pileup[, dssp_pileup]
    """
    AA         = '-ACDEFGHIKLMNPQRSTVWY'
    DSSP_CODES = {ord(c): c for c in 'LHEGIBTSCL '}

    flat_segs = segments.merge(
        master_metadata[["alignment_id", "start_idx"]],
        on="alignment_id", how="inner"
    ).dropna(subset=["start_idx"])

    if flat_segs.empty:
        return pd.DataFrame(columns=["alignment_id", "sequ_pileup"])

    # --- build flat index arrays ---
    q_base           = flat_segs["q_start"].values.astype(np.int64) - 1
    s_base           = (flat_segs["start_idx"].astype(np.int64).values +
                        flat_segs["s_start"].values.astype(np.int64)) -1
    lengths_vector   = flat_segs["length"].values.astype(np.int64)
    aln_ids_flat     = flat_segs["alignment_id"].values

    repeats          = np.repeat(np.arange(len(lengths_vector)), lengths_vector)
    offsets          = np.concatenate([np.arange(l) for l in lengths_vector])

    global_q_indices = (q_base[repeats] + offsets).astype(np.int64)
    global_s_indices = (s_base[repeats] + offsets).astype(np.int64)

    valid = (
        (global_q_indices >= 0) & (global_q_indices < query_length) &
        (global_s_indices >= 0) & (global_s_indices < len(target_seqs))
    )
    global_q_indices       = global_q_indices[valid]
    global_s_indices       = global_s_indices[valid]
    alignment_ids_expanded = aln_ids_flat[repeats][valid]

    # --- per-alignment pileup dicts ---
    unique_aln_ids  = flat_segs["alignment_id"].unique()
    sequ_dict = {aln_id: ['.'] * query_length for aln_id in unique_aln_ids}
    dssp_dict = {aln_id: ['.'] * query_length for aln_id in unique_aln_ids} \
                if target_dssp is not None else None

    s_chars = target_seqs[global_s_indices]
    for q_pos, s_char, aln_id in zip(global_q_indices, s_chars,
                                      alignment_ids_expanded):
        sequ_dict[aln_id][q_pos] = (AA[s_char] if s_char < len(AA)
                                    else chr(int(s_char)))

    if target_dssp is not None:
        d_chars = target_dssp[global_s_indices]
        for q_pos, d_char, aln_id in zip(global_q_indices, d_chars,
                                          alignment_ids_expanded):
            dssp_dict[aln_id][q_pos] = DSSP_CODES.get(int(d_char), '.')

    # --- join to strings ---
    rows = [{"alignment_id": aln_id,
             "sequ_pileup":  ''.join(sequ_dict[aln_id]),
             **({"dssp_pileup": ''.join(dssp_dict[aln_id])}
                if dssp_dict else {})}
            for aln_id in unique_aln_ids]

    return pd.DataFrame(rows)

import numpy as np
import pandas as pd


def best_pfam_per_alignment(
    pfam_df: pd.DataFrame,
    segments_df: pd.DataFrame,
    master_metadata: pd.DataFrame,
    pfam_col: str = "query_accession",
    min_overlap_len=30,  # don't assign based on very short alignment
    min_coverage=0.5
    ) -> pd.DataFrame:
    """For each alignment_id, find the best Pfam hit that overlaps the aligned

    range of that specific alignment in the target. Alignments with no
    qualifying hits return 'Unassigned'.
    """
    # Quick guard rail: if segments is empty, nothing to assign
    if segments_df is None or segments_df.empty:
        return pd.DataFrame(columns=["alignment_id", "pfam", "clan"])

    # 1. compute aligned range per alignment_id from segments
    segs = segments_df.copy()
    segs["s_end"] = segs["s_start"] + segs["length"]
    aln_bounds = (
        segs.groupby("alignment_id")
        .agg(
            aln_start=("s_start", "min"),
            aln_end=("s_end", "max"),
        )
        .reset_index()
    )
    aln_bounds["aln_len"] = aln_bounds["aln_end"] - aln_bounds["aln_start"]

    # Keep a record of ALL distinct alignment_ids we intend to report on
    all_alignments = pd.DataFrame(
        {"alignment_id": aln_bounds["alignment_id"].unique()}
    )

    if pfam_df is None or pfam_df.empty:
        # Fast-track return all as Unassigned if Pfam data is completely missing
        all_alignments["pfam"] = "Unassigned"
        all_alignments["clan"] = "Unassigned"
        return all_alignments

    # 2. link alignment_id -> target_id via master_metadata
    aln_bounds = aln_bounds.merge(
        master_metadata[["alignment_id", "target_id"]],
        on="alignment_id",
        how="inner",
    )

    # 3. merge with pfam on target_id
    pfam = pfam_df.copy()
    pfam = pfam.rename(columns={"target_name": "target_id"})

    pfam["pf_start"] = pfam["env_from"] - 1
    pfam["pf_end"] = pfam["env_to"]
    pfam["pf_len"] = pfam["pf_end"] - pfam["pf_start"]

    df = pfam.merge(aln_bounds, on="target_id", how="inner")

    # 4. compute overlap and apply coverage constraints
    if not df.empty:
        overlap_start = np.maximum(df["pf_start"], df["aln_start"])
        overlap_end = np.minimum(df["pf_end"], df["aln_end"])
        df["overlap_len"] = (overlap_end - overlap_start).clip(lower=0)

        df["query_coverage"] = df["overlap_len"] / df["aln_len"]
        df["sbjct_coverage"] = df["overlap_len"] / df["pf_len"]

        # Filter for valid hits
        df = df[
            (df["overlap_len"] >= min_overlap_len)
            | (df["query_coverage"] > min_coverage)
            | (df["sbjct_coverage"] > min_coverage)
        ]

    # 5. Extract the best hit per alignment (if any valid hits exist)
    if not df.empty:
        df = df.sort_values(
            by=["overlap_len", "e_value"], ascending=[False, True]
        )
        best_hits = (
            df.groupby("alignment_id")
            .first()
            .reset_index()[["alignment_id", pfam_col, "clan"]]
            .rename(columns={pfam_col: "pfam"})
        )
    else:
        best_hits = pd.DataFrame(columns=["alignment_id", "pfam", "clan"])

    # 6. Map back onto ALL alignments to capture the missing ones as 'Unassigned'
    final_df = all_alignments.merge(best_hits, on="alignment_id", how="left")
    final_df["pfam"] = final_df["pfam"].fillna("Unassigned")
    final_df["clan"] = final_df["clan"].fillna("Unassigned")

    return final_df  # columns: alignment_id, pfam, clan

def best_pfam_clan_per_target( ##>> NOT USED?
    pfam_df: pd.DataFrame,
    id_col: str = "target_name",
    pfam_col: str = "query_accession",
    clan_col: str = "clan",
    score_col: Optional[str] = "score",
) -> pd.DataFrame:
    """
    Reduce a many-to-one Pfam dataframe to one best Pfam/clan per target.

    Parameters
    ----------
    pfam_df   : DataFrame with at least id_col, pfam_col, and clan_col.
    id_col    : Column identifying the protein (default 'target_name').
    pfam_col  : Column containing the Pfam accession.
    clan_col  : Column containing the associated clan.
    score_col : If present, choose the highest-scoring hit; otherwise choose
                the first hit encountered.

    Returns
    -------
    pd.DataFrame
        Indexed by id_col with columns:
            - 'pfam'
            - 'clan'
    """
    if pfam_df is None or pfam_df.empty:
        return pd.DataFrame(columns=["pfam", "clan"])

    df = pfam_df.copy()

    if score_col and score_col in df.columns:
        # Highest score first so groupby().first() picks the best hit
        df = df.sort_values(score_col, ascending=False)

    best = (
        df.groupby(id_col, sort=False)[[pfam_col, clan_col]]
          .first()
          .rename(columns={
              pfam_col: "pfam",
              clan_col: "clan",
          })
    )

    return best

import pandas as pd
import numpy as np
from typing import Optional

import numpy as np
import pandas as pd


import numpy as np
import pandas as pd



def add_pfam_rank(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort by z_score descending, assign pfam_rank = order of first
    appearance of each unique pfam value.

    pfam_rank = 1 for the pfam of the highest z_score hit,
                2 for the next new pfam encountered, etc.
    Targets with the same pfam share the same pfam_rank.
    NaN pfam gets rank = n_unique + 1 (last).

    Parameters
    ----------
    df : dataframe with columns target_id, z_score, pfam

    Returns
    -------
    df sorted by z_score descending with new column pfam_rank (int)
    """
    df = df.sort_values("z_score", ascending=False).copy()

    rank      = 1
    seen      = {}   # pfam → rank assigned at first appearance

    ranks = []
    for pfam in df["pfam"]:
        if pd.isna(pfam):
            ranks.append(None)   # handle at end
            continue
        if pfam not in seen:
            seen[pfam] = rank
            rank += 1
        ranks.append(seen[pfam])

    df["pfam_rank"] = ranks

    # NaN pfam gets rank after all named pfams
    nan_rank = rank
    df["pfam_rank"] = df["pfam_rank"].fillna(nan_rank).astype(int)

    return df

# =====================================================================
# 1. KABSCH #KERNEL (TARGET -> QUERY FRAME)
# =====================================================================

def _batch_kabsch_target_to_query(
    all_q_coords: np.ndarray,
    all_t_coords: np.ndarray,
    offsets: np.ndarray,
    weights: np.ndarray = None,
):
    """
    Vectorized (optionally weighted) Kabsch superposition matching
    `svd_superimpose` math. Superimposes target slices (all_t_coords)
    onto query slices (all_q_coords).

    Parameters
    ----------
    all_q_coords : (total_points, 3) float32
    all_t_coords : (total_points, 3) float32
    offsets : (n_alignments + 1,) int64
    weights : (total_points,) float32, optional
        Per-point weight, same concatenated layout as all_q_coords.
        None (default) reproduces the original unweighted behavior
        exactly -- every point weight 1.0. Zero weights are valid (e.g.
        masking out flexible/non-core residues), but every alignment
        segment must have at least one nonzero weight, or its weighted
        centroid/RMSD normalization divides by zero.

    Returns
    -------
    rmsds : (n_alignments,) float32
        Per-alignment WEIGHTED RMSD (normalized by sum of weights, not
        point count, when weights are given).
    R, t : rotation matrices and translation vectors per alignment.
    """
    n_aligns = len(offsets) - 1
    lengths = np.diff(offsets)  # POINT COUNTS -- stays count-based; used only
                                  # for broadcasting per-alignment values back
                                  # out to per-point arrays via np.repeat.

    if weights is None:
        weights = np.ones(len(all_q_coords), dtype=np.float32)
    weights = weights.astype(np.float32)

    # sum of WEIGHTS per alignment -- replaces `lengths` as the normalizer
    # for means and RMSD (lengths itself is still needed separately, see above)
    w_sums = np.add.reduceat(weights, offsets[:-1])
    if np.any(w_sums <= 0):
        raise ValueError(
            "at least one alignment segment has zero total weight -- "
            "every segment needs at least one nonzero-weight point"
        )

    repeat_counts = lengths

    # 1. Weighted centroids: sum(w_i * X_i) / sum(w_i)
    q_sums = np.add.reduceat(all_q_coords * weights[:, None], offsets[:-1], axis=0)
    t_sums = np.add.reduceat(all_t_coords * weights[:, None], offsets[:-1], axis=0)
    q_means = (q_sums / w_sums[:, None]).astype(np.float32)
    t_means = (t_sums / w_sums[:, None]).astype(np.float32)

    q_means_rep = np.repeat(q_means, repeat_counts, axis=0)
    t_means_rep = np.repeat(t_means, repeat_counts, axis=0)

    # 2. Mean-center coordinates (unchanged -- uses weighted means)
    P_c = all_t_coords - t_means_rep
    Q_c = all_q_coords - q_means_rep

    # 3. WEIGHTED covariance matrices: H = sum_i w_i * outer(P_c_i, Q_c_i)
    P_x, P_y, P_z = P_c[:, 0], P_c[:, 1], P_c[:, 2]
    Q_x, Q_y, Q_z = Q_c[:, 0], Q_c[:, 1], Q_c[:, 2]

    H = np.zeros((n_aligns, 3, 3), dtype=np.float32)
    H[:, 0, 0] = np.add.reduceat(weights * P_x * Q_x, offsets[:-1])
    H[:, 0, 1] = np.add.reduceat(weights * P_x * Q_y, offsets[:-1])
    H[:, 0, 2] = np.add.reduceat(weights * P_x * Q_z, offsets[:-1])
    H[:, 1, 0] = np.add.reduceat(weights * P_y * Q_x, offsets[:-1])
    H[:, 1, 1] = np.add.reduceat(weights * P_y * Q_y, offsets[:-1])
    H[:, 1, 2] = np.add.reduceat(weights * P_y * Q_z, offsets[:-1])
    H[:, 2, 0] = np.add.reduceat(weights * P_z * Q_x, offsets[:-1])
    H[:, 2, 1] = np.add.reduceat(weights * P_z * Q_y, offsets[:-1])
    H[:, 2, 2] = np.add.reduceat(weights * P_z * Q_z, offsets[:-1])

    # 4-6. SVD, reflection correction, rotation matrix -- UNCHANGED. These
    # act purely on H, which already encodes the weighting.
    U, S, Vt = np.linalg.svd(H)
    Vt_T = np.swapaxes(Vt, -1, -2)
    U_T = np.swapaxes(U, -1, -2)
    dets = np.linalg.det(Vt_T @ U_T)
    D = np.repeat(np.eye(3, dtype=np.float32)[None, :, :], n_aligns, axis=0)
    D[:, 2, 2] = dets
    R = Vt_T @ D @ U_T

    # translation -- unchanged formula, uses weighted means
    t = q_means - np.matmul(R, t_means[:, :, None]).squeeze(-1)

    # 7. Apply transformation to ALL target points (unweighted -- every
    # point still gets transformed, weighting only affected the FIT)
    R_rep = np.repeat(R, repeat_counts, axis=0)
    t_rep = np.repeat(t, repeat_counts, axis=0)
    P_transformed = np.matmul(R_rep, all_t_coords[:, :, None]).squeeze(-1) + t_rep

    # 8. WEIGHTED RMSD: sqrt( sum(w_i * ||diff_i||^2) / sum(w_i) )
    diff_sq = np.sum((P_transformed - all_q_coords) ** 2, axis=1)
    msd = np.add.reduceat(weights * diff_sq, offsets[:-1]) / w_sums

    return np.sqrt(msd).astype(np.float32), R, t

def old_batch_kabsch_target_to_query(
    all_q_coords: np.ndarray,
    all_t_coords: np.ndarray,
    all_weights: np.ndarray,
    offsets: np.ndarray
) -> np.ndarray:
    """
    Vectorized Kabsch superposition matching `svd_superimpose` math.
    Superimposes target slices (all_t_coords) onto query slices (all_q_coords).

    Parameters
    ----------
    all_q_coords : (total_points, 3) float32
        Concatenated reference (query) coordinates.
    all_t_coords : (total_points, 3) float32
        Concatenated subject (target) coordinates.
    all_weights : (total_points, 1) float32
        Concatenated position-specific weights - T.B.I.
    offsets : (n_alignments + 1,) int64
        Boundary offsets for contiguous slices in all_q_coords / all_t_coords.

    Returns
    -------
    rmsds : (n_alignments,) float32
        Per-alignment RMSD values after optimal Kabsch superposition.
    """
    n_aligns = len(offsets) - 1
    rmsds = np.zeros(n_aligns, dtype=np.float32)

    # Calculate lengths per alignment
    lengths = np.diff(offsets)

    # 1. Compute Centroids per alignment
    q_sums = np.add.reduceat(all_q_coords, offsets[:-1], axis=0)
    t_sums = np.add.reduceat(all_t_coords, offsets[:-1], axis=0)

    q_means = (q_sums / lengths[:, None]).astype(np.float32)
    t_means = (t_sums / lengths[:, None]).astype(np.float32)

    # Repeat centroids for broadcasted element-wise subtraction
    repeat_counts = lengths
    q_means_rep = np.repeat(q_means, repeat_counts, axis=0)
    t_means_rep = np.repeat(t_means, repeat_counts, axis=0)

    # 2. Mean-center coordinates
    P_c = all_t_coords - t_means_rep  # Subject (Target) centered
    Q_c = all_q_coords - q_means_rep  # Reference (Query) centered

    # 3. Compute Covariance Matrices H = P_c.T @ Q_c for each alignment batch
    # Outer product equivalent across flattened segment blocks
    P_x, P_y, P_z = P_c[:, 0], P_c[:, 1], P_c[:, 2]
    Q_x, Q_y, Q_z = Q_c[:, 0], Q_c[:, 1], Q_c[:, 2]

    # Pre-allocate H matrices: shape (n_aligns, 3, 3)
    H = np.zeros((n_aligns, 3, 3), dtype=np.float32)

    H[:, 0, 0] = np.add.reduceat(P_x * Q_x, offsets[:-1])
    H[:, 0, 1] = np.add.reduceat(P_x * Q_y, offsets[:-1])
    H[:, 0, 2] = np.add.reduceat(P_x * Q_z, offsets[:-1])

    H[:, 1, 0] = np.add.reduceat(P_y * Q_x, offsets[:-1])
    H[:, 1, 1] = np.add.reduceat(P_y * Q_y, offsets[:-1])
    H[:, 1, 2] = np.add.reduceat(P_y * Q_z, offsets[:-1])

    H[:, 2, 0] = np.add.reduceat(P_z * Q_x, offsets[:-1])
    H[:, 2, 1] = np.add.reduceat(P_z * Q_y, offsets[:-1])
    H[:, 2, 2] = np.add.reduceat(P_z * Q_z, offsets[:-1])

    # 4. Vectorized SVD
    U, S, Vt = np.linalg.svd(H)

    # 5. Reflection correction (determinants of Vt.T @ U.T)
    # Vt is (N, 3, 3), U is (N, 3, 3)
    Vt_T = np.swapaxes(Vt, -1, -2)
    U_T = np.swapaxes(U, -1, -2)

    dets = np.linalg.det(Vt_T @ U_T)

    # Construct D matrices (N, 3, 3)
    D = np.repeat(np.eye(3, dtype=np.float32)[None, :, :], n_aligns, axis=0)
    D[:, 2, 2] = dets

    # Rotation matrices: R = Vt.T @ D @ U.T -> shape (n_aligns, 3, 3)
    R = Vt_T @ D @ U_T

    # 6. Compute Translation vectors t = Q_mean - R @ P_mean
    # R: (N, 3, 3), t_means: (N, 3, 1)
    t = q_means - np.matmul(R, t_means[:, :, None]).squeeze(-1)

    # 7. Apply Transformation to All Target Points
    R_rep = np.repeat(R, repeat_counts, axis=0)
    t_rep = np.repeat(t, repeat_counts, axis=0)

    # Superimposed Target coordinates: (R @ P.T).T + t
    P_transformed = np.matmul(R_rep, all_t_coords[:, :, None]).squeeze(-1) + t_rep

    # 8. Compute RMSDs
    diff_sq = np.sum((P_transformed - all_q_coords) ** 2, axis=1)
    msd = np.add.reduceat(diff_sq, offsets[:-1]) / lengths

    return np.sqrt(msd).astype(np.float32), R, t


# =====================================================================
# 2. MAIN METRICS COMPUTATION (USE THIS FUNCTION)
# =====================================================================

import numpy as np
import pandas as pd

def compute_metrics(
    clipped_segs: pd.DataFrame,
    master_metadata: pd.DataFrame,
    query_ascii: np.ndarray,
    target_seqs: np.ndarray,
    query_length: int,
    query_ca_coords: np.ndarray,   # Shape (Q, 3)
    target_ca_coords: np.ndarray,  # Shape (T, 3)
    query_weights: np.ndarray = None,  # Shape (Q,) — position-specific weights
    show_histogram = False
) -> pd.DataFrame:
    """
    Vectorized calculation engine. Computes metrics and superimposes target 
    coordinates into the query frame using clipped segments with optional weighting.
    """
    if clipped_segs.empty:
        return pd.DataFrame()

    # Default to uniform weights of 1.0 if not provided
    if query_weights is None:
        query_weights = np.ones(len(query_ascii), dtype=np.float64)

    clipped_segs = clipped_segs.copy()
    clipped_segs["s_end"] = clipped_segs["s_start"] + clipped_segs["length"] - 1

    # 1. Aggregations per alignment
    seg_agg = clipped_segs.groupby("alignment_id").agg(
        alignment_length=("length", "sum"),
        target_left=("s_start", "min"),
        target_right=("s_end", "max")
    ).reset_index()

    stats_df = pd.merge(seg_agg, master_metadata, on="alignment_id", how="left")

    # 2. Coverages
    stats_df["query_coverage"] = stats_df["alignment_length"] / query_length
    stats_df["target_coverage"] = stats_df["alignment_length"] / stats_df["target_length"]

    # 3. Flat index mapping
    flat_segs = clipped_segs.merge(master_metadata[["alignment_id", "start_idx"]], on="alignment_id", how="inner")
    flat_segs = flat_segs.dropna(subset=["start_idx"])

    # CRITICAL: guarantee each alignment_id's rows are physically contiguous
    # before building the flat index arrays below. clip_segments() processes
    # domain_ranges independently and concatenates per-range results, so an
    # alignment with segments in multiple (discontiguous) domain ranges ends
    # up split across two non-adjacent blocks in flat_segs -- np.unique()
    # further down correctly COUNTS each alignment's total points but has no
    # way to know they aren't contiguous, and np.add.reduceat() silently
    # slices whatever consecutive rows happen to fall at the computed offset,
    # which can also corrupt the RMSD of OTHER, unrelated alignments whose
    # true boundaries no longer line up once ANY prior alignment is split.
    flat_segs = flat_segs.sort_values("alignment_id", kind="stable").reset_index(drop=True)

    final_cols = [
        "alignment_id", "target_id", "z_score", "query_coverage", 
        "target_coverage", "alignment_length", "sequence_identity", 
        "rmsd_new", "rotation", "translation", "target_left", "target_right"
    ]

    if flat_segs.empty:
        stats_df["sequence_identity"] = 0.0
        stats_df["rmsd_new"] = np.nan
        stats_df["rotation"] = None
        stats_df["translation"] = None
        return stats_df[[c for c in final_cols if c in stats_df.columns]]

    q_base = flat_segs["q_start"].values 
    start_idx_values = flat_segs["start_idx"].astype(np.int64).values
    s_base = start_idx_values + flat_segs["s_start"].values
    lengths_vector = flat_segs["length"].values

    repeats = np.repeat(np.arange(len(lengths_vector)), lengths_vector)
    offsets = np.concatenate([np.arange(l) for l in lengths_vector])

    global_q_indices = (q_base[repeats] + offsets).astype(np.int64)
    global_s_indices = (s_base[repeats] + offsets).astype(np.int64)

    valid_bounds_mask = (
        (global_q_indices >= 0) & (global_q_indices < len(query_ascii)) &
        (global_s_indices >= 0) & (global_s_indices < len(target_seqs)) &
        (global_q_indices < len(query_ca_coords)) &
        (global_s_indices < len(target_ca_coords)) &
        (global_q_indices < len(query_weights))
    )

    global_q_indices = global_q_indices[valid_bounds_mask]
    global_s_indices = global_s_indices[valid_bounds_mask]
    alignment_ids_expanded = flat_segs["alignment_id"].values[repeats][valid_bounds_mask]

    if len(global_q_indices) == 0:
        stats_df["sequence_identity"] = 0.0
        stats_df["rmsd_new"] = np.nan
        stats_df["rotation"] = None
        stats_df["translation"] = None
        return stats_df[[c for c in final_cols if c in stats_df.columns]]

    # 4. Sequence Identity Aggregation
    is_identical = (query_ascii[global_q_indices] == target_seqs[global_s_indices]).astype(int)
    identity_counts = pd.DataFrame({"alignment_id": alignment_ids_expanded, "match": is_identical})
    agg = identity_counts.groupby("alignment_id").agg(matches=("match", "sum"), total=("match", "count"))
    agg["sequence_identity"] = agg["matches"] / agg["total"]

    stats_df = stats_df.merge(agg["sequence_identity"].reset_index(), on="alignment_id", how="left")
    stats_df["sequence_identity"] = stats_df["sequence_identity"].fillna(0.0)

    # 5. Target-to-Query Structural Alignment
    unique_align_ids, align_counts = np.unique(alignment_ids_expanded, return_counts=True)

    offsets_numba = np.zeros(len(unique_align_ids) + 1, dtype=np.int64)
    offsets_numba[1:] = np.cumsum(align_counts)

    q_coords_flat = query_ca_coords[global_q_indices].astype(np.float64)
    t_coords_flat = target_ca_coords[global_s_indices].astype(np.float64)
    w_flat = query_weights[global_q_indices].astype(np.float64)  # Extract per-residue weights

    rmsds, rotations, translations = _batch_kabsch_target_to_query(
        q_coords_flat, t_coords_flat, offsets_numba, w_flat
    )
    if show_histogram:
        daliscope.analysis.occupancy.plot_score_distribution(rmsds, cutoff=5.0, title='rmsds')

    kabsch_df = pd.DataFrame({
        "alignment_id": unique_align_ids,
        "rmsd_new": rmsds,
        "rotation": list(rotations),
        "translation": list(translations)
    })

    stats_df = stats_df.merge(kabsch_df, on="alignment_id", how="left")

    return stats_df[[c for c in final_cols if c in stats_df.columns]]

def unweighted_compute_metrics(
    clipped_segs: pd.DataFrame,
    master_metadata: pd.DataFrame,
    query_ascii: np.ndarray,
    target_seqs: np.ndarray,
    query_length: int,
    query_ca_coords: np.ndarray,  # Shape (Q, 3)
    target_ca_coords: np.ndarray   # Shape (T, 3)
) -> pd.DataFrame:
    """
    Vectorized calculation engine. Computes metrics and superimposes target 
    coordinates into the query frame using clipped segments.
    """
    if clipped_segs.empty:
        return pd.DataFrame()

    clipped_segs = clipped_segs.copy()
    clipped_segs["s_end"] = clipped_segs["s_start"] + clipped_segs["length"] - 1

    # 1. Aggregations per alignment
    seg_agg = clipped_segs.groupby("alignment_id").agg(
        alignment_length=("length", "sum"),
        target_left=("s_start", "min"),
        target_right=("s_end", "max")
    ).reset_index()

    stats_df = pd.merge(seg_agg, master_metadata, on="alignment_id", how="left")

    # 2. Coverages
    stats_df["query_coverage"] = stats_df["alignment_length"] / query_length
    stats_df["target_coverage"] = stats_df["alignment_length"] / stats_df["target_length"]

    # 3. Flat index mapping
    flat_segs = clipped_segs.merge(master_metadata[["alignment_id", "start_idx"]], on="alignment_id", how="inner")
    flat_segs = flat_segs.dropna(subset=["start_idx"])

    final_cols = [
        "alignment_id", "target_id", "z_score", "query_coverage", 
        "target_coverage", "alignment_length", "sequence_identity", 
        "rmsd_new", "rotation", "translation", "target_left", "target_right"
    ]

    if flat_segs.empty:
        stats_df["sequence_identity"] = 0.0
        stats_df["rmsd_new"] = np.nan
        stats_df["rotation"] = None
        stats_df["translation"] = None
        return stats_df[[c for c in final_cols if c in stats_df.columns]]

    q_base = flat_segs["q_start"].values 
    start_idx_values = flat_segs["start_idx"].astype(np.int64).values
    s_base = start_idx_values + flat_segs["s_start"].values
    lengths_vector = flat_segs["length"].values

    repeats = np.repeat(np.arange(len(lengths_vector)), lengths_vector)
    offsets = np.concatenate([np.arange(l) for l in lengths_vector])

    global_q_indices = (q_base[repeats] + offsets).astype(np.int64)
    global_s_indices = (s_base[repeats] + offsets).astype(np.int64)

    valid_bounds_mask = (
        (global_q_indices >= 0) & (global_q_indices < len(query_ascii)) &
        (global_s_indices >= 0) & (global_s_indices < len(target_seqs)) &
        (global_q_indices < len(query_ca_coords)) &
        (global_s_indices < len(target_ca_coords))
    )

    global_q_indices = global_q_indices[valid_bounds_mask]
    global_s_indices = global_s_indices[valid_bounds_mask]
    alignment_ids_expanded = flat_segs["alignment_id"].values[repeats][valid_bounds_mask]

    if len(global_q_indices) == 0:
        stats_df["sequence_identity"] = 0.0
        stats_df["rmsd_new"] = np.nan
        stats_df["rotation"] = None
        stats_df["translation"] = None
        return stats_df[[c for c in final_cols if c in stats_df.columns]]

    # 4. Sequence Identity Aggregation
    is_identical = (query_ascii[global_q_indices] == target_seqs[global_s_indices]).astype(int)
    identity_counts = pd.DataFrame({"alignment_id": alignment_ids_expanded, "match": is_identical})
    agg = identity_counts.groupby("alignment_id").agg(matches=("match", "sum"), total=("match", "count"))
    agg["sequence_identity"] = agg["matches"] / agg["total"]

    stats_df = stats_df.merge(agg["sequence_identity"].reset_index(), on="alignment_id", how="left")
    stats_df["sequence_identity"] = stats_df["sequence_identity"].fillna(0.0)

    # 5. Target-to-Query Structural Alignment
    # Guarantee offsets match the EXACT order of global_q_indices
    unique_align_ids, align_counts = np.unique(alignment_ids_expanded, return_counts=True)

    offsets_numba = np.zeros(len(unique_align_ids) + 1, dtype=np.int64)
    offsets_numba[1:] = np.cumsum(align_counts)

    q_coords_flat = query_ca_coords[global_q_indices].astype(np.float64)
    t_coords_flat = target_ca_coords[global_s_indices].astype(np.float64)

    print(len(q_coords_flat), len(t_coords_flat), len(offsets_numba))

    rmsds, rotations, translations = _batch_kabsch_target_to_query(
        q_coords_flat, t_coords_flat, offsets_numba
    )

    kabsch_df = pd.DataFrame({
        "alignment_id": unique_align_ids,
        "rmsd_new": rmsds,
        "rotation": list(rotations),
        "translation": list(translations)
    })

    stats_df = stats_df.merge(kabsch_df, on="alignment_id", how="left")

    return stats_df[[c for c in final_cols if c in stats_df.columns]]


def oldish_compute_metrics(clipped_segs: pd.DataFrame, master_metadata: pd.DataFrame,
                    query_ascii: np.ndarray, target_seqs: np.ndarray, query_length: int) -> pd.DataFrame:
    """
    Vectorized calculation engine. Computes alignment metrics including target_left 
    and target_right boundary coordinates of aligned residues.
    """
    if clipped_segs.empty:
        return pd.DataFrame()

    # Calculate end coordinate for each segment on target
    clipped_segs = clipped_segs.copy()
    clipped_segs["s_end"] = clipped_segs["s_start"] + clipped_segs["length"] - 1

    # 1. Compute total alignment length, target_left, and target_right per alignment
    seg_agg = clipped_segs.groupby("alignment_id").agg(
        alignment_length=("length", "sum"),
        target_left=("s_start", "min"),
        target_right=("s_end", "max")
    ).reset_index()

    stats_df = pd.merge(seg_agg, master_metadata, on="alignment_id", how="left")

    # 2. Compute coverages
    stats_df["query_coverage"] = stats_df["alignment_length"] / query_length
    stats_df["target_coverage"] = stats_df["alignment_length"] / stats_df["target_length"]

    # 3. Vectorized Sequence Identity
    flat_segs = clipped_segs.merge(master_metadata[["alignment_id", "start_idx"]], on="alignment_id", how="inner")
    flat_segs = flat_segs.dropna(subset=["start_idx"])

    if flat_segs.empty:
        stats_df["sequence_identity"] = 0.0
        final_cols = ["alignment_id", "target_id", "z_score", "query_coverage", "target_coverage", 
                      "alignment_length", "sequence_identity", "target_left", "target_right"]
        return stats_df[[c for c in final_cols if c in stats_df.columns]]

    # Secure execution data primitives
    q_base = flat_segs["q_start"].values 
    start_idx_values = flat_segs["start_idx"].astype(np.int64).values
    s_base = start_idx_values + flat_segs["s_start"].values
    lengths_vector = flat_segs["length"].values

    # Expand arrays to generate flat index coordinates
    repeats = np.repeat(np.arange(len(lengths_vector)), lengths_vector)
    offsets = np.concatenate([np.arange(l) for l in lengths_vector])

    global_q_indices = (q_base[repeats] + offsets).astype(np.int64)
    global_s_indices = (s_base[repeats] + offsets).astype(np.int64)

    # Final boundary array check to shield from bad coordinate values
    valid_bounds_mask = (
        (global_q_indices >= 0) & (global_q_indices < len(query_ascii)) &
        (global_s_indices >= 0) & (global_s_indices < len(target_seqs))
    )

    global_q_indices = global_q_indices[valid_bounds_mask]
    global_s_indices = global_s_indices[valid_bounds_mask]
    alignment_ids_expanded = flat_segs["alignment_id"].values[repeats][valid_bounds_mask]

    if len(global_q_indices) == 0:
        stats_df["sequence_identity"] = 0.0
        final_cols = ["alignment_id", "target_id", "z_score", "query_coverage", "target_coverage", 
                      "alignment_length", "sequence_identity", "target_left", "target_right"]
        return stats_df[[c for c in final_cols if c in stats_df.columns]]

    # Vectorized comparison matching raw ASCII integers
    is_identical = (query_ascii[global_q_indices] == target_seqs[global_s_indices]).astype(int)

    # Aggregate matches
    identity_counts = pd.DataFrame({"alignment_id": alignment_ids_expanded, "match": is_identical})
    agg = identity_counts.groupby("alignment_id").agg(matches=("match", "sum"), total=("match", "count"))
    agg["sequence_identity"] = agg["matches"] / agg["total"]

    # Combine back to master statistics compilation frame
    stats_df = stats_df.merge(agg["sequence_identity"].reset_index(), on="alignment_id", how="left")
    stats_df["sequence_identity"] = stats_df["sequence_identity"].fillna(0.0)

    final_cols = [
        "alignment_id", "target_id", "z_score", "query_coverage", 
        "target_coverage", "alignment_length", "sequence_identity", 
        "target_left", "target_right"
    ]
    return stats_df[[c for c in final_cols if c in stats_df.columns]]

def old_compute_metrics(clipped_segs: pd.DataFrame, master_metadata: pd.DataFrame,
                    query_ascii: np.ndarray, target_seqs: np.ndarray, query_length: int) -> pd.DataFrame:
    """
    Vectorized calculation engine. Fully robust against residual NaN integer overflows.
    """
    if clipped_segs.empty:
        return pd.DataFrame()

    # 1. Compute total equivalent positions per alignment
    lengths = clipped_segs.groupby("alignment_id")["length"].sum().reset_index(name="alignment_length")
    stats_df = pd.merge(lengths, master_metadata, on="alignment_id", how="left")

    # 2. Compute coverages (safely works with 0 lengths now)
    stats_df["query_coverage"] = stats_df["alignment_length"] / query_length
    stats_df["target_coverage"] = stats_df["alignment_length"] / stats_df["target_length"]

    # 2. Compute coverages
    stats_df["query_coverage"] = stats_df["alignment_length"] / query_length
    stats_df["target_coverage"] = stats_df["alignment_length"] / stats_df["target_length"]

    # 3. Vectorized Sequence Identity
    # Pull across metadata mapping
    flat_segs = clipped_segs.merge(master_metadata[["alignment_id", "start_idx"]], on="alignment_id", how="inner")

    # HARD FILTER: Drop rows where start_idx is missing before any conversions happen
    flat_segs = flat_segs.dropna(subset=["start_idx"])

    if flat_segs.empty:
        stats_df["sequence_identity"] = 0.0
        return stats_df

    # Secure execution data primitives
    q_base = flat_segs["q_start"].values 
    # Forcing float values safely to integers now that we guaranteed NaNs are gone
    start_idx_values = flat_segs["start_idx"].astype(np.int64).values
    s_base = start_idx_values + flat_segs["s_start"].values
    lengths_vector = flat_segs["length"].values

    # Expand arrays to generate flat index coordinates
    repeats = np.repeat(np.arange(len(lengths_vector)), lengths_vector)
    offsets = np.concatenate([np.arange(l) for l in lengths_vector])

    global_q_indices = (q_base[repeats] + offsets).astype(np.int64)
    global_s_indices = (s_base[repeats] + offsets).astype(np.int64)

    # Final boundary array check to shield from bad coordinate values
    valid_bounds_mask = (
        (global_q_indices >= 0) & (global_q_indices < len(query_ascii)) &
        (global_s_indices >= 0) & (global_s_indices < len(target_seqs))
    )

    global_q_indices = global_q_indices[valid_bounds_mask]
    global_s_indices = global_s_indices[valid_bounds_mask]
    alignment_ids_expanded = flat_segs["alignment_id"].values[repeats][valid_bounds_mask]

    if len(global_q_indices) == 0:
        stats_df["sequence_identity"] = 0.0
        return stats_df

    # Vectorized comparison matching raw ASCII integers
    is_identical = (query_ascii[global_q_indices] == target_seqs[global_s_indices]).astype(int)

    # Aggregate matches
    identity_counts = pd.DataFrame({"alignment_id": alignment_ids_expanded, "match": is_identical})
    agg = identity_counts.groupby("alignment_id").agg(matches=("match", "sum"), total=("match", "count"))
    agg["sequence_identity"] = agg["matches"] / agg["total"]

    # Combine back to master statistics compilation frame
    stats_df = stats_df.merge(agg["sequence_identity"].reset_index(), on="alignment_id", how="left")
    stats_df["sequence_identity"] = stats_df["sequence_identity"].fillna(0.0)

    final_cols = ["alignment_id", "target_id", "z_score", "query_coverage", "target_coverage", "alignment_length", "sequence_identity"]
    return stats_df[final_cols]

# ------------------------------------------------------------------ #
# Parser for domain range string                                      #
# ------------------------------------------------------------------ #

def parse_domain_ranges(range_str: str) -> list:
    tokens = [t.strip() for t in re.split(r'[\s,]+', range_str.strip()) if t.strip()]
    if not tokens:
        raise ValueError("No domain ranges provided")
    ranges = []
    for token in tokens:
        m = re.fullmatch(r'(\d+)\s*-\s*(\d+)', token)
        if not m:
            raise ValueError(
                f"Cannot parse '{token}' — expected format: start-end (e.g. 10-150)"
            )
        start, end = int(m.group(1)), int(m.group(2))
        if start >= end:
            raise ValueError(f"Range {token}: start must be less than end")
        ranges.append((start, end))
    return ranges

def clip_segments(segments_df: pd.DataFrame, domain_ranges: list) -> pd.DataFrame:
    if segments_df.empty or not domain_ranges:
        return pd.DataFrame()

    clipped_list = []

    for q_dom_start, q_dom_end in domain_ranges:
        df = segments_df.copy()
        df['q_end'] = df['q_start'] + df['length'] - 1

        # Keep segments that overlap with the current domain range
        df = df[(df['q_end'] >= q_dom_start) & (df['q_start'] <= q_dom_end)].copy()
        if df.empty:
            continue

        # Left offset trimming: MUST apply identical offset to BOTH q_start and s_start
        left_clip = np.maximum(0, q_dom_start - df['q_start'])
        right_clip = np.maximum(0, df['q_end'] - q_dom_end)

        df['q_start'] += left_clip
        df['s_start'] += left_clip  # <--- CRITICAL: Shift target sequence start in tandem
        df['length'] -= (left_clip + right_clip)

        clipped_list.append(df.drop(columns=['q_end']))

    if not clipped_list:
        return pd.DataFrame(columns=segments_df.columns)

    return pd.concat(clipped_list, ignore_index=True)


def old_clip_segments(segments_df: pd.DataFrame, domain_ranges: list[tuple[int, int]]) -> pd.DataFrame:
    """
    Vectorized matrix clipping of structural segments against multiple domain intervals.
    Accepts raw dataframes; completely decoupled from the Project object.
    """
    if segments_df.empty or not domain_ranges:
        return segments_df.copy()

    df = segments_df.copy()
    df["q_end"] = df["q_start"] + df["length"] - 1

    # Vectorized computation across all domain boundaries simultaneously
    clipped_blocks = []
    for d_start, d_end in domain_ranges:
        # Evaluate overlaps for all alignment fragments in parallel
        overlap_start = np.maximum(df["q_start"].values, d_start)
        overlap_end = np.minimum(df["q_end"].values, d_end)
        overlap_len = overlap_end - overlap_start + 1

        valid_mask = overlap_len > 0
        if not np.any(valid_mask):
            continue

        sub_df = df[valid_mask].copy()
        # Vectorized proportional coordinate shifts for target starts
        sub_df["s_start"] = sub_df["s_start"] + (overlap_start[valid_mask] - sub_df["q_start"])
        sub_df["q_start"] = overlap_start[valid_mask]
        sub_df["length"] = overlap_len[valid_mask]
        clipped_blocks.append(sub_df)

    if not clipped_blocks:
        return pd.DataFrame(columns=["alignment_id", "q_start", "s_start", "length"])

    result_df = pd.concat(clipped_blocks, ignore_index=True)
    return result_df.drop_duplicates(subset=["alignment_id", "q_start", "s_start"]).drop(columns=["q_end"])
