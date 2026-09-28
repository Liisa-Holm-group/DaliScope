"""Mutate segment boundaries"""

import numpy as np
from scipy import sparse
from collections import defaultdict, deque
from typing import List, Tuple, Dict, Any
from dataclasses import replace, is_dataclass

def OBSOLETE_split_multi_sheet_strands(
    segments: List[Tuple[str, int, int]],
    strand_ranges: List[Tuple[int, int]],
    adj: sparse.coo_matrix,
    ss_str: str,
    contact_cutoff: float = 8.0,
    min_strand_len: int = 2
) -> Tuple[List[Tuple[str, int, int]], List[Tuple[int, int]], str]:
    """
    Detects long strands that participate in two disjoint contacts and cuts 
    them at the point of zero internal/local contact continuity.
    
    Args:
        segments: List of SSE tuples [('E', start, end), ...]
        strand_ranges: List of (start, end) strand tuples
        adj: Sparse distance matrix between residues (COO matrix)
        ss_str: Secondary structure string (e.g., "--EEEE---HHHH--")
        contact_cutoff: Max distance (Å) to consider two residues in contact
        min_strand_len: Minimum residue length to retain a split sub-strand
        
    Returns:
        new_segments, new_strand_ranges, new_ss_str
    """
    # 1. Binarize contact matrix within cutoff
    adj_csr = adj.tocsr()

    # 2. Track new segments and strands
    new_segments = []
    new_strand_ranges = []
    ss_chars = list(ss_str)

    # Convert strand ranges to a set for quick lookup
    strand_set = set(strand_ranges)

    for seg in segments:
        sse_type, start_res, end_res = seg[0], seg[1], seg[2]

        # Only inspect Beta strands ('E')
        if sse_type != 'E' or (start_res, end_res) not in strand_set:
            new_segments.append(seg)
            continue

        length = end_res - start_res + 1

        # Short strands (< 5 residues) rarely cross distinct sheets without contact
        if length < 5:
            new_segments.append(seg)
            new_strand_ranges.append((start_res, end_res))
            continue

        # Extract local contact matrix for residues in this strand
        local_dist = adj_csr[start_res:end_res + 1, start_res:end_res + 1].toarray()
        # Non-zero entries within threshold count as contacts
        contact_map = (local_dist > 0) & (local_dist <= contact_cutoff)

        # 3. Check for sequential connectivity along the strand length
        # A cut point occurs if no residue in section [0..k] contacts any residue in [k+1..length-1]
        cut_indices = []
        for k in range(0, length - 1):
            left_to_right_contacts = contact_map[:k + 1, k + 1:]
            if not np.any(left_to_right_contacts):
                cut_indices.append(k)

        # If no internal disconnection was found, keep strand intact
        if not cut_indices:
            new_segments.append(seg)
            new_strand_ranges.append((start_res, end_res))
            continue

        # 4. Perform the cuts
        # Choose cut points and create sub-strands
        sub_ranges = []
        curr_start = start_res

        for cut_k in cut_indices:
            cut_res = start_res + cut_k
            if (cut_res - curr_start + 1) >= min_strand_len:
                sub_ranges.append((curr_start, cut_res))
                curr_start = cut_res + 1
            else:
                # Residues before cut point are too short -> revert to coil/loop
                for r in range(curr_start, cut_res + 1):
                    ss_chars[r] = '-'
                curr_start = cut_res + 1

        # Final sub-segment after last cut
        if (end_res - curr_start + 1) >= min_strand_len:
            sub_ranges.append((curr_start, end_res))
        else:
            for r in range(curr_start, end_res + 1):
                ss_chars[r] = '-'

        # Add validated sub-strands to output
        for s_start, s_end in sub_ranges:
            new_segments.append(('E', s_start, s_end))
            new_strand_ranges.append((s_start, s_end))

    new_ss_str = "".join(ss_chars)
    return new_segments, new_strand_ranges, new_ss_str

from dataclasses import is_dataclass, replace
from typing import Any, Tuple
import numpy as np


def split_strands_by_saba_hbond_transitions(
    segments: np.ndarray,
    strand_ranges: np.ndarray,
    sheet_partners: np.ndarray,
    pairs: np.ndarray,
    ss_arr: np.ndarray,
    min_strand_len: int = 2,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Detects long/chimeric beta strands that H-bond with two distinct sheet regions,

    splits them, and re-maps existing frozen sheet_partners pairings to the
    new sub-strands.

    All inputs and return values are NumPy arrays.
    """
    n_res = len(ss_arr)
    new_ss_arr = ss_arr.copy()

    # Convert strand ranges into a lookup set/mapping
    strand_tuples = [tuple(r) for r in strand_ranges]
    strand_set = set(strand_tuples)

    # 1. Map each residue index to its parent strand ID
    res_to_strand_idx = np.full(n_res, -1, dtype=int)
    for st_idx, (s_start, s_end) in enumerate(strand_tuples):
        res_to_strand_idx[s_start : s_end + 1] = st_idx

    # 2. Map each residue index to all strands it forms SABA H-bonds with
    res_hbond_partners = {}
    for p in pairs:
        i, j = p.i, p.j
        st_i = res_to_strand_idx[i] if 0 <= i < n_res else -1
        st_j = res_to_strand_idx[j] if 0 <= j < n_res else -1

        if st_i != -1 and st_j != -1 and st_i != st_j:
            res_hbond_partners.setdefault(i, set()).add(st_j)
            res_hbond_partners.setdefault(j, set()).add(st_i)

    new_segments = []
    new_strand_ranges = []
    old_to_new_strands = {}

    # 3. Analyze each segment for H-bond partner transitions
    for seg in segments:
        sse_type, start_res, end_res = seg[0], int(seg[1]), int(seg[2])
        st_tuple = (start_res, end_res)

        if sse_type != "E" or st_tuple not in strand_set:
            new_segments.append((sse_type, start_res, end_res))
            continue

        length = end_res - start_res + 1
        if length < 5:  # Skip short strands
            new_segments.append((sse_type, start_res, end_res))
            new_strand_ranges.append(st_tuple)
            old_to_new_strands[st_tuple] = [st_tuple]
            continue

        strand_profile = [
            res_hbond_partners.get(r_idx, set())
            for r_idx in range(start_res, end_res + 1)
        ]

        cut_indices = []
        for k in range(0, length - 1):
            left_partners = set().union(*strand_profile[: k + 1])
            right_partners = set().union(*strand_profile[k + 1 :])

            if (
                left_partners
                and right_partners
                and left_partners.isdisjoint(right_partners)
            ):
                if not (
                    len(strand_profile[k]) == 0 and len(strand_profile[k + 1]) == 0
                ):
                    cut_indices.append(k)

        if not cut_indices:
            new_segments.append((sse_type, start_res, end_res))
            new_strand_ranges.append(st_tuple)
            old_to_new_strands[st_tuple] = [st_tuple]
            continue

        # 4. Perform strand splitting
        sub_ranges = []
        curr_start = start_res

        filtered_cuts = []
        for c in cut_indices:
            if not filtered_cuts or (c - filtered_cuts[-1]) >= min_strand_len:
                filtered_cuts.append(c)

        for cut_k in filtered_cuts:
            cut_res = start_res + cut_k
            if (cut_res - curr_start + 1) >= min_strand_len:
                sub_ranges.append((curr_start, cut_res))
            else:
                new_ss_arr[curr_start : cut_res + 1] = "-"
            curr_start = cut_res + 1

        if (end_res - curr_start + 1) >= min_strand_len:
            sub_ranges.append((curr_start, end_res))
        else:
            new_ss_arr[curr_start : end_res + 1] = "-"

        old_to_new_strands[st_tuple] = sub_ranges

        for s_start, s_end in sub_ranges:
            new_segments.append(("E", s_start, s_end))
            new_strand_ranges.append((s_start, s_end))

    # 5. Re-map sheet_partners pairings (handling frozen/immutable objects)
    def map_to_sub_strand(
        orig_strand: Tuple[int, int], ref_res: int
    ) -> Tuple[int, int]:
        sub_strands = old_to_new_strands.get(orig_strand, [orig_strand])
        for s_start, s_end in sub_strands:
            if s_start <= ref_res <= s_end:
                return (s_start, s_end)
        return sub_strands[0]

    updated_sheet_partners = []
    for partner in sheet_partners:
        s_a = (
            tuple(partner.strand_a)
            if isinstance(partner.strand_a, np.ndarray)
            else partner.strand_a
        )
        s_b = (
            tuple(partner.strand_b)
            if isinstance(partner.strand_b, np.ndarray)
            else partner.strand_b
        )

        res_a = getattr(partner, "res_a", s_a[0])
        res_b = getattr(partner, "res_b", s_b[0])

        new_sa = map_to_sub_strand(s_a, res_a)
        new_sb = map_to_sub_strand(s_b, res_b)

        if is_dataclass(partner):
            new_partner = replace(partner, strand_a=new_sa, strand_b=new_sb)
        elif hasattr(partner, "_replace"):
            new_partner = partner._replace(strand_a=new_sa, strand_b=new_sb)
        else:
            partner_cls = type(partner)
            fields = {k: v for k, v in partner.__dict__.items()}
            fields["strand_a"] = new_sa
            fields["strand_b"] = new_sb
            new_partner = partner_cls(**fields)

        updated_sheet_partners.append(new_partner)

    # Convert all output items to NumPy arrays
    segments_arr = np.array(new_segments, dtype=object)
    strand_ranges_arr = np.asarray(new_strand_ranges, dtype=int)
    sheet_partners_arr = np.asarray(updated_sheet_partners, dtype=object)

    return segments_arr, strand_ranges_arr, sheet_partners_arr, new_ss_arr


def old_split_strands_by_saba_hbond_transitions(
    segments: List[Tuple[str, int, int]],
    strand_ranges: List[Tuple[int, int]],
    sheet_partners: List[Any],
    pairs: List[Any],
    ss_str: str,
    min_strand_len: int = 2
) -> Tuple[List[Tuple[str, int, int]], List[Tuple[int, int]], List[Any], str]:
    """
    Detects long/chimeric beta strands that H-bond with two distinct sheet regions,
    splits them, and re-maps existing frozen sheet_partners pairings to the new sub-strands.
    """
    # 1. Map each residue index to its parent strand ID
    res_to_strand_idx: Dict[int, int] = {}
    for st_idx, (s_start, s_end) in enumerate(strand_ranges):
        for r_idx in range(s_start, s_end + 1):
            res_to_strand_idx[r_idx] = st_idx

    # 2. Map each residue index to all strands it forms SABA H-bonds with
    res_hbond_partners: Dict[int, set] = defaultdict(set)
    for p in pairs:
        i, j = p.i, p.j
        st_i = res_to_strand_idx.get(i)
        st_j = res_to_strand_idx.get(j)

        if st_i is not None and st_j is not None and st_i != st_j:
            res_hbond_partners[i].add(st_j)
            res_hbond_partners[j].add(st_i)

    new_segments = []
    new_strand_ranges = []
    ss_chars = list(ss_str)
    strand_set = set(strand_ranges)

    old_to_new_strands: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}

    # 3. Analyze each segment for H-bond partner transitions
    for seg in segments:
        sse_type, start_res, end_res = seg[0], seg[1], seg[2]

        if sse_type != 'E' or (start_res, end_res) not in strand_set:
            new_segments.append(seg)
            continue

        length = end_res - start_res + 1
        if length < 5:  # Skip short strands
            new_segments.append(seg)
            new_strand_ranges.append((start_res, end_res))
            old_to_new_strands[(start_res, end_res)] = [(start_res, end_res)]
            continue

        strand_profile = [
            res_hbond_partners.get(r_idx, set())
            for r_idx in range(start_res, end_res + 1)
        ]

        cut_indices = []
        for k in range(0, length - 1):
            left_partners = set().union(*strand_profile[:k + 1])
            right_partners = set().union(*strand_profile[k + 1:])

            if left_partners and right_partners and left_partners.isdisjoint(right_partners):
                if not (len(strand_profile[k]) == 0 and len(strand_profile[k+1]) == 0):
                    cut_indices.append(k)

        if not cut_indices:
            new_segments.append(seg)
            new_strand_ranges.append((start_res, end_res))
            old_to_new_strands[(start_res, end_res)] = [(start_res, end_res)]
            continue

        # 4. Perform strand splitting
        sub_ranges = []
        curr_start = start_res

        filtered_cuts = []
        for c in cut_indices:
            if not filtered_cuts or (c - filtered_cuts[-1]) >= min_strand_len:
                filtered_cuts.append(c)

        for cut_k in filtered_cuts:
            cut_res = start_res + cut_k
            if (cut_res - curr_start + 1) >= min_strand_len:
                sub_ranges.append((curr_start, cut_res))
                curr_start = cut_res + 1
            else:
                for r in range(curr_start, cut_res + 1):
                    ss_chars[r] = '-'
                curr_start = cut_res + 1

        if (end_res - curr_start + 1) >= min_strand_len:
            sub_ranges.append((curr_start, end_res))
        else:
            for r in range(curr_start, end_res + 1):
                ss_chars[r] = '-'

        old_to_new_strands[(start_res, end_res)] = sub_ranges

        for s_start, s_end in sub_ranges:
            new_segments.append(('E', s_start, s_end))
            new_strand_ranges.append((s_start, s_end))

    # 5. Re-map sheet_partners pairings (handling frozen/immutable objects)
    def map_to_sub_strand(orig_strand: Tuple[int, int], ref_res: int) -> Tuple[int, int]:
        sub_strands = old_to_new_strands.get(orig_strand, [orig_strand])
        for s_start, s_end in sub_strands:
            if s_start <= ref_res <= s_end:
                return (s_start, s_end)
        return sub_strands[0]

    updated_sheet_partners = []
    for partner in sheet_partners:
        s_a = partner.strand_a
        s_b = partner.strand_b

        res_a = getattr(partner, 'res_a', s_a[0])
        res_b = getattr(partner, 'res_b', s_b[0])

        new_sa = map_to_sub_strand(s_a, res_a)
        new_sb = map_to_sub_strand(s_b, res_b)

        if is_dataclass(partner):
            new_partner = replace(partner, strand_a=new_sa, strand_b=new_sb)
        elif hasattr(partner, '_replace'):
            new_partner = partner._replace(strand_a=new_sa, strand_b=new_sb)
        else:
            partner_cls = type(partner)
            fields = {k: v for k, v in partner.__dict__.items()}
            fields['strand_a'] = new_sa
            fields['strand_b'] = new_sb
            new_partner = partner_cls(**fields)

        updated_sheet_partners.append(new_partner)

    new_ss_str = "".join(ss_chars)
    return new_segments, new_strand_ranges, updated_sheet_partners, new_ss_str
