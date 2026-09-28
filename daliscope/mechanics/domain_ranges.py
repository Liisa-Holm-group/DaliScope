"""Convert domain boundaries without changing the selected query residues.

PDB/PUU domain strings and community masks use one-based inclusive ranges.
Clipping and aligned-array slices use zero-based half-open ranges.
"""
from numbers import Integral
import re


def _validated_ranges(ranges, *, pdb):
    result = []
    for pair in ranges:
        try:
            start, end = pair
        except (TypeError, ValueError) as exc:
            raise ValueError("Each domain range must contain two integer boundaries") from exc
        if any(isinstance(value, bool) or not isinstance(value, Integral) for value in (start, end)):
            raise ValueError("Domain boundaries must be integers")
        start, end = int(start), int(end)
        if (pdb and (start < 1 or end < start)) or (not pdb and (start < 0 or end <= start)):
            convention = "one-based inclusive" if pdb else "zero-based half-open"
            raise ValueError(f"Invalid {convention} domain range: {start}-{end}")
        result.append((start, end))
    return result


def pdb_ranges_to_clipping(ranges):
    """Convert one-based inclusive ranges to zero-based half-open ranges."""
    return [(start - 1, end) for start, end in _validated_ranges(ranges, pdb=True)]


def clipping_ranges_to_pdb(ranges):
    """Convert zero-based half-open ranges to one-based inclusive ranges."""
    return [(start + 1, end) for start, end in _validated_ranges(ranges, pdb=False)]


def _convert_domain_string(text, converter):
    if not isinstance(text, str):
        raise ValueError("Domain ranges must be supplied as a string")
    if not text.strip():
        return ""
    domains = []
    for domain in text.split(","):
        ranges = []
        for segment in domain.split("_"):
            match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", segment)
            if match is None:
                raise ValueError(f"Invalid domain segment: {segment!r}; expected start-end")
            ranges.append(tuple(map(int, match.groups())))
        domains.append("_".join(f"{start}-{end}" for start, end in converter(ranges)))
    return ", ".join(domains)


def pdb_domain_string_to_clipping(text):
    """Convert PDB/PUU text, preserving comma-separated and discontinuous domains."""
    return _convert_domain_string(text, pdb_ranges_to_clipping)


def clipping_domain_string_to_pdb(text):
    """Convert clipping text back to the domain viewer's PDB residue convention."""
    return _convert_domain_string(text, clipping_ranges_to_pdb)
