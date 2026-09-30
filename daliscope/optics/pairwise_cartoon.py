"""Pairwise all-atom cartoons using Daliviewer target structures.

The DALI pack stores each target's C-alpha coordinates after subtracting that
chain's centroid. Its per-hit R and t act on those centered pack coordinates.
The readpdb.cgi endpoint returns uncentered coordinates, so this module centers
the requested chain before applying R and t.
"""

import re
import urllib.request

import numpy as np
import py3Dmol


_DALIVIEWER_URL = "http://ekhidna2.biocenter.helsinki.fi/cgi-bin/daliviewer/readpdb.cgi"
_IDENTITY = "1.000000,0.000000,0.000000,0.000000,1.000000,0.000000,0.000000,0.000000,1.000000"
_REQUEST_TIMEOUT_SECONDS = 20


def _target_parts(target_id):
    target_id = str(target_id).strip()
    if not re.fullmatch(r"[A-Za-z0-9]{5}", target_id):
        raise ValueError(
            f"Cannot fetch target {target_id!r} from Daliviewer: expected a "
            "four-character DALI structure code followed by one chain ID."
        )
    return target_id[:4].lower(), target_id[4]


def _atom_record(line):
    return line.startswith(("ATOM  ", "HETATM"))


def _ca_coordinates(lines):
    coords = []
    for line in lines:
        if not _atom_record(line) or len(line) < 54:
            continue
        if line[12:16].strip() != "CA" or line[16] not in (" ", "A"):
            continue
        # Calcium ions also have atom name CA; only peptide C-alpha atoms
        # belong in the centroid used by the DALI pack.
        if line[76:78].strip().upper() == "CA":
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords, dtype=float)


def extract_pdb_from_html(html_text):
    """Extract PDB records from Daliviewer's HTML-wrapped or plain response."""
    match = re.search(r"<pre(?:\s[^>]*)?>(.*?)</pre\s*>", html_text, re.DOTALL | re.IGNORECASE)
    if match:
        return match.group(1).strip()
    valid_prefixes = ("ATOM", "HETATM", "REMARK", "HEADER", "TITLE", "MODEL", "ENDMDL", "TER", "END")
    return "\n".join(line for line in html_text.splitlines() if line.startswith(valid_prefixes))


def generate_daliviewer_url(target_id, base_url=_DALIVIEWER_URL):
    """Request raw coordinates for the structure code in a DALI target ID."""
    pdbid, _ = _target_parts(target_id)
    return f"{base_url}?html=1&jobid=pdb&pdbid={pdbid}&u={_IDENTITY}&t=0,0,0"


def transform_pdb(pdb_text, R, t, chain=None, expected_ca_count=None):
    """Center one target chain and transform it into the query's coordinates.

    Daliviewer returns original coordinates. The corresponding DALI-pack
    coordinates are centered on that chain's C-alpha centroid; R and t were
    computed from those centered coordinates. For standalone use, omitting
    ``chain`` transforms the whole PDB as in the author-supplied function.
    """
    lines = pdb_text.splitlines()
    atom_lines = [line for line in lines if _atom_record(line)]
    if chain is not None:
        if len(chain) != 1:
            raise ValueError("chain must be a single-character chain ID")
        atom_lines = [line for line in atom_lines if len(line) > 21 and line[21] == chain]
    if not atom_lines:
        detail = f"chain {chain!r}" if chain is not None else "any atom records"
        raise ValueError(f"Daliviewer returned no {detail}; the target structure may be unavailable.")

    ca_coords = _ca_coordinates(atom_lines)
    if not len(ca_coords):
        raise ValueError(f"Daliviewer returned no C-alpha coordinates for chain {chain!r}.")
    if expected_ca_count is not None and len(ca_coords) != int(expected_ca_count):
        raise ValueError(
            f"Daliviewer returned {len(ca_coords)} C-alpha residues for chain {chain!r}, "
            f"but this DALI pack expects {int(expected_ca_count)}. "
            "The remote structure may differ from the pack; use a different target."
        )

    R = np.asarray(R, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    if R.shape != (3, 3) or t.shape != (3,) or not (np.isfinite(R).all() and np.isfinite(t).all()):
        raise ValueError("This alignment has no valid 3D rotation and translation (R, t).")
    centroid = ca_coords.mean(axis=0)

    transformed = []
    for line in atom_lines:
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        point = (xyz - centroid) @ R.T + t
        transformed.append(
            line[:30] + f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}" + line[54:]
        )
    # Keep non-coordinate metadata, then emit only the selected chain's atoms.
    header = [line for line in lines if not _atom_record(line) and not line.startswith(("TER", "END"))]
    return "\n".join(header + transformed + ["TER", "END"])


def plot_superimposed_cartoons(target_pdb_str, query_pdb_str):
    view = py3Dmol.view(width=700, height=500)
    view.addModel(target_pdb_str, "pdb")
    view.addModel(query_pdb_str, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"colorscheme": "cyanCarbon"}})
    view.setStyle({"model": 1}, {"cartoon": {"colorscheme": "orangeCarbon"}})
    view.zoomTo()
    view.show()


def pairwise_cartoons(df, target_ix, query_pdb_str):
    """Display a full-atom target/query cartoon for one row of a DALI view.

    This optional plot requires the public Daliviewer service. A failed fetch
    or coordinate mismatch raises an actionable error instead of showing a
    misleading blank or misaligned cartoon.
    """
    row = df.iloc[target_ix]
    target_id = row["target_id"]
    _, chain = _target_parts(target_id)
    url = generate_daliviewer_url(target_id)
    try:
        with urllib.request.urlopen(url, timeout=_REQUEST_TIMEOUT_SECONDS) as response:
            html_text = response.read().decode("utf-8", errors="replace")
    except (OSError, TimeoutError) as exc:
        raise RuntimeError(
            f"Could not fetch target {target_id!r} from Daliviewer. "
            "Check the network connection or try another target."
        ) from exc

    raw_target_pdb = extract_pdb_from_html(html_text)
    try:
        transformed_target_pdb = transform_pdb(
            raw_target_pdb, row["R"], row["t"], chain=chain,
            expected_ca_count=row.get("target_length"),
        )
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Could not render target {target_id!r}: {exc} "
            "Try another target or use the pack-based C-alpha superimposition viewer."
        ) from exc
    plot_superimposed_cartoons(transformed_target_pdb, query_pdb_str)
