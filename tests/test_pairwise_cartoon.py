"""Scientific alignment and failure-mode checks for the optional full-atom cartoon."""

from urllib.error import URLError
from urllib.parse import parse_qs, urlsplit

import numpy as np
import pandas as pd
import pytest

from daliscope.optics import pairwise_cartoon


def _atom(serial, chain, residue, point, atom="CA"):
    x, y, z = point
    return (
        f"ATOM  {serial:5d} {atom:^4s} ALA {chain}{residue:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
    )


def _xyz(pdb_text, chain):
    return np.array([
        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        for line in pdb_text.splitlines()
        if line.startswith("ATOM  ") and line[21] == chain and line[12:16].strip() == "CA"
    ])


def _two_chain_pdb():
    lines = [
        _atom(1, "A", 1, (-100, 10, 0)),
        _atom(2, "A", 2, (-99, 10, 0)),
        _atom(3, "B", 1, (20, 30, 40)),
        _atom(4, "B", 2, (21, 30, 40)),
        _atom(5, "B", 3, (20, 32, 40)),
        _atom(6, "B", 3, (20.5, 32.5, 40), atom="CB"),
        "TER", "END",
    ]
    return "\n".join(lines)


def test_transform_uses_selected_chain_centroid_and_pack_rotation():
    remote = _two_chain_pdb()
    raw_b = _xyz(remote, "B")
    centered_pack_b = raw_b - raw_b.mean(axis=0)
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    t = np.array([4.0, 5.0, 6.0])

    result = pairwise_cartoon.transform_pdb(remote, R, t, chain="B", expected_ca_count=3)

    assert not any(line.startswith("ATOM  ") and line[21] == "A" for line in result.splitlines())
    np.testing.assert_allclose(_xyz(result, "B"), centered_pack_b @ R.T + t, atol=0.001)
    assert any(line[12:16].strip() == "CB" for line in result.splitlines() if line.startswith("ATOM  "))


def test_url_supports_afdb_pseudo_id_and_preserves_chain_in_caller():
    parsed = urlsplit(pairwise_cartoon.generate_daliviewer_url("v9vaA"))
    params = parse_qs(parsed.query)
    assert params["jobid"] == ["pdb"]
    assert params["pdbid"] == ["v9va"]
    assert params["html"] == ["1"]
    assert pairwise_cartoon._target_parts("5nnlB") == ("5nnl", "B")
    with pytest.raises(ValueError, match="four-character DALI structure code"):
        pairwise_cartoon.generate_daliviewer_url("AF-Q9H3P7-F1")


def test_pairwise_cartoons_selects_target_chain_and_displays(monkeypatch):
    remote = _two_chain_pdb()
    seen = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return ("<PRE>\n" + remote + "\n</PRE>").encode()

    def fake_urlopen(url, *, timeout):
        seen.append((url, timeout))
        return Response()

    monkeypatch.setattr(pairwise_cartoon.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(pairwise_cartoon, "plot_superimposed_cartoons",
                        lambda target, query: seen.append((target, query)))
    view = pd.DataFrame([{"target_id": "5nnlB", "target_length": 3,
                          "R": np.eye(3), "t": np.zeros(3)}])

    pairwise_cartoon.pairwise_cartoons(view, 0, "QUERY PDB")

    assert seen[0][1] == pairwise_cartoon._REQUEST_TIMEOUT_SECONDS
    assert parse_qs(urlsplit(seen[0][0]).query)["pdbid"] == ["5nnl"]
    assert seen[1][1] == "QUERY PDB"
    assert _xyz(seen[1][0], "A").size == 0
    assert len(_xyz(seen[1][0], "B")) == 3


def test_unavailable_or_wrong_structure_fails_clearly(monkeypatch):
    view = pd.DataFrame([{"target_id": "zzzzA", "target_length": 3,
                          "R": np.eye(3), "t": np.zeros(3)}])

    class EmptyResponse:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            # readpdb.cgi replies HTTP 200 with only REMARK lines for unknown IDs.
            return b"<PRE>REMARK  None\nTER\nEND</PRE>"

    monkeypatch.setattr(pairwise_cartoon.urllib.request, "urlopen",
                        lambda url, *, timeout: EmptyResponse())
    with pytest.raises(RuntimeError, match="no chain 'A'.*pack-based C-alpha"):
        pairwise_cartoon.pairwise_cartoons(view, 0, "QUERY PDB")

    def fail_request(url, *, timeout):
        raise URLError("network unavailable")

    monkeypatch.setattr(pairwise_cartoon.urllib.request, "urlopen", fail_request)
    with pytest.raises(RuntimeError, match="Could not fetch target 'zzzzA'.*network connection"):
        pairwise_cartoon.pairwise_cartoons(view, 0, "QUERY PDB")


def test_ca_count_mismatch_prevents_misleading_superimposition():
    with pytest.raises(ValueError, match="DALI pack expects 4"):
        pairwise_cartoon.transform_pdb(_two_chain_pdb(), np.eye(3), np.zeros(3),
                                       chain="B", expected_ca_count=4)
