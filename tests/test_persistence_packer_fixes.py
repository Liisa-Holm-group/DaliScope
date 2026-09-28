"""Regression checks for session file collisions and invalid query pack output."""
from pathlib import Path
import json
import os
import pickle
import sqlite3
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from daliscope.core.project import Project

ROOT = Path(__file__).resolve().parents[1]


def test_session_roundtrip_preserves_colliding_display_names(tmp_path):
    project = Project()
    first = pd.DataFrame({"value": [1], "coords": [np.array([1, 2, 3])]})
    second = pd.DataFrame({"value": [2], "coords": [np.array([4, 5, 6])]})
    assert project.add_view("domain 1", first)
    assert project.add_view("domain_1", second)
    project.add_model("model 1", {"value": 1})
    project.add_model("model_1", {"value": 2})

    session = project.save(tmp_path / "session")
    restored = Project.open(session)
    pd.testing.assert_frame_equal(restored.views["domain 1"], first)
    pd.testing.assert_frame_equal(restored.views["domain_1"], second)
    assert restored.models == {"model 1": {"value": 1}, "model_1": {"value": 2}}


def test_open_accepts_legacy_010_storage_filenames(tmp_path):
    session = tmp_path / "legacy"
    (session / "views").mkdir(parents=True)
    (session / "models").mkdir()
    expected = pd.DataFrame({"value": [1], "coords": [np.array([1, 2, 3])]})
    expected[["value"]].to_parquet(session / "views/domain_1.parquet", index=True)
    with (session / "views/domain_1_arrays.pkl").open("wb") as handle:
        pickle.dump({"coords": expected["coords"].tolist()}, handle)
    with (session / "models/model_1.pkl").open("wb") as handle:
        pickle.dump({"value": 1}, handle)
    manifest = {
        "saved_at": "2026-09-28T12:00:00",
        "daliscope_version": "0.1.0",
        "pack_path": None,
        "views": {"domain 1": {
            "parquet": "domain_1.parquet", "arrays": "domain_1_arrays.pkl",
            "array_cols": ["coords"], "row_count": 1, "provenance": {},
        }},
        "models": {"model 1": {"pkl": "model_1.pkl", "provenance": {}}},
        "plots": {},
    }
    (session / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    restored = Project.open(session)
    pd.testing.assert_frame_equal(restored.views["domain 1"], expected)
    assert restored.models["model 1"] == {"value": 1}


@pytest.fixture
def generator_inputs(tmp_path):
    shards = tmp_path / "shards"
    shards.mkdir()
    coords = np.array([[0, 0, 0], [3, 0, 0], [6, 1, 0], [9, 1, 1], [12, 0, 1]], dtype=np.float32)
    np.savez(shards / "tiny.npz", **{
        "1.1abcA_coords": coords,
        "1.1abcA_plddt": np.full(5, 90, dtype=np.uint8),
        "1.1abcA_seq": np.frombuffer(b"AAAAA", dtype=np.uint8),
        "1.1abcA_dssp": np.frombuffer(b"HHHHH", dtype=np.uint8),
        "1.1abcA_compnd": np.array("Synthetic target"),
    })
    ledger = tmp_path / "ledger.db"
    conn = sqlite3.connect(ledger)
    try:
        conn.execute("CREATE TABLE proteins (dali_id TEXT, n_res INT, umid TEXT, shard_path TEXT, cx REAL, cy REAL, cz REAL, ns_id INT)")
        conn.execute("INSERT INTO proteins VALUES ('1abcA', 5, '1.1abcA', 'tiny.npz', 0, 0, 0, 1)")
        conn.commit()
    finally:
        conn.close()
    pfam = tmp_path / "pfam.db"
    conn = sqlite3.connect(pfam)
    try:
        conn.execute("CREATE TABLE hmmer_hits (target_name TEXT, query_accession TEXT, e_value REAL, env_from INT, env_to INT, database TEXT)")
        conn.execute("CREATE TABLE pfam_entries (AC TEXT, CL TEXT)")
        conn.commit()
    finally:
        conn.close()
    names = tmp_path / "pfam_names.tsv"
    names.write_text("pfam\tclan\tclan_short\tshort\tname\n", encoding="utf-8")
    dali = tmp_path / "results.tsv"
    dali.write_text("target_id\tnamespace\tz_score\tqstarts\tsstarts\tlengths\n1abcA\tPDB\t10\t[1]\t[1]\t[5]\n", encoding="utf-8")
    pdb = tmp_path / "query.pdb"
    lines = []
    serial = 0
    for residue, point in enumerate(coords, 1):
        for atom in ("N", "CA", "C", "O"):
            serial += 1
            x, y, z = point
            lines.append(f"ATOM  {serial:5d} {atom:^4s} ALA A{residue:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 90.00          {atom[0]:>2s}\n")
    pdb.write_text("".join(lines) + "END\n", encoding="utf-8")
    dat = tmp_path / "query.dat"
    dat.write_text('>>>> testA 5\n-sequence "AAAAA"\n-dssp "HHHHH"\n-compnd "Synthetic query"\n', encoding="utf-8")
    return {
        "tsv": dali, "pdb": pdb, "dat": dat,
        "ledger": ledger, "shards": shards, "pfam": pfam, "names": names,
    }


@pytest.mark.parametrize("chain,metadata_length", [("B", 5), ("A", 6)])
def test_pack_cli_rejects_invalid_query_before_overwriting(generator_inputs, tmp_path, chain, metadata_length):
    files = generator_inputs
    if metadata_length != 5:
        files["dat"].write_text('>>>> testA 6\n-sequence "AAAAAA"\n-dssp "HHHHHH"\n-compnd "Synthetic query"\n', encoding="utf-8")
    output = tmp_path / "output"
    output.mkdir()
    destination = output / "test.tar.gz"
    previous = b"previous output must remain unchanged"
    destination.write_bytes(previous)
    command = [
        sys.executable, "-m", "daliscope.packer.make_colab_pack",
        str(files["tsv"]), "PDB", "test", str(files["pdb"]), chain, str(files["dat"]),
        "--ledger-db", str(files["ledger"]), "--shard-dir", str(files["shards"]),
        "--pfam-db", str(files["pfam"]), "--pfam-names", str(files["names"]),
        "--output-dir", str(output), "--overwrite",
    ]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                            env={**os.environ, "PYTHONUTF8": "1"})
    assert result.returncode != 0, result.stdout
    assert "Query.pdb" in result.stderr
    assert destination.read_bytes() == previous
    assert not list(output.glob("daliscope-pack-*"))
