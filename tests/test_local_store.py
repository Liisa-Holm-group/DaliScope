from contextlib import closing
from pathlib import Path
import sqlite3
import zipfile

import numpy as np
import pytest

from daliscope.core.project import Project
from daliscope.packer.build_store import build_store
from daliscope.packer.make_colab_pack import main as pack


@pytest.fixture
def example(tmp_path):
    archive = Path(__file__).resolve().parents[1] / "examples" / "local-pack" / "inputs.zip"
    with zipfile.ZipFile(archive) as inputs:
        inputs.extractall(tmp_path)
    return tmp_path


def test_rebuild_public_example(example):
    output = example / "rebuilt"
    assert build_store(example / "results.tsv", 2, example / "targets.tsv", output) == 4
    for store in (example / "store", output):
        with closing(sqlite3.connect(store / "master_ledger.db")) as conn:
            assert conn.execute("SELECT count(*), sum(n_res) FROM proteins WHERE ns_id=2").fetchone() == (4, 1965)
    name = "shards/shard-00000.npz"
    with np.load(output / name) as rebuilt, np.load(example / "store" / name) as bundled:
        assert set(rebuilt.files) == set(bundled.files)
        for field in bundled.files:
            if field.endswith("_coords"):
                np.testing.assert_allclose(rebuilt[field], bundled[field], atol=0.0001)
            else:
                np.testing.assert_array_equal(rebuilt[field], bundled[field])
    assert pack([str(example / "results.tsv"), "AFDB2", "example", str(example / "query.pdb.gz"),
                 "A", str(example / "query.dat"), "--ledger-db", str(output / "master_ledger.db"),
                 "--shard-dir", str(output / "shards"), "--pfam-db", str(example / "pfam-39.0.db"),
                 "--output-dir", str(example)]) == 0
    project = Project.load_pack(str(example / "example.tar.gz"))
    assert project.query_length == 435
    assert (len(project.coords), len(project.segments), len(project.master_metadata)) == (1965, 63, 4)
    assert len(project._pfam) == 7 and project._pfam["target_name"].nunique() == 4


def test_custom_ids_without_annotations(example):
    results = example / "local.tsv"
    results.write_text("target_id\tnamespace\tz_score\tqstarts\tsstarts\tlengths\n"
                       "001\tLOCAL\t80\t[1]\t[1]\t[435]\nNA\tLOCAL\t80\t[1]\t[1]\t[435]\n",
                       encoding="utf-8")
    manifest = example / "local-targets.tsv"
    manifest.write_text("target_id\tdat\n001\tquery.dat\nNA\tquery.dat\nother\tunused.dat\n",
                        encoding="utf-8")
    store = example / "local"
    assert build_store(results, 5, manifest, store) == 2
    assert pack([str(results), "LOCAL", "private", str(example / "query.pdb.gz"), "A",
                 str(example / "query.dat"), "--ledger-db", str(store / "master_ledger.db"),
                 "--shard-dir", str(store / "shards"), "--output-dir", str(example)]) == 0
    project = Project.load_pack(str(example / "private.tar.gz"))
    view = next(iter(project.views.values()))
    assert set(view["target_id"]) == {"001", "NA"}
    assert view["pfam"].eq("Unassigned").all()


def test_incomplete_or_mismatched_inputs_leave_no_store(example):
    results = example / "results.tsv"
    manifest = example / "targets.tsv"
    original = manifest.read_text(encoding="utf-8")
    output = example / "new"
    with pytest.raises(FileExistsError):
        build_store(results, 2, manifest, example / "store")
    manifest.write_text(original + original.splitlines()[1] + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate manifest entry"):
        build_store(results, 2, manifest, output)
    manifest.write_text("\n".join(original.splitlines()[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing 1 manifest"):
        build_store(results, 2, manifest, output)
    manifest.write_text(original.replace("lpa4A.pdb.gz", "eszmA.pdb.gz"), encoding="utf-8")
    with pytest.raises(ValueError, match="lpa4A: model and .dat"):
        build_store(results, 2, manifest, output)
    assert not output.exists()
