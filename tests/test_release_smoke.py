from pathlib import Path
import io
import sqlite3
import tarfile

import numpy as np
import pandas as pd
import pytest

from daliscope.core.project import Project, initialize_project
from daliscope.packer.pfam_query import extract_hits
from daliscope.mechanics.partitioning import dp_optimal_partition

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "notebooks" / "data" / "3ubpC_PDB25.tar.gz"


@pytest.fixture(scope="module")
def project():
    project = Project.load_pack(str(DEMO))
    name = next(iter(project.views))
    project.add_subset("z8", project.views[name]["z_score"] >= 8, name)
    return project


def test_default_view_keeps_all_aligned_query_positions(project):
    name = next(iter(project.views))
    expected = project.segments.groupby("alignment_id")["length"].sum()
    actual = project.views[name].set_index("alignment_id")["alignment_length"]
    pd.testing.assert_series_equal(actual.sort_index(), expected.loc[actual.index].sort_index(), check_names=False, check_dtype=False)
    assert project.query_length == 570
    assert project.query_ca_coords.shape == (570, 3)
    assert len(project.views[name]) > 1000


def test_invalid_initialize_path_has_consistent_return_shape(tmp_path):
    assert initialize_project(str(tmp_path / "missing.tar.gz")) == (None, None, None)


def test_missing_archive_fails_clearly(tmp_path):
    with pytest.raises(FileNotFoundError):
        Project.load_pack(str(tmp_path / "missing.tar.gz"))


def test_session_roundtrip_preserves_views_and_pack_reference(project, tmp_path):
    destination = project.save(str(tmp_path / "session"))
    restored = Project.open(str(destination))
    assert restored.pack_path == project.pack_path
    assert restored.query_length == project.query_length
    assert restored.views.keys() == project.views.keys()
    for name, df in project.views.items():
        pd.testing.assert_frame_equal(restored.views[name], df, check_dtype=False)
    second = restored.save(str(tmp_path / "session-again"))
    assert Project.open(str(second)).pack_path == project.pack_path


def test_session_overwrite_requires_explicit_option(project, tmp_path):
    destination = project.save(str(tmp_path / "session"))
    with pytest.raises(FileExistsError):
        project.save(str(destination))


def test_pack_without_optional_annotations_loads(tmp_path):
    destination = tmp_path / "no-annotations.tar.gz"
    with tarfile.open(DEMO) as original, tarfile.open(destination, "w:gz") as output:
        for member in original.getmembers():
            if member.isfile() and member.name not in {"pfam.tsv", "mini_pfam_names.tsv", "mini_ligand.tsv", "mini_ligand_names.tsv"}:
                output.addfile(member, original.extractfile(member))
    loaded = Project.load_pack(str(destination))
    df = next(iter(loaded.views.values()))
    assert len(df) > 1000
    assert df["pfam"].eq("Unassigned").all()


@pytest.mark.parametrize("namespace,label", [("1", "PDB"), (4, "VIRO3D")])
def test_pfam_database_filtering_uses_numeric_namespace(tmp_path, namespace, label):
    database = tmp_path / "annotations.db"
    with sqlite3.connect(database) as conn:
        conn.execute("CREATE TABLE hmmer_hits (target_name TEXT, query_accession TEXT, e_value REAL, env_from INT, env_to INT, database TEXT)")
        conn.execute("CREATE TABLE pfam_entries (AC TEXT, CL TEXT)")
        conn.execute("INSERT INTO pfam_entries VALUES ('PF00001', 'CL0001')")
        conn.executemany("INSERT INTO hmmer_hits VALUES (?, 'PF00001', 1e-8, 1, 10, ?)", [("hitA", label), ("unwantedA", label), ("hitA", "OTHER")])
    output = tmp_path / "pfam.tsv"
    extract_hits(database, ["hitA", "hitA"], namespace, output)
    result = pd.read_csv(output, sep="\t")
    assert len(result) == 1
    assert result.loc[0, "database"] == label
    assert result.loc[0, "clan"] == "CL0001"


def test_missing_pfam_database_does_not_create_an_empty_database(tmp_path):
    database = tmp_path / "missing.db"
    with pytest.raises(FileNotFoundError):
        extract_hits(database, ["hitA"], 1, tmp_path / "pfam.tsv")
    assert not database.exists()


def test_partition_accepts_weighted_edges_without_external_module():
    hierarchy = {0: (1, 1), 1: (1, 2), 2: (2, 1), 3: (2, 2)}
    edges = [(0, 1, 4.0), (2, 3, 4.0), (1, 2, 0.1)]
    matrix = np.zeros((4, 4))
    for u, v, weight in edges:
        matrix[u, v] = matrix[v, u] = weight
    expected = dp_optimal_partition(hierarchy, matrix)
    actual = dp_optimal_partition(hierarchy, edges)
    assert actual == expected

import subprocess
import sys


def test_generator_end_to_end_with_local_database_fixture(tmp_path):
    from daliscope.packer.make_colab_pack import main
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    coords = np.array([[0, 0, 0], [3, 0, 0], [6, 1, 0], [9, 1, 1], [12, 0, 1]], dtype=np.float32)
    np.savez(shard_dir / "tiny.npz", **{
        "1.1abcA_coords": coords, "1.1abcA_plddt": np.full(5, 90, dtype=np.uint8),
        "1.1abcA_seq": np.frombuffer(b"AAAAA", dtype=np.uint8),
        "1.1abcA_dssp": np.frombuffer(b"HHHHH", dtype=np.uint8),
        "1.1abcA_compnd": np.array("Synthetic test target"),
    })
    ledger = tmp_path / "ledger.db"
    with sqlite3.connect(ledger) as conn:
        conn.execute("CREATE TABLE proteins (dali_id TEXT, n_res INT, umid TEXT, shard_path TEXT, cx REAL, cy REAL, cz REAL, ns_id INT)")
        conn.execute("INSERT INTO proteins VALUES ('1abcA', 5, '1.1abcA', 'tiny.npz', 0, 0, 0, 1)")
    pfam_db = tmp_path / "pfam.db"
    with sqlite3.connect(pfam_db) as conn:
        conn.execute("CREATE TABLE hmmer_hits (target_name TEXT, query_accession TEXT, e_value REAL, env_from INT, env_to INT, database TEXT)")
        conn.execute("CREATE TABLE pfam_entries (AC TEXT, CL TEXT)")
        conn.execute("INSERT INTO hmmer_hits VALUES ('1abcA', 'PF00001', 1e-10, 1, 5, 'PDB')")
        conn.execute("INSERT INTO pfam_entries VALUES ('PF00001', 'CL0001')")
    names = tmp_path / "pfam_names.tsv"
    names.write_text("pfam\tclan\tclan_short\tshort\tname\nPF00001\tCL0001\tExample\tExample\tSynthetic family\n", encoding="utf-8")
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
    output = tmp_path / "output"
    arguments = [str(dali), "PDB", "test", str(pdb), "A", str(dat),
                 "--ledger-db", str(ledger), "--shard-dir", str(shard_dir),
                 "--pfam-db", str(pfam_db), "--pfam-names", str(names),
                 "--output-dir", str(output)]
    assert main(arguments) == 0
    archive = output / "test.tar.gz"
    with tarfile.open(archive) as pack:
        assert all("/" not in m.name for m in pack.getmembers())
        assert "mini_pfam_names.tsv" in pack.getnames()
    loaded = Project.load_pack(str(archive))
    view = next(iter(loaded.views.values()))
    assert loaded.query_length == 5
    assert len(view) == 1
    assert view.iloc[0]["alignment_length"] == 5
    assert view.iloc[0]["pfam"] == "PF00001"
    assert not list(output.glob("daliscope-pack-*"))
    with pytest.raises(SystemExit):
        main(arguments)


def test_pack_cli_help_is_available_without_database_access():
    result = subprocess.run([sys.executable, "-m", "daliscope.packer.make_colab_pack", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "--ledger-db" in result.stdout
    assert "--upload-to" in result.stdout


def test_bounded_structural_community_pipeline(project):
    from daliscope.analysis.communities import run_community_detection
    full = next(iter(project.views))
    selected = project.views[full].sort_values("z_score", ascending=False).head(20).copy()
    project.add_view("community_smoke", selected, parent=full)
    modules, count = run_community_detection(
        "community_smoke", (1, 150), project, k_neighbors=5,
        infomap_args="--two-level --num-trials 2 --seed 42 --silent", verbose=False,
    )
    assert count >= 1
    assert set(modules["target_id"]) == set(selected["target_id"])
    assert modules["module"].notna().all()


def test_missing_query_structure_has_a_clear_error(tmp_path):
    destination = tmp_path / "no-query.tar.gz"
    with tarfile.open(DEMO) as original, tarfile.open(destination, "w:gz") as output:
        for member in original.getmembers():
            if member.isfile() and member.name != "Query.pdb":
                output.addfile(member, original.extractfile(member))
    with pytest.raises(FileNotFoundError, match="Query.pdb"):
        Project.load_pack(str(destination))


def test_empty_pfam_annotation_table_is_supported(tmp_path):
    destination = tmp_path / "empty-pfam.tar.gz"
    with tarfile.open(DEMO) as original, tarfile.open(destination, "w:gz") as output:
        for member in original.getmembers():
            if not member.isfile():
                continue
            if member.name == "pfam.tsv":
                header = original.extractfile(member).readline()
                info = tarfile.TarInfo("pfam.tsv")
                info.size = len(header)
                output.addfile(info, io.BytesIO(header))
            else:
                output.addfile(member, original.extractfile(member))
    loaded = Project.load_pack(str(destination))
    assert next(iter(loaded.views.values()))["pfam"].eq("Unassigned").all()
