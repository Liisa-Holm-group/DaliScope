"""Guard the hosted bootstrap against dependency changes and mixed checkouts."""
import importlib.util
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_setup():
    path = ROOT / "scripts/colab_setup.py"
    assert path.is_file(), "Colab bootstrap is missing"
    spec = importlib.util.spec_from_file_location("colab_setup", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_constraints_preserve_managed_and_loaded_packages(monkeypatch):
    setup = load_setup()
    versions = {"numpy": "2.0.2", "pandas": "2.2.3", "scipy": "1.15.3", "matplotlib": "3.10.0"}
    monkeypatch.setattr(setup.metadata, "requires", lambda name: ["ipython==7.34.0", "ipykernel==6.17.1", "pandas==2.2.3"])
    monkeypatch.setattr(setup.metadata, "version", lambda name: versions[name])
    constraints = setup.runtime_constraints()
    assert "ipython==7.34.0" in constraints
    assert "ipykernel==6.17.1" in constraints
    for name, version in versions.items():
        assert f"{name}=={version}" in constraints
    assert len(constraints) == len(set(constraints))


def checkout(tmp_path):
    for args in (["init"], ["config", "user.name", "Bootstrap test"],
                 ["config", "user.email", "bootstrap@example.invalid"]):
        subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)
    (tmp_path / "tracked.txt").write_text("release")
    subprocess.run(["git", "add", "."], cwd=tmp_path, check=True)
    subprocess.run(["git", "commit", "-m", "release"], cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(["git", "tag", "v0.1.2"], cwd=tmp_path, check=True)
    return tmp_path


def test_refuses_reused_checkout_at_other_commit(tmp_path):
    setup = load_setup()
    repo = checkout(tmp_path)
    setup.verify_checkout(repo, "v0.1.2")
    (repo / "tracked.txt").write_text("different version")
    subprocess.run(["git", "commit", "-am", "different"], cwd=repo, check=True, capture_output=True)
    with pytest.raises(RuntimeError, match="different DaliScope version"):
        setup.verify_checkout(repo, "v0.1.2")


def test_refuses_modified_release_files(tmp_path):
    setup = load_setup()
    repo = checkout(tmp_path)
    (repo / "tracked.txt").write_text("local edits")
    with pytest.raises(RuntimeError, match="modified tracked files"):
        setup.verify_checkout(repo, "v0.1.2")


def test_comprehensive_notebook_preserves_colab_renderer():
    import json
    import plotly.io as pio

    notebook = json.loads((ROOT / "notebooks/RunDaliScope-1.ipynb").read_text(encoding="utf-8"))
    sources = ["".join(cell["source"]) for cell in notebook["cells"]
               if cell["cell_type"] == "code" and "pio.renderers.default" in "".join(cell["source"])]
    previous = pio.renderers.default
    try:
        for source in sources:
            exec("\n".join(line for line in source.splitlines() if not line.startswith("%")), {"IN_COLAB": True})
        assert pio.renderers.default == "colab"
    finally:
        pio.renderers.default = previous


@pytest.fixture
def colab_runtime(tmp_path, monkeypatch):
    """Keep Git real while isolating runtime state and external Colab services."""
    import sys
    from types import ModuleType, SimpleNamespace

    setup = load_setup()
    repo = tmp_path / "DaliScope"
    repo.mkdir()
    checkout(repo)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "path", sys.path.copy())
    # setenv records even an initially absent key, so the helper's write is undone.
    monkeypatch.setenv("DALISCOPE_COLAB_READY", "not-ready-for-this-test")
    for name in list(sys.modules):
        if name == "daliscope" or name.startswith("daliscope."):
            monkeypatch.delitem(sys.modules, name)

    manager_calls = []
    google = ModuleType("google")
    google.__path__ = []
    colab = ModuleType("google.colab")
    colab.__path__ = []
    output = ModuleType("google.colab.output")
    output.enable_custom_widget_manager = lambda: manager_calls.append(Path.cwd())
    colab.output = output
    google.colab = colab
    for module in (google, colab, output):
        monkeypatch.setitem(sys.modules, module.__name__, module)

    real_requires = setup.metadata.requires

    def requirements(name):
        if name == "google-colab":
            return ["ipython==7.34.0", "ipykernel==6.17.1", "pandas==2.2.3"]
        return real_requires(name)

    monkeypatch.setattr(setup.metadata, "requires", requirements)
    real_run = subprocess.run
    installations = []

    def run_with_fake_pip(command, *args, **kwargs):
        if tuple(command[1:4]) == ("-m", "pip", "install"):
            installations.append(tuple(command))
            return subprocess.CompletedProcess(command, 0)
        return real_run(command, *args, **kwargs)

    monkeypatch.setattr(setup.subprocess, "run", run_with_fake_pip)
    return SimpleNamespace(setup=setup, repo=repo, installations=installations,
                           manager_calls=manager_calls)


def imported_module(monkeypatch, name, path, version="0.1.2"):
    """Model a module already present in a running notebook kernel."""
    import sys
    from types import ModuleType

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f'__version__ = "{version}"\n', encoding="utf-8")
    module = ModuleType(name)
    module.__file__ = str(path)
    module.__version__ = version
    monkeypatch.setitem(sys.modules, name, module)
    return module


def test_repeated_setup_installs_once_but_reinitializes_notebook(colab_runtime, monkeypatch):
    runtime = colab_runtime
    root = runtime.repo.resolve()
    setup = runtime.setup
    setup.setup_colab(root, "v0.1.2")
    assert Path.cwd() == root
    assert setup.sys.path[0] == str(root)
    imported_module(monkeypatch, "daliscope", root / "daliscope" / "__init__.py")

    # A user may change directories or remove the path between tutorial runs.
    setup.os.chdir(root.parent)
    setup.sys.path.remove(str(root))
    setup.setup_colab(root, "v0.1.2")

    assert len(runtime.installations) == 1
    assert f"{root}[colab]" in runtime.installations[0]
    assert runtime.manager_calls == [root, root]
    assert Path.cwd() == root
    assert setup.sys.path[0] == str(root)
    assert setup.sys.path.count(str(root)) == 1


def test_refuses_already_imported_different_version_before_install(colab_runtime, monkeypatch):
    runtime = colab_runtime
    initial_cwd = Path.cwd()
    imported_module(monkeypatch, "daliscope", runtime.repo / "daliscope" / "__init__.py",
                    version="0.1.1")

    with pytest.raises(RuntimeError, match="already imported"):
        runtime.setup.setup_colab(runtime.repo, "v0.1.2")

    assert runtime.installations == []
    assert runtime.manager_calls == []
    assert Path.cwd() == initial_cwd
    assert str(runtime.repo.resolve()) not in runtime.setup.sys.path


@pytest.mark.parametrize("loaded_name", ["daliscope", "daliscope.core.project"])
def test_refuses_imported_modules_from_another_checkout(colab_runtime, monkeypatch, loaded_name):
    runtime = colab_runtime
    initial_cwd = Path.cwd()
    imported_module(monkeypatch, "daliscope", runtime.repo / "daliscope" / "__init__.py")
    other_root = runtime.repo.parent / "another-checkout"
    relative_file = "__init__.py" if loaded_name == "daliscope" else "core/project.py"
    imported_module(monkeypatch, loaded_name, other_root / "daliscope" / relative_file)

    with pytest.raises(RuntimeError):
        runtime.setup.setup_colab(runtime.repo, "v0.1.2")

    assert runtime.installations == []
    assert runtime.manager_calls == []
    assert Path.cwd() == initial_cwd
    assert str(runtime.repo.resolve()) not in runtime.setup.sys.path
