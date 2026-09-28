"""Set up a versioned DaliScope checkout inside a Google Colab runtime."""
from importlib import metadata
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def runtime_constraints():
    """Keep Colab's declared dependencies and already-loaded numerical libraries."""
    requirements = list(metadata.requires("google-colab") or [])
    for name in ("numpy", "pandas", "scipy", "matplotlib"):
        try:
            requirements.append(f"{name}=={metadata.version(name)}")
        except metadata.PackageNotFoundError:
            pass
    return list(dict.fromkeys(requirements))


def verify_checkout(root, release):
    """Fail safely if a previous tutorial left a different or edited checkout."""
    try:
        head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        expected = subprocess.check_output(["git", "rev-parse", f"{release}^{{commit}}"],
                                           cwd=root, text=True, stderr=subprocess.DEVNULL).strip()
    except subprocess.CalledProcessError as error:
        raise RuntimeError("This runtime contains a different DaliScope version. Start a fresh runtime and run setup again.") from error
    if head != expected:
        raise RuntimeError("This runtime contains a different DaliScope version. Start a fresh runtime and run setup again.")
    if subprocess.run(["git", "diff", "--quiet", "HEAD", "--"], cwd=root).returncode:
        raise RuntimeError("DaliScope has modified tracked files. Save your edits and start a fresh runtime.")


def setup_colab(root, release):
    from google.colab import output

    root = Path(root).resolve()
    verify_checkout(root, release)
    loaded = sys.modules.get("daliscope")
    if loaded is not None and loaded.__version__ != release.removeprefix("v"):
        raise RuntimeError("Another DaliScope version is already imported. Restart the runtime and run setup again.")
    for name, module in tuple(sys.modules.items()):
        if name == "daliscope" or name.startswith("daliscope."):
            source = getattr(module, "__file__", None)
            if source is None or not Path(source).resolve().is_relative_to(root / "daliscope"):
                raise RuntimeError("Another DaliScope checkout is already imported. Restart the runtime and run setup again.")
    marker = f"{os.getpid()}:{sys.executable}:{root}:{release}"
    if os.environ.get("DALISCOPE_COLAB_READY") != marker:
        print(f"Installing DaliScope {release} in the Colab runtime...", flush=True)
        with tempfile.TemporaryDirectory(prefix="daliscope-colab-") as work:
            constraints = Path(work) / "constraints.txt"
            constraints.write_text("\n".join(runtime_constraints()) + "\n", encoding="utf-8")
            subprocess.run([sys.executable, "-m", "pip", "install", "--constraint", str(constraints),
                            "-e", f"{root}[colab]"], check=True)
        os.environ["DALISCOPE_COLAB_READY"] = marker
    # Editable-install .pth files are only loaded when Python starts; this kernel is already running.
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    os.chdir(root)
    output.enable_custom_widget_manager()
    print(f"Ready: DaliScope {release}. Continue with the next cell.", flush=True)
