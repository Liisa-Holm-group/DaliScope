"""The CLI must honor cell timeouts and clean up its kernel."""
from pathlib import Path
import subprocess
import sys

import nbformat

ROOT = Path(__file__).resolve().parents[1]


def test_cell_timeout_is_enforced(tmp_path):
    path = tmp_path / f"timeout-guard-{tmp_path.name}.ipynb"
    notebook = nbformat.v4.new_notebook(cells=[
        nbformat.v4.new_code_cell("import time; time.sleep(3)")
    ])
    nbformat.write(notebook, path)
    output = ROOT / ".local-work/notebooks-executed" / path.name
    try:
        result = subprocess.run([sys.executable, str(ROOT / "scripts/validate_notebooks.py"),
                                 str(path), "--timeout", "1"],
                                capture_output=True, text=True, timeout=30)
        assert result.returncode != 0, "A cell exceeded its timeout without failing"
        assert "CellTimeoutError" in result.stderr
    finally:
        output.unlink(missing_ok=True)
