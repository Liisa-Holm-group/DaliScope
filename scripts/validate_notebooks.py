"""Execute selected notebooks in a fresh kernel and save local validation outputs."""
import argparse
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebooks", nargs="+", type=Path)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = root / ".local-work/notebooks-executed"
    output.mkdir(parents=True, exist_ok=True)
    for path in args.notebooks:
        path = path.resolve()
        notebook = nbformat.read(path, as_version=4)
        manager = KernelManager(kernel_name="python3")
        # Use this interpreter, independently of any user-installed kernel spec.
        manager.kernel_spec.argv = [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"]
        client = NotebookClient(notebook, km=manager, timeout=args.timeout,
                                resources={"metadata": {"path": str(path.parent)}})
        try:
            client.execute()
        finally:
            nbformat.write(notebook, output / path.name)
            if manager.has_kernel:
                manager.shutdown_kernel(now=True)
        print(f"PASS {path.name}: {sum(c.cell_type == 'code' for c in notebook.cells)} code cells", flush=True)


if __name__ == "__main__":
    main()
