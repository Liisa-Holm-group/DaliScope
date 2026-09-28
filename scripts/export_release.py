"""Export the staged release files without Git history, caches, or local outputs."""
import argparse
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    destination = args.destination.resolve()
    if destination == root or root in destination.parents:
        parser.error("destination must be outside the source repository")
    if destination.exists() and any(destination.iterdir()):
        parser.error("destination must be absent or empty")
    paths = subprocess.check_output(["git", "ls-files", "--cached", "-z"], cwd=root).decode("utf-8").split("\0")
    destination.mkdir(parents=True, exist_ok=True)
    count = 0
    for name in filter(None, paths):
        source = (root / name).resolve()
        target = (destination / name).resolve()
        if root not in source.parents or destination not in target.parents:
            raise ValueError(f"unexpected path in Git index: {name}")
        if not source.is_file():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        count += 1
    print(f"Exported {count} files to {destination}")


if __name__ == "__main__":
    main()
