"""Download the optional GOLD example and verify its recorded checksum."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=("GOLD",))
    parser.add_argument("--tag", default="v0.1.1", help="GitHub release tag")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    item = json.loads((root / "notebooks/data/datasets.json").read_text(encoding="utf-8"))[args.dataset]
    output = root / "notebooks/data" / item["file"]
    if output.exists():
        if hashlib.sha256(output.read_bytes()).hexdigest() == item["sha256"]:
            print(f"Already present and verified: {output}")
            return
        parser.error(f"existing file has a different checksum: {output}; move it aside before downloading")
    url = f"https://github.com/Liisa-Holm-group/DaliScope/releases/download/{args.tag}/{item['file']}"
    temp_path = None
    try:
        with urlopen(url, timeout=60) as response, tempfile.NamedTemporaryFile(dir=output.parent, delete=False) as temp:
            temp_path = Path(temp.name)
            digest = hashlib.sha256()
            while chunk := response.read(1024 * 1024):
                digest.update(chunk)
                temp.write(chunk)
        if digest.hexdigest() != item["sha256"]:
            parser.error("download checksum does not match the frozen example")
        temp_path.replace(output)
    except HTTPError as error:
        parser.error(f"release asset is unavailable ({error.code}); publish the release first or use the supplied local pack")
    except (URLError, TimeoutError) as error:
        parser.error(f"download failed: {error}")
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()
    print(f"Downloaded and verified: {output}")


if __name__ == "__main__":
    main()
