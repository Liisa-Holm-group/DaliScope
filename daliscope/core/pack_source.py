"""Resolve a local pack path or download an HTTP(S) pack without extracting it."""
import gzip
import hashlib
from http.client import IncompleteRead
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
import zlib
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import urlopen


_REQUIRED_MEMBERS = frozenset({
    "query_meta.tsv", "summary.tsv", "id_mapping.tsv", "segments.tsv",
    "colab_pack_index.tsv", "colab_pack.npz", "Query.pdb",
})
_CHUNK_SIZE = 1024 * 1024


def _validate_download(path: Path) -> None:
    """Check archive structure and gzip completeness; Project validates its data."""
    try:
        with gzip.open(path, "rb") as compressed:
            with tarfile.open(fileobj=compressed, mode="r|") as archive:
                present = {
                    PurePosixPath(member.name).name
                    for member in archive if member.isfile()
                }
            # Tar reading can stop before the gzip trailer; consume it to check CRC
            # and detect a truncated download even when its tar headers are intact.
            while compressed.read(_CHUNK_SIZE):
                pass
    except (OSError, EOFError, tarfile.TarError, zlib.error) as error:
        raise ValueError(
            "Downloaded file is not a complete DaliScope pack (gzip-tar). "
            "The link may have expired or returned an HTML page."
        ) from error
    missing = _REQUIRED_MEMBERS - present
    if missing:
        raise ValueError(
            "Downloaded archive is not a DaliScope pack; missing required files: "
            + ", ".join(sorted(missing))
        )


def resolve_pack_source(
    source: str | Path,
    download_dir: str | Path = "data/downloaded",
    timeout: float = 120,
) -> Path:
    """Return a local pack path, downloading HTTP(S) sources when needed.

    Local paths keep their original meaning, including relative and Windows paths.
    Remote archives are checked for required contents and cached by their SHA-256;
    no members are extracted. The printed digest records the downloaded bytes and
    does not imply comparison with an upstream checksum.
    """
    source_text = str(source)
    parsed = urlsplit(source_text)
    if parsed.scheme.lower() not in {"http", "https"}:
        return Path(source)
    if not parsed.netloc:
        raise ValueError("A pack URL must include an HTTP or HTTPS server name.")

    cache = Path(download_dir)
    cache.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        # Keep the original URL, including its query string, for job-specific URLs.
        with urlopen(source_text, timeout=timeout) as response:
            content_length = response.headers.get("Content-Length")
            expected = int(content_length) if content_length and content_length.isdigit() else None
            digest = hashlib.sha256()
            received = 0
            with tempfile.NamedTemporaryFile(
                dir=cache, prefix=".daliscope-", suffix=".part", delete=False
            ) as output:
                temporary = Path(output.name)
                while chunk := response.read(_CHUNK_SIZE):
                    output.write(chunk)
                    digest.update(chunk)
                    received += len(chunk)
        if expected is not None and received != expected:
            raise ValueError(
                f"Incomplete DaliScope pack download: received {received} bytes, expected {expected}."
            )
        _validate_download(temporary)
        checksum = digest.hexdigest()
        destination = cache / f"{checksum}.tar.gz"
        # Publish only a complete, validated download. Identical content has the
        # same name; different jobs with the same URL basename cannot collide.
        temporary.replace(destination)
        print(f"Pack file: {destination}\nSHA-256 (computed): {checksum}")
        return destination
    except HTTPError as error:
        raise ValueError(
            f"Could not download DaliScope pack (HTTP {error.code}). "
            "The result link may have expired; check the URL or use a local pack."
        ) from error
    except (URLError, TimeoutError, OSError, IncompleteRead) as error:
        raise ValueError(
            "Could not download or save a complete DaliScope pack. "
            "Check the URL or download the file and use its local path."
        ) from error
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
