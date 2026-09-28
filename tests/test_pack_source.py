"""Exercise local and remote pack sources without depending on a public server."""
from contextlib import contextmanager
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib
import importlib.util
import io
from pathlib import Path
import tarfile
import threading
from urllib.parse import urlsplit

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "notebooks/data/3ubpC_PDB25.tar.gz"
REQUIRED = (
    "query_meta.tsv", "summary.tsv", "id_mapping.tsv", "segments.tsv",
    "colab_pack_index.tsv", "colab_pack.npz", "Query.pdb",
)


def resolver():
    assert importlib.util.find_spec("daliscope.core.pack_source") is not None, "Pack URL support is missing"
    return importlib.import_module("daliscope.core.pack_source").resolve_pack_source


def archive_bytes(*, omit=(), extra=()):
    """A small archive for transport tests, not a scientific input fixture."""
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name, payload in [(name, b"") for name in REQUIRED if name not in omit] + list(extra):
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    return output.getvalue()


@pytest.fixture
def serve():
    @contextmanager
    def start(routes):
        requests = []

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                requests.append(self.path)
                status, payload, headers = routes.get(self.path, routes.get(urlsplit(self.path).path, (404, b"Expired job", {})))
                self.send_response(status)
                self.send_header("Content-Length", str(headers.get("Content-Length", len(payload))))
                for name, value in headers.items():
                    if name != "Content-Length":
                        self.send_header(name, str(value))
                self.end_headers()
                self.wfile.write(payload)
                self.wfile.flush()
                self.close_connection = True

            def log_message(self, *args):
                pass

        server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        worker = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
        worker.start()
        try:
            yield f"http://127.0.0.1:{server.server_port}", requests
        finally:
            server.shutdown()
            server.server_close()
            worker.join(timeout=5)

    return start


@pytest.mark.parametrize("source", ["data/my query.tar.gz", Path("data/my query.tar.gz"), r"C:\data\my-query.tar.gz"])
def test_local_paths_keep_their_original_semantics(source, tmp_path):
    assert resolver()(source, tmp_path / "unused") == Path(source)
    assert not (tmp_path / "unused").exists()


def test_http_query_download_has_recorded_checksum_and_loads_real_pack(serve, tmp_path, capsys):
    source_bytes = DEMO.read_bytes()
    digest = hashlib.sha256(source_bytes).hexdigest()
    with serve({"/download": (200, source_bytes, {"Content-Type": "application/octet-stream"})}) as (base, requests):
        source = base + "/download?jobid=CNRC3YTJrfM&method=COLAB"
        downloaded = resolver()(source, tmp_path / "cache")
    assert requests == ["/download?jobid=CNRC3YTJrfM&method=COLAB"]
    assert downloaded.read_bytes() == source_bytes
    assert digest in downloaded.name
    printed = capsys.readouterr().out
    assert str(downloaded) in printed
    assert digest in printed
    assert set(downloaded.parent.iterdir()) == {downloaded}

    from daliscope.core.project import Project
    loaded = Project.load_pack(str(downloaded))
    assert loaded.query_length == 570
    assert len(next(iter(loaded.views.values()))) > 1000


def test_same_url_basename_keeps_distinct_job_content(serve, tmp_path):
    first = archive_bytes(extra=[("marker.txt", b"first")])
    second = archive_bytes(extra=[("marker.txt", b"second")])
    with serve({"/job.tar.gz?id=1": (200, first, {}), "/job.tar.gz?id=2": (200, second, {})}) as (base, _):
        path1 = resolver()(base + "/job.tar.gz?id=1", tmp_path)
        path2 = resolver()(base + "/job.tar.gz?id=2", tmp_path)
    assert path1 != path2
    assert path1.read_bytes() == first
    assert path2.read_bytes() == second
    assert set(tmp_path.iterdir()) == {path1, path2}


@pytest.mark.parametrize("problem,expected", [
    ("html", "gzip-tar"), ("expired", "HTTP 404"), ("truncated", "complete"), ("corrupted", "complete"),
    ("missing", "Query.pdb"), ("disconnected", "complete"),
])
def test_failed_download_preserves_existing_cache_and_removes_partial(problem, expected, serve, tmp_path):
    previous = archive_bytes()
    old_path = tmp_path / (hashlib.sha256(previous).hexdigest() + ".tar.gz")
    old_path.write_bytes(previous)
    payload = archive_bytes(extra=[("marker.txt", b"new")])
    status, headers = 200, {}
    if problem == "html":
        payload = b"<html><body>This job has expired</body></html>"
        headers = {"Content-Type": "text/html"}
    elif problem == "expired":
        status, payload = 404, b"Expired job"
    elif problem == "truncated":
        payload = payload[:-8]  # Keep tar headers but remove the gzip trailer.
    elif problem == "corrupted":
        # Valid first gzip/tar member followed by an invalid DEFLATE member.
        payload += bytes.fromhex("1f8b0800000000000003ff0000000000000000")
    elif problem == "missing":
        payload = archive_bytes(omit=("Query.pdb",))
    elif problem == "disconnected":
        headers = {"Content-Length": len(payload) + 100}
    with serve({"/job.tar.gz": (status, payload, headers)}) as (base, _):
        with pytest.raises(ValueError, match=expected):
            resolver()(base + "/job.tar.gz", tmp_path)
    assert old_path.read_bytes() == previous
    assert set(tmp_path.iterdir()) == {old_path}


def test_archive_is_never_extracted_to_the_filesystem(serve, tmp_path, monkeypatch):
    payload = archive_bytes(extra=[("../outside.txt", b"must stay in the archive")])

    def forbid_extraction(*args, **kwargs):
        pytest.fail("Pack source resolution must not extract archive members")

    monkeypatch.setattr(tarfile.TarFile, "extract", forbid_extraction)
    monkeypatch.setattr(tarfile.TarFile, "extractall", forbid_extraction)
    with serve({"/job.tar.gz": (200, payload, {})}) as (base, _):
        path = resolver()(base + "/job.tar.gz", tmp_path / "cache")
    assert path.is_file()
    assert not (tmp_path / "outside.txt").exists()
    assert set(tmp_path.iterdir()) == {tmp_path / "cache"}


def test_https_request_preserves_query_and_uses_requested_timeout(tmp_path, monkeypatch):
    resolve = resolver()
    module = importlib.import_module("daliscope.core.pack_source")
    payload = archive_bytes()
    seen = []

    class Response(io.BytesIO):
        headers = {"Content-Length": str(len(payload))}

    def open_response(url, *, timeout):
        seen.append((url, timeout))
        return Response(payload)

    monkeypatch.setattr(module, "urlopen", open_response)
    source = "https://example.invalid/download?job=a%2Fb&method=COLAB"
    path = resolve(source, tmp_path, timeout=17)
    assert seen == [(source, 17)]
    assert path.read_bytes() == payload


def test_run_notebook_input_cell_resolves_a_supplied_pack_url(serve, tmp_path, monkeypatch):
    import ast
    import json

    notebook = json.loads((ROOT / "notebooks/RunDaliScope.ipynb").read_text(encoding="utf-8"))
    cell = next(cell for cell in notebook["cells"] if cell.get("id") == "5bde54e0-8511-4272-aafe-1fe96b07dfea")
    payload = DEMO.read_bytes()
    monkeypatch.chdir(tmp_path)
    with serve({"/job.tar.gz": (200, payload, {})}) as (base, requests):
        source = base + "/job.tar.gz?jobid=CNRC3YTJrfM"
        tree = ast.parse("".join(cell["source"]))
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "pack_source" for target in node.targets):
                node.value = ast.Constant(source)
        runtime = {"IN_COLAB": False, "pack_source": source}
        exec(compile(ast.fix_missing_locations(tree), "RunDaliScope input cell", "exec"), runtime)
        assert requests == ["/job.tar.gz?jobid=CNRC3YTJrfM"], "Run input cell did not download the supplied pack_source"
    path = Path(runtime["pack_path"])
    assert path.read_bytes() == payload
