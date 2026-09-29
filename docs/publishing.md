# Publishing DaliScope

The public repository is [Liisa-Holm-group/DaliScope](https://github.com/Liisa-Holm-group/DaliScope). Maintain its clone in `DaliScope`; retain the original project backup and laboratory databases separately.

## Prepare a patch release

1. Work on a branch and review the complete diff.
2. Update `daliscope.__version__`, `CITATION.cff`, the download script's default tag, data instructions, notebook setup tags, current Colab/Pages links, release notes, and validation record together. Preserve version references describing historical tests or archived scientific output. If only the added navigation header in an archived HTML file changes, update its saved checksum in `docs/example-outputs/outputs.json`; keep the original download record unchanged.
3. Install `.[notebook,packer,dev]`, run `python -m pip check` and `python -m pytest -q`, and execute the affected notebooks in fresh kernels with `scripts/validate_notebooks.py`. GOLD is required for the advanced example.
4. Build with `python -m build` from the final clean checkout. Verify that the source distribution includes the reference requirements and three small packs, and that neither package includes GOLD or local files. If a custom source ZIP is needed, export tracked files with `scripts/export_release.py` to an empty directory outside the repository; new helpers must already be tracked.
5. Merge only after the Windows and Ubuntu GitHub Actions checks pass. Create the release from the tested commit, with its corresponding version tag.

The software license is MIT, copyright Liisa Holm group. Contact: `hao.liu@helsinki.fi`. A manuscript describing DaliScope is being prepared for submission; update this status after submission and add the paper citation when available.

## Attach files

For future releases, attach the wheel, source distribution, and a checksum manifest. Attach the unchanged `GOLD.tar.gz` as well: it is about 50 MiB and excluded from ordinary Git commits. Its dataset checksum is in `notebooks/data/datasets.json`.

GitHub provides source archives for each release, so an additional custom source ZIP is optional. Keep existing release attachments and their checksum manifests available to preserve download links. Colab tutorials install the package and obtain their example data automatically; users do not need to download these attachments manually.

Use the Release attachment basenames in `checksums.json`, with `bytes` and `sha256` for each file. Keep every attachment in one directory for verification; do not encode the maintainer's local folder layout in the manifest. The manifest itself is not included in its own hash list.

For example, with a validated tag and authenticated GitHub CLI:

```bash
gh release create v0.1.7 --repo Liisa-Holm-group/DaliScope --verify-tag --draft --title "DaliScope 0.1.6" --notes-file docs/release-notes.md
gh release upload v0.1.7 /path/to/assets/GOLD.tar.gz /path/to/assets/daliscope-0.1.6-py3-none-any.whl /path/to/assets/daliscope-0.1.6.tar.gz /path/to/assets/checksums.json --repo Liisa-Holm-group/DaliScope
```

Verify the draft attachments against `checksums.json`, then publish with `gh release edit v0.1.7 --draft=false --repo Liisa-Holm-group/DaliScope` and check their public downloads.

Once the release exists, `python scripts/download_example.py GOLD` downloads and verifies the optional pack. For earlier versions, pass the corresponding tag explicitly, such as `--tag v0.1.0`.

## Preserve reproducibility

Keep earlier version tags and released artifacts available. Document any coordinate or parameter corrections that change historical results, and avoid silently replacing earlier scientific outputs. Record the DaliScope version, data checksum, domains, filters, and clustering/motif parameters with each analysis.

Inspect current annotation summaries and structures when interpreting community modules. Module numbers can change between runs. Browser widget interactions and production laboratory databases require their own checks; the validation record describes what was executed locally and on GitHub.