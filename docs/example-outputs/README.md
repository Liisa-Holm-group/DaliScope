# HTML example outputs

These author-supplied notebook exports let users browse saved figures, tables and text **without installing software or running code**.

| Current example | Browse HTML | Original export |
| --- | --- | --- |
| WorkedExample1 | [HTML example output](https://liisa-holm-group.github.io/DaliScope/example-outputs/WorkedExample1.html) | `SignatureTest.html` |
| WorkedExample2 | [HTML example output](https://liisa-holm-group.github.io/DaliScope/example-outputs/WorkedExample2.html) | `WorkedExample-6.html` |

The [example-output landing page](https://liisa-holm-group.github.io/DaliScope/) links to both exports and their current runnable notebooks.

## What the HTML provides

Figures, tables and text are saved in the HTML; all 35 images in WorkedExample1 and all 11 images in WorkedExample2 are embedded. Python-backed notebook controls cannot execute from these saved pages. Embedded 3D views require their external 3Dmol viewer and browser WebGL support; exported RequireJS and MathJax resources also load from a CDN.

To change analysis inputs or run the interactive workflow, use the corresponding notebook in [Colab](../colab.md) or local Jupyter.

## Source and preservation

Liisa Holm supplied these temporary URLs. Both files were downloaded on **28 September 2026**:

- WorkedExample1: [SignatureTest.html](http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp/SignatureTest.html).
- WorkedExample2: [WorkedExample-6.html](http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp/WorkedExample-6.html).

The temporary source URLs may expire. These archived copies remain in this repository and are served by GitHub Pages.

Their browser titles and a small header identify them as HTML example outputs and map the original names to the current notebook names. Scientific code, saved outputs, embedded figures and viewer data are otherwise unchanged. Reversing those two presentation changes reproduces the downloaded source bytes exactly.

The exports do not record their execution version, commit or data-pack checksum. They are preserved author-supplied snapshots, not results regenerated with current releases. Links to runnable notebooks use v0.1.7 independently of these snapshots.

Original download SHA-256:

- `SignatureTest.html`: `32938ce14f0edc7ff0c306366d13b8b7ef37b9f5bb50d2a037311814cd5dfdd3`
- `WorkedExample-6.html`: `333b222d85c22f4bad6db9f7a21f6911476919674f9cc9bf16f98341ab15ce81`

[outputs.json](outputs.json) records the original URLs, retrieval timestamps, original and saved file checksums, and the two presentation changes.

## Site publishing

GitHub Pages serves the `main` branch's `/docs` directory. `docs/.nojekyll` keeps the exported HTML as static files. Earlier release tags and downloadable assets retain their original contents.
