# Release notes

## DaliScope 0.1.8

Integrated Liisa Holm group's revised WorkedExample1 and WorkedExample2 notebooks. The examples add domain and motif diagrams, GOLD community illustrations, and a pairwise structure cartoon. Their GitHub sources omit saved execution output so the examples run cleanly from the first cell. Both tutorials keep the release-pinned Colab setup, their example data loading, and the corrected conversion from one-based PDB residues to zero-based clipping ranges. Claims about occupancy-weighted rigid-body fitting were corrected to match the implemented uniform least-squares fit.

Applying the stated PDB residue ranges through the existing conversion changes the case-study populations relative to the authors' saved notebook outputs: ZN `dom_5_FOLD` is 1,173 rather than 1,176 targets, and GOLD `dom_FOLD` is 6,465 rather than 6,708. In the corrected GOLD run, the two inspected annotation-selected modules each remained one child on recursion; earlier claims of additional splits were not retained. Figures and manuscript numbers should be taken from the corrected run.

The new pairwise cartoon helper retrieves a target structure from the DALI viewer endpoint and displays it with the query. It depends on that endpoint when the cartoon cell runs. Example HTML pages remain author-supplied earlier snapshots with unknown execution versions; they were not regenerated from these revised notebooks.

The bundled datasets, DALI/Pfam annotations and core analysis algorithms are unchanged from 0.1.7. Current notebook links, source setup and GOLD release downloads use 0.1.8. See the [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.8/docs/validation.md) for tested runtimes and limitations.

## DaliScope 0.1.7

The plane-bisector selector embeds its existing ipympl canvas in an ipywidgets container, bypassing the initial direct-display PNG preview involved in the reported `AttributeError: 'NoneType' object has no attribute 'dpi'`. The live canvas, controls and event callbacks are retained, with the existing display fallback for non-widget backends. Interactive output requires a working widget frontend; the direct static PNG preview is no longer emitted.

`RunDaliScope.ipynb` explains the actual Z-score/query-coverage axes, selected side and visible crop bounds. An empty bisector selection now gives instructions for adjusting the existing controls and rerunning the harvest cell, without registering an empty view. The Section 1 heading typo is corrected.

Scientific selection calculations, default analysis parameters and frozen input packs are unchanged from 0.1.6. Current notebook links, source setup and release downloads use 0.1.7. Earlier tags, assets and author-supplied HTML scientific outputs remain available. See the [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.7/docs/validation.md) for exact tested source and frontend limits.

## DaliScope 0.1.6

`RunDaliScope.ipynb` accepts a direct HTTP/HTTPS DALI data-pack download URL or a local file path through its native Colab `pack_source` field. The default remains the bundled 3ubpC example. Users can open the standard versioned Colab link, paste their generated pack's download URL, connect to a CPU runtime, and run the notebook.

The pack-source helper downloads URL inputs, checks the complete gzip-tar stream and required pack members, computes and reports SHA-256, and atomically saves the archive in a runtime cache. Local paths keep their existing meaning. The DALI results HTML page is not a data-pack input. Generated server packs are temporary; keep the downloaded input when recording or restoring an analysis.

The Colab guide supplies the canonical v0.1.6 `YOUR_COLAB_URL` and an HTML launch-link snippet for the DALI CGI template. Applying that snippet is a server-side deployment step; this release does not modify the hosted CGI. The fixed GitHub notebook link opens the tutorial, where the user supplies the pack URL.

Two author-supplied worked-example HTML exports are archived for browsing saved figures, tables and text without installation or Python execution. Their original filenames map to the current WorkedExample1 and WorkedExample2 notebooks. Download URLs, timestamps and original checksums are recorded; the browser title and added navigation header are the only presentation changes. Their execution version is unknown, and the saved scientific results were not regenerated for this release.

Current README, Colab/Pages links, notebook setup tags, and release downloads use 0.1.6. Archived HTML navigation links use the same version; their original scientific content remains unchanged. Scientific algorithms, default analysis parameters, and frozen input packs are unchanged from 0.1.5. Earlier tags and release assets remain available. See the [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.6/docs/validation.md) for the checks performed and their limits.

## DaliScope 0.1.5

Bundled the complete Pfam descriptions reference for local pack generation and clarified that DALI data packs are self-contained. Tutorial computations and example datasets are unchanged from 0.1.4. The standardized notebook names remain `RunDaliScope.ipynb`, `WorkedExample1.ipynb`, and `WorkedExample2.ipynb`.

## DaliScope 0.1.3

Refreshed Matplotlib backend discovery after live installation and preserved existing widget figures while later display/profile plots are drawn. The notebooks use standard Python 3 kernel metadata for Colab. Scientific algorithms and example data are unchanged from 0.1.2. Hosted Colab execution and frontend checks are recorded in the [validation record](validation.md).
