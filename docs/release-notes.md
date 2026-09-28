# DaliScope 0.1.6

`RunDaliScope.ipynb` accepts a direct HTTP/HTTPS DALI data-pack download URL or a local file path through its native Colab `pack_source` field. The default remains the bundled 3ubpC example. Users can open the standard versioned Colab link, paste their generated pack's download URL, connect to a CPU runtime, and run the notebook.

The pack-source helper downloads URL inputs, checks the complete gzip-tar stream and required pack members, computes and reports SHA-256, and atomically saves the archive in a runtime cache. Local paths keep their existing meaning. The DALI results HTML page is not a data-pack input. Generated server packs are temporary; keep the downloaded input when recording or restoring an analysis.

The Colab guide supplies the canonical v0.1.6 `YOUR_COLAB_URL` and an HTML launch-link snippet for the DALI CGI template. Applying that snippet is a server-side deployment step; this release does not modify the hosted CGI. The fixed GitHub notebook link opens the tutorial, where the user supplies the pack URL.

Two author-supplied worked-example HTML exports are archived for browsing saved figures, tables and text without installation or Python execution. Their original filenames map to the current WorkedExample1 and WorkedExample2 notebooks. Download URLs, timestamps and original checksums are recorded; the browser title and added navigation header are the only presentation changes. Their execution version is unknown, and the saved scientific results were not regenerated for this release.

Current README, Colab/Pages links, notebook setup tags, and release downloads use 0.1.6. Archived HTML navigation links use the same version; their original scientific content remains unchanged. Scientific algorithms, default analysis parameters, and frozen input packs are unchanged from 0.1.5. Earlier tags and release assets remain available. See the [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.6/docs/validation.md) for the checks performed and their limits.
