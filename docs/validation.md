# Validation record

## 0.1.6 URL data-pack input

`RunDaliScope.ipynb` now accepts a direct HTTP/HTTPS data-pack download URL or a local path in its native Colab `pack_source` field. URL downloads are streamed to a temporary file, checked for a complete gzip-tar archive and the seven required regular input files, then published under a SHA-256 cache filename. The helper does not extract archive members to the filesystem. Failed downloads stop without loading the bundled example or replacing a valid cached input. The printed checksum records downloaded bytes rather than verification against a separate upstream checksum.

The final Windows suite passed **63 tests** in **38.70 seconds**, with the same 13 existing Infomap/ipympl/traitlets warnings. Fourteen new checks cover local path compatibility, HTTP query strings, HTTPS URL handling, complete and malformed gzip streams, HTML/404 responses, missing members, interrupted downloads, preservation of cached files, absence of filesystem extraction, and the notebook input cell's use of a supplied URL. Dependency checking passed.

On **28 September 2026**, the DALI server generated a fresh pack for JOBID **CNRC3YTJrfM**. Its archive was downloaded through the new URL input and the complete comprehensive notebook ran in a fresh **Windows reference kernel**: **29 non-empty code cells** passed in **39.52 seconds**, with query **3ubpC**, query length **570**, and **1,048 targets**. The only validation override was setting `pack_source` to the server's generated download URL. The input SHA-256 was `075cd1a35da92514517b810954b860e92119efd5430caa6303d9a0b9289efc9f`. The generated pack's page stated a one-hour lifetime, independently of the original search results' one-week lifetime.

All four notebook schemas passed. After normalizing version references, quickstart and WorkedExample1/2 match 0.1.5 exactly; RunDaliScope changes only its input instructions and input cell. Scientific analysis cells, default parameters, cell IDs and metadata are preserved, and source notebooks contain no executed outputs. The archived example HTML scientific content also remains byte-for-byte reproducible from the supplied exports after reversing the added browser title/header; only current notebook navigation links changed for 0.1.6.

This is a local computational test of the URL path, not a new hosted Colab frontend test. The fixed GitHub/Colab launch URL opens the notebook; users paste the generated archive URL before running cells. Arbitrary URL query parameters are not consumed. The DALI CGI maintainer must apply the documented button URL to the hosted template. Earlier hosted widget evidence and remaining frontend limitations below remain historical.

## 0.1.5 bundled Pfam descriptions

The local packer now defaults to the installed `daliscope/packer/data/pfam_names.tsv` instead of requiring a descriptions table outside the repository. The explicit `--pfam-names` option remains available. The README links to the official self-contained DALI data-pack service and clarifies that analysis reads Pfam data from the pack itself.

The fixed official Pfam 38.2 table contains **30,134 unique families**. Independent comparison with its upstream archive verified the only changes: adding a five-column header and filling empty `clan` fields with the Pfam accession. Other fields, including empty `clan_short`, whitespace, and row order, are preserved. All **760 rows / 3,800 fields** of the originally supplied mini descriptions match exactly. This does not establish the contents of the unavailable private full table. Source URLs, hashes, conversion, and CC0-1.0 terms are recorded in the [bundled reference notes](../daliscope/packer/data/README.md).

The complete Windows unit suite passed **49 tests** in **36.41 seconds**, with the same 13 existing Infomap/ipympl/traitlets warnings. The two new checks exercise default pack generation from an unrelated directory without an external names table, and standalone mini extraction preserving unassigned clan labels. Existing explicit names-table generation and invalid-query protection remain covered.

The wheel and source distribution include the full reference, its source/license notes, and the compatibility Bash wrapper; the reference bytes match the recorded SHA-256 in both archives. A wheel installed into a separate directory with `--no-deps --no-index --target` passed default mini extraction and generated a loadable synthetic pack from an unrelated working directory, without `--pfam-names`. Module paths were checked to ensure the test used the installed wheel rather than the source checkout. The wrapper was checked for inclusion; it was not executed on Windows.

All four notebook schemas passed validation. Parsed notebooks are identical to 0.1.4 after normalizing only `v0.1.4` to `v0.1.5`; scientific code, parameters, metadata, cell IDs, and example data are unchanged, and source notebooks contain no executed output. The computational and hosted/frontend evidence below remains historical; this patch does not claim another full GOLD run or new hosted mouse checks.

## 0.1.4 notebook naming update

The 0.1.4 tutorials are `RunDaliScope.ipynb`, `WorkedExample1.ipynb`, and `WorkedExample2.ipynb`; `01_quickstart.ipynb` retains its name. README and notebook references, Colab links/setup tags, test paths, and release downloads use the same 0.1.4 version.

All four notebook schemas were validated and their parsed contents compared with the 0.1.3 originals. The only notebook changes are the requested filenames/display metadata and filename/version references; computational code, analysis parameters, cell IDs, and frozen inputs are preserved. Source notebooks contain no executed outputs. The Windows Colab-setup and figure-lifecycle tests passed **14 tests**, with 11 existing upstream warnings, using the renamed comprehensive notebook. GitHub Actions executes the quickstart, `RunDaliScope.ipynb`, and `WorkedExample1.ipynb` using their new paths.

The full Windows/Linux and hosted results below were recorded on 0.1.3, when the tutorials had their original filenames. They are historical computational/frontend evidence. This naming patch does not claim new hosted-Colab mouse checks or a repeated complete GOLD computation. Earlier release tags and their notebook names remain available.

## 0.1.3 validation history

Status on **28 September 2026** for DaliScope **0.1.3**, source commit [e8cc5cd](https://github.com/Liisa-Holm-group/DaliScope/commit/e8cc5cdb996fc5d3a03e57e45dc236e180461cc5). Completed checks and checks awaiting final results are listed separately.

## Completed checks

| Check | Environment | Result |
| --- | --- | --- |
| Full unit suite | Windows | **47 passed** in 32.59 seconds; 13 warnings: 2 existing Infomap warnings and 11 upstream ipympl/traitlets warnings |
| Full unit suite in an isolated copy | Linux, Python 3.12.3 | **47 passed** in 25.52 seconds; 13 warnings |
| Targeted Colab setup and figure-lifecycle suite | Windows and Linux, Matplotlib 3.10 | **14 passed on each platform** |
| All four complete tutorials | Fresh Windows kernels | Quickstart: 8 code cells; comprehensive workflow: 30; SignatureTest: 21; GOLD: 15; **74 total**, no notebook errors |
| All four complete tutorials | Fresh Linux kernels in the isolated copy | Quickstart: 8 code cells; comprehensive workflow: 30; SignatureTest: 21; GOLD: 15; **74 total**, no notebook errors |
| GitHub Actions | Windows, Ubuntu, Colab dependency environment | Installation/dependency/test checks passed; included-data Linux tutorials completed |
| Complete comprehensive workflow | Cold hosted Colab runtime, source commit e8cc5cd | **29 non-empty code cells passed**, plus one empty cell; automatic setup, no extra backend patch |
| Hit-plane view redraw after all later analysis | Same hosted runtime | Switching to `dom_1` visibly updated the points and y-axis to approximately 0.25; original canvas comm remained open |
| Static plot controls after all later analysis | Same hosted runtime | Selecting `pfam` then Refresh in the Pfam/Clan violin panel regenerated all three distributions with `by pfam` titles |

The Linux full-suite environment used NumPy 2.0.2, pandas 2.2.3, Matplotlib 3.11.2, IPython 7.34.0, ipykernel 6.17.1, ipympl 0.10.0, and Plotly 7.1.0.

An initial isolated Linux comprehensive run stalled in the client waiting for I/O after the profile cell. A diagnostic confirmed that the kernel had completed that cell and was idle; the interrupted attempt is not counted as passing. A fresh run without diagnostic intervention completed all 29 non-empty code cells in 26.53 seconds. No production source or analysis input was changed for the retry.

The comprehensive workflow keeps its community pipeline off by default. Notebook completion checks the configured default path, rather than every optional parameter combination.

The Linux GOLD tutorial completed in 276.98 seconds. Across each platform, the notebooks contain 74 code cells: 71 non-empty cells execute, and three empty cells are skipped.

## Hosted test provenance

The cold hosted test uses the actual public 0.1.3 source commit in a fresh Google Colab runtime. Because the release tag was not yet published, a visible first cell clones the feature branch and creates a local `v0.1.3` test tag for setup verification. The remaining original tutorial cells are unchanged. This procedure is distinct from opening the final public release link and downloading its assets; it must not be reported as that final onboarding check.

The hosted runtime diagnostic reported Python 3.13.15, google-colab 1.0.0, IPython 7.34.0, ipykernel 6.17.1, NumPy 2.1.3, pandas 2.2.3, SciPy 1.16.3, Matplotlib 3.10.0, Plotly 7.1.0, anywidget 0.9.21, and ipympl 0.10.0. DaliScope 0.1.3 imported from `/content/DaliScope/daliscope/__init__.py`, at source commit e8cc5cd. After the complete workflow and browser view switch, its hit-plane canvas comm was still open and its y-axis spanned approximately 0.00175–0.25175. No scientific parameters were changed during execution; the view switch was performed afterwards.

## Scope and practical limits

The 0.1.3 changes refresh cached Matplotlib backend entries after live installation, preserve figures owned by existing widgets during later drawing, and use generic Python 3 notebook metadata. Scientific algorithms, analysis parameters, domain conventions, and frozen example data are unchanged from 0.1.2. Local users still explicitly select **Python (DaliScope)**.

Backend refresh uses Matplotlib's private `backend_registry._validate_and_store_entry_points` method to register only new entries and retain existing backends such as `inline`. Future Matplotlib changes may require compatibility adjustments; the targeted checks above exercise Matplotlib 3.10.

Direct mouse rotation/zoom in py3Dmol remains unverified because browser automation refused fractional input coordinates inside the viewer iframe. That tool limitation is separate from DaliScope computation and widget emission. Hosted SignatureTest and the full GOLD analysis, Plotly mouse controls and direct Matplotlib canvas dragging, and browser file downloads have not yet been recorded as completed. See the [Colab guide](colab.md) for the manual checks.

The local reference environment uses Python 3.12; [Windows reference requirements](requirements-windows-py312.txt) remain available. Cloud runtime libraries can change, so setup preserves the runtime's installed numerical dependencies. Production protein/Pfam/ligand databases are not included; pack generation is exercised with a synthetic fixture. Record version, pack checksum, domains, filters, and clustering/motif parameters with scientific outputs.

## Release verification

The release workflow builds from the final clean merged commit, records that commit and each attachment's SHA-256/size in `checksums.json`, and checks the wheel, source distribution, and clean source ZIP for matching versions and unintended local files. The source ZIP is compared with the recorded Git tree. The separate GOLD attachment retains SHA-256 `d3fc732890400446bd128fa178e9a3e13c3c28b1dac5c3656c76459185f660ed`. Public attachment downloads and the versioned tutorial entry points are checked during publication. The hosted code test above records its own exact pre-release source and procedure.

## Earlier validation

The immutable [0.1.2 release-time record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.2/docs/validation.md) documents the previous Windows/Linux checks. After that release, its eight-cell quickstart completed on hosted Colab with output creation, four-view session restoration, and visible py3Dmol rendering. Its comprehensive tutorial completed after applying the backend refresh in a temporary visible cell; domain Apply/Reset and centering were also checked. Those historical results do not replace the cold 0.1.3 checks above.
