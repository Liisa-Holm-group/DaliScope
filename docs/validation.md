# Validation record

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
