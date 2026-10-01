# Using DaliScope in Google Colab

## Open and run a tutorial

Open one of the versioned links in the [README](../README.md#run-in-your-browser-with-colab), connect to a runtime, and select **Runtime → Run all**. A standard **CPU** runtime is sufficient for the tutorial code; GPU acceleration is not used by these workflows.

The first code cell handles Colab setup. It obtains **DaliScope v0.1.9** in `/content/DaliScope`, installs the package with `.[colab]` in the runtime's Python environment, and enables Google's custom widget manager. After live installation it refreshes Matplotlib's cached backend discovery so the running kernel can request the newly installed widget backend. Installing Colab dependencies is done once per runtime. A GitHub notebook link supplies only the notebook, so the source clone and data files are still necessary; the setup cell provides them automatically.

Use the Colab runtime's kernel. Local venv creation, Conda commands, ipykernel registration, and the local **Python (DaliScope)** kernel selection are not part of this procedure. In a local Jupyter session, the first Colab setup cell is a no-op and the README's local installation applies.

The v0.1.9 notebooks use generic **Python 3** kernel metadata to avoid Colab's unrecognized local-runtime notice. This does not select your local environment: when using local Jupyter, still choose the registered **Python (DaliScope)** kernel described in the README.

Start with `01_quickstart.ipynb`. It should report query `3ubpC`, length **570**, and **1,048** rows in the default view, then show plots, a structure viewer, and a session save/restore confirmation. Longer examples run their default analysis choices first. After completion, explore the display controls. For revised data, domains, filters, or motifs, start again from the cell that loads a fresh `Project`, then execute the following sections in order and apply your new control choices before registering views. Reusing an existing view name does not replace its earlier analysis; use a fresh project or a new unique view name.

## Data and working directory

The cloned repository includes the 3ubpC, dotp, and ZN example packs. Quickstart sets the working directory to `/content/DaliScope`; the longer notebooks use `/content/DaliScope/notebooks`, where paths such as `data/ZN_full.tar.gz` resolve. Keep the setup and path cells before the loading cells.

`WorkedExample2.ipynb` uses the optional **GOLD** pack, approximately 50 MiB. In Colab it downloads the release asset when absent and verifies its SHA-256 against `notebooks/data/datasets.json` before analysis. To download or verify it manually in a Colab code cell:

```python
!python /content/DaliScope/scripts/download_example.py GOLD --tag v0.1.9
```

## Use your own DALI data pack

1. Generate a self-contained pack using the [DALI data-pack service](http://ekhidna2.biocenter.helsinki.fi/dali/colab.html). On the generated page, copy the actual archive download URL. Use the pack link, not the DALI search-results HTML page or the page containing the download button.
2. Open [RunDaliScope v0.1.9 in Colab](https://colab.research.google.com/github/Liisa-Holm-group/DaliScope/blob/v0.1.9/notebooks/RunDaliScope.ipynb). Before executing the analysis, paste the download URL into its **`pack_source`** string field. HTTP and HTTPS download URLs are supported. The default field value loads the bundled example instead.
3. Connect to a standard CPU runtime and select **Runtime → Run all**. The loading cell downloads the archive, checks its complete gzip-tar stream and required pack files, and prints the cached path and SHA-256 before loading the project. Record that checksum with your analysis.

For an uploaded pack, use Colab's Files panel and enter its absolute path in `pack_source`, such as `/content/my-query.tar.gz`. In local Jupyter, edit the same `pack_source` string in the code cell to a local path or URL. The archive must follow the [data-pack contract](data-format.md); a lone PDB is not a DaliScope input pack.

The generated server download is temporary: the current download page states that the pack is available for **one hour**, separately from the search results' one-week lifetime. Download and keep the input pack before the link expires. If necessary, regenerate it from the original DALI job. After the data-loading cell completes, save the runtime copy to your computer with:

```python
from google.colab import files
files.download(str(pack_path))
```

The runtime cache is also temporary. A computed SHA-256 identifies the downloaded bytes; it is not a comparison against a separately supplied server checksum. Keep that input file alongside outputs and sessions. If you change `pack_source` after analysis has started, reload a fresh project and execute the subsequent cells in order.

The comprehensive `RunDaliScope.ipynb` leaves its expensive community pipeline off by default. A **Pipeline is OFF** message for that step is expected. The signature step is controlled by its IC condition. GOLD executes community detection, recursive examples, and STRUCTAL cohesion evaluation; the cohesion stage can take several minutes on the frozen example. Let these computation cells finish before running dependent cells. Runtime depends on the population and cloud hardware.

## Use the browser controls

The setup enables Colab's custom widget manager for the domain viewer, ipywidgets, Plotly controls, and Matplotlib widget backend. Version 0.1.3 keeps existing widget figures available when later display/profile plots are drawn and cleaned up. After later sections finish, return to the comprehensive tutorial's hit-plane plot, switch its selected view, and check that the plotted points update along with the controls. Some viewers also load JavaScript from external servers. An emitted widget without working mouse controls is not a completed frontend check.

- Rotate and zoom the py3Dmol structure; use available target controls to change the selected structure.
- In the domain viewer, edit the PDB domain text and click **Apply**. Verify that the coloring and domain-centering buttons respond. **Reset** restores the initial ranges. If testing an edit before running a case study, reset it before registering the domains.
- Hover over Plotly points, zoom, and switch the selected view. Changes should update the displayed data.
- In panels with **↺ Refresh**, choose the view or grouping, then click that button to regenerate the plot and refresh available views. For example, the Pfam/Clan violin panel changes its chart titles and distributions after selecting `pfam` and clicking Refresh.
- Apply selection/profile choices before the cells that lock their values and register a view. If that view was already registered, reload a fresh project and work through the sections again, or choose a new unique view name. Display changes do not automatically rebuild subsequent analysis.

Domain text uses **one-based inclusive PDB residues**; the notebooks explicitly convert it to zero-based half-open clipping ranges. GOLD's stated core is PDB residues **23–123**. See [usage](usage.md) for the API conventions.

If a canvas is blank, confirm the setup cell completed, reconnect if needed, and rerun its display cell. If installation requests a runtime restart, restart once and rerun the setup and path cells before analysis. A newly allocated runtime has to repeat setup. Report persistent failures with the full traceback or visible browser error, notebook/cell, DaliScope version, and runtime Python version.

## Keep and download results

The Files panel shows files on the cloud runtime rather than on your computer. Quickstart writes `hit-landscape.png`, `hits.tsv`, and its saved session under `/content/DaliScope/results/quickstart/`. Download the figure or table with a code cell:

```python
from google.colab import files
files.download("/content/DaliScope/results/quickstart/hit-landscape.png")
files.download("/content/DaliScope/results/quickstart/hits.tsv")
```

To download the complete results directory:

```python
from pathlib import Path
import shutil
from google.colab import files

results = Path("/content/DaliScope/results")
if results.is_dir():
    archive = shutil.make_archive("/content/daliscope-results", "zip", results)
    files.download(archive)
else:
    print("No results directory yet; run the export cells first.")
```

Download outputs before deleting or losing the runtime. Saving the notebook to Drive preserves the notebook document; it does not preserve every file under `/content`. A DaliScope session records the original pack path. A session referencing `/content/...` is not automatically portable to a Windows path; keep the pack and its recorded path available when restoring structures. Figures and TSV exports can be used independently.

## Source version and validation limits

The v0.1.9 Colab links, setup source checkout, and GOLD release selection use the same version. The notebooks include Liisa Holm group's revised worked examples and a structural-cartoon helper. Version 0.1.9 clarifies worked-example interpretation and plot labels. The notebooks retain the corrected PDB-residue clipping conversion and their version-pinned setup. The existing URL/local-path `pack_source` input remains available. Keep the software version fixed when recording an analysis. Record the data checksum, domains, filters, and clustering/motif parameters alongside outputs; selecting a different pack or changing those choices can change results. The [validation record](validation.md) separates current checks from earlier hosted frontend evidence.

The 0.1.3 comprehensive workflow completed all 29 non-empty tutorial code cells in a fresh hosted Colab runtime on 28 September 2026. Setup installed the dependencies and refreshed backend discovery automatically. After the final analysis, switching the existing hit-plane widget to `dom_1` visibly updated the axes and plotted points. This pre-release check used the public source commit, with only the initial clone adapted to its branch because the release tag was not yet available; see [validation](validation.md) for exact provenance.

The earlier published v0.1.2 quickstart completed its eight original code cells on hosted Colab, including PNG/TSV/session creation, restoration, and visible py3Dmol rendering. Domain Apply/Reset and centering were checked there. Direct mouse rotation/zoom in py3Dmol, Plotly mouse controls, and browser file downloads remain unverified. Hosted WorkedExample1 and the complete hosted GOLD analysis have not been recorded. The immutable v0.1.2 documentation and release assets retain their earlier snapshot.

## Server launcher configuration

For the DALI server's Launch Colab button, set `YOUR_COLAB_URL` to this permanent notebook URL:

```text
YOUR_COLAB_URL = "https://colab.research.google.com/github/Liisa-Holm-group/DaliScope/blob/latest/notebooks/RunDaliScope.ipynb"
```

The CGI can render the following HTML. Replace `ACTUAL_PACK_DOWNLOAD_URL` with the generated archive URL in the server template:

```html
<a href="https://colab.research.google.com/github/Liisa-Holm-group/DaliScope/blob/latest/notebooks/RunDaliScope.ipynb" target="_blank" rel="noopener">Launch DaliScope in Colab</a>
<p>Copy the <a href="ACTUAL_PACK_DOWNLOAD_URL">DALI data-pack download URL</a> into the notebook's <code>pack_source</code> field, connect to a CPU runtime, and choose Runtime &rarr; Run all.</p>
```

The GitHub/Colab link opens the notebook; users then paste the pack download URL into its native field. Server maintainers must apply this replacement to the hosted CGI template. Changes to this repository do not deploy that server template.

The `latest` branch is advanced by the DaliScope release maintainer only after a formal release is published and its public downloads and entry points have been verified. Each notebook keeps its explicit release tag internally, so its installed source and release data match the opened notebook. The DALI server maintainer applies the permanent button URL once; subsequent DaliScope releases do not require another CGI link change. See the [publishing procedure](publishing.md#update-the-permanent-colab-entry).

For reproducible analyses, use a versioned notebook link, such as [RunDaliScope v0.1.9](https://colab.research.google.com/github/Liisa-Holm-group/DaliScope/blob/v0.1.9/notebooks/RunDaliScope.ipynb), and record the package version and input checksum. A saved notebook copy and an existing runtime do not update automatically when `latest` advances; reopen the permanent link in a fresh runtime to use a new release.
