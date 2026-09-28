# Using DaliScope in Google Colab

## Open and run a tutorial

Open one of the versioned links in the [README](../README.md#run-in-your-browser-with-colab), connect to a runtime, and select **Runtime → Run all**. A standard **CPU** runtime is sufficient for the tutorial code; GPU acceleration is not used by these workflows.

The first code cell handles Colab setup. It obtains **DaliScope v0.1.2** in `/content/DaliScope`, installs the package with `.[colab]` in the runtime's Python environment, and enables Google's custom widget manager. Installing Colab dependencies is done once per runtime. A GitHub notebook link supplies only the notebook, so the source clone and data files are still necessary; the setup cell provides them automatically.

Use the Colab runtime's kernel. Local venv creation, Conda commands, ipykernel registration, and the local **Python (DaliScope)** kernel selection are not part of this procedure. In a local Jupyter session, the first Colab setup cell is a no-op and the README's local installation applies.

Start with `01_quickstart.ipynb`. It should report query `3ubpC`, length **570**, and **1,048** rows in the default view, then show plots, a structure viewer, and a session save/restore confirmation. Longer examples run their default analysis choices first. After completion, adjust controls or parameters and rerun the affected downstream cells to apply the changes.

## Data and working directory

The cloned repository includes the 3ubpC, dotp, and ZN example packs. Quickstart sets the working directory to `/content/DaliScope`; the longer notebooks use `/content/DaliScope/notebooks`, where paths such as `data/ZN_full.tar.gz` resolve. Keep the setup and path cells before the loading cells.

`WorkedExample-6.ipynb` uses the optional **GOLD** pack, approximately 50 MiB. In Colab it downloads the release asset when absent and verifies its SHA-256 against `notebooks/data/datasets.json` before analysis. To download or verify it manually in a Colab code cell:

```python
!python /content/DaliScope/scripts/download_example.py GOLD --tag v0.1.2
```

For your own pack, upload it through Colab's Files panel and set `pack_path` to its absolute runtime path, such as `/content/my-query.tar.gz`. It must follow the [data-pack contract](data-format.md); a lone PDB is not a DaliScope input pack.

The comprehensive `RunDaliScope-1.ipynb` leaves its expensive community pipeline off by default. A **Pipeline is OFF** message for that step is expected. The signature step is controlled by its IC condition. GOLD executes community detection, recursive examples, and cohesion evaluation; let these computation cells finish before running dependent cells. Runtime depends on the population and cloud hardware.

## Use the browser controls

The setup enables Colab's custom widget manager for the domain viewer, ipywidgets, Plotly controls, and Matplotlib widget backend. Some viewers also load JavaScript from external servers. An emitted widget without working mouse controls is not a completed frontend check.

- Rotate and zoom the py3Dmol structure; use available target controls to change the selected structure.
- In the domain viewer, edit the PDB domain text and click **Apply**. Verify that the coloring and domain-centering buttons respond. **Reset** restores the initial ranges. If testing an edit before running a case study, reset it before registering the domains.
- Hover over Plotly points, zoom, and switch the selected view. Changes should update the displayed data.
- Adjust selection/profile controls, then rerun the later cells that lock the selected values and register a view. Changing a displayed control does not automatically rebuild every subsequent analysis.

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

The v0.1.2 Colab links, setup source checkout, and GOLD release selection use the same version. Keep that version fixed when recording an analysis. Record the data checksum, domains, filters, and clustering/motif parameters alongside outputs; upgrading source or correcting domain ranges can change results.

Actual Colab frontend interaction has not yet been verified. Local notebook execution and GitHub Actions checks establish computations and widget emission for the recorded environments, rather than proving that every Colab runtime or browser widget works. Consult [validation](validation.md) for the checks that were executed. The procedures here describe how to perform that additional cloud/browser check.
