# DaliScope

DaliScope is a notebook-based toolkit for exploring DALI structural search results. It combines domain-level filtering, Pfam annotation, sequence and secondary-structure profiles, interactive 3D views, motif analysis, and structural community detection.

DaliScope analyzes an existing DALI data pack. Running a DALI search and provisioning the protein/Pfam databases are separate steps.

![Example DALI hit landscape](docs/images/overview.png)

## Quick start

Use Python 3.12 for the reference environment. The package requires Python 3.10 or newer; platform and version coverage are recorded in [validation](docs/validation.md).

```bash
git clone https://github.com/Liisa-Holm-group/DaliScope.git
cd DaliScope
python -m venv .venv
```

Activate the environment on Linux/macOS:

```bash
source .venv/bin/activate
```

Or on Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

If PowerShell does not allow activation, use `.\.venv\Scripts\python.exe` in place of `python` in the following commands.

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[notebook]"
python -m ipykernel install --user --name daliscope --display-name "Python (DaliScope)"
python -m jupyterlab
```

Open **`notebooks/01_quickstart.ipynb`**, select the **Python (DaliScope)** kernel, and run the cells from top to bottom. The small `3ubpC_PDB25.tar.gz` example is included in the repository. The notebook loads it, plots the hit landscape and domain architecture, creates a filtered view, shows a 3D superimposition, and saves a session.

## Documentation

- [Usage guide](docs/usage.md): choosing data, views, parameters, plots, and saved sessions.
- [Data-pack format](docs/data-format.md): required files, table columns, and coordinate conventions.
- [Preparing your own data](docs/prepare-data.md): a portable pack generator and its database requirements.
- [Examples and data provenance](notebooks/data/README.md): included datasets and the optional GOLD download.
- [Validation](docs/validation.md): checks performed and practical limitations.
- [Publishing guide](docs/publishing.md): creating the group repository and distributing larger data.

## Tutorials

| Notebook | Purpose | Input |
| --- | --- | --- |
| `01_quickstart.ipynb` | Short introduction with static and interactive output | Included 3ubpC pack |
| `RunDaliScope-1.ipynb` | Comprehensive workflow with user-controlled domain and motif steps | Included 3ubpC pack |
| `SignatureTest.ipynb` | Multidomain analysis and sequence signatures | Included ZN pack |
| `WorkedExample-6.ipynb` | GOLD-like fold communities and STRUCTAL/Infomap analysis | Optional GOLD pack |

The comprehensive workflow contains interactive and optional steps. Execute it section by section, adjusting the parameters for your query. Case-study parameters and biological interpretations apply to their frozen input datasets.

## Minimal Python example

Run this from the repository root after installation:

```python
from daliscope.core.project import Project

project = Project.load_pack("notebooks/data/3ubpC_PDB25.tar.gz")
full_view = next(iter(project.views))
hits = project.views[full_view]
print(project.query_identifier, project.query_length, len(hits))
print(hits[["target_id", "z_score", "query_coverage", "rmsd"]].head())
```

The default view is named `dom_0-<query_length>`. Discover the view name or use the name returned by the notebook loading helper.

## Development

```bash
python -m pip install -e ".[notebook,packer,dev]"
python -m pytest -q
python -m build
```

Tests exercise the bundled example, session persistence, optional annotations, database namespace handling, and graph partition input. GitHub Actions also runs the quickstart notebook in a fresh kernel.

## Contact and citation

Maintained by the **Liisa Holm group**. Contact: **hao.liu@helsinki.fi**.

A manuscript describing DaliScope is being prepared for submission. Publication details will be added when available. When reporting an analysis, record the DaliScope version, data-pack checksum, selected domains, filters, and clustering/motif parameters, and cite DALI and the annotation sources used in your work.

## Licensing

DaliScope code and its original documentation are distributed under the [MIT License](LICENSE), copyright 2026 Liisa Holm group. Third-party structures and annotations retain their source providers' terms; see the [data provenance notes](notebooks/data/README.md). A software license does not replace the licenses of source data or accompanying publications.
