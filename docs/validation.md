# Validation record

Checks performed on 28 September 2026 for DaliScope 0.1.2. This record separates local/Linux compatibility checks from execution on Google's hosted Colab service.

## Completed checks

- 41 tests passed on Windows/Python 3.12.14 and Ubuntu 24.04/Python 3.12.3. Coverage includes frozen example loading, persistence, pack generation, domain conventions, stale/mixed Colab checkouts, repeated setup, renderer selection, and enforcement of notebook cell timeouts.
- All four tutorials completed in fresh Linux kernels: quickstart (8 code cells), comprehensive workflow (30), SignatureTest (21), and GOLD (15), totaling 74 code cells. The environment used IPython 7.34.0, ipykernel 6.17.1, pandas 2.2.3, and NumPy 2.0.2. No JupyterLab installation was needed for this lightweight environment.
- Quickstart, the comprehensive workflow, and SignatureTest also completed in fresh Windows kernels. The GOLD full community, recursive, cohesion, and transitive-superimposition analyses completed on Linux.
- An independent Linux environment installed the official `google-colab` package from [Google's colabtools source](https://github.com/googlecolab/colabtools/blob/58990611c039377dbb9a39eaa10f47e327f3f4e1/setup.py), including its nine declared dependencies. Running the real bootstrap twice in an IPython kernel performed one real pip installation and called the real custom widget manager twice, without stubs. Google-managed and already-loaded NumPy/pandas/SciPy/Matplotlib versions remained unchanged. Dependency checks reported no broken requirements.
- The optional GOLD asset downloaded from the public 0.1.1 release into a clean Linux validation directory and matched the recorded SHA-256. Its data is unchanged in 0.1.2.
- Source notebooks have no saved outputs or execution counts. Original frozen data-pack checksums are unchanged.

GitHub Actions checks Windows and Ubuntu/Python 3.12 installation, dependencies, the full test suite, quickstart, and builds. A separate Ubuntu job installs the Colab-compatible IPython/kernel/pandas versions plus NumPy 2.0.2 and executes the three included-data tutorials. Consult the [Actions run](https://github.com/Liisa-Holm-group/DaliScope/actions) for its exact tested commit.

Release verification checks the wheel, source distribution, clean source ZIP, and frozen GOLD pack against the checksum manifest. The source distribution contains the Colab setup script, tutorials, guides, Windows reference requirements, and three small packs. GOLD remains a separate Release attachment. The source ZIP is compared with the release commit's tracked files; archives exclude caches, local outputs, and laboratory archives.

## Reference environments

The Windows/Python 3.12 reference requirements are in [requirements-windows-py312.txt](requirements-windows-py312.txt). Install those requirements, then install the repository with `python -m pip install -e ".[notebook,packer,dev]"`. This is a recorded Windows environment, rather than a cross-platform lock file.

The Linux compatibility environment used NumPy 2.0.2, pandas 2.2.3, SciPy 1.18.1, Matplotlib 3.11.2, Plotly 7.1.0, Numba 0.67.0, Infomap 2.15.1, and PyRoaring 1.1.0. IPython 7.34.0 and ipykernel 6.17.1 match Google's declared requirements. Colab's actual runtime image can change; the setup obtains constraints from the runtime's installed `google-colab` metadata and keeps its existing numerical libraries.

## Practical limits

Hosted Colab execution and browser mouse interaction have not yet been verified; the available test browser requires a Google login. The real widget-manager calls above establish Python-side integration and output emission, not working browser controls. Follow the [Colab guide](colab.md) to check py3Dmol rotation/zoom, domain Apply/Reset/centering, Plotly controls, and Matplotlib widgets. CDN access and Colab resource availability depend on the current runtime and browser.

The complete production protein/Pfam/ligand databases are not included. Production pack generation needs a group-server check; a synthetic fixture exercises the portable command's integration contract.

Clipping uses zero-based half-open ranges; UI/PDB and community masks use one-based inclusive ranges. Use the conversion helpers in [usage](usage.md). The algorithms, parameters, and data in 0.1.2 are unchanged from 0.1.1. Boundary corrections introduced in 0.1.1 can change results relative to 0.1.0; earlier releases and their [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.1/docs/validation.md) remain available.

Record version, pack checksum, domains, filters, and clustering/motif parameters with analysis outputs. Module identifiers are run-specific; inspect current annotations and structures before drawing biological conclusions. Two existing Infomap pending-deprecation warnings remain, and Infomap is constrained below version 3.
