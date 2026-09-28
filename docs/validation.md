# Validation record

Local validation on 28 September 2026 with Windows and Python 3.12.14, for DaliScope 0.1.1.

## Completed checks

- Editable installation of `.[notebook,packer,dev]` completed successfully; `python -m pip check` found no broken requirements.
- 32 regression tests passed. The existing release tests cover pack loading, view/session round trips, missing annotations, weighted graph input, a bounded structural community pipeline, and synthetic pack generation. New tests cover colliding view/model names, legacy session manifests, invalid pack inputs preserving previous outputs, and conversion between PDB/UI and clipping coordinates.
- All four notebooks completed in fresh kernels: `01_quickstart.ipynb` (7 code cells), `RunDaliScope-1.ipynb` (29), `SignatureTest.ipynb` (20), and `WorkedExample-6.ipynb` (14), totaling 70 code cells.
- The GOLD example completed full community detection, both annotation-selected recursive examples, module-cohesion evaluation, and the transitive-superimposition example. Representatives are selected from current annotation counts rather than historical module identifiers.
- The custom-filter examples were exercised against the bundled 3ubp pack using the documented column names.
- An independently installed wheel passed dependency checks, imports of all 44 published modules, example loading from outside the source directory, domain conversion, and collision-free session restoration.
- The source notebooks contain no saved outputs or execution counts.

Release verification checks the wheel, source distribution, clean source ZIP, and frozen GOLD pack against the checksum manifest. The source distribution includes `requirements-windows-py312.txt` and the three small example packs. Package/source archives exclude GOLD, caches, local outputs, and laboratory archives; GOLD is a separate attachment. The source ZIP is compared with the release commit's tracked files. Checksum entries use the public attachment basenames.

GitHub Actions runs installation, dependency checks, all regression tests, quickstart execution, and package builds on Ubuntu and Windows with Python 3.12. Results are available on [GitHub Actions](https://github.com/Liisa-Holm-group/DaliScope/actions); consult the individual run for its exact tested commit.

## Reference environment

Major installed versions: NumPy 2.5.3, pandas 2.3.3, SciPy 1.18.1, Matplotlib 3.11.2, Plotly 7.1.0, Numba 0.67.0, Infomap 2.15.1, and PyRoaring 1.1.0. The complete Windows/Python 3.12 environment is recorded in [requirements-windows-py312.txt](requirements-windows-py312.txt).

For the Windows reference environment, install those pinned requirements and then install the repository with `python -m pip install -e ".[notebook,packer,dev]"`. This pin file records a tested environment; it is not a cross-platform lock file.

## Fixes checked in 0.1.1

- Distinct session names such as `domain 1` and `domain_1` retain their own view/model files. Existing 0.1.0 manifests remain readable.
- Pack generation rejects missing query CA coordinates and lengths that disagree with metadata, before replacing any existing output.
- The GOLD tutorial's stated PDB residues 23–123 map to clipping range `[22, 123)` and to the same physical residues in community masks. Multi-domain and discontinuous ranges are also tested.
- Custom-filter instructions use the actual column names. Coverage is described relative to the full query length, without inferring protein/domain classification from coverage alone.
- Tutorials describe domain-clipped superimposition followed by occupancy-based filtering. Filtering inherits the fit; it does not perform an occupancy-weighted refit.
- Release metadata, version references, reference requirements, and attachment checksum names are kept consistent.

## Practical limits

Headless execution checks Python computations and widget emission. It does not verify browser mouse interaction or CDN availability. Check the 3D/domain widgets in JupyterLab before a demonstration. The complete production protein/Pfam/ligand databases are not included, so production database generation must be checked on the group server. The synthetic fixture verifies the portable command's integration contract.

Clipping uses zero-based half-open ranges; UI/PDB and fingerprint/community masks use one-based inclusive ranges. Use the conversion helpers documented in [usage](usage.md). Correcting domain boundaries can change selected populations, motifs, and community assignments relative to 0.1.0, although the underlying fitting, motif, and clustering algorithms are unchanged. Record the software version, data checksum, ranges, filters, and parameters with analysis results.

Module identifiers are run-specific. Interpret current annotations and structures for the chosen inputs; executing a notebook does not validate its biological conclusions. Two pending-deprecation warnings from Infomap's current result-access APIs remain; the dependency is constrained below version 3.
