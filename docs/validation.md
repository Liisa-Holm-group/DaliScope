# Validation record

Local validation on Windows with Python 3.12.14, for DaliScope 0.1.0.

## Completed checks

- Fresh virtual environment installation of `.[notebook,packer,dev]` completed successfully.
- `python -m pip check`: no broken requirements.
- 15 regression tests passed. They cover the bundled full view, filtered-view/session round trips (including nested fields and row indices), missing inputs, absent/empty Pfam annotations, numeric database namespaces, weighted graph input, a bounded structural community pipeline, and an end-to-end synthetic pack generation fixture.
- `01_quickstart.ipynb`: all 7 code cells completed in a fresh kernel.
- `RunDaliScope-1.ipynb`: all 29 code cells completed with its default optional-pipeline settings.
- `SignatureTest.ipynb`: all 20 code cells completed in a fresh kernel.
- All published analysis, mechanics, optics, and widgets modules imported.
- Four frozen example packs passed coordinate/sequence array, index extent, query-length checks, and checksum verification.
- The README overview image was generated from the bundled quickstart and visually inspected.

- `WorkedExample-6.ipynb`: all 14 code cells completed, including full GOLD community detection, recursive examples, and module-cohesion evaluation.
- All four notebooks completed 70 code cells in fresh kernels. Interactive browser actions remain a separate check.
- Wheel and source distribution built successfully; wheel installation and example loading passed from outside the source directory.
- Source distribution includes the three small example packs and excludes GOLD, caches, and local archives.
- CITATION.cff passed the official schema; local documentation links and staged credential-pattern checks passed.
- The clean release snapshot contains 92 files and is approximately 14.9 MiB; history/local outputs are kept separately.

## Reference environment

Major installed versions: NumPy 2.5.3, pandas 2.3.3, SciPy 1.18.1, Matplotlib 3.11.2, Plotly 7.1.0, Numba 0.67.0, Infomap 2.15.1, and PyRoaring 1.1.0. The complete Windows/Python 3.12 environment is recorded in `requirements-windows-py312.txt`.

For the Windows reference environment, install those pinned requirements and then install the repository with `python -m pip install -e ".[notebook,packer,dev]"`. This pin file records a tested environment; it is not a cross-platform lock file.

## Fixes found during release preparation

- The default full query view clipped away the first position; it now spans `[0, query_length)`.
- Invalid input paths returned two values from a helper otherwise returning three.
- Packs without Pfam annotations were accepted by the archive reader but failed during metadata joining.
- Session restore dropped the pack reference; view saves dropped row indices, complex-field types, and column order.
- Weighted-edge partitioning attempted to import an absent external module.
- Pfam extraction passed string namespaces to an integer-indexed namespace mapping.
- Pack generation depended on fixed laboratory paths and always attempted a server upload.

## Practical limits

Headless execution checks Python computations and widget emission. It does not verify browser mouse interaction or CDN availability. Check the 3D/domain widgets in JupyterLab before a demonstration. The complete production protein/Pfam/ligand databases are not included, so production database generation must be checked on the group server. The synthetic fixture verifies the portable command's integration contract.

Clipping uses zero-based half-open ranges; fingerprint/community masks use one-based inclusive ranges. Convert boundaries explicitly between these APIs. Legacy case-study interpretations and user thresholds need scientific review for the chosen inputs; executing a notebook does not validate its biological conclusions.

GitHub Actions is configured for Ubuntu and Windows on Python 3.12. Current hosted results are available on [GitHub Actions](https://github.com/Liisa-Holm-group/DaliScope/actions). Consult the individual runs for the OS/Python combinations that passed; the local checks recorded here establish the Windows reference environment.
