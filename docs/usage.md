# Usage guide

## Environment and input

Use the [Colab entry](colab.md) for browser setup, or follow the root README to install locally and select the DaliScope Jupyter kernel. Install before launching the kernel. If you update packages in an existing session, restart that kernel once and run from the beginning.

Start with `notebooks/01_quickstart.ipynb`. It locates the repository whether the kernel starts in the root or `notebooks`. Longer tutorials set their working directory to `notebooks`, where paths such as `data/ZN_full.tar.gz` apply.

The input is a `.tar.gz` data pack containing structures, sequences, alignment blocks, and metadata. A lone PDB or raw result table is insufficient; see [data format](data-format.md) and [data preparation](prepare-data.md).

```python
from daliscope.analysis.data_preview import load_and_inspect_project
project, FULL_DF, FULL_VIEW = load_and_inspect_project(pack_path)
```

`project` owns the pack, derived views, models, and plot recipes. `FULL_DF` is the default view's DataFrame, and `FULL_VIEW` is its name. Inspect registered views with `project.list_inventory()`.

## Read the overview

- `z_score` is supplied by the DALI search. Clipping does not recompute it.
- `query_coverage` and `target_coverage` divide aligned length by the full query and target sequence lengths, respectively, including in clipped views. A small clipped domain can have low query coverage despite being almost completely aligned; these fractions alone do not establish a single-domain or multidomain architecture.
- `sequence_identity` is an aligned-residue fraction between 0 and 1.
- `rmsd` is the current view's superimposition metric. Domain clipping recomputes a uniform least-squares fit over the retained aligned residues. Row-filtered subsets, including `_FOLD` views, inherit the parent fit and metrics.
- `pfam` and `clan` describe the aligned target region. `Unassigned` indicates an annotation gap, rather than an established new function.

Automatic `FILTERED_0` and `FILTERED_1` populations are exploratory selections. Their names do not establish a biological classification; inspect the structures and annotations.

## Select a population

```python
project.add_subset(
    "z8", FULL_DF["z_score"] >= 8, FULL_VIEW,
    function="z_score_filter", parameters={"z_score_min": 8},
)
```

The score cutoff is an example, not a universal confidence threshold.

The clipping API takes **zero-based, half-open query ranges**: the start is included and the end is excluded. Displayed PDB residues are generally one-based. For displayed residues 23 through 123, pass `22-123` to clipping.

```python
from daliscope.mechanics.metrics import register_multiple_domains
from daliscope.mechanics.domain_ranges import pdb_domain_string_to_clipping
from daliscope.analysis.occupancy import create_fold_subset
pdb_domains = "23-123"  # One-based, inclusive displayed residues.
register_multiple_domains(pdb_domain_string_to_clipping(pdb_domains), project, custom_name="domain")
create_fold_subset(project, "domain", cutoff=0.8)
```

Separate domains with commas, and discontinuous parts of one domain with underscores. A direct clipping string can be `0-50_100-150, 200-300`. The equivalent PDB string is `1-50_101-150, 201-300`.

`DomNetViewer.domain_string` and its text box use **one-based, inclusive PDB/PUU ranges**. Pass `viewer.clipping_domain_string` to `register_multiple_domains` to preserve the highlighted residues. `pdb_domain_string_to_clipping` and its inverse, `clipping_domain_string_to_pdb`, also support multiple and discontinuous domains. Do not pass the displayed string directly to clipping. Motif and plotting helpers have their own position arguments; inspect the aligned residue and the helper's docstring before defining a signature.

The fold selection registers `<view>_FOLD` using an occupancy-weighted coverage score. A cutoff of 0.8 is more restrictive than 0.5. This measures recurring structural coverage in the chosen population rather than an unweighted percentage of domain residues. The selection does not apply occupancy weights to the coordinate fit or recompute the retained rows' metrics.

## Plot and explore

```python
project.viz.architecture(FULL_VIEW, renderer="matplotlib", max_targets=10)
project.viz.plot_msa(FULL_VIEW, 0, project.query_length,
                     plot_type="heatmap", data_col="dssp_pileup")
```

The structural fingerprint/community API uses **one-based, inclusive** ranges (for example `(1, 150)`), while clipping uses zero-based, half-open ranges. Convert the stored clipping provenance to select exactly the same query residues, including discontinuous domains:

```python
from daliscope.mechanics.domain_ranges import clipping_ranges_to_pdb
community_ranges = clipping_ranges_to_pdb(
    project.provenance["domain"].parameters["domain_ranges"]
)
# Use community_ranges as domain_range when calling run_community_detection.
```

For example, clipping `(22, 123)` converts to community/PDB `(23, 123)`. A literal `(23, 123)` has different meanings in these two APIs.

The longer tutorials show domain widgets, sequence signatures, and STRUCTAL/Infomap communities. The comprehensive `RunDaliScope-1.ipynb` controls expensive optional steps with execution flags; the GOLD worked example executes its full community and recursive analyses. Try a small population first; record domain ranges, neighbor counts, Infomap arguments, and motif parameters.

Drag in the 3D viewer to rotate, scroll to zoom, and use its target controls to change structures. Plotly supports hover labels and zoom. These browser interactions need a manual check beyond headless execution.

## Export and resume

```python
from pathlib import Path
output = Path("results/my-analysis")
output.mkdir(parents=True, exist_ok=True)
FULL_DF[["target_id", "z_score", "query_coverage", "rmsd"]].to_csv(
    output / "hits.tsv", sep="\t", index=False,
)
project.save(str(output / "session"))
restored = type(project).open(str(output / "session"))
```

`results/` is ignored by Git. Replacing an existing session requires `overwrite=True`. Sessions record the original pack path; keep that pack accessible at the recorded path to restore structures. Sessions contain Python pickle data and should only be opened from trusted sources.

Use `fig.savefig(...)` for Matplotlib PNG/PDF export and `fig.write_html(...)` for interactive Plotly HTML. Optional Plotly static image export requires Kaleido and a suitable browser installation.

## Troubleshooting

| Symptom | Action |
| --- | --- |
| Missing Python module | Install in the selected kernel environment; check `sys.executable` |
| Pack not found | Check the working directory and input path |
| View not found | Run preceding load/domain cells; inspect `project.list_inventory()` |
| Duplicate view name | Choose another name or create a fresh project |
| Blank widget/canvas | Check the kernel and notebook extras, then reload the page and rerun its display cell |
| Interactive plots fail offline | Some viewers use a JavaScript CDN; use static Matplotlib plots or provide network access |
| Missing GOLD input | Follow `notebooks/data/README.md` |
| Restored session lacks structures | Restore its original pack at the recorded path |

For a bug report, include the traceback, OS, Python/DaliScope versions, notebook and cell name, and pack checksum.

