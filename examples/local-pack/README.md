# Build a DaliScope pack from four real hits

This example starts with four alkaline phosphatase hits from an AFDB2 search with PDB chain **1ajaA**. The download contains the six-column alignment TSV, query PDB and DALI metadata, a four-target structure store, target input files for rebuilding it, and the seven matching Pfam 39.0 hits. All example inputs work offline after installing the packer.

### Create and open the pack

Use Python 3.12. Get the current checkout with `git clone https://github.com/Liisa-Holm-group/DaliScope.git` and enter it with `cd DaliScope`. Run these commands from the repository root in the same environment. Each command also runs as written in PowerShell.

```bash
python -m pip install -e ".[packer]"
python -m zipfile -e examples/local-pack/inputs.zip local-pack-inputs
cd local-pack-inputs
daliscope-pack results.tsv AFDB2 phosphatase query.pdb.gz A query.dat --ledger-db store/master_ledger.db --shard-dir store/shards --pfam-db pfam-39.0.db --output-dir output
python -c "from daliscope.core.project import Project; p = Project.load_pack('output/phosphatase.tar.gz'); print('Query length:', p.query_length)"
```

Expect **4 alignments, 1,965 coordinate rows, 63 alignment segments, 4/4 successful superimpositions and a query length of 435**. All four targets have Pfam assignments. To explore the result in Colab, upload `output/phosphatase.tar.gz` and use its runtime path as `pack_source`.

### Rebuild the structure store

Still in `local-pack-inputs`, run:

```bash
daliscope-build-store results.tsv AFDB2 targets.tsv --output-dir rebuilt
daliscope-pack results.tsv AFDB2 rebuilt-phosphatase query.pdb.gz A query.dat --ledger-db rebuilt/master_ledger.db --shard-dir rebuilt/shards --pfam-db pfam-39.0.db --output-dir output
```

The builder selects the targets in `results.tsv` from the manifest `targets.tsv`. The example's complete reference store has only those four targets. For a different search, use your own manifest and target files and the appropriate full [Pfam snapshot](https://doi.org/10.5281/zenodo.24465684). See [Preparing your own data](../../docs/prepare-data.md).

### Data sources

The query coordinates are from [PDB 1AJA](https://www.rcsb.org/structure/1AJA). PDB archive data are available under [CC0](https://www.wwpdb.org/about/usage-policies). DALI alignment segments and DaliLite metadata are supplied for this example. The target files and their derived coordinate arrays use the AlphaFold Database v2 models used for the search. AlphaFold DB data are attributed to Google DeepMind and EMBL-EBI under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See the [AlphaFold Database publication](https://doi.org/10.1093/nar/gkad1011).

| DALI target | UniProt accession | Organism |
| --- | --- | --- |
| eszmA | P00634 | *Escherichia coli* K-12 |
| lpa4A | P19405 | *Bacillus subtilis* 168 |
| eivhA | Q7K3X8 | *Drosophila melanogaster* |
| eh9sA | Q9VRM9 | *Drosophila melanogaster* |

The example Pfam SQLite file contains only the AFDB2 hit rows for these four IDs. The compiled annotations are from [Liu and Holm's DaliScope Pfam snapshots](https://doi.org/10.5281/zenodo.24465684), licensed under CC BY 4.0. The bundled Pfam family descriptions retain CC0.
