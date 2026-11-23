# Preparing a DaliScope data pack

A DaliScope pack is a `.tar.gz` archive containing a query, the aligned targets and their annotations. If you have a Dali or ISS web JOBID, enter it in the [data-pack service](http://ekhidna2.biocenter.helsinki.fi/dali/colab.html) and download the archive. If you already have the archive, open [RunDaliScope](https://colab.research.google.com/github/Liisa-Holm-group/DaliScope/blob/latest/notebooks/RunDaliScope.ipynb) and set `pack_source` to the archive URL or its local path in the runtime.

## Try local packing

Start from a Python 3.12 environment and get the current source:

```bash
git clone https://github.com/Liisa-Holm-group/DaliScope.git
cd DaliScope
```

The [four-hit tutorial](../examples/local-pack/README.md) includes a real alignment, query, miniature target structure store and the matching Pfam annotations. It runs without access to the group servers. Rebuilding its structure store is optional.

## Pack your own search

1. **Install the tools.** Clone the repository and run `python -m pip install -e ".[packer]"` from its root. Use a checkout containing `daliscope-build-store` to build a structure store. Both commands accept `--help`.
2. **Collect your search inputs.** Provide the query PDB (optionally gzip-compressed), its chain ID and the matching DALI `.dat` file. Export a tab-separated result with `target_id`, `namespace`, `z_score`, `qstarts`, `sstarts` and `lengths`. The last three columns give lists of one-based query starts, target starts and segment lengths for each alignment. A legacy DALI text report with structural equivalences also works. A ranking containing only target IDs and scores lacks these segments.
3. **Prepare the target structures.** If you have a compatible SQLite ledger and NPZ shards, use them directly. Otherwise, [build a store for your hits](#build-a-structure-store-for-your-hits) from the target files used in your search.
4. **Choose annotations.** Download and extract one ZIP from the [DaliScope Pfam dataset](https://doi.org/10.5281/zenodo.24465684). Pass both files from the same ZIP. The paper and release `v0.1.9` use **38.2**. For new analyses with **39.0**, use the packer from commit `32a97cd` or a checkout containing that change. These downloads supply Pfam hits and names. The target coordinates come from step 3.
5. **Create and load the archive.** Run the example command below with your paths. Check the reported packed-target count against the number of distinct hits. The packer reports and skips missing structures. Open the resulting `output/my_query.tar.gz` in DaliScope, or upload it to Colab.

```bash
daliscope-pack results.tsv AFDB2 my_query query.pdb A query.dat --ledger-db structures/master_ledger.db --shard-dir structures/shards --pfam-db pfam/pfam_data.db --pfam-names pfam/pfam_names.tsv --output-dir output
```

The second argument identifies the **target** database. For PDB searches, use `PDB` and optionally add `--ligand-db /path/to/PDB_ligands.db`. The Python module `python -m daliscope.packer.make_colab_pack` accepts the same arguments. Existing archives are preserved unless you supply `--overwrite`.

## Build a structure store for your hits

Use the DaliLite chain files from the search. If you only retained the original PDB structures, regenerate the chain files with the same preprocessing and residue order using the DaliLite import utility. See the [DaliLite v5 manual](https://github.com/liisaholm/DaliLite.v5/blob/main/MANUAL.html).

Prepare a tab-separated manifest. Use the exact target IDs in your alignment. File paths are relative to the manifest, or may be absolute. This example shows an AFDB2 model with its original DALI metadata and chain.

```tsv
target_id	dat	model	chain
eszmA	targets/eszmA.dat	targets/eszmA.pdb.gz	A
```

```bash
daliscope-build-store results.tsv AFDB2 targets.tsv --output-dir structures
```

The output is `structures/master_ledger.db` and `structures/shards/`. The tool selects only the distinct hits in the alignment, so a manifest may cover a larger collection. It reports missing targets or inconsistent coordinates and creates the final directory only after the build succeeds. Choose a new output directory for each build.

- **Experimental PDB structures.** Use `PDB`. The target DaliLite `.dat` files already contain ordered Cα coordinates, sequences and secondary structure. The manifest needs only `target_id` and `dat`. Experimental structures use the packer's conventional confidence placeholder of 100.
- **Predicted structures.** For `AFDB2`, include the matching model in the `model` column. The builder reads the first model, defaults to chain A, and takes Cα confidence from its B-factor field. PDB and mmCIF files, including gzipped files, are supported. The coordinate order and residue set must match the DALI file used by the alignment. In particular, use the corresponding clipped model if the search uses clipped target numbering.
- **Your own identifiers.** Use `LOCAL` for both commands. Supply a target DALI file for each hit and include `model` for predicted targets with pLDDT. Omit `--pfam-db` when you have no annotations. To include your own annotations, use a compatible database with `database=LOCAL` and `target_name` values matching the search IDs. The Zenodo snapshots use the Dali reference IDs and namespaces. A UniProt accession or lab identifier must first be mapped to its matching reference record before using those annotations.

For example, to pack hits from your own target collection without Pfam, run these commands from the directory containing your manifest, query and results.

```bash
daliscope-build-store results.tsv LOCAL targets.tsv --output-dir my-structures
daliscope-pack results.tsv LOCAL my_query query.pdb A query.dat --ledger-db my-structures/master_ledger.db --shard-dir my-structures/shards --output-dir output
```

For custom IDs containing leading zeros or values such as `NA`, use the checkout containing this guide to load the pack in local Jupyter. The published v0.1.9 notebooks predate preservation of these IDs. Core coordinate-based views work with the packed structures. The optional full-atom cartoon fetches coordinates using Dali server identifiers.

One invocation handles one target namespace. Use separate packs for searches against different identifier systems.

## Required external schemas

The protein SQLite ledger has a `proteins` table with `dali_id`, `n_res`, `umid`, `shard_path`, `cx`, `cy`, `cz`, and `ns_id`. Each referenced NPZ shard provides `<umid>_coords`, `<umid>_plddt`, `<umid>_seq`, `<umid>_dssp`, and `<umid>_compnd`.

The Pfam annotation database has either `hmmer_hits_with_CL` (Pfam 39.0 databases) with `target_name`, `query_accession`, `e_value`, `env_from`, `env_to`, `database`, and `CL`, or `hmmer_hits` (Pfam 38.2 databases) with the first six of these columns plus `pfam_entries` with `AC` and `CL`. When both are present, `hmmer_hits_with_CL` is used. Namespace labels in the `database` column are `PDB`, `AFDB2`, `BFVD`, and `Viro3D`; custom databases can also use `LOCAL`. Labels are matched without regard to case, so namespace 4 (`VIRO3D`) selects the `Viro3D` hits.

The optional ligand database has `ligand` (`pdb_id`, `chain_id`, `res_name`, `hetatm`) and `ligand_names` (`res_name`, `compound`). The complete structure and ligand resources are maintained separately. The Pfam snapshots are linked above.
