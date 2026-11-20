# Preparing a DaliScope data pack

## Starting from an existing pack

A [DALI data pack](http://ekhidna2.biocenter.helsinki.fi/dali/colab.html) is self-contained: it includes the structures, sequences, Pfam annotations, and family/clan descriptions needed for DaliScope analysis. Use the linked DALI service to generate a pack from a completed DALI search job, then download it and keep a permanent copy. Use its `.tar.gz` path in the notebook. No separate protein ledger, Pfam database, or description-table download is needed to analyze the pack.

## Local pack generation

The local pack generator consumes DALI alignment output and query files plus pre-existing protein and Pfam annotation databases. It does not download structures or create those databases. Ask the group maintainers for their current database paths/schema when working on the group server. Pfam descriptions are included in the installed DaliScope package.

Install `python -m pip install -e ".[notebook,packer]"`, then inspect `daliscope-pack --help`.

```bash
daliscope-pack results.tsv PDB my_query query.pdb A query.dat \
  --ledger-db /path/to/master_ledger.db \
  --shard-dir /path/to/shards \
  --pfam-db /path/to/pfam_data.db \
  --output-dir /path/to/output
```

On PowerShell, put the command on one line or use PowerShell's line-continuation syntax. Alternatively run `python -m daliscope.packer.make_colab_pack` with the same arguments. The historical `make_colab_pack.csh` is a Bash wrapper for this Python command.

Arguments:

- `results.tsv`: a supported DALI TSV/text report with hit scores and alignment equivalences. Modern TSV accepts `target_id`, `namespace`, `z_score`, `qstarts`, `sstarts`, and `lengths`; the starts/lengths contain lists such as `[1, 20]`.
- `PDB`: target namespace, also accepting `AFDB2`, `BFVD`, `VIRO3D`, or numeric IDs 1-4.
- `my_query`: output basename; the result is `my_query.tar.gz`.
- `query.pdb A`: query structure and selected chain. A gzipped PDB is supported.
- `query.dat`: DALI metadata containing the matching query sequence, secondary structure, length, and domain hierarchy.
- `--ledger-db` / `--shard-dir`: target structures and sequence records.
- `--pfam-db`: target Pfam annotations in a local SQLite database.
- `--pfam-names`: optional replacement for the bundled complete Pfam 39.0 description table. For another annotation snapshot, provide a descriptions TSV with at least the five columns listed below. For a Pfam 38.2 annotation database, use the 38.2 table of release 0.1.9 (see the reference README).

The default descriptions file is `daliscope/packer/data/pfam_names.tsv`, resolved relative to the installed package, so it works from any current directory. DaliScope reads its `pfam`, `clan`, `clan_short`, `short`, and `name` columns; the bundled file also has `type`, `clan_name`, and `clan_entries`, as in the DALI server table. The generator selects the families used by the current search and includes their descriptions as `mini_pfam_names.tsv` inside the new pack. See the [reference source and license](../daliscope/packer/data/README.md) for the fixed upstream version, conversion rule, and checksums. Descriptions and retired accessions can differ between Pfam releases; record the annotation and description versions used to generate a new pack.

For PDB targets, add `--ligand-db /path/to/PDB_ligands.db` to include ligand records. Ligand extraction is optional and currently supports PDB only.

The generator runs intermediate steps in a temporary directory, packages the required files at the archive root, and preserves existing output unless `--overwrite` is given. It uses the same Python interpreter as the calling command. It creates a local pack by default; upload requires the explicit `--upload-to user@host:/directory/` option and a configured `scp`/SSH connection.

## Required external schemas

The protein SQLite ledger has a `proteins` table with `dali_id`, `n_res`, `umid`, `shard_path`, `cx`, `cy`, `cz`, and `ns_id`. Each referenced NPZ shard provides `<umid>_coords`, `<umid>_plddt`, `<umid>_seq`, `<umid>_dssp`, and `<umid>_compnd`.

The Pfam annotation database has either `hmmer_hits_with_CL` (Pfam 39.0 databases) with `target_name`, `query_accession`, `e_value`, `env_from`, `env_to`, `database`, and `CL`, or `hmmer_hits` (Pfam 38.2 databases) with the first six of these columns plus `pfam_entries` with `AC` and `CL`. When both are present, `hmmer_hits_with_CL` is used. Namespace labels in the `database` column are `PDB`, `AFDB2`, `BFVD`, and `Viro3D`; they are matched without regard to case, so namespace 4 (`VIRO3D`) selects the `Viro3D` hits.

The optional ligand database has `ligand` (`pdb_id`, `chain_id`, `res_name`, `hetatm`) and `ligand_names` (`res_name`, `compound`). These per-protein annotation, structure, and ligand databases are server resources and are not bundled with this repository; the Pfam descriptions reference is bundled.

Missing proteins/shards are reported and skipped by the coordinate producer. Inspect the reported counts and archive contents before interpreting incomplete populations. A synthetic database fixture is exercised in the release tests; production database provisioning remains a group-server task.
