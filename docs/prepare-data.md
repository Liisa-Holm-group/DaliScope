# Preparing a DaliScope data pack

## Starting from an existing pack

If your DALI workflow already provides a DaliScope/Colab pack, download it and keep a permanent local copy. Use its `.tar.gz` path in the notebook. No protein ledger or annotation database is needed to analyze the included example packs.

## Local pack generation

The pack generator consumes DALI alignment output and query files plus pre-existing local protein and Pfam databases. It does not download structures or create those databases. Ask the group maintainers for their current database paths/schema when working on the group server.

Install `python -m pip install -e ".[notebook,packer]"`, then inspect `daliscope-pack --help`.

```bash
daliscope-pack results.tsv PDB my_query query.pdb A query.dat \
  --ledger-db /path/to/master_ledger.db \
  --shard-dir /path/to/shards \
  --pfam-db /path/to/pfam_data.db \
  --pfam-names /path/to/pfam_names.tsv \
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
- `--pfam-db` / `--pfam-names`: target Pfam annotations and descriptions.

For PDB targets, add `--ligand-db /path/to/PDB_ligands.db` to include ligand records. Ligand extraction is optional and currently supports PDB only.

The generator runs intermediate steps in a temporary directory, packages the required files at the archive root, and preserves existing output unless `--overwrite` is given. It uses the same Python interpreter as the calling command. It creates a local pack by default; upload requires the explicit `--upload-to user@host:/directory/` option and a configured `scp`/SSH connection.

## Required external schemas

The protein SQLite ledger has a `proteins` table with `dali_id`, `n_res`, `umid`, `shard_path`, `cx`, `cy`, `cz`, and `ns_id`. Each referenced NPZ shard provides `<umid>_coords`, `<umid>_plddt`, `<umid>_seq`, `<umid>_dssp`, and `<umid>_compnd`.

The Pfam database has `hmmer_hits` with `target_name`, `query_accession`, `e_value`, `env_from`, `env_to`, and `database`, plus `pfam_entries` with `AC` and `CL`. Namespace labels in `hmmer_hits.database` are `PDB`, `AFDB2`, `BFVD`, and `VIRO3D`. The Pfam descriptions TSV has `pfam`, `clan`, `clan_short`, `short`, and `name`.

The optional ligand database has `ligand` (`pdb_id`, `chain_id`, `res_name`, `hetatm`) and `ligand_names` (`res_name`, `compound`). The group databases and their release metadata are not bundled with this repository.

Missing proteins/shards are reported and skipped by the coordinate producer. Inspect the reported counts and archive contents before interpreting incomplete populations. A synthetic database fixture is exercised in the release tests; production database provisioning remains a group-server task.

## Related tools

The group's [ISS_tools](https://github.com/Liisa-Holm-group/ISS_tools) provides DALI output conversion and Pfam annotation workflows. Confirm that its output contains the equivalence fields expected by this generator. See [data format](data-format.md) for the final archive contract.
