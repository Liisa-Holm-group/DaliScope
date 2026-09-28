# Data-pack format

A pack is a gzip-compressed tar archive. Files are read in memory and identified by basename. Put one copy of each file at the archive root.

## Required contents

| File | Required fields or arrays | Purpose |
| --- | --- | --- |
| `query_meta.tsv` | `id`, `length`, `sequence`, `dssp`, `compnd`, `domain_string` | One query metadata row |
| `summary.tsv` | `alignment_id`, `protein_id`, `namespace_id`, `z_score` | Alignments |
| `id_mapping.tsv` | `protein_id`, `target_id`, `namespace_raw`, `namespace_id` | Target mapping |
| `segments.tsv` | `alignment_id`, `q_start`, `s_start`, `length` | Aligned blocks |
| `colab_pack_index.tsv` | `umid`, `dali_id`, `length`, `start_idx`, `end_idx`, `compnd`, `centroid` | Target offsets and descriptions |
| `colab_pack.npz` | `coords`, `seq`, `dssp`; optional `plddt` | Target coordinate and sequence arrays |
| `Query.pdb` | Query C-alpha coordinates in metadata sequence order | Query visualization and superimposition |

`coords` has shape `(N, 3)`. `seq` and `dssp` have length `N` and store character codes. Target offsets are zero-based, half-open: `[start_idx, end_idx)`, with `end_idx - start_idx == length`. Coordinate units are Angstroms.

The query metadata length must match the sequence length and number of query C-alpha atoms. The query PDB must use the same residue order as the DALI metadata; the preparation command renumbers the selected chain.

Alignment blocks in `segments.tsv` use **one-based** DALI start positions and positive lengths. The loader subtracts one from the starts. This differs from the half-open ranges accepted by the clipping API.

## Optional annotations

| File | Fields | Purpose |
| --- | --- | --- |
| `pfam.tsv` | `target_name`, `query_accession`, `e_value`, `env_from`, `env_to`, `database`, `clan` | Pfam hits; envelope positions are one-based |
| `mini_pfam_names.tsv` | `pfam`, `clan`, `clan_short`, `short`, `name` | Family and clan descriptions |
| `mini_ligand.tsv` | `pdb_id`, `chain_id`, `res_name`, `hetatm` | Ligand records for PDB targets |
| `mini_ligand_names.tsv` | `res_name`, `compound` | Ligand names |

Without Pfam annotations, retained hits are labeled `Unassigned`. The standard preparation command requires a Pfam database; a custom producer can omit the annotation files.

## Namespaces and provenance

Namespaces are `1=PDB`, `2=AFDB2`, `3=BFVD`, and `4=VIRO3D`. Alignments whose targets were not packed can be absent from derived views; compare summary and packed-target counts.

Record source/database versions, query and chain, search settings, generation date, and SHA-256 checksum. The frozen examples' identifiers and checksums are in `notebooks/data/datasets.json`; upstream database release dates were not recorded in the supplied packs.
