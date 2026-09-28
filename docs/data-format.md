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

`query_meta.tsv`'s `domain_string` records the PUU suggestion using **one-based, inclusive PDB/sequence residue numbers**, with commas separating domains and underscores separating discontinuous segments. For example, `1-50_101-150, 151-200` represents two domains. The domain viewer displays this convention. Convert it with `pdb_domain_string_to_clipping` before calling the clipping API; the example becomes `0-50_100-150, 150-200`.

The fingerprint/community API also uses one-based inclusive ranges. Convert clipping provenance with `clipping_ranges_to_pdb` instead of reusing its numeric boundaries unchanged. This preserves both endpoints and supports multiple segments. See the [usage guide](usage.md) for the corresponding calls.

## Annotation contents

Official [DALI data packs](http://ekhidna2.biocenter.helsinki.fi/dali/colab.html) include Pfam hits and their descriptions. DaliScope reads them from the pack itself; no external Pfam database is needed for analysis.

| File | Fields | Purpose |
| --- | --- | --- |
| `pfam.tsv` | `target_name`, `query_accession`, `e_value`, `env_from`, `env_to`, `database`, `clan` | Pfam hits; envelope positions are one-based |
| `mini_pfam_names.tsv` | `pfam`, `clan`, `clan_short`, `short`, `name` | Family and clan descriptions |
| `mini_ligand.tsv` | `pdb_id`, `chain_id`, `res_name`, `hetatm` | Ligand records for PDB targets |
| `mini_ligand_names.tsv` | `res_name`, `compound` | Ligand names |

The local preparation command requires a Pfam annotation database and bundles descriptions selected from its packaged reference (or an explicit `--pfam-names` replacement). The loader also accepts custom packs that omit these annotation files; in that compatibility mode, retained hits are labeled `Unassigned`. Ligand files are optional.

## Namespaces and provenance

Namespaces are `1=PDB`, `2=AFDB2`, `3=BFVD`, and `4=VIRO3D`. Alignments whose targets were not packed can be absent from derived views; compare summary and packed-target counts.

Record source/database versions, query and chain, search settings, generation date, and SHA-256 checksum. The frozen examples' identifiers and checksums are in `notebooks/data/datasets.json`; upstream database release dates were not recorded in the supplied packs.
