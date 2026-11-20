# Bundled Pfam descriptions

`pfam_names.tsv` is the complete family and clan description table of
**Pfam 39.0**. It is the table that the DALI data-pack service uses for new
packs since 17 November 2026, so locally generated packs carry the same
descriptions. The packer uses this local reference by default; it needs no
separate download of description data.

The table contains **53,275 Pfam families** and eight tab-separated columns:

| Column | Content |
|---|---|
| `pfam` | Pfam accession |
| `clan` | Clan accession; the family's own accession when it has no clan |
| `name` | Family description (Pfam `DE`) |
| `short` | Family identifier (Pfam `ID`) |
| `type` | Family type (Pfam `TP`) |
| `clan_name` | Clan description; empty when the family has no clan |
| `clan_short` | Clan identifier; empty when the family has no clan |
| `clan_entries` | Number of families in the clan; 0 when the family has no clan |

DaliScope reads `pfam`, `clan`, `clan_short`, `short`, and `name`; the other
columns are kept so that the file matches the DALI server table. The families,
clans, identifiers, descriptions, and types are identical to the `pfam_entries`
table of the DALI server's Pfam 39.0 annotation database, which was parsed from
the Pfam 39.0 release files `Pfam-A.hmm.dat.gz` and `Pfam-C.gz`. Of the 1,074
clans, every `clan_entries` value equals the number of the clan's families in
this table. File SHA-256 (UTF-8, LF line endings):

```text
ef2d3ebd84a0fc4b0d1473db16931b87ca49fc1a9c864584ae0786540c40fb23
```

Compared with the Pfam 38.2 table bundled up to release 0.1.9, 23,179 families
are new, 38 families were retired, and 1,724 families have a changed clan,
identifier, or description. Packs made from a Pfam 38.2 annotation database,
including the example packs and the analyses in the DaliScope paper, should use
the [Pfam 38.2 table of release 0.1.9](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.9/daliscope/packer/data/pfam_names.tsv)
through `--pfam-names`. An existing DALI data pack already includes its
annotations and descriptions; loading it does not use this file.

Pfam data are distributed under **CC0-1.0**, as stated in the
[Pfam documentation](https://pfam-docs.readthedocs.io/en/latest/pfam.html).
The DaliScope code's MIT license does not replace those data terms. Cite Pfam
when using its annotations in an analysis.

For a different annotation snapshot, supply `--pfam-names /path/to/pfam_names.tsv`
with at least the columns `pfam`, `clan`, `clan_short`, `short`, and `name`. The
bundled descriptions do not contain per-protein annotations: generating a new
pack still requires the protein ledger/shards and the corresponding Pfam
annotation database (`--pfam-db`).
