# Bundled Pfam descriptions

`pfam_names.tsv` is the complete family and clan description table from the
fixed [Pfam 38.2 release](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/).
The packer uses this local reference by default; it needs no separate download
of description data. The table contains **30,134 distinct Pfam families** and
five columns: `pfam`, `clan`, `clan_short`, `short`, and `name`.

Source: [Pfam-A.clans.tsv.gz](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.clans.tsv.gz)
(553,777 bytes). Source SHA-256:

```text
86062b7ef1a0e0caee0c28cef479ac0d294c80789c51cbf75e513b368ce3a6f6
```

To produce the bundled file, decompress the source as UTF-8, add the five-column
header above, and replace an empty `clan` field with that row's Pfam accession.
Keep `clan_short` empty when it is empty upstream, and preserve every other
field, whitespace, and row order. Write the result as tab-separated UTF-8 with
LF line endings. Bundled file SHA-256:

```text
79bc7f952518b525f08f64394616e068a2380e9afa3a9934c30e843f35b29f8e
```

The derived values exactly reproduce all five fields of the 760 rows in the
original DaliScope `mini_pfam_names.tsv` supplied with the project. The unavailable
private full `pfam_names.tsv` is not claimed to be an identical copy of this file.
The existing example packs are unchanged.

Pfam data are distributed under **CC0-1.0**, as stated in the release's
[copyright notice](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/relnotes.txt)
and [Pfam documentation](https://pfam-docs.readthedocs.io/en/latest/pfam.html).
The DaliScope code's MIT license does not replace those data terms. Cite Pfam
when using its annotations in an analysis.

For a different annotation snapshot, supply `--pfam-names /path/to/pfam_names.tsv`
with the same five-column header. The bundled descriptions do not contain
per-protein annotations: generating a new pack still requires the protein
ledger/shards and the corresponding Pfam annotation database (`--pfam-db`).
An existing DALI data pack already includes its annotations and descriptions.
