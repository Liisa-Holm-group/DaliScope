# Example data

These are frozen DaliScope input packs supplied with the project. They contain search-result metadata, packed target C-alpha structures, sequences, secondary-structure codes, Pfam annotations, and query structures.

| Pack | Query | Query residues | Packed targets | Raw alignments | Distribution |
| --- | --- | ---: | ---: | ---: | --- |
| `3ubpC_PDB25.tar.gz` | 3ubpC | 570 | 1,048 | 1,120 | Included; quickstart |
| `dotp_ISS.tar.gz` | 9b8eA | 222 | 814 | 816 | Included; additional example |
| `ZN_full.tar.gz` | c9kdA | 1,724 | 1,573 | 1,584 | Included; motif case study |
| `GOLD.tar.gz` | 5tdqA | 123 | 12,245 | 12,316 | Optional release asset |

Raw alignments and packed targets differ because multiple alignments can map to a target and some targets were not retrieved during pack generation. See `datasets.json` for SHA-256 checksums and namespace IDs.

## GOLD download

The GOLD example is distributed with the `v0.1.1` GitHub release. Run from the repository root:

```bash
python scripts/download_example.py GOLD
```

The script checks SHA-256 against `datasets.json` before accepting the download. Before that release exists, use the original supplied GOLD pack in this directory. The large pack is intentionally excluded from regular Git commits.

## Provenance and terms

Structures originate from the [Protein Data Bank](https://www.wwpdb.org/about/usage) and the target databases identified by the pack namespaces (including AFDB2). Pfam family/clan annotations come from [Pfam/InterPro](https://www.ebi.ac.uk/interpro/). Query identifiers for case studies are recorded exactly as supplied; do not assume that every local identifier is a public PDB accession.

The upstream database release dates, original DALI run settings, and a complete generation ledger were not supplied. These examples support reproducing the packaged analysis rather than reconstructing the original searches. Add those records when producing new packs.

The repository's MIT license applies to DaliScope's own code and documentation. Third-party structures and annotations retain the source providers' terms and attribution requirements. Cite DALI and the relevant data sources when reporting results.
