# DaliScope 0.1.1

This patch corrects session persistence, invalid data-pack input handling, and tutorial instructions discovered during independent review.

- Store views and models with unique internal file identifiers, preserving distinct display names such as `domain 1` and `domain_1`. Existing 0.1.0 session manifests remain readable.
- Validate generated query CA coordinates against query metadata before publishing a pack. Invalid chains or lengths fail without replacing an existing output.
- Convert one-based inclusive UI/PDB domain ranges to zero-based half-open clipping ranges. The GOLD tutorial now uses the stated PDB residues 23–123 consistently for clipping and community analysis.
- Select GOLD tutorial representatives from current annotation counts instead of historical module numbers, and show current case-study outputs.
- Correct custom-filter column names and describe query coverage using the full query length. Coverage alone does not establish domain classification.
- Describe the executed workflow as domain-clipped superimposition followed by occupancy-based filtering. Filtering does not perform an occupancy-weighted refit.
- Include the Windows reference requirements in the source distribution and use flat Release attachment names in the checksum manifest.

The corrected range conversion can change selected populations, motifs, and community assignments relative to 0.1.0. Record the software version and ranges with analysis results. Module identifiers are run-specific; inspect current annotations and structures before interpreting them. The underlying superimposition, motif, and clustering algorithms are unchanged.

See [validation](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.1/docs/validation.md) for executed checks and their limits. The original example data packs retain their recorded checksums. [DaliScope 0.1.0](https://github.com/Liisa-Holm-group/DaliScope/releases/tag/v0.1.0) remains available for reproducing earlier runs.