# DaliScope 0.1.5

DALI data packs include the Pfam annotations and descriptions needed for analysis. The README now links directly to the official self-contained pack service and removes the misleading separate-database requirement for pack users.

For local pack generation, include the complete Pfam 38.2 descriptions reference in `daliscope/packer/data/pfam_names.tsv`, with its fixed upstream source, conversion rule, checksums, and CC0-1.0 terms. The table is included in the wheel, source distribution, and source ZIP. The packer resolves it relative to its installed package instead of an external `pfamdata` directory. `--pfam-names` remains available for a different annotation snapshot; protein ledgers/shards and per-protein annotation databases remain explicit server inputs.

The reference reproduces all five fields of the 760 rows in the supplied original mini table. It is an official full reference, not a claimed byte-for-byte restoration of the unavailable private full file.

New regression checks exercise pack generation without an external description table from an unrelated directory, and standalone extraction of bundled family/clan descriptions. Existing explicit-reference generation remains supported.

README, Colab links/setup tags, and release downloads use the same 0.1.5 version. Scientific algorithms, tutorial computations, parameters, standardized notebook names, and frozen example packs are unchanged from 0.1.4. Earlier tags and release assets remain available. See the [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.5/docs/validation.md) for current packaging checks and historical computational/frontend evidence.
