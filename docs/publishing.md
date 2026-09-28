# Publishing DaliScope to the group organization

## Create the repository

Create an **empty** repository named `DaliScope` under `Liisa-Holm-group`. Leave automatic README, license, and .gitignore creation disabled: these files are supplied here. The agreed software license is MIT with the copyright holder Liisa Holm group. Review the original scientific case-study interpretations before public release.

The existing working directory is associated with a different origin. Keep that history and remote intact. Publish the clean exported release directory as the new repository.

## Publish the prepared snapshot

From the exported release directory:

```bash
git init -b main
git add .
git commit -m "Prepare DaliScope with installation, tutorials, and validation"
git remote add origin https://github.com/Liisa-Holm-group/DaliScope.git
git push -u origin main
```

The contact is `hao.liu@helsinki.fi`. The repository includes an MIT LICENSE and a CITATION.cff software citation. A manuscript describing DaliScope is being prepared for submission. Update this status after submission and add the paper citation when its details become available.

## Distribute the optional GOLD data

`GOLD.tar.gz` is about 50 MiB and is excluded from regular Git commits. Once the code is published, attach the original pack as a GitHub Release asset named `GOLD.tar.gz` on tag `v0.1.0`. Its checksum is in `notebooks/data/datasets.json`.

With authenticated GitHub CLI, from the prepared repository:

```bash
gh release create v0.1.0 --repo Liisa-Holm-group/DaliScope --title "DaliScope 0.1.0" --notes-file docs/release-notes.md
gh release upload v0.1.0 /path/to/original/GOLD.tar.gz --repo Liisa-Holm-group/DaliScope
```

After that release exists, `python scripts/download_example.py GOLD` downloads the verified file for the advanced notebook. Until then, retain the supplied local GOLD file and copy it into `notebooks/data/` when running that example.

## Final checks

- Confirm the repository owner/name and chosen visibility.
- Confirm MIT and group copyright, contact, and scientific descriptions.
- Verify bundled structures/annotations and their provenance with the group; third-party data keeps its original terms.
- Check GitHub Actions results after the push; local Windows checks do not establish Linux or other Python-version support.
- Open the quickstart in a browser and check 3D/widget interactions.
- Keep the original project backup and laboratory databases outside the new repository.
