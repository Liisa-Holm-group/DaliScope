# Worked-example outputs

The saved pages show the paper's DaliScope v0.1.9 analyses with Pfam 38.2. Both come from the notebook outputs used for the supplementary PDFs.

| Worked example | Data pack | Result shown |
| --- | --- | --- |
| [Rat KDM3B (ST1)](https://liisa-holm-group.github.io/DaliScope/example-outputs/WorkedExample1.html) | `ZN_full.tar.gz` | 1,173 targets in the Domain 5 subset. The final diagram shows the structural and sequence features identified in the query. |
| [GOLD domain (ST2)](https://liisa-holm-group.github.io/DaliScope/example-outputs/WorkedExample2.html) | `GOLD.tar.gz` | 6,465 targets and the 23-community partition used for Figure 2. |

For ST2 and Figure 2b, the clan and family heatmaps group rows by the module with the highest count, then by decreasing peak count. Ties between peak modules use the lowest module number. This display uses the saved community assignments and the same family cutoff (more than 20 targets per family in a module). The tagged v0.1.9 notebook displays the rows in its original order. The plotting function on `main` uses the order in ST2.

[Open the interactive notebooks](../colab.md) to change inputs and operate Python-backed controls. The HTML pages include saved plots and text. Embedded 3D views depend on WebGL and their external viewer.

[outputs.json](outputs.json) records data-pack and HTML file checksums. The earlier [September HTML exports](https://github.com/Liisa-Holm-group/DaliScope/tree/dae7eaf72f058eecdf9f04488c9782f1b16fd530/docs/example-outputs) remain available with their provenance.
