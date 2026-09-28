# DaliScope 0.1.3

This patch addresses setup and figure-lifecycle failures found in the comprehensive Google Colab tutorial. Matplotlib 3.10 could cache its available backends before setup installed ipympl, causing `%matplotlib widget` to reject the backend even though ipympl was present. Later display/profile drawing and cleanup could also close an earlier widget's figure, leaving its controls responsive without redrawing the plot.

- Refresh Matplotlib's cached backend discovery after live pip installation so the tutorial can request ipympl without recreating the kernel.
- Keep the widget backend active through the comprehensive workflow and limit later display/profile cleanup to newly created figures, preserving figures already owned by interactive widgets.
- Use generic `python3` kernel metadata in all four notebooks to avoid Colab's unrecognized `daliscope` runtime notice. Local instructions still register and explicitly select **Python (DaliScope)**.
- Keep source checkout, tutorial links, and optional GOLD release selection on the same v0.1.3 tag.

The scientific algorithms, example data, domain conventions, and analysis parameters are unchanged from 0.1.2. Windows and Linux each passed 47 tests. The 0.1.3 comprehensive workflow completed in a fresh hosted Colab runtime with automatic setup; switching the hit-plane view after all later analysis visibly redrew its axes and data. See the [validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.3/docs/validation.md) for the tested source and exact scope.

Version [0.1.2](https://github.com/Liisa-Holm-group/DaliScope/releases/tag/v0.1.2) introduced the Colab entry point, automatic setup and GOLD download, preserved cloud numerical dependencies, and consistent local Python 3.12 setup. Its immutable [release-time validation](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.2/docs/validation.md) remains available. After publication on 28 September 2026, its quickstart completed all eight original cells on hosted Colab, including output creation, session restoration, and visible py3Dmol rendering; its comprehensive workflow exposed the backend cache failure fixed by this patch. Domain Apply/Reset and centering controls were subsequently verified; direct py3Dmol mouse rotation/zoom remains unverified. The [current validation record](https://github.com/Liisa-Holm-group/DaliScope/blob/main/docs/validation.md) records post-release checks without altering earlier tags or release assets.

Earlier releases remain available: [0.1.1](https://github.com/Liisa-Holm-group/DaliScope/releases/tag/v0.1.1) and [0.1.0](https://github.com/Liisa-Holm-group/DaliScope/releases/tag/v0.1.0).
