# DaliScope 0.1.2

This release adds a browser entry through Google Colab and makes the local Python setup consistent across the README and tutorials.

- Open any of the four tutorials directly in Colab. Their first cell clones the matching version, installs the lightweight `colab` extra once per runtime, and enables custom widgets.
- Allow IPython 7.34, matching the version declared by Google's Colab package. Preserve Colab-managed dependencies and already-installed numerical libraries during setup instead of installing a local Jupyter server/kernel.
- Select Colab's Plotly renderer in the cloud and retain the local renderer outside Colab.
- Download the optional GOLD release attachment automatically in its Colab tutorial, checking the unchanged frozen-data SHA-256 before analysis.
- Document CPU runtimes, output downloads, runtime resets, and the Python 3.12 local setup. Conda is an optional alternative to venv.
- Reject stale, modified, or mixed source checkouts, and add dependency-compatibility coverage to Ubuntu CI.
- Declare Jinja2 for table previews and use asynchronous notebook validation so cell timeouts are enforced.

The scientific algorithms, example data, domain conventions, and analysis parameters are unchanged from 0.1.1. See [validation](https://github.com/Liisa-Holm-group/DaliScope/blob/v0.1.2/docs/validation.md) for the executed checks. Hosted Colab execution and browser mouse interaction require a Google login and are not yet verified; Linux compatibility checks are recorded separately.

Earlier releases remain available: [0.1.1](https://github.com/Liisa-Holm-group/DaliScope/releases/tag/v0.1.1) and [0.1.0](https://github.com/Liisa-Holm-group/DaliScope/releases/tag/v0.1.0).
