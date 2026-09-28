#!/usr/bin/env bash
# Compatibility wrapper. This file uses Bash despite its historical .csh suffix.
# Run --help to see the required local database configuration.
set -euo pipefail
exec python -m daliscope.packer.make_colab_pack "$@"
