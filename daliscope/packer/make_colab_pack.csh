#!/usr/bin/env bash
# Compatibility wrapper. This file uses Bash despite its historical .csh suffix.
# Uses bundled Pfam descriptions; run --help for local annotation/protein database inputs.
set -euo pipefail
exec python -m daliscope.packer.make_colab_pack "$@"
