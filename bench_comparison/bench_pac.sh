#!/bin/bash
# Compatibility wrapper for the renamed PACK* benchmark entry point.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/bench_packstar.sh" "$@"
