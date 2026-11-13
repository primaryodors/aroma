#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: ./dock.sh RECEPTOR ODORANT [options]"
    echo "Example: ./dock.sh OR1A1 linalool"
    exit 2
fi

PROT="$1"
LIG="$2"

make bin/aromadock bin/ic bin/cavity_search test/moiety_test bin/phew bin/ringflip

# Batch ligands follow the receptor page list, with the top ligand first.
exec python3 run.py "$@"
