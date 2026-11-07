#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: ./dock.sh RECEPTOR ODORANT [options]"
    echo "Example: ./dock.sh OR1A1 linalool"
    exit 2
fi

PROT="$1"
LIG="$2"

make bin/aromadock bin/ic bin/cavity_search test/moiety_test bin/phew

python3 run.py "$@"
