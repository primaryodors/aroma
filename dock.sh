#!/bin/bash

PROT="$1"
LIG="$2"

make bin/aromadock bin/ic bin/cavity_search test/moiety_test bin/phew

python3 run.py "$@"
