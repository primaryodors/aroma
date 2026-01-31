#!/bin/bash

PROT="$1"
LIG="$2"

make bin/aromadock bin/ic bin/cavity_search

python3 run.py "$@"
