#!/bin/bash

PROT="$1"
LIG="$2"

make bin/aromadock

python3 run.py "$PROT" "$LIG"
