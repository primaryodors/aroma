#!/bin/bash

PROT="$1"
LIG="$2"

make bin/reec

python3 run.py "$PROT" "$LIG"
