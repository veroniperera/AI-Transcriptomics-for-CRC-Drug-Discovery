#!/usr/bin/env bash
set -euo pipefail

for f in R/0[1-6]*.R; do Rscript "$f"; done
for f in python/0[1-6]*.py; do python "$f"; done

python python/07_dude_enrichment.py analyze
