#!/usr/bin/env bash
# Compile the extended abstract. latexmk runs pdflatex -> bibtex -> pdflatex
# as many times as the citations and references need, so this is the whole
# build. Run from anywhere:  bash docs/extended_abstract/build.sh
set -euo pipefail
cd "$(dirname "$0")"
latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error extended_abstract.tex
echo "built: $(pwd)/extended_abstract.pdf"
