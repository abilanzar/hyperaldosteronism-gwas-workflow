#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

Rscript scripts/fixed_effect_meta.R \
  --manifest config/cohorts.example.tsv \
  --analysis european \
  --out results/european_meta.tsv \
  --qc-out results/european_qc.tsv

Rscript scripts/fixed_effect_meta.R \
  --manifest config/cohorts.example.tsv \
  --analysis multi \
  --out results/multi_ancestry_meta.tsv \
  --qc-out results/multi_ancestry_qc.tsv

Rscript -e '
  eur <- read.delim("results/european_meta.tsv")
  multi <- read.delim("results/multi_ancestry_meta.tsv")
  eur_qc <- read.delim("results/european_qc.tsv")
  stopifnot(nrow(eur) == 4L)
  stopifnot(nrow(multi) == 4L)
  stopifnot(all(eur$K == 4L))
  stopifnot(all(multi$K == 6L))
  stopifnot(any(eur_qc$ACTION == "swapped"))
  stopifnot(any(eur_qc$ACTION == "strand_complement"))
  stopifnot(all(is.finite(eur$BETA)))
  stopifnot(all(eur$P >= 0 & eur$P <= 1))
  message("Synthetic workflow smoke test passed.")
'

