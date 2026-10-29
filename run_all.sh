#!/usr/bin/env bash
# Full pipeline: genotype prep -> GWAS (FarmCPU, BLINK, MLM) -> consensus -> candidate genes -> Manhattan plots
set -euo pipefail
cd "$(dirname "$0")"

Rscript scripts/01_prepare_genotypes.R
Rscript scripts/02_run_gwas.R
Rscript scripts/03_consensus_markers.R
Rscript scripts/04_candidate_genes.R
Rscript scripts/05_manhattan_plots.R
