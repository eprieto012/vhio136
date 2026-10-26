#!/usr/bin/env bash
# =============================================================================
# sbatch-run-m01.sh — submit the m01 analyses (rootpath and vhio) to Slurm
# =============================================================================
#
# Each .Rmd is run on a compute node through interactive-R.sh (same conda R + renv
# library as interactive sessions). The two batches are independent, so both
# jobs can run at the same time.
#
# Usage (from any directory):
#   scripts/sbatch-run-m01.sh
#
# Follow the jobs:
#   squeue -u $USER
#   tail -f results/vhio136/<batch>/m01_<batch>_<jobid>.log
#
# Outputs per batch (results/vhio136/<batch>/):
#   seurat/  plots/          objects and PDFs written by the .Rmd itself
#   knit/                    knitr report (.md) and chunk figures
#   m01_<batch>_<jobid>.log  console output and errors

set -euo pipefail

# Run from the repository root, wherever this script is called from.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

# Submit one .Rmd for one batch.
#   $1  batch name (rootpath / vhio): sets job name and output folder
#   $2  .Rmd file in scripts/
#
# knitr::knit is used instead of rmarkdown::render because the compute nodes
# have no pandoc. Paths inside R start with "../" because interactive-R.sh starts R
# from renv_vhio136/.
submit_rmd() {
  local batch="$1" rmd="$2"
  local out="results/vhio136/$batch"
  mkdir -p "$out/knit"

  sbatch -p short -c 8 --mem=120G -t 10:00:00 \
    -J "m01_$batch" -o "$out/m01_${batch}_%j.log" \
    --wrap "scripts/interactive-R.sh -s -e 'knitr::opts_chunk\$set(fig.path = \"../$out/knit/\"); knitr::knit(\"../scripts/$rmd\", output = \"../$out/knit/m01_$batch.md\")'"
}

submit_rmd rootpath m01_data_load_rootpath_latest.Rmd
submit_rmd vhio     m01_data_load_vhio_latest.Rmd
