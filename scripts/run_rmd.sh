#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=vhio136_rmd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=endikaprieto@vhio.net
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --output=logs/%x_%j.log

# Usage (submit from the project root, vhio136/):
# cd /mnt/bioinfnas/immuno/eprieto/analyses/vhio136
# sbatch -J m01_vhio scripts/run_rmd.sh scripts/single-cell/m01_data_load_vhio_latest.Rmd
# sbatch -J m01_rootpath scripts/run_rmd.sh scripts/single-cell/m01_data_load_rootpath_latest.Rmd
# sbatch -J m01_yossef scripts/run_rmd.sh scripts/single-cell/m01_data_load_yossef_latest.Rmd

set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ ! -x "$REPO/scripts/interactive-R.sh" ]]; then
  echo "Submit from the project root (vhio136/): $REPO" >&2
  exit 1
fi

if [[ $# -ne 1 || ! -f "$1" ]]; then
  echo "Usage: sbatch -J <name> scripts/run_rmd.sh scripts/<analysis>/<file>.Rmd" >&2
  exit 1
fi
RMD="$(realpath "$1")"

exec "$REPO/scripts/interactive-R.sh" --no-save --no-restore -s -e "
  code <- tempfile(fileext = '.R')
  invisible(knitr::purl('$RMD', output = code, quiet = TRUE, documentation = 0))
  pdf(NULL)
  source(code, echo = TRUE, max.deparse.length = Inf)
"