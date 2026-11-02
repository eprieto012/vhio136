#!/bin/bash
#SBATCH --job-name=cellranger_yossef
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log

# Cell Ranger multi for one library of the Yossef et al. cohort (12 libraries, dbGaP FASTQs),
# with the settings that produced data/single-cell/yossef-cohort/:
#   Cell Ranger 9.0.0
#   GEX reference refdata-gex-GRCh38-2024-A, VDJ reference refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
#   GEX + VDJ libraries (chain auto-detected); create-bam,false
#
# Usage:
#   1. Copy multiconfig.csv (one per library) and fill in the [placeholders].
#   2. sbatch run-cellranger-multi.sh <library_id> <multiconfig.csv> [<data_dir>]
#      e.g. sbatch run-cellranger-multi.sh RY4246 RY4246.csv /path/to/vhio136/data/single-cell/yossef-cohort/RY4246
# Cell Ranger writes <library_id>/ in the current directory. If <data_dir> is given, the files read by
# m01 are copied there (raw and filtered matrices, filtered contigs, clonotypes, web summary).

set -euo pipefail

# Cell Ranger 9.0.0 on PATH, or through the container used for data/ (Docker Hub litd/docker-cellranger:v9.0.0):
# singularity build cellranger-9.0.0.sif docker://litd/docker-cellranger@sha256:91fe8317edca1ae48f73738e521b5baf872d2b3be0fd1b4db2707f50ebac22b2
# export CELLRANGER="singularity exec -B /path/to/data cellranger-9.0.0.sif cellranger"
CELLRANGER="${CELLRANGER:-cellranger}"

if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch $0 <library_id> <multiconfig.csv> [<data_dir>]" >&2
  exit 1
fi
id="$1"
config="$2"
data_dir="${3:-}"

$CELLRANGER multi --id="$id" --csv="$config" --localcores="${SLURM_CPUS_PER_TASK:-16}" --localmem=120

if [[ -n "$data_dir" ]]; then
  outs="$id/outs"
  mkdir -p "$data_dir"
  cp -r "$outs/multi/count/raw_feature_bc_matrix" "$data_dir/"
  cp -r "$outs/per_sample_outs/$id/count/sample_filtered_feature_bc_matrix" "$data_dir/"
  cp "$outs/per_sample_outs/$id/vdj_t/filtered_contig_annotations.csv" "$outs/per_sample_outs/$id/vdj_t/clonotypes.csv" "$data_dir/"
  cp "$outs/per_sample_outs/$id/web_summary.html" "$data_dir/"
fi
