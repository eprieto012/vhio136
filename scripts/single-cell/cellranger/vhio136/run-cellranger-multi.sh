#!/bin/bash
#SBATCH --job-name=cellranger_vhio136
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log

# Cell Ranger multi for one library of vhio136, with the settings that produced data/single-cell/rootpath/ and data/single-cell/vhio/:
#   Cell Ranger 7.1.0
#   GEX reference refdata-gex-GRCh38-2020-A, VDJ reference refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0
#   multiconfig_rootpath.csv: 8 libraries (pre_*, post2_d42_*), GEX + VDJ-T
#   multiconfig_vhio.csv:     4 libraries (post1_d14_*, post6_d140_*), GEX + Antibody Capture (TotalSeq-C) + VDJ-T
# Both declare the [feature] reference. RootPath has no antibody library, so its matrices carry the 138 antibody
# features with zero counts; keep the reference anyway: m01 reads the matrices as a list (toc$`Gene Expression`),
# which Read10X only returns when there is more than one feature type.
#
# Usage:
#   1. Copy the config of the batch (one per library) and fill in the [placeholders].
#   2. sbatch run-cellranger-multi.sh <library_id> <multiconfig.csv> [<data_dir>]
#      e.g. sbatch run-cellranger-multi.sh pre_bulk1 pre_bulk1.csv /path/to/vhio136/data/single-cell/rootpath/pre_bulk1
#           sbatch run-cellranger-multi.sh post1_pd1 post1_pd1.csv /path/to/vhio136/data/single-cell/vhio/post1_d14_pd1
# Cell Ranger writes <library_id>/ in the current directory. If <data_dir> is given, the files read by
# m01 are copied there (raw and filtered matrices, filtered contigs, clonotypes, web summary).

set -euo pipefail

# Cell Ranger 7.1.0 on PATH, or through the container used for data/ (Docker Hub litd/docker-cellranger:v7.1.0):
# singularity build cellranger-7.1.0.sif docker://litd/docker-cellranger@sha256:5b27ba5674db54cb9820d912cef3ba841dc582f137fba7c8353a91881ae1f469
# export CELLRANGER="singularity exec -B /path/to/data cellranger-7.1.0.sif cellranger"
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
