#!/bin/bash
#SBATCH -p long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%A_%a.log

# Align the TCR library screen FASTQs to the library and count reads per construct (README, section 5.1).
#
# 1. Create the tcrlib_hits environment and point ENV_HITS to it (README, section 1).
# 2. Place the FASTQs as $RAW_DIR/<screen>/<sample>_{1,2}.fastq.gz, with <screen> and <sample> as in
#    data/aux/tcr-library/tcr_library_screen_samples.tsv (the original NAS layout, $RAW_DIR/<fastq_dir>/, also works).
# 3. From the project root (vhio136/), submit one SLURM array task per sample of a screen:
#      export RAW_DIR=/path/to/fastq ENV_HITS=/path/to/envs/tcrlib_hits
#      scripts/tcr-library/hits/run-hits.sh screen-1    # 14 samples, sequenced 2022
#      scripts/tcr-library/hits/run-hits.sh screen-2    # 12 samples, sequenced 2024
#    SLURM resources can be changed with PARTITION, CPUS, MEM and WALLTIME (defaults: long, 16, 32gb, 2-00:00:00).
#
# Output: results/tcr-library/<screen>/<sample>/<sample>.counts.txt, plus the SAM, the fastp report, the log and
# tmp/ (intermediate FASTQs). Only the counts are used downstream; the SAM and tmp/ are kept and can be deleted by hand.
# Log of each task: logs/hits_<screen>_<job>_<task>.log. Next step: the two Rmd of README section 5.2.

set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
HITS="$REPO/scripts/tcr-library/hits"
ENV_HITS="${ENV_HITS:-/mnt/bioinfnas/immuno/eprieto/envs/tcrlib_hits}"
RAW_DIR="${RAW_DIR:-/mnt/petasan_immuno/raw_data/itag/VHIO/360RIO136}"
PARTITION="${PARTITION:-long}"
CPUS="${CPUS:-16}"
MEM="${MEM:-32gb}"
WALLTIME="${WALLTIME:-2-00:00:00}"
SAMPLES="$REPO/data/aux/tcr-library/tcr_library_screen_samples.tsv"
LIB="$REPO/data/aux/tcr-library/TCR_lib_RIO136_AAV.fasta"
export PATH="$ENV_HITS/bin:$PATH"

if [[ $# -ne 1 ]]; then
  echo "Usage: scripts/tcr-library/hits/run-hits.sh <screen> (screen-1 | screen-2)" >&2
  exit 1
fi
SCREEN="$1"
if [[ ! -f "$SAMPLES" || ! -x "$ENV_HITS/bin/bowtie2" ]]; then
  echo "Submit from the project root (vhio136/); needs $SAMPLES and $ENV_HITS" >&2
  exit 1
fi
rows() { awk -F'\t' -v s="$SCREEN" 'NR > 1 && $1 == s' "$SAMPLES"; }
n=$(rows | wc -l)
if [[ "$n" -eq 0 ]]; then
  echo "No samples for screen '$SCREEN' in $SAMPLES" >&2
  exit 1
fi

# --- Submission mode (outside a SLURM job): index the library once, then one array task per sample ---
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  python3 "$HITS/bowtie2index.py" "$LIB"
  mkdir -p "$REPO/logs" "$REPO/results/tcr-library/$SCREEN"
  exec sbatch -J "hits_$SCREEN" -p "$PARTITION" --cpus-per-task="$CPUS" --mem="$MEM" --time="$WALLTIME" --array=1-"$n" "$0" "$SCREEN"
fi

# --- Task mode: sample number SLURM_ARRAY_TASK_ID of this screen ---
read -r sample fastq_dir < <(rows | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f2,3)
fastq_path="$RAW_DIR/$SCREEN"
[[ -f "$fastq_path/${sample}_1.fastq.gz" ]] || fastq_path="$RAW_DIR/$fastq_dir"
r1="$fastq_path/${sample}_1.fastq.gz"
r2="$fastq_path/${sample}_2.fastq.gz"
if [[ ! -f "$r1" || ! -f "$r2" ]]; then
  echo "FASTQs of $SCREEN / $sample not found in $RAW_DIR/$SCREEN/ or $RAW_DIR/$fastq_dir/" >&2
  exit 1
fi
echo "$SCREEN / $sample: $r1 $r2"
cd "$REPO/results/tcr-library/$SCREEN"
python3 "$HITS/align_sample.py" "$r1" "$r2" "$LIB" -o "$sample" -t "$SLURM_CPUS_PER_TASK"
test -s "$sample/$sample.counts.txt"
echo "Counts written: results/tcr-library/$SCREEN/$sample/$sample.counts.txt"
