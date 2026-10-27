#!/usr/bin/env bash
# Start R with the project environment: conda R + renv library (renv_vhio136/).
# Usage (not on the front end; get a compute node first):
#   srun --partition=highmem --nodes=1 --ntasks=1 --cpus-per-task=25 --mem=100G \
#        --nodelist=bioinf3.vhio.org --time=0-8:00:00 --pty /bin/bash -i
#   scripts/interactive-R.sh
# R starts in renv_vhio136/, so project paths are "../setup", "../results", ...

ENV_REPRO="${ENV_REPRO:-/mnt/bioinfnas/immuno/eprieto/Rmagic/miniforge3/envs/renv_vhio136}"
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../renv_vhio136" && pwd)" || exit 1

export PATH="$ENV_REPRO/bin:$PATH"
export PKG_CONFIG_LIBDIR="$ENV_REPRO/lib/pkgconfig:$ENV_REPRO/share/pkgconfig"
export PKG_CONFIG_PATH="$PKG_CONFIG_LIBDIR"
# en_US.UTF-8 is not installed on the nodes; R already fell back to "C"
export LC_ALL=C.UTF-8

# renv_vhio136/.Rprofile activates renv; R_LIBS* unset so only the renv library is used.
cd "$PROJECT"
exec env -u R_LIBS -u R_LIBS_SITE -u R_LIBS_USER "$ENV_REPRO/bin/R" "$@"