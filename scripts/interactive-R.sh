#!/usr/bin/env bash
# =============================================================================
# interactive-R.sh — start an interactive R session in the locked environment
# =============================================================================
#
# The environment has two layers, and this script connects them:
#
#   1. R itself (R 4.4.0 + system libraries)  -> conda env "renv_vhio136"
#   2. R packages (Seurat, ProjecTILs, ...)    -> renv library in renv_vhio136/
#
# Usage
# -----
#   Terminal:
#     scripts/interactive-R.sh
#
#   VS Code: set "r.rterm.linux" to this file, then "R: Create R Terminal".
#
#   Another machine: point to your own environments, e.g.
#     ENV_REPRO=/path/to/conda_env MAGIC_PYTHON=/path/to/magic_env/bin/python scripts/interactive-R.sh

# --- 1. Where is R? ------------------------------------------------------------
# Conda environment that provides R. The value after ":-" is the default;
# an ENV_REPRO variable already set in your shell takes precedence.
ENV_REPRO="${ENV_REPRO:-/mnt/bioinfnas/immuno/eprieto/Rmagic/miniforge3/envs/renv_vhio136}"

# --- 2. Where is the renv project? --------------------------------------------
# renv_vhio136/ sits next to scripts/. It is located relative to this file
# (BASH_SOURCE), so the script works from any directory. This only works when
# the file is executed, not when its lines are pasted into a terminal.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="${PROJECT:-$(cd "$SCRIPT_DIR/../renv_vhio136" && pwd)}"

# --- 3. Use the conda R and its libraries --------------------------------------
# Conda binaries go first on PATH, and pkg-config looks at the conda libraries
# in case a package ever has to be compiled.
export PATH="$ENV_REPRO/bin:$PATH"
export PKG_CONFIG_LIBDIR="$ENV_REPRO/lib/pkgconfig:$ENV_REPRO/share/pkgconfig"
export PKG_CONFIG_PATH="$PKG_CONFIG_LIBDIR"

# --- 4. Activate renv -----------------------------------------------------------
# R reads the .Rprofile of the directory it starts in: renv_vhio136/.Rprofile
# runs renv/activate.R, which switches R to the project's package library.
# R's working directory is therefore renv_vhio136/.
cd "$PROJECT" || { echo "renv project not found: $PROJECT" >&2; exit 1; }

# --- 5. Start R -----------------------------------------------------------------
# Unset R_LIBS* so that no package library from outside the project can be
# loaded by mistake: only the renv library (and base R) remain visible.
# "exec" replaces this shell with R; "$@" forwards any arguments to R
# (VS Code adds --no-save --no-restore; sbatch-run-m01.sh passes -e '...').
exec env -u R_LIBS -u R_LIBS_SITE -u R_LIBS_USER "$ENV_REPRO/bin/R" "$@"