# --- Project paths ---
if (!exists("setup_dir", inherits = TRUE) || !dir.exists(setup_dir)) {
  stop("setup_dir must point to the project's setup directory before sourcing setup.R.")
}
setup_dir <- normalizePath(setup_dir, mustWork = TRUE)
project_dir <- normalizePath(file.path(setup_dir, ".."), mustWork = TRUE)
datapath_load <- file.path(project_dir, "data")
aux_dir <- file.path(datapath_load, "aux")  # auxiliary input tables (TCR library, repertoires, databases, signatures)
results_dir <- file.path(project_dir, "results")
datapath <- results_dir
if (!dir.exists(datapath_load)) {
  stop("Data directory not found: ", datapath_load)
}

# --- R Markdown options ---
knitr::opts_chunk$set(echo = FALSE)

# --- Thread limits ---
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1
)

# --- Future globals ---
options(future.globals.maxSize = 10 * 1024^3)

# --- Reproducibility ---
SEED <- 123
PDFCLUSTER_CD4_SEED <- 123
PDFCLUSTER_CD8A_SEED <- 532
REPRO_NCORES <- 40
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

# --- Configure Python environment and verify MAGIC availability ---
library(reticulate)
reticulate::use_python(
  Sys.getenv("MAGIC_PYTHON", "/mnt/bioinfnas/immuno/eprieto/Rmagic/miniforge3/envs/magic_env/bin/python"),
  required = TRUE
)
reticulate::py_config()
stopifnot(reticulate::py_module_available("magic"))
reticulate::import("magic")

# --- Set random seeds in Python for reproducibility ---
reticulate::py_run_string(sprintf("
import numpy as np
import random
np.random.seed(%d)
random.seed(%d)
", SEED, SEED))

# --- Shared libraries, plotting styles, and functions ---
source(file.path(setup_dir, "libraries.R"))
source(file.path(setup_dir, "colors.R"))
source(file.path(setup_dir, "functions.R"))