source("renv/activate.R")
# Never write Rplots.pdf: in non-interactive sessions (Rscript, run-rmd.sh, srun) and in interactive sessions
# without a display (interactive-R.sh on a compute node) the default graphics device is a null PDF
if (!interactive() || !nzchar(Sys.getenv("DISPLAY"))) options(device = function(...) grDevices::pdf(file = NULL, ...))
