source("renv/activate.R")
# Never write Rplots.pdf: in non-interactive sessions (Rscript, run_rmd.sh, srun) the default graphics device is a null PDF
if (!interactive()) options(device = function(...) grDevices::pdf(file = NULL, ...))
