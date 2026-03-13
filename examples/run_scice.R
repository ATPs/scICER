#!/usr/bin/env Rscript

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush(stdout())
}

args <- commandArgs(trailingOnly = TRUE)
input_qs <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/data1/xlab/researches/20250709_scICE/20260209_reviews/20260305_pancreatic/pancreatic_harmony1.qs"
}
output_rds <- if (length(args) >= 2) args[[2]] else "scice_results.rds"

qread_threads <- as.integer(Sys.getenv("SCICE_QREAD_THREADS", "40"))

log_step("Loading required packages...")
library(Seurat)
library(scICER)
library(qs)
log_step("Package loading complete.")

log_step(sprintf("Starting qs::qread from: %s", input_qs))
data2 <- qs::qread(input_qs, nthread = qread_threads)
log_step(sprintf("Finished qs::qread. Cells: %d, Features: %d", ncol(data2), nrow(data2)))

log_step("Starting scICE_clustering()")
scice_results <- scICE_clustering(
  object = data2,
  cluster_range = 2:10,
  remove_threshold = Inf,
  n_workers = 10,
  n_trials = 10,
  n_bootstrap = 10,
  seed = 123,
  verbose = TRUE,
  graph_name = "RNA_snn"
)
log_step("Finished scICE_clustering()")

saveRDS(scice_results, output_rds)
log_step(sprintf("Saved results to: %s", output_rds))
log_step("Run complete.")
