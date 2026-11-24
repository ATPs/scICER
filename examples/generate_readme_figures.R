#!/usr/bin/env Rscript
# From the repository root, pass the archived ZhengMix_4757 Seurat .qs file.
# With no argument, redraw figures from publication/readme/zhengmix4757-example.rds.
pkgload::load_all('.', quiet = TRUE)
library(Seurat)
library(ggplot2)

output_dir <- file.path('publication', 'readme')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
example_file <- file.path(output_dir, 'zhengmix4757-example.rds')
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop('Supply at most one archived Seurat .qs file.')

if (length(args) == 1L) {
  archived <- qs::qread(args[[1]], nthreads = 2)
  # Keep the original graph, UMAP and annotations; retain just two count features.
  counts <- GetAssayData(archived, assay = 'RNA', layer = 'counts')
  object <- CreateSeuratObject(CreateAssayObject(counts[1:2, , drop = FALSE]),
                             meta.data = archived[[]][, 'celltype', drop = FALSE])
  object[['umap_graph']] <- archived[['umap_graph']]
  object[['umap']] <- archived[['umap2']]
  rm(archived, counts)
  gc()
  parameters <- list(
    graph_name = 'umap_graph', cluster_range = 2:15, n_workers = 4,
    n_trials = 15, n_bootstrap = 100, min_cluster_size = 2,
    remove_threshold = Inf, seed = 123, verbose = TRUE
  )
  results <- do.call(scICE_clustering, c(list(object = object), parameters))
  saveRDS(list(object = object, results = results, parameters = parameters),
          example_file, compress = 'xz')
}

example <- readRDS(example_file)
object <- example$object
results <- example$results
stopifnot(ncol(object) == 4688L, 6L %in% results$n_cluster[results$ic < 1.005],
          9L %in% results$n_cluster)
set.seed(123)  # Fix jitter placement in plot_ic().
ic_plot <- plot_ic(results, threshold = 1.005, title = 'ZhengMix_4757: clustering consistency')
stability_plot <- plot_stability(results, cluster_number = 9)
object <- get_robust_labels(results, threshold = 1.005, return_seurat = TRUE, object = object)
cluster_plot <- DimPlot(object, reduction = 'umap', group.by = 'clusters_6', label = TRUE)

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(output_dir, name), plot = plot, device = ragg::agg_png,
         width = width, height = height, units = 'in', dpi = 160, bg = 'white')
}
save_figure(ic_plot, 'plot-ic.png', 9, 5.5)
save_figure(stability_plot, 'plot-stability.png', 7, 4.5)
save_figure(cluster_plot, 'cluster-labels.png', 7, 5.5)
labels <- get_robust_labels(results, threshold = 1.005)
stopifnot(setequal(labels$cell_id, Cells(example$object)), 'clusters_6' %in% names(labels))
print(extract_consistent_clusters(results, threshold = 1.005))
cat('Figure generation completed\n')
