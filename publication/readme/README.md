# README plotting example: ZhengMix_4757

These figures are generated with the current scICER source using a real benchmark
dataset from the paper. The example contains 4,688 retained cells; `ZhengMix_4757` is the original dataset name.

The source Seurat object is archived at:

```
/data1/xlab/researches/20250709_scICE/20250709_data_download/20250806_R_scICER_sclensJuliaUMAPinR/ZhengMix_4757/seurat.scICER.qs
```

The root-level analysis script `20250806_run_R_sclensJuliaUMAPinR.R` rebuilt a
UMAP fuzzy graph from the saved scLENS-derived embedding using `uwot::umap()`
with cosine distance, 15 neighbors and `min_dist = 0.1`.
Its `umap_graph` graph, `umap2` coordinates and `celltype` annotations are preserved.
For the example file, the coordinates are named `umap` and only two original
count features are retained. The object is intended for graph-based clustering
and plotting; use the original full expression data for marker analysis.

The example RDS contains `object` (the lightweight Seurat object), `results` (the
scICER result), and `parameters`. The analysis settings are:

```r
graph_name = "umap_graph"
cluster_range = 2:15
n_workers = 4
n_trials = 15
n_bootstrap = 100
min_cluster_size = 2
remove_threshold = Inf
seed = 123
```

This is a new seeded run on the archived graph, so the figures illustrate current
package usage. They do not replay the original paper run, which used different
worker settings and did not fix a seed for this workflow. For repeatability, keep the
configuration, including `n_workers`, and software versions fixed.

From the repository root, rebuild the example and figures with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript examples/generate_readme_figures.R \
  /data1/xlab/researches/20250709_scICE/20250709_data_download/20250806_R_scICER_sclensJuliaUMAPinR/ZhengMix_4757/seurat.scICER.qs
```

Or redraw the figures from the bundled result without rerunning clustering:

```bash
Rscript examples/generate_readme_figures.R
```

The script uses existing installations of pkgload, Seurat, qs, ggplot2 and ragg.
All files in this directory are excluded from R package builds by `.Rbuildignore`.

Other root-level scripts used as graph-method references were
`20250729_run_R_scICER.R` (RNA SNN), `20250730_run_R_knn.R` (RNA KNN),
`20250804_run_R_scICER.R` (SCT SNN), `20250804_run_R_scICER.SCT.knn.R`
(SCT KNN), `20250716_run_R_sclens.R` (scLENS SNN) and
`20250722_run_R_sclens_knn.R` (scLENS KNN). Each selects a graph independently
of the reduction used for visualization.

In this run, the returned cluster numbers were 3–7 and 9–15; requested targets
2 and 8 were not recovered by the resolution search. At threshold 1.005,
K = 3, 4, 5 and 6 were stable. The histogram illustrates K = 9. The UMAP illustrates
K = 6 as one stable candidate, with biological interpretation left to the user.

Generated with R 4.5.2, scICER 1.1.2, Seurat 6.0.0 and ggplot2 4.0.3.
