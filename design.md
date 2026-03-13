# scICE_clustering Design Document

## 1. Scope

This document explains the current implementation of `scICE_clustering()` in **scICER**.

It covers:

- every public parameter (meaning, where it is used, and what happens if you change it),
- the full runtime workflow from Seurat graph input to final outputs,
- the internal function call chain and behavior of each major helper,
- runtime, parallelism, memory, and failure modes.

This description matches the code on `main` at commit `0d07f1d`.

## 2. Quick Mental Model

`scICE_clustering()` is a multi-stage pipeline:

1. Validate inputs and initialize runtime context.
2. Extract a Seurat graph and convert it to `igraph`.
3. For each target cluster number `k`, search a gamma (resolution) range likely to produce that `k`.
4. Optionally filter unstable `k` values.
5. For each remaining `k`, run intensive optimization:
   - evaluate many gamma values,
   - keep only gammas whose median cluster count matches target `k`,
   - compute IC/MEI stability metrics,
   - select `best_labels`.
6. Return an `scICE` object and annotate which `k` values are consistent under `ic_threshold`.

Important: current implementation uses **effective cluster counting** when `min_cluster_size > 1`.

- Effective cluster count = number of clusters with size `>= min_cluster_size`.
- If all clusters are smaller than threshold, effective count is `0`.
- Final small-cluster merge is applied **once** and only to `best_labels`.

## 3. Public API: Parameter-by-Parameter

Function signature:

```r
scICE_clustering(
  object,
  graph_name = NULL,
  cluster_range = 1:20,
  n_workers = 10,
  n_trials = 15,
  n_bootstrap = 100,
  seed = NULL,
  beta = 0.1,
  n_iterations = 10,
  max_iterations = 150,
  ic_threshold = Inf,
  objective_function = "CPM",
  remove_threshold = 1.15,
  min_cluster_size = 2,
  resolution_tolerance = 1e-8,
  verbose = TRUE
)
```

### 3.1 `object`

- Meaning: input Seurat object.
- Validation: must inherit from `Seurat`.
- Used by:
  - graph lookup (`object@graphs[[graph_name]]`),
  - final `cell_names` storage (`Cells(object)`).
- If invalid: immediate `stop()`.

### 3.2 `graph_name`

- Meaning: which Seurat graph slot to cluster (for example `"RNA_snn"`, `"harmony_snn"`).
- Default: `<DefaultAssay(object)>_snn` when `NULL`.
- Validation: must exist in `names(object@graphs)`.
- Consequences of change:
  - changes graph topology and weights,
  - can drastically change runtime and discovered cluster structure.

### 3.3 `cluster_range`

- Meaning: target cluster numbers `k` to evaluate.
- Used by:
  - resolution search (one range per `k`),
  - optimization (one worker task per valid `k`),
  - final consistency reporting.
- Consequences of change:
  - larger range -> more total work,
  - each additional `k` adds another full search + optimization branch.
- Practical note: this argument is not heavily validated (for positivity/ordering) in current code, so pass a clean integer vector.

### 3.4 `n_workers`

- Meaning: global parallel worker budget.
- Effective workers:
  - Unix-like: `min(requested, detectCores() - 1)`,
  - Windows: forced to `1` (no fork backend).
- Used by multiple nested stages (search/filter/optimize/bootstrap).
- Consequences of change:
  - usually faster up to a point,
  - too high can increase memory pressure and scheduling overhead,
  - nested stages may internally downscale worker counts.

### 3.5 `n_trials`

- Meaning: number of Leiden trials per gamma in optimization Phase 1.
- Used in `optimize_clustering()`.
- Consequences of change:
  - linear increase in compute cost,
  - higher robustness of median cluster count and IC estimates.

### 3.6 `n_bootstrap`

- Meaning: number of bootstrap resamples for final IC uncertainty.
- Used in optimization Phase 5.
- Consequences of change:
  - linear increase in bootstrap compute,
  - no extra Leiden runs (bootstrap resamples existing best clustering matrix).

### 3.7 `seed`

- Meaning: deterministic seed root for reproducibility.
- Behavior:
  - if `NULL`, results are stochastic.
  - if set, code derives stage-specific seeds for cluster search, gamma trials, and bootstrap.
- Consequence: improves repeatability, but does not eliminate all cross-platform runtime differences.

### 3.8 `beta`

- Meaning: Leiden randomness/noise parameter.
- Used in `leiden_clustering()` calls during optimization.
- Consequence:
  - affects exploration and partition variability,
  - can alter IC landscape and best gamma selection.

### 3.9 `n_iterations`

- Meaning: Leiden iteration count for optimization Phase 1.
- Does **not** control preliminary search iterations (those use lighter fixed settings).
- Consequence:
  - larger values can refine partitions but increase runtime.

### 3.10 `max_iterations`

- Meaning: upper limit for iterative refinement in optimization Phase 4.
- Consequence:
  - larger values allow more refinement when IC is not perfect,
  - can significantly increase runtime in hard cases.

### 3.11 `ic_threshold`

- Meaning: post-hoc threshold defining which tested `k` values are called “consistent”.
- Important: does not affect resolution search or optimization mechanics; it affects final reporting (`consistent_clusters`).
- Consequence:
  - stricter threshold -> fewer reported consistent cluster numbers,
  - looser threshold -> more reported candidates.

### 3.12 `objective_function`

- Meaning: Leiden objective (`"modularity"` or `"CPM"`).
- Current behavior:
  - exact string `"modularity"` uses modularity branch,
  - any other value goes to CPM branch in wrapper logic.
- Consequence:
  - changes gamma interpretation and optimization landscape.

### 3.13 `remove_threshold`

- Meaning: pre-optimization exclusion threshold for problematic cluster numbers.
- Behavior:
  - `Inf`: skip filtering entirely,
  - finite: run coarse IC screening over sampled gammas and exclude unstable `k`.
- Consequence:
  - lower threshold can exclude more `k` early (faster, but may remove viable targets),
  - `Inf` maximizes recall but may spend more compute downstream.

### 3.14 `min_cluster_size`

- Meaning: minimum cells for a cluster to count as effective.
- Default: `2`.
- Current semantics:
  - search/optimization target matching uses effective count only,
  - IC/MEI are computed on raw (unmerged) trial labels,
  - only final `best_labels` gets one deterministic merge pass.
- Edge rule: all-small trial => effective count `0`.
- Consequence:
  - increasing this value can make target `k` harder to hit,
  - can change “has output vs no output” even when raw clusters look reasonable.

### 3.15 `resolution_tolerance`

- Meaning: tolerance for binary search range convergence.
- Used in:
  - CPM lower search bound (`log(resolution_tolerance)`),
  - stopping condition for bound search loops,
  - narrow-range gamma sequence behavior.
- Consequence:
  - smaller tolerance can increase search iterations,
  - too large may produce coarse bounds.

### 3.16 `verbose`

- Meaning: enable detailed timestamped logging.
- Includes stage summaries and heartbeat messages.
- Consequence:
  - better observability,
  - minor overhead from log I/O.

## 4. Output Object (`class = "scICE"`)

Returned list fields (main ones):

- `gamma`: selected gamma per returned `k`.
- `labels`: raw extracted clustering arrays used for IC/MEI.
- `ic`: median bootstrap IC per `k`.
- `ic_vec`: bootstrap IC vectors.
- `n_cluster`: cluster numbers represented in returned entries.
- `best_labels`: final labels (with one final merge if `min_cluster_size > 1`).
- `n_iter`: final iteration count per `k`.
- `mei`: MEI scores.
- `k`: alias/metadata for iteration count.
- `excluded`, `exclusion_reason`: filtering metadata.
- `cell_names`: cells from input Seurat object.
- `cluster_range_tested`: original requested `cluster_range`.
- `consistent_clusters`: subset of `n_cluster` with `ic < ic_threshold`.
- `graph_name`, `min_cluster_size`: run metadata.

## 5. Full Workflow and Internal Call Graph

### 5.1 Entry Point: `scICE_clustering()`

Main responsibilities:

1. clear cache,
2. validate `object`, `graph_name`, `n_workers`, `min_cluster_size`,
3. compute effective worker count,
4. create runtime context (memory/spill settings),
5. extract graph from Seurat,
6. call `graph_to_igraph()`,
7. call `clustering_main()`,
8. post-process consistency summary and attach metadata.

### 5.2 Graph Conversion: `graph_to_igraph()`

Key behavior:

- accepts sparse graph (`dgCMatrix`/`Matrix`), converts if needed,
- reads sparse slots directly (`@i`, `@p`, `@x`) for memory efficiency,
- creates undirected weighted `igraph`.

Why this matters:

- large datasets spend non-trivial time here,
- conversion preserves the graph structure fed to Leiden.

### 5.3 Core Orchestration: `clustering_main()`

Stages:

1. clear cache + set run seed,
2. compute search bounds (`start_g`, `end_g`),
3. call `find_resolution_ranges()` for all target `k`,
4. optional filtering stage (skipped when `remove_threshold = Inf`),
5. call `optimize_clustering()` per valid `k` in parallel,
6. compute MEI for successful branches,
7. merge successful and excluded entries and return compatibility list.

### 5.4 Resolution Search: `find_resolution_ranges()`

Goal: find `[left_bound, right_bound]` gamma range for each target `k`.

#### 5.4.1 Preliminary settings

Based on graph size (`vcount`):

- `n_preliminary_trials`: 3 (>=200k), 5 (>=100k), else 15,
- `n_iter_preliminary`: 3 (>=200k), else 5,
- `max_search_iterations`: 30 (>=200k), else 50.

#### 5.4.2 Trial execution helper chain

- `run_preliminary_trials()`
  - runs up to `n_preliminary_trials` per gamma,
  - supports parallel trial launch and early stop,
  - emits heartbeat with running trials.
- Each trial uses `run_single_trial_count()`:
  - calls `cached_leiden_clustering()`,
  - computes:
    - `effective_count = count_effective_clusters(labels, min_cluster_size)`,
    - `raw_count = length(unique(labels))`.

#### 5.4.3 Lower-bound search logic

At each mid gamma:

- if `effective >= target`: move `right = mid`,
- else if `raw < target`: move `left = mid`,
- else (over-fragmented: raw high but effective low): move `right = mid`.

Interpretation:

- last branch tries to avoid drifting into high-gamma all-small regimes.

#### 5.4.4 Upper-bound search logic

At each mid gamma:

- if `effective >= target`: move `left = mid`,
- else if `raw < target`: move `left = mid`,
- else (over-fragmented): move `right = mid`.

#### 5.4.5 Post-processing

- produce cluster-specific ranges,
- adjust near-overlapping adjacent ranges,
- return named list `{k -> c(left, right)}`.

### 5.5 Optional Filtering Stage (inside `clustering_main()`)

- skipped entirely when `remove_threshold = Inf`.
- otherwise:
  - sample a few gamma points in each range,
  - compute IC scores,
  - exclude a target `k` if `min(IC)` is above threshold.

Result:

- excluded targets are retained as metadata entries (`excluded = TRUE`) but not optimized.

### 5.6 Intensive Optimization: `optimize_clustering()`

This is the most expensive part.

#### 5.6.1 Phase 1: evaluate gamma sequence

- build 11-step gamma sequence (5-step in narrow large-graph case),
- for each gamma:
  - run `n_trials` Leiden calls (`run_leiden_trial()` -> `leiden_clustering()`),
  - compute median effective cluster count across trials,
  - if median count != target `k`, mark invalid and skip IC,
  - if median count == target, compute IC and keep matrix reference.

Critical point:

- **exact equality** is required (`mean_clusters == target_clusters`).
- If no gamma passes this exact criterion, optimization returns `NULL` for that `k`.

#### 5.6.2 Phase 2: target-count filtering

- summarize how many gamma values produced each effective count,
- keep only exact matches to target.

#### 5.6.3 Phase 3: best gamma selection

- choose first perfect IC (`==1`) if present,
- else choose minimum IC.

#### 5.6.4 Phase 4: iterative refinement (optional)

Triggered when best IC is not perfect and multiple gammas remain:

- increase iteration budget by `delta_n = 2`,
- warm-start Leiden from previous memberships,
- keep stable/better gammas,
- stop on convergence criteria or `max_iterations`.

#### 5.6.5 Phase 5: bootstrap

- resample columns from best clustering matrix,
- recompute IC distribution (`ic_bootstrap`),
- use median as `ic_median`.

#### 5.6.6 Final label extraction and merge

- extract raw representative labels via `extract_clustering_array()` + `get_best_clustering()`,
- apply one final merge pass: `merge_small_clusters_to_neighbors()`.

Returned fields include:

- raw labels array (`labels`),
- final merged `best_labels`,
- IC statistics and iteration metadata.

## 6. Internal Helpers and Their Roles

### 6.1 `leiden_clustering()`

- wrapper over `igraph::cluster_leiden()`.
- passes edge weights when present.
- returns **0-based integer labels**.
- CPM/modularity branch selected by `objective_function`.

### 6.2 `cached_leiden_clustering()`

- process-local cache keyed by resolution/objective/iterations/beta/suffix.
- avoids repeated identical Leiden computations.

### 6.3 `count_effective_clusters()`

- counts clusters with size `>= min_cluster_size`.
- all-small scenario returns `0`.

### 6.4 `merge_small_clusters_to_neighbors()`

One-pass deterministic merge logic:

1. split clusters into large (`>= min_cluster_size`) and small.
2. If no large cluster exists:
   - collapse all cells into one cluster (largest size; tie -> smallest ID).
3. Else for each small cluster:
   - compute mean connectivity to each large cluster using SNN submatrix,
   - pick maximum connectivity target,
   - tie-break by smallest target ID.
4. relabel to contiguous 0-based IDs.

### 6.5 ECS/IC/MEI stack

- `calculate_ecs()` uses ClustAssess element-level similarity,
- `calculate_ic_from_extracted()` computes probability-weighted pairwise similarity aggregate,
- `calculate_mei_from_array()` computes per-cell stability from pairwise ECS.

### 6.6 Parallel and runtime utilities

- `cross_platform_mclapply()`: `mclapply` on Unix, `lapply` on Windows.
- heartbeat:
  - `get_heartbeat_interval_seconds()` reads `options("scICER.internal_heartbeat_seconds")`, default 60s,
  - `create_heartbeat_logger()` throttles periodic liveness logs.
- memory/spill:
  - `cap_workers_by_memory()`,
  - `create_runtime_context()`,
  - `activate_runtime_spill()` with `qs`,
  - `store/load/release_cluster_matrix()`.

## 7. Runtime and Scaling Characteristics

### 7.1 Major cost drivers

- number of cells/vertices,
- number of target `k` values,
- resolution search iterations,
- `n_trials` and gamma steps in Phase 1,
- `max_iterations` impact on Phase 4.

### 7.2 Approximate Leiden call budget

Per target `k` (worst-case rough upper bound):

- Resolution search:
  - up to `(lower + upper iterations) * n_preliminary_trials`
  - on large graphs: `60 * 3 = 180` preliminary Leiden runs.
- Optimization Phase 1:
  - `n_steps * n_trials` (typically `11 * n_trials`).
- Phase 4:
  - additional rounds, data-dependent.

Large `cluster_range` multiplies this cost.

## 8. Why "No Output" Can Happen

The pipeline can return no valid cluster solutions even without exceptions.

Common patterns:

1. No gamma exactly matches target effective count in Phase 2.
2. All targets are excluded by filtering (`remove_threshold` finite).
3. Effective-count semantics with higher `min_cluster_size` makes target hard to hit.
4. Strong stochastic/non-monotonic behavior near gamma boundaries.

Symptoms in logs:

- many lines like `median effective clusters = 0` with very high raw clusters,
- final `ERROR - No gammas produced target effective cluster count X`,
- final result vectors length `0`.

## 9. Relationship to Seurat Resolution/Gamma

In Leiden context (`algorithm = 4` in Seurat language):

- scICER `gamma` and Seurat `FindClusters(..., resolution=...)` represent the same conceptual resolution control.
- However, numeric equivalence is not guaranteed across pipelines because of:
  - graph construction differences,
  - objective function settings,
  - iteration/randomness choices,
  - post-processing semantics (especially effective counting and singleton handling).

## 10. Parameter Tuning Guidance (Practical)

For large datasets:

1. Start with small `cluster_range` (for example 1-3 nearby `k` values).
2. Keep `remove_threshold = Inf` while debugging search behavior.
3. Use moderate `n_trials` first, then increase if needed.
4. Keep `verbose = TRUE` and heartbeat enabled for long runs.
5. Treat `min_cluster_size` as a strong semantic switch, not a cosmetic post-processing option.

## 11. Internal Debug Options (Advanced)

The code reads several internal R options:

- `scICER.internal_heartbeat_seconds`
- `scICER.internal_memory_budget_bytes`
- `scICER.internal_worker_memory_overhead`
- `scICER.internal_force_spill`
- `scICER.internal_spill_threshold_bytes`
- `scICER.internal_spill_dir`

These are not part of the public stable API, but are useful for profiling and large-run operations.

## 12. End-to-End Call Chain (Reference)

```text
scICE_clustering
  -> graph_to_igraph
  -> clustering_main
       -> find_resolution_ranges
            -> run_preliminary_trials
                 -> run_single_trial_count
                      -> cached_leiden_clustering
                           -> leiden_clustering
                      -> count_effective_clusters
       -> (optional) filtering stage
       -> optimize_clustering (per k)
            -> evaluate_gamma
                 -> run_leiden_trial
                      -> leiden_clustering
                 -> count_effective_clusters
                 -> extract_clustering_array
                 -> calculate_ic_from_extracted
                      -> calculate_ecs
            -> iterative refinement (optional)
            -> bootstrap IC
            -> get_best_clustering
            -> merge_small_clusters_to_neighbors
       -> calculate_mei_from_array
            -> calculate_ecs
  -> consistency post-processing (ic_threshold)
```

