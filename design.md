# scICE_clustering Design Document

## 1. Scope

This document explains the current implementation of `scICE_clustering()` in **scICER**.

It covers:

- every public parameter (meaning, where it is used, and what happens if you change it),
- the full runtime workflow from Seurat graph input to final outputs,
- the internal function call chain and behavior of each major helper,
- runtime, parallelism, memory, and failure modes.

This description matches the code on `main` after the raw-cluster-aware
resolution search and gamma-admission updates on 2026-03-17, plus the
`plot_ic()` gamma-label placement update on 2026-03-17.

## 1.1 Visualization Note

- `plot_ic(..., show_gamma = TRUE)` now renders gamma values directly under each
  x-axis cluster label as a two-line tick label and rotates the x-axis text by
  45 degrees.
- The plot subtitle remains focused on threshold and exclusion context instead
  of listing all selected gamma values inline.

## 2. Quick Mental Model

`scICE_clustering()` is a multi-stage pipeline:

1. Validate inputs and initialize runtime context.
2. Extract a Seurat graph and convert it to `igraph`.
3. For each target cluster number `k`, run a raw-cluster-aware binary search to
   find a coarse gamma range, then optionally clamp that range to a raw-cluster
   plateau/bracket near the target.
4. Optionally filter unstable `k` values.
5. For each remaining `k`, run intensive optimization:
   - evaluate many gamma values,
   - collect both effective-cluster and raw-cluster medians across trials,
   - compute IC using all trials for any gamma admitted by the effective/raw
     candidate rules,
   - choose one candidate family via the ordered admission ladder:
     `raw_strict_soft -> strict_soft -> relaxed_soft -> strict_hard -> relaxed_hard -> relaxed_unguarded -> raw_relaxed_soft -> raw_relaxed_hard -> raw_relaxed_unguarded`,
   - refine multi-gamma candidate sets by minimum raw-cluster median gap to `k`,
   - compute IC/MEI stability metrics,
   - select `best_labels`.
6. Return an `scICE` object and annotate which `k` values are consistent under `ic_threshold`.

Important: current implementation uses **effective cluster counting** when `min_cluster_size > 1`.

- Effective cluster count = number of clusters with size `>= min_cluster_size`.
- Raw cluster count = number of unique cluster labels before any final merge.
- If all clusters are smaller than threshold, effective count is `0`.
- Resolution search now reuses **one representative preliminary clustering per
  gamma step** and summarizes nominal median effective/raw counts from that
  single run.
- Raw cluster medians are used as guardrails during resolution search and as a
  bounded fallback family during gamma admission.
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
  - optimization now splits this budget into:
    - active outer cluster workers (`k`-level queue),
    - per-cluster nested workers for gamma/bootstrap,
  - outer optimization uses dynamic queue scheduling (`mc.preschedule = FALSE`) to reduce long-tail imbalance.

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
  - effective count still drives the main target-matching objective,
  - resolution search and optimization also inspect raw cluster medians to
    reject pathological over-fragmented matches,
  - resolution search additionally clamps coarse gamma bounds to a raw-cluster
    plateau/bracket/near-target region when `min_cluster_size > 1`,
  - optimization can fall back to bounded raw-count families when effective
    admission families are empty,
  - IC/MEI are computed on raw (unmerged) trial labels,
  - only final `best_labels` gets one deterministic merge pass.
- Edge rule: all-small trial => effective count `0`.
- Consequence:
  - increasing this value can make target `k` harder to hit,
  - can change “has output vs no output” even when raw clusters look
    reasonable, because raw-count guards/fallbacks are also activated.

### 3.15 `resolution_tolerance`

- Meaning: tolerance for binary search range convergence.
- Used in:
  - CPM lower search bound (`log(resolution_tolerance)`),
  - stopping condition for bound search loops,
  - narrow-range gamma sequence behavior,
  - plateau-clamp probe grid generation.
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
- `effective_cluster_median`: median effective cluster count at the selected gamma.
- `raw_cluster_median`: median raw cluster count at the selected gamma.
- `admission_mode`: winning Phase-2 admission family for the selected gamma.
- `best_labels_raw_cluster_count`: raw cluster count of the selected best trial
  before the final merge.
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
   - worker layout is computed as outer `k` workers + per-`k` nested budget,
   - outer `k` tasks run with dynamic queue scheduling to reduce stragglers,
6. compute MEI for successful branches,
7. merge successful and excluded entries and return compatibility list,
8. emit per-`k` selection diagnostics:
   - selected gamma,
   - selected effective/raw medians,
   - winning admission family,
   - raw cluster count before final merge.

### 5.4 Resolution Search: `find_resolution_ranges()`

Goal: find `[left_bound, right_bound]` gamma range for each target `k` while
avoiding high-gamma regimes where raw clusters explode but effective clusters
collapse after `min_cluster_size` filtering.

#### 5.4.1 Preliminary settings

Based on graph size (`vcount`):

- `n_preliminary_trials`: 3 (>=200k), 5 (>=100k), else 15,
- `n_iter_preliminary`: 3 (>=200k), else 5,
- `max_search_iterations`: 30 (>=200k), else 50.

#### 5.4.2 Trial execution helper chain

- `run_preliminary_trials()`
  - now runs **one representative** preliminary Leiden call per gamma step via
    `run_single_trial_count(..., use_cache = TRUE)`,
  - then fills nominal vectors of length `n_preliminary_trials` with that single
    result so downstream median-based logic stays unchanged,
  - still computes/logs worker mode and capacity, but no longer launches a full
    preliminary burst per gamma step.
- Each trial uses `run_single_trial_count()`:
  - calls `cached_leiden_clustering()`,
  - computes:
    - `effective_count = count_effective_clusters(labels, min_cluster_size)`,
    - `raw_count = length(unique(labels))`.

Important implication:

- Actual preliminary Leiden call count is now approximately the number of gamma
  steps evaluated, not `gamma steps * n_preliminary_trials`.
- The logged “trials per step” remains the nominal setting (`3`, `5`, or `15`).

#### 5.4.3 Search-state classification

At each gamma probe, `classify_resolution_search_state()` evaluates:

- median raw cluster count,
- median effective cluster count,
- whether raw count is below target, in an acceptable search band, or above it,
- whether the step is over-fragmented (`raw` not below target but `effective < target`).

The search-band ceiling is defined by `raw_cluster_search_upper(target)`:

- `target + 1` when `target <= 10`,
- otherwise `max(target + 2, ceiling(target * 1.05))`.

This is stricter than the later optimization guards and is used only to steer
the binary search away from obviously over-fragmented regions.

#### 5.4.4 Lower-bound and upper-bound search logic

Lower-bound search action:

- if raw median is below target, increase gamma,
- otherwise decrease gamma.

Upper-bound search action:

- if raw median is below target, increase gamma,
- if raw median is in-band and effective median meets target, increase gamma,
- otherwise decrease gamma.

Interpretation:

- raw-below means gamma is too small,
- raw-in-band plus effective-at-target supports widening upward,
- over-fragmented or raw-above-band states push the search downward so the
  range does not drift into all-small-cluster regimes.

#### 5.4.5 Plateau clamp and post-processing

After coarse lower/upper bounds are found:

- if `min_cluster_size > 1`, build a probe gamma grid across the coarse range
  using the same helper as optimization (`build_gamma_sequence_for_range()`),
- evaluate representative raw medians on that grid,
- stabilize the raw-median curve with `cummax()` to suppress non-monotone noise,
- clamp the final range by priority:
  - exact raw plateau (`raw_exact`),
  - nearest raw crossing/bracket (`raw_bracket`),
  - near-target raw region with `abs(raw_median - target) <= 1` (`raw_near_target`),
  - otherwise keep the coarse range (`coarse`).

Additional details:

- CPM degenerate bounds are widened on the **gamma scale** before being mapped
  back to log space.
- Adjacent cluster ranges are still nudged apart when they nearly overlap.
- Final return value is a named list `{k -> c(left, right)}` on the actual gamma scale.

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

- build gamma sequence with the shared helper `build_gamma_sequence_for_range()`
  (typically 11 steps, or 5 for narrow large-graph cases),
- for each gamma:
  - run `n_trials` Leiden calls (`run_leiden_trial()` -> `leiden_clustering()`),
  - compute:
    - `median_effective_clusters`,
    - `median_clusters_int = as.integer(median_effective_clusters)`,
    - `raw_cluster_median`,
    - `raw_cluster_median_int = as.integer(raw_cluster_median)`,
    - effective hit count (`effective_count == k`),
    - raw hit count (`raw_count == k`),
    - effective/raw median gaps to target,
    - soft/hard raw-cluster guard flags.

Admission flags:

- effective strict: `as.integer(median_effective_clusters) == k`,
- effective relaxed: `effective_hit_count >= 1 && abs(median_effective_clusters - k) <= 1`,
- raw strict: `as.integer(raw_cluster_median) == k`,
- raw relaxed: `raw_hit_count >= 1 && abs(raw_cluster_median - k) <= 1`.

Phase-1 IC rule:

- if **none** of the four admission flags pass, mark the gamma invalid and skip IC,
- otherwise compute IC from **all trials** at that gamma and keep the full trial matrix.

#### 5.6.2 Phase 2: target-count filtering

- summarize effective-count and raw-count diagnostics across the gamma grid,
- choose the first non-empty candidate family from this ordered ladder:
  1. `raw_strict_soft`
  2. `strict_soft`
  3. `relaxed_soft`
  4. `strict_hard`
  5. `relaxed_hard`
  6. `relaxed_unguarded`
  7. `raw_relaxed_soft`
  8. `raw_relaxed_hard`
  9. `raw_relaxed_unguarded`

Where:

- `soft` and `hard` refer to raw-cluster ceilings from `raw_cluster_guard_limits(target)`:
  - soft = `max(target + 3, ceiling(target * 1.1))`
  - hard = `max(target + 5, ceiling(target * 1.5))`
- `raw_strict_*` intentionally outranks effective families when there is an
  exact bounded raw-count plateau.

Tie refinement:

- if the winning family contains multiple gammas and `min_cluster_size > 1`,
  keep only the gammas with the smallest raw-median gap to the target.

If all families are empty, optimization returns `NULL`.

#### 5.6.3 Phase 3: best gamma selection

- choose first perfect IC (`==1`) if present,
- else choose minimum IC,
- retain diagnostics for the selected gamma:
  - `effective_cluster_median`,
  - `raw_cluster_median`,
  - `admission_mode`.

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
- record `best_labels_raw_cluster_count`,
- apply one final merge pass: `merge_small_clusters_to_neighbors()`,
- log final selected diagnostics including:
  - selected gamma,
  - selected effective/raw medians,
  - winning admission family,
  - raw and final cluster counts for `best_labels`.

Returned fields include:

- raw labels array (`labels`),
- final merged `best_labels`,
- selected effective/raw medians,
- admission family metadata,
- raw cluster count before final merge,
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

### 6.4 Raw-count guard and search helpers

- `raw_cluster_guard_limits()`:
  - defines optimization-stage raw-cluster ceilings:
    - soft = `max(target + 3, ceiling(target * 1.1))`,
    - hard = `max(target + 5, ceiling(target * 1.5))`.
- `raw_cluster_search_upper()`:
  - defines the stricter resolution-search raw-count ceiling.
- `build_gamma_sequence_for_range()`:
  - centralizes gamma-grid generation for both optimization and plateau clamp.
- `classify_resolution_search_state()`:
  - converts raw/effective medians into lower/upper search actions.
- `stabilize_probe_raw_medians()`:
  - applies `cummax()` to raw-median probe curves before plateau clamp.
- `clamp_gamma_range_to_raw_plateau()`:
  - shrinks coarse ranges to exact raw plateaus, brackets, or near-target raw regions.
- `select_gamma_admission()`:
  - implements the ordered candidate-family ladder used in optimization Phase 2.
- `refine_gamma_candidates_by_raw_gap()`:
  - breaks ties by keeping gammas with the minimum raw-median gap to target.

### 6.5 `merge_small_clusters_to_neighbors()`

One-pass deterministic merge logic:

1. split clusters into large (`>= min_cluster_size`) and small.
2. If no large cluster exists:
   - collapse all cells into one cluster (largest size; tie -> smallest ID).
3. Else for each small cluster:
   - compute mean connectivity to each large cluster using SNN submatrix,
   - pick maximum connectivity target,
   - tie-break by smallest target ID.
4. relabel to contiguous 0-based IDs.

### 6.6 ECS/IC/MEI stack

- `calculate_ecs()` uses ClustAssess element-level similarity,
- `calculate_ic_from_extracted()` computes probability-weighted pairwise similarity aggregate,
- `calculate_mei_from_array()` computes per-cell stability from pairwise ECS.

### 6.7 Parallel and runtime utilities

- `cross_platform_mclapply()`:
  - `mclapply` on Unix and `lapply` on Windows,
  - supports `mc.preschedule` on Unix callers.
- heartbeat:
  - `get_heartbeat_interval_seconds()` reads `options("scICER.internal_heartbeat_seconds")`, default 60s,
  - `create_heartbeat_logger()` throttles periodic liveness logs.
- memory/spill:
  - `cap_workers_by_memory()`,
  - `create_runtime_context()`,
  - `activate_runtime_spill()` with `qs`,
  - `store/load/release_cluster_matrix()`.

Optimization-stage worker allocation details:

- starts from memory-capped optimization worker budget,
- computes active outer `k` workers,
- passes the full optimization budget into each outer worker,
- each optimizer then derives its nested gamma/bootstrap workers from that
  budget and its current gamma-grid size,
- runs outer `k` loop with dynamic queue (`mc.preschedule = FALSE`).

Resolution-search preliminary trial allocation details:

- computes per-target worker capacity from global search budget and active outer `k` workers,
- still logs serial/parallel preliminary mode based on worker capacity,
- but actual preliminary evaluation now reuses one representative cached Leiden
  result per gamma step rather than launching a full per-step trial burst.

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
  - coarse lower + upper binary search now costs about one representative Leiden
    call per search step, not `search steps * n_preliminary_trials`,
  - large graphs: up to about `30 + 30 = 60` representative preliminary runs,
  - smaller graphs: up to about `50 + 50 = 100` representative preliminary runs,
  - if `min_cluster_size > 1`, plateau clamp adds another probe pass:
    typically 11 representative runs (or 5 for narrow large-graph cases).
- Optimization Phase 1:
  - `n_steps * n_trials` (typically `11 * n_trials`).
- Phase 4:
  - additional rounds, data-dependent.

Large `cluster_range` multiplies this cost.

## 8. Why "No Output" Can Happen

The pipeline can return no valid cluster solutions even without exceptions.

Common patterns:

1. No gamma survives the Phase-2 admission ladder after raw guards and tie refinement.
2. Raw cluster counts reach over-fragmented regimes, so the search/optimizer avoids them even if a few effective hits appear.
3. All targets are excluded by filtering (`remove_threshold` finite).
4. Effective-count semantics with higher `min_cluster_size` makes target hard to hit.
5. Strong stochastic/non-monotonic behavior near gamma boundaries.

Symptoms in logs:

- many lines with `over_fragmented = TRUE` or very large `median_raw`,
- plateau clamp modes such as `coarse` when no clean raw plateau/bracket is found,
- final errors about no gamma satisfying raw-guarded strict/relaxed admission or bounded raw-count fallback admission for effective cluster count `X`,
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
- `scICER.internal_preliminary_parallel_min_workers`
- `scICER.internal_preliminary_parallel_max_workers`
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
            -> classify_resolution_search_state
            -> build_gamma_sequence_for_range
            -> clamp_gamma_range_to_raw_plateau
       -> (optional) filtering stage
       -> optimize_clustering (per k)
            -> evaluate_gamma
                 -> run_leiden_trial
                      -> leiden_clustering
                 -> count_effective_clusters
                 -> passes_raw_cluster_guard
                 -> extract_clustering_array
                 -> calculate_ic_from_extracted
                      -> calculate_ecs
            -> select_gamma_admission
            -> refine_gamma_candidates_by_raw_gap
            -> iterative refinement (optional)
            -> bootstrap IC
            -> get_best_clustering
            -> merge_small_clusters_to_neighbors
       -> calculate_mei_from_array
            -> calculate_ecs
  -> consistency post-processing (ic_threshold)
```
