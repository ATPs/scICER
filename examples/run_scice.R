#!/usr/bin/env Rscript

script_name <- "Rscript examples/run_scice.R"
default_input_qs <- "/data1/xlab/researches/20250709_scICE/20260209_reviews/20260305_pancreatic/pancreatic_harmony1.qs"

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush(stdout())
}

rss_gib <- function() {
  status_file <- "/proc/self/status"
  if (!file.exists(status_file)) {
    return(NA_real_)
  }

  status_lines <- readLines(status_file, warn = FALSE)
  rss_line <- status_lines[grepl("^VmRSS:\\s+[0-9]+\\s+kB$", status_lines)]
  if (length(rss_line) == 0L) {
    return(NA_real_)
  }

  as.numeric(sub("^VmRSS:\\s+([0-9]+)\\s+kB$", "\\1", rss_line[[1]])) / 1024^2
}

format_gib <- function(bytes) {
  sprintf("%.2f GiB", as.numeric(bytes) / 1024^3)
}

log_memory <- function(label, object = NULL) {
  rss_text <- if (is.na(rss_gib())) "NA" else sprintf("%.2f GiB", rss_gib())
  size_text <- if (is.null(object)) {
    NULL
  } else {
    format_gib(object.size(object))
  }

  parts <- c(paste("RSS =", rss_text))
  if (!is.null(size_text)) {
    parts <- c(parts, paste("object.size =", size_text))
  }

  log_step(sprintf("%s: %s", label, paste(parts, collapse = " | ")))
}

read_env_string <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  value
}

parse_int_vector <- function(value, default = NULL, option_label = "integer vector") {
  if (is.null(value) || !nzchar(value)) {
    return(default)
  }

  value <- gsub("\\s+", "", value)
  if (!nzchar(value)) {
    return(default)
  }

  parsed <- if (grepl("^[0-9]+:[0-9]+$", value)) {
    bounds <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    seq.int(bounds[[1]], bounds[[2]])
  } else {
    as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  }

  if (length(parsed) == 0L || anyNA(parsed)) {
    stop("Invalid value for ", option_label, ": ", value)
  }

  parsed
}

parse_numeric_vector <- function(value, default = NULL, option_label = "numeric vector") {
  if (is.null(value) || !nzchar(value)) {
    return(default)
  }

  tokens <- strsplit(gsub("\\s+", "", value), ",", fixed = TRUE)[[1]]
  parsed <- as.numeric(tokens)
  if (length(parsed) == 0L || anyNA(parsed) || any(!is.finite(parsed))) {
    stop("Invalid value for ", option_label, ": ", value)
  }

  parsed
}

parse_numeric_scalar <- function(value, option_label) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || is.na(parsed)) {
    stop("Invalid numeric value for ", option_label, ": ", value)
  }
  parsed
}

parse_integer_scalar <- function(value, option_label, min_value = NULL) {
  parsed <- parse_numeric_scalar(value, option_label)
  if (!is.finite(parsed) || abs(parsed - round(parsed)) > .Machine$double.eps^0.5) {
    stop(option_label, " must be an integer value. Received: ", value)
  }

  parsed <- as.integer(round(parsed))
  if (!is.null(min_value) && parsed < min_value) {
    stop(option_label, " must be >= ", min_value, ". Received: ", parsed)
  }

  parsed
}

is_bool_literal <- function(value) {
  tolower(value) %in% c("1", "0", "true", "false", "yes", "no", "y", "n", "on", "off")
}

parse_bool_string <- function(value, option_label) {
  value <- tolower(as.character(value)[[1]])
  if (value %in% c("1", "true", "yes", "y", "on")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "no", "n", "off")) {
    return(FALSE)
  }
  stop(
    "Invalid boolean value for ", option_label, ": ", value,
    ". Accepted values: true/false, 1/0, yes/no, on/off."
  )
}

build_light_scice_object <- function(seurat_obj, graph_name) {
  assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  keep_feature <- rownames(seurat_obj[[assay_name]])[[1]]

  if (is.null(keep_feature) || !nzchar(keep_feature)) {
    stop("Failed to determine one feature to retain in the lightweight Seurat object.")
  }

  Seurat::DietSeurat(
    object = seurat_obj,
    assays = assay_name,
    graphs = graph_name,
    dimreducs = NULL,
    features = keep_feature,
    misc = FALSE
  )
}

print_help <- function() {
  help_lines <- c(
    paste("Usage:", script_name, "--input_qs PATH [--output_qs PATH]"),
    "",
    "Purpose:",
    "  Run scICER on a very large Seurat object stored in qs format while reducing",
    "  clustering-phase memory usage. The script reads the full object, keeps only",
    "  one feature plus the target graph in a lightweight Seurat object, removes",
    "  the full object before scICE_clustering(), then reloads the original object",
    "  only to merge metadata back in, and finally saves the augmented Seurat object",
    "  back to qs format.",
    "",
    "Precedence:",
    "  Named CLI option > environment variable > built-in default.",
    "  This script no longer uses positional arguments.",
    "",
    "Options:",
    "  --help, -h",
      "      Print this help text and exit.",
    "",
    "  --input_qs PATH",
      "      Path to the input qs file containing the original Seurat object.",
      paste("      Env: SCICE_INPUT_QS | Default:", default_input_qs),
    "",
    "  --output_qs PATH",
    "      Path for qs::qsave(full_object_with_metadata, ...).",
    "      If omitted, it defaults to --input_qs and overwrites the input file.",
    "      Env: SCICE_OUTPUT_QS | Default: same as --input_qs",
    "",
    "  --graph_name NAME",
    "      Graph slot to cluster, for example RNA_snn or harmony_snn.",
    "      Env: SCICE_GRAPH_NAME | Default: RNA_snn",
    "",
    "  --qread_threads INT",
    "      Number of threads passed to qs::qread().",
    "      Env: SCICE_QREAD_THREADS | Default: 40",
    "",
    "  --n_workers INT",
    "      scICER worker budget used inside scICE_clustering().",
    "      Env: SCICE_N_WORKERS | Default: 40",
    "",
    "  --n_trials INT",
    "      Number of repeated Leiden trials per resolution/gamma.",
    "      Env: SCICE_N_TRIALS | Default: 4",
    "",
    "  --n_bootstrap INT",
    "      Number of bootstrap iterations used for IC estimation.",
    "      Env: SCICE_N_BOOTSTRAP | Default: 20",
    "",
    "  --seed INT",
    "      Random seed passed to scICE_clustering().",
    "      Env: SCICE_SEED | Default: 123",
    "",
    "  --plot_threshold NUM",
    "      Threshold passed to plot_ic(). Ignored when plotting is skipped.",
    "      Env: SCICE_PLOT_THRESHOLD | Default: 1.005",
    "",
    "  --label_threshold NUM",
    "      Threshold passed to get_robust_labels() when merging metadata back.",
    "      Use Inf to merge all returned cluster solutions.",
    "      Env: SCICE_LABEL_THRESHOLD | Default: 1.01",
    "",
    "  --cluster_range SPEC",
    "      Cluster-number targets for cluster-range mode.",
    "      Accepted forms: 4:5 or 4,5 or 4,5,6.",
    "      Ignored when --resolution is provided.",
    "      Env: SCICE_CLUSTER_RANGE | Default: 4:5",
    "",
    "  --resolution LIST",
    "      Comma-separated manual gamma values, for example 0.01,0.02.",
    "      When set, the script uses manual-resolution mode and skips cluster_range search.",
    "      Env: SCICE_RESOLUTION | Default: empty",
    "",
    "  --remove_threshold NUM",
    "      remove_threshold passed to scICE_clustering().",
    "      Env: SCICE_REMOVE_THRESHOLD | Default: Inf",
    "",
    "  --skip_plot[=BOOL]",
    "      Skip plot_ic() rendering. If specified without a value, it is treated as true.",
    "      Accepted BOOL values: true/false, 1/0, yes/no, on/off.",
    "      Env: SCICE_SKIP_PLOT | Default: FALSE",
    "",
    "What the script does:",
    "  1. qs::qread() the original Seurat object.",
    "  2. Build a lightweight Seurat object with DietSeurat() that keeps one feature",
    "     plus the requested graph.",
    "  3. rm(full_object); gc() before clustering to reduce RSS.",
    "  4. Run scICE_clustering(..., verbose = TRUE).",
    "  5. Optionally render plot_ic().",
    "  6. Reload the original Seurat object with qs::qread().",
    "  7. Add clustering metadata back with get_robust_labels(..., return_seurat = TRUE).",
    "  8. Save the augmented Seurat object with qs::qsave().",
    "",
    "Examples:",
    paste("  ", script_name, "--help"),
    paste(
      "  ", script_name,
      "--input_qs input.qs --output_qs augmented_with_scice.qs --cluster_range 4:5",
      "--qread_threads 40 --n_workers 40 --n_trials 4 --n_bootstrap 20 --seed 123"
    ),
    paste(
      "  SCICE_RESOLUTION=0.01,0.02 SCICE_SKIP_PLOT=1",
      script_name,
      "--input_qs input.qs"
    )
  )

  cat(paste(help_lines, collapse = "\n"), "\n")
  flush(stdout())
}

parse_cli_args <- function(args) {
  known_options <- c(
    "input_qs",
    "output_qs",
    "graph_name",
    "qread_threads",
    "n_workers",
    "n_trials",
    "n_bootstrap",
    "seed",
    "plot_threshold",
    "label_threshold",
    "cluster_range",
    "resolution",
    "remove_threshold",
    "skip_plot"
  )

  named <- list()
  positional <- character(0)
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("--help", "-h")) {
      return(list(help = TRUE, named = named, positional = positional))
    }

    if (!startsWith(arg, "--")) {
      positional <- c(positional, arg)
      i <- i + 1L
      next
    }

    raw_option <- substring(arg, 3L)
    if (!nzchar(raw_option)) {
      stop("Encountered an empty option name. Run --help for usage.")
    }

    option_name <- raw_option
    option_value <- NULL
    if (grepl("=", raw_option, fixed = TRUE)) {
      parts <- strsplit(raw_option, "=", fixed = TRUE)[[1]]
      option_name <- parts[[1]]
      option_value <- paste(parts[-1], collapse = "=")
    }

    if (!(option_name %in% known_options)) {
      stop("Unknown option: --", option_name, ". Run --help for usage.")
    }

    if (identical(option_name, "skip_plot")) {
      if (is.null(option_value)) {
        next_arg_is_bool <- i < length(args) && is_bool_literal(args[[i + 1L]])
        if (next_arg_is_bool) {
          option_value <- args[[i + 1L]]
          i <- i + 1L
        } else {
          option_value <- "true"
        }
      }
    } else if (is.null(option_value)) {
      if (i >= length(args) || startsWith(args[[i + 1L]], "--")) {
        stop("Option --", option_name, " requires a value. Run --help for usage.")
      }
      option_value <- args[[i + 1L]]
      i <- i + 1L
    }

    named[[option_name]] <- option_value
    i <- i + 1L
  }

  if (length(positional) > 0L) {
    stop("Positional arguments are no longer supported. Use --input_qs and --output_qs. Run --help for usage.")
  }

  list(help = FALSE, named = named, positional = positional)
}

resolve_path_option <- function(named, option_name, env_name, default) {
  if (!is.null(named[[option_name]])) {
    return(named[[option_name]])
  }
  read_env_string(env_name, default)
}

resolve_scalar_option <- function(named, option_name, env_name, default) {
  if (!is.null(named[[option_name]])) {
    return(named[[option_name]])
  }
  read_env_string(env_name, as.character(default))
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(cli$help)) {
  print_help()
  quit(save = "no", status = 0)
}

input_qs <- resolve_path_option(
  named = cli$named,
  option_name = "input_qs",
  env_name = "SCICE_INPUT_QS",
  default = default_input_qs
)
output_qs <- resolve_path_option(
  named = cli$named,
  option_name = "output_qs",
  env_name = "SCICE_OUTPUT_QS",
  default = ""
)

graph_name <- resolve_scalar_option(cli$named, "graph_name", "SCICE_GRAPH_NAME", "RNA_snn")
qread_threads <- parse_integer_scalar(
  resolve_scalar_option(cli$named, "qread_threads", "SCICE_QREAD_THREADS", 40),
  option_label = "--qread_threads",
  min_value = 1L
)
n_workers <- parse_integer_scalar(
  resolve_scalar_option(cli$named, "n_workers", "SCICE_N_WORKERS", 40),
  option_label = "--n_workers",
  min_value = 1L
)
n_trials <- parse_integer_scalar(
  resolve_scalar_option(cli$named, "n_trials", "SCICE_N_TRIALS", 4),
  option_label = "--n_trials",
  min_value = 1L
)
n_bootstrap <- parse_integer_scalar(
  resolve_scalar_option(cli$named, "n_bootstrap", "SCICE_N_BOOTSTRAP", 20),
  option_label = "--n_bootstrap",
  min_value = 1L
)
seed <- parse_integer_scalar(
  resolve_scalar_option(cli$named, "seed", "SCICE_SEED", 123),
  option_label = "--seed"
)
plot_threshold <- parse_numeric_scalar(
  resolve_scalar_option(cli$named, "plot_threshold", "SCICE_PLOT_THRESHOLD", 1.005),
  option_label = "--plot_threshold"
)
label_threshold <- parse_numeric_scalar(
  resolve_scalar_option(cli$named, "label_threshold", "SCICE_LABEL_THRESHOLD", 1.01),
  option_label = "--label_threshold"
)
cluster_range <- parse_int_vector(
  resolve_scalar_option(cli$named, "cluster_range", "SCICE_CLUSTER_RANGE", "4:5"),
  default = 4:5,
  option_label = "--cluster_range"
)
resolution_values <- parse_numeric_vector(
  resolve_scalar_option(cli$named, "resolution", "SCICE_RESOLUTION", ""),
  default = NULL,
  option_label = "--resolution"
)
remove_threshold <- parse_numeric_scalar(
  resolve_scalar_option(cli$named, "remove_threshold", "SCICE_REMOVE_THRESHOLD", Inf),
  option_label = "--remove_threshold"
)
skip_plot <- parse_bool_string(
  resolve_scalar_option(cli$named, "skip_plot", "SCICE_SKIP_PLOT", FALSE),
  option_label = "--skip_plot"
)

if (!nzchar(input_qs)) {
  stop("Input qs path is empty. Use --input_qs or SCICE_INPUT_QS.")
}
if (!nzchar(output_qs)) {
  output_qs <- input_qs
}

log_step("Loading required packages...")
library(Seurat)
library(scICER)
library(qs)
log_step("Package loading complete.")
log_step(sprintf("Input qs: %s", input_qs))
log_step(sprintf("Output qs: %s", output_qs))
if (identical(normalizePath(output_qs, winslash = "/", mustWork = FALSE), normalizePath(input_qs, winslash = "/", mustWork = FALSE))) {
  log_step("Output qs matches input_qs; the input file will be overwritten after metadata merge.")
}
log_step(sprintf("Graph name: %s", graph_name))
log_step(sprintf("qread threads: %d", qread_threads))
log_step(sprintf("n_workers: %d | n_trials: %d | n_bootstrap: %d | seed: %d", n_workers, n_trials, n_bootstrap, seed))
if (!is.null(resolution_values)) {
  log_step(sprintf("Manual resolutions: %s", paste(signif(resolution_values, 6), collapse = ", ")))
} else {
  log_step(sprintf("Cluster range: %s", paste(cluster_range, collapse = ", ")))
}

log_step(sprintf("Starting qs::qread with %d thread(s)", qread_threads))
full_object <- qs::qread(input_qs, nthread = qread_threads)
log_memory("After qs::qread", full_object)

log_step("Building lightweight Seurat object for scICER...")
light_object <- build_light_scice_object(full_object, graph_name = graph_name)
log_memory("After building lightweight object", light_object)

log_step("Dropping the full Seurat object before scICE_clustering()")
rm(full_object)
invisible(gc())
log_memory("After rm(full_object) + gc()", light_object)

scice_args <- list(
  object = light_object,
  graph_name = graph_name,
  remove_threshold = remove_threshold,
  n_workers = n_workers,
  n_trials = n_trials,
  n_bootstrap = n_bootstrap,
  seed = seed,
  verbose = TRUE
)

if (!is.null(resolution_values)) {
  scice_args$resolution <- resolution_values
  log_step(sprintf(
    "Running manual-resolution mode with gamma value(s): %s",
    paste(signif(resolution_values, 6), collapse = ", ")
  ))
} else {
  scice_args$cluster_range <- cluster_range
  log_step(sprintf(
    "Running cluster_range mode with cluster_range = %s",
    paste(cluster_range, collapse = ", ")
  ))
}

scice_results <- do.call(scICE_clustering, scice_args)
log_memory("After scICE_clustering()", scice_results)

if (!skip_plot) {
  log_step(sprintf("Rendering plot_ic(threshold = %.4f)", plot_threshold))
  print(plot_ic(scice_results, threshold = plot_threshold))
} else {
  log_step("Skipping plot_ic() because --skip_plot is enabled.")
}

log_step("Dropping the lightweight Seurat object before reloading the original object...")
rm(light_object)
invisible(gc())
log_memory("After rm(light_object) + gc()")

log_step("Reloading the original Seurat object for metadata merge...")
full_object <- qs::qread(input_qs, nthread = qread_threads)
log_memory("After reloading the original object", full_object)

full_object <- get_robust_labels(
  scice_results,
  return_seurat = TRUE,
  object = full_object,
  threshold = label_threshold
)

if (is.null(full_object)) {
  stop("No metadata columns were added back to the original Seurat object.")
}

cluster_columns <- grep("^clusters_", colnames(full_object@meta.data), value = TRUE)
log_step(sprintf(
  "Metadata merge complete. Added %d column(s): %s",
  length(cluster_columns),
  paste(cluster_columns, collapse = ", ")
))
log_memory("After metadata merge", full_object)

rm(scice_results)
invisible(gc())
log_memory("After rm(scice_results) + gc()", full_object)

qs::qsave(full_object, output_qs, preset = "balanced")
log_step(sprintf("Saved augmented Seurat object to: %s", output_qs))

log_step("Run complete.")
