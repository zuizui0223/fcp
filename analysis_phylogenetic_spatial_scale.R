#!/usr/bin/env Rscript

# Topology-based phylogenetic logistic sensitivity analysis for the JBI manuscript.
# Exact, unique Open Tree Taxonomy matches are used to induce a synthetic-tree
# topology. Polytomies are resolved repeatedly and assigned Grafen branch lengths.
# This is a topology sensitivity analysis, not a time-calibrated species phylogeny.

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phylolm)
  library(rotl)
})

parse_args <- function(args) {
  out <- list(
    dataset = NULL,
    outdir = NULL,
    analysis_label = "baseline_original_occurrences",
    replicates = 100L,
    seed = 20260724L,
    min_cells = 20L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    if (!key %in% names(out)) stop("Unknown argument: --", key)
    value <- args[[i + 1L]]
    if (key %in% c("replicates", "seed", "min_cells")) value <- as.integer(value)
    out[[key]] <- value
    i <- i + 2L
  }
  if (is.null(out$dataset) || is.null(out$outdir)) {
    stop("--dataset and --outdir are required")
  }
  out
}

retry_call <- function(fun, attempts = 4L, label = "API call") {
  last <- NULL
  for (i in seq_len(attempts)) {
    result <- try(fun(), silent = TRUE)
    if (!inherits(result, "try-error")) return(result)
    last <- result
    if (i < attempts) Sys.sleep(2^(i - 1L))
  }
  stop(label, " failed after ", attempts, " attempts: ", as.character(last))
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

prepare_tree <- function(tree, seed) {
  set.seed(seed)
  resolved <- if (is.binary(tree)) tree else multi2di(tree, random = TRUE)
  resolved$edge.length <- NULL
  resolved <- compute.brlen(resolved, method = "Grafen", power = 1)
  positive <- resolved$edge.length[
    is.finite(resolved$edge.length) & resolved$edge.length > 0
  ]
  if (!length(positive)) stop("Could not assign positive branch lengths")
  resolved$edge.length[
    !is.finite(resolved$edge.length) | resolved$edge.length <= 0
  ] <- min(positive) * 1e-6
  resolved
}

failed_fit_row <- function(label, replicate_id, seed, status, n_species) {
  data.frame(
    analysis_label = label,
    replicate = replicate_id,
    seed = seed,
    status = status,
    n_species = n_species,
    estimate = NA_real_,
    std_error = NA_real_,
    odds_ratio = NA_real_,
    odds_ratio_ci_low = NA_real_,
    odds_ratio_ci_high = NA_real_,
    p_value = NA_real_,
    alpha = NA_real_,
    log_likelihood = NA_real_,
    stringsAsFactors = FALSE
  )
}

fit_one <- function(data, tree, replicate_id, seed, analysis_label) {
  tr <- prepare_tree(tree, seed)
  kept <- intersect(tr$tip.label, rownames(data))
  to_drop <- setdiff(tr$tip.label, kept)
  if (length(to_drop)) tr <- drop.tip(tr, to_drop)
  d <- data[tr$tip.label, , drop = FALSE]
  if (nrow(d) < 20L || length(unique(d$among)) < 2L) {
    return(failed_fit_row(
      analysis_label, replicate_id, seed, "not_estimable", nrow(d)
    ))
  }

  fit <- try(
    phyloglm(
      among ~ moisture_z + effort_z,
      data = d,
      phy = tr,
      method = "logistic_MPLE",
      boot = 0
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(failed_fit_row(
      analysis_label,
      replicate_id,
      seed,
      paste0("failed: ", substr(as.character(fit), 1, 300)),
      nrow(d)
    ))
  }

  beta <- unname(coef(fit)[["moisture_z"]])
  covariance <- try(vcov(fit), silent = TRUE)
  se <- if (!inherits(covariance, "try-error")) {
    sqrt(diag(covariance))[["moisture_z"]]
  } else {
    NA_real_
  }
  z <- beta / se
  p <- if (is.finite(z)) 2 * pnorm(-abs(z)) else NA_real_
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  ll <- try(as.numeric(logLik(fit)), silent = TRUE)

  data.frame(
    analysis_label = analysis_label,
    replicate = replicate_id,
    seed = seed,
    status = "complete",
    n_species = nrow(d),
    estimate = beta,
    std_error = se,
    odds_ratio = exp(beta),
    odds_ratio_ci_low = exp(ci_low),
    odds_ratio_ci_high = exp(ci_high),
    p_value = p,
    alpha = if (!is.null(fit$alpha)) unname(fit$alpha) else NA_real_,
    log_likelihood = if (!inherits(ll, "try-error")) ll else NA_real_,
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
error_log <- file.path(args$outdir, "phylogenetic_error.log")
options(error = function() {
  lines <- c(
    paste0("ERROR: ", geterrmessage()),
    "TRACEBACK:",
    capture.output(traceback(20))
  )
  writeLines(lines, error_log)
  message(paste(lines, collapse = "\n"))
  quit(save = "no", status = 1, runLast = FALSE)
})

raw <- read.csv(args$dataset, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "canonical_name", "family", "spatial_scale", "classification_source",
  "n_climate_cells", "moisture_breadth"
)
missing <- setdiff(required, names(raw))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

data <- raw[
  raw$classification_source == "baseline_unambiguous" &
    raw$n_climate_cells >= args$min_cells &
    raw$spatial_scale %in% c("within_population", "among_population"),
  required,
  drop = FALSE
]
data <- data[!duplicated(data$canonical_name), , drop = FALSE]
data$among <- as.integer(data$spatial_scale == "among_population")
data$moisture_z <- zscore(data$moisture_breadth)
data$effort_z <- zscore(log1p(data$n_climate_cells))
data <- data[complete.cases(data[, c("among", "moisture_z", "effort_z")]), , drop = FALSE]
if (nrow(data) < 20L) stop("Fewer than 20 complete species before tree matching")

resolved <- retry_call(
  function() tnrs_match_names(
    data$canonical_name,
    context_name = "Land plants",
    do_approximate_matching = FALSE,
    include_suppressed = FALSE
  ),
  label = "Open Tree TNRS"
)
if (!nrow(resolved)) stop("Open Tree TNRS returned no rows")

resolved$query_name <- data$canonical_name[
  match(tolower(resolved$search_string), tolower(data$canonical_name))
]
resolved$exact_unique_match <- (
  !is.na(resolved$query_name) &
    !is.na(resolved$ott_id) &
    !is.na(resolved$approximate_match) &
    !resolved$approximate_match &
    is.finite(resolved$score) &
    resolved$score == 1 &
    !is.na(resolved$number_matches) &
    resolved$number_matches == 1
)
resolved$in_synthetic_tree <- FALSE

eligible_index <- which(resolved$exact_unique_match)
if (length(eligible_index)) {
  ids <- as.integer(resolved$ott_id[eligible_index])
  membership <- retry_call(
    function() is_in_tree(ids),
    label = "Open Tree synthetic-tree membership"
  )
  if (length(membership) != length(ids)) {
    stop("Open Tree membership result length did not match the requested OTT ids")
  }
  resolved$in_synthetic_tree[eligible_index] <- as.logical(unname(membership))
}

resolved$duplicate_ott_id <- FALSE
non_missing_ids <- which(!is.na(resolved$ott_id))
resolved$duplicate_ott_id[non_missing_ids] <- (
  duplicated(resolved$ott_id[non_missing_ids]) |
    duplicated(resolved$ott_id[non_missing_ids], fromLast = TRUE)
)
resolved$model_eligible <- (
  resolved$exact_unique_match &
    resolved$in_synthetic_tree &
    !resolved$duplicate_ott_id
)

audit_columns <- intersect(
  c(
    "query_name", "search_string", "unique_name", "approximate_match", "score",
    "ott_id", "is_synonym", "flags", "number_matches", "exact_unique_match",
    "in_synthetic_tree", "duplicate_ott_id", "model_eligible"
  ),
  names(resolved)
)
write.csv(
  resolved[, audit_columns, drop = FALSE],
  file.path(args$outdir, "phylogenetic_name_resolution.csv"),
  row.names = FALSE
)

eligible <- resolved[resolved$model_eligible, , drop = FALSE]
if (nrow(eligible) < 20L) {
  stop(
    "Fewer than 20 species had exact, unique Open Tree matches in the synthetic tree: ",
    nrow(eligible)
  )
}

ott_ids <- as.integer(eligible$ott_id)
tree <- retry_call(
  function() tol_induced_subtree(ott_ids = ott_ids, label_format = "id"),
  label = "Open Tree induced subtree"
)

tip_ids <- suppressWarnings(as.integer(gsub("[^0-9]", "", tree$tip.label)))
tip_names <- eligible$query_name[match(tip_ids, as.integer(eligible$ott_id))]
if (anyNA(tip_names)) {
  stop("Could not map all induced-subtree tip labels back to input species")
}
tree$tip.label <- tip_names
tips_to_drop <- setdiff(tree$tip.label, eligible$query_name)
if (length(tips_to_drop)) tree <- drop.tip(tree, tips_to_drop)

data <- data[data$canonical_name %in% tree$tip.label, , drop = FALSE]
rownames(data) <- data$canonical_name
write.tree(tree, file = file.path(args$outdir, "opentree_induced_topology.tre"))

fits <- vector("list", args$replicates)
for (i in seq_len(args$replicates)) {
  fits[[i]] <- fit_one(
    data,
    tree,
    replicate_id = i,
    seed = args$seed + i - 1L,
    analysis_label = args$analysis_label
  )
}
results <- do.call(rbind, fits)
write.csv(
  results,
  file.path(args$outdir, "phylogenetic_logistic_replicates.csv"),
  row.names = FALSE
)

complete <- results[results$status == "complete", , drop = FALSE]
if (nrow(complete)) {
  summary <- data.frame(
    analysis_label = args$analysis_label,
    n_input_species = nrow(raw),
    n_model_species = unique(complete$n_species)[1],
    n_replicates_requested = args$replicates,
    n_replicates_complete = nrow(complete),
    median_estimate = median(complete$estimate),
    min_estimate = min(complete$estimate),
    max_estimate = max(complete$estimate),
    median_odds_ratio = median(complete$odds_ratio),
    min_odds_ratio = min(complete$odds_ratio),
    max_odds_ratio = max(complete$odds_ratio),
    median_ci_low = median(complete$odds_ratio_ci_low, na.rm = TRUE),
    median_ci_high = median(complete$odds_ratio_ci_high, na.rm = TRUE),
    median_p_value = median(complete$p_value, na.rm = TRUE),
    fraction_negative = mean(complete$estimate < 0),
    median_alpha = median(complete$alpha, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
} else {
  summary <- data.frame(
    analysis_label = args$analysis_label,
    n_input_species = nrow(raw),
    n_model_species = nrow(data),
    n_replicates_requested = args$replicates,
    n_replicates_complete = 0,
    median_estimate = NA_real_,
    min_estimate = NA_real_,
    max_estimate = NA_real_,
    median_odds_ratio = NA_real_,
    min_odds_ratio = NA_real_,
    max_odds_ratio = NA_real_,
    median_ci_low = NA_real_,
    median_ci_high = NA_real_,
    median_p_value = NA_real_,
    fraction_negative = NA_real_,
    median_alpha = NA_real_,
    stringsAsFactors = FALSE
  )
}
write.csv(
  summary,
  file.path(args$outdir, "phylogenetic_logistic_summary.csv"),
  row.names = FALSE
)

manifest <- list(
  analysis_label = args$analysis_label,
  input_dataset = args$dataset,
  min_cells = args$min_cells,
  formula = "among ~ moisture_z + effort_z",
  method = "phylolm::phyloglm(method='logistic_MPLE')",
  phylogeny = paste(
    "Open Tree of Life induced synthetic-tree topology; exact, unique,",
    "non-approximate name matches only; duplicate OTT IDs excluded"
  ),
  polytomy_treatment = paste(
    args$replicates,
    "random ape::multi2di resolutions; results summarized across resolutions"
  ),
  branch_lengths = "Grafen branch lengths assigned separately to every resolved topology",
  seed = args$seed,
  package_versions = list(
    R = R.version.string,
    ape = as.character(packageVersion("ape")),
    jsonlite = as.character(packageVersion("jsonlite")),
    phylolm = as.character(packageVersion("phylolm")),
    rotl = as.character(packageVersion("rotl"))
  ),
  n_name_queries = nrow(resolved),
  n_model_eligible_names = sum(resolved$model_eligible),
  results = summary,
  interpretation_guard = paste(
    "Topology-based phylogenetic sensitivity analysis. The Open Tree topology and",
    "Grafen branch lengths are not a time-calibrated species phylogeny. Results do not",
    "establish causation, local adaptation or morph-specific climatic tolerance."
  )
)
write_json(
  manifest,
  file.path(args$outdir, "phylogenetic_logistic_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)
print(manifest)
