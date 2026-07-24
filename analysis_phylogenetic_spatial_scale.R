#!/usr/bin/env Rscript

# Topology-based phylogenetic logistic sensitivity analysis for the JBI manuscript.
# The analysis uses exact Open Tree of Life name matches, an induced synthetic-tree
# topology, random resolutions of polytomies, Grafen branch lengths, and phyloglm
# with the logistic_MPLE method. It is a sensitivity analysis, not a time-calibrated
# species phylogeny.

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
    value <- args[[i + 1L]]
    if (!key %in% names(out)) stop("Unknown argument: --", key)
    if (key %in% c("replicates", "seed", "min_cells")) value <- as.integer(value)
    out[[key]] <- value
    i <- i + 2L
  }
  if (is.null(out$dataset) || is.null(out$outdir)) {
    stop("--dataset and --outdir are required")
  }
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_name_resolution <- function(names) {
  result <- try(
    tnrs_match_names(
      names,
      context_name = "Land plants",
      do_approximate_matching = FALSE,
      include_suppressed = FALSE
    ),
    silent = TRUE
  )
  if (inherits(result, "try-error")) {
    result <- tnrs_match_names(
      names,
      context_name = "All life",
      do_approximate_matching = FALSE,
      include_suppressed = FALSE
    )
  }
  result
}

prepare_tree <- function(tree, seed) {
  set.seed(seed)
  resolved <- if (is.binary.tree(tree)) tree else multi2di(tree, random = TRUE)
  resolved$edge.length <- NULL
  resolved <- compute.brlen(resolved, method = "Grafen", power = 1)
  positive <- resolved$edge.length[is.finite(resolved$edge.length) & resolved$edge.length > 0]
  if (!length(positive)) stop("Could not assign positive branch lengths")
  resolved$edge.length[!is.finite(resolved$edge.length) | resolved$edge.length <= 0] <-
    min(positive) * 1e-6
  resolved
}

fit_one <- function(data, tree, replicate_id, seed, analysis_label) {
  tr <- prepare_tree(tree, seed)
  kept <- intersect(tr$tip.label, rownames(data))
  tr <- drop.tip(tr, setdiff(tr$tip.label, kept))
  d <- data[tr$tip.label, , drop = FALSE]
  if (nrow(d) < 20L || length(unique(d$among)) < 2L) {
    return(data.frame(
      analysis_label = analysis_label,
      replicate = replicate_id,
      seed = seed,
      status = "not_estimable",
      n_species = nrow(d),
      estimate = NA_real_,
      std_error = NA_real_,
      odds_ratio = NA_real_,
      odds_ratio_ci_low = NA_real_,
      odds_ratio_ci_high = NA_real_,
      p_value = NA_real_,
      alpha = NA_real_,
      log_likelihood = NA_real_,
      stringsAsFactors = FALSE
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
    return(data.frame(
      analysis_label = analysis_label,
      replicate = replicate_id,
      seed = seed,
      status = paste0("failed: ", substr(as.character(fit), 1, 200)),
      n_species = nrow(d),
      estimate = NA_real_,
      std_error = NA_real_,
      odds_ratio = NA_real_,
      odds_ratio_ci_low = NA_real_,
      odds_ratio_ci_high = NA_real_,
      p_value = NA_real_,
      alpha = NA_real_,
      log_likelihood = NA_real_,
      stringsAsFactors = FALSE
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

resolved <- safe_name_resolution(data$canonical_name)
resolved$query_name <- data$canonical_name[
  match(tolower(resolved$search_string), tolower(data$canonical_name))
]
resolved$exact_eligible <- (
  !is.na(resolved$query_name) &
    !is.na(resolved$ott_id) &
    !resolved$approximate_match &
    resolved$score == 1
)

tree_membership <- rep(FALSE, nrow(resolved))
if (any(resolved$exact_eligible)) {
  eligible_ids <- resolved$ott_id[resolved$exact_eligible]
  in_tree <- is_in_tree(eligible_ids)
  tree_membership[resolved$exact_eligible] <- as.logical(in_tree)
}
resolved$in_synthetic_tree <- tree_membership
resolved$duplicate_ott_id <- duplicated(resolved$ott_id) | duplicated(resolved$ott_id, fromLast = TRUE)
resolved$model_eligible <- resolved$exact_eligible & resolved$in_synthetic_tree & !resolved$duplicate_ott_id

audit_columns <- intersect(
  c(
    "query_name", "search_string", "unique_name", "approximate_match", "score",
    "ott_id", "is_synonym", "flags", "number_matches", "exact_eligible",
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
  stop("Fewer than 20 species had exact, unique Open Tree matches in the synthetic tree")
}

ott_ids <- eligible$ott_id
names(ott_ids) <- eligible$query_name
tree <- tol_induced_subtree(ott_ids = ott_ids, label_format = "id")

tip_ids <- suppressWarnings(as.numeric(gsub("[^0-9]", "", tree$tip.label)))
tip_names <- eligible$query_name[match(tip_ids, eligible$ott_id)]
if (anyNA(tip_names)) {
  stop("Could not map all induced-subtree tip labels back to input species")
}
tree$tip.label <- tip_names
tree <- drop.tip(tree, setdiff(tree$tip.label, eligible$query_name))

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
    median_ci_low = median(complete$odds_ratio_ci_low),
    median_ci_high = median(complete$odds_ratio_ci_high),
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
    "Open Tree of Life induced synthetic-tree topology; exact non-approximate",
    "name matches only; duplicate OTT IDs excluded"
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
