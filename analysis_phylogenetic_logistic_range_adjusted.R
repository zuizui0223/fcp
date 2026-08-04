#!/usr/bin/env Rscript

# Topology-aware sensitivity analysis for the automated, range-adjusted
# 107-species dataset. Outputs are isolated from the historical 34-species
# analysis and retain a complete OpenTree matching audit.

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phylolm)
  library(rotl)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop(sprintf("Missing value after %s", flag))
  args[[hit + 1]]
}

input <- get_arg("--dataset")
outdir <- get_arg("--outdir")
if (is.null(input) || is.null(outdir)) stop("--dataset and --outdir are required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

components <- c(
  "bio1_q95q05", "bio5_q95q05", "bio6_q95q05", "bio7_q95q05",
  "bio12_q95q05", "bio14_q95q05", "bio15_q95q05", "bio17_q95q05"
)
base_covariates <- c(
  "geographic_hull_area_km2",
  "n_occurrence_records",
  "n_supporting_p1_records"
)

d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("canonical_name", "spatial_scale", components, base_covariates)
missing <- setdiff(required, names(d))
if (length(missing)) stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
d <- d[d$spatial_scale %in% c("within_population", "among_population"), , drop = FALSE]
d <- d[!duplicated(d$canonical_name), , drop = FALSE]
d$among <- as.integer(d$spatial_scale == "among_population")
if (nrow(d) < 20 || length(unique(d$among)) < 2) stop("Insufficient binary model dataset")

# OpenTree name resolution. Keep the raw return unchanged for audit because
# column availability can vary across rotl/OpenTree versions.
tnrs <- tnrs_match_names(unique(d$canonical_name), context_name = "Land plants")
write.csv(tnrs, file.path(outdir, "opentree_tnrs_audit.csv"), row.names = FALSE)
if (!all(c("search_string", "ott_id") %in% names(tnrs))) {
  stop(sprintf("Unexpected TNRS schema: %s", paste(names(tnrs), collapse = ", ")))
}
matched <- tnrs[!is.na(tnrs$ott_id) & nzchar(tnrs$search_string), , drop = FALSE]
matched <- matched[!duplicated(matched$search_string), , drop = FALSE]
matched <- matched[!duplicated(matched$ott_id), , drop = FALSE]
write.csv(matched, file.path(outdir, "opentree_tnrs_matched_unique.csv"), row.names = FALSE)
if (nrow(matched) < 20) stop(sprintf("Only %d species matched OpenTree", nrow(matched)))

# rotl::tol_induced_subtree has returned a phylo object directly in current
# releases, while some wrappers expose it under $phylo. Accept both formats.
subtree_result <- tol_induced_subtree(ott_ids = as.numeric(matched$ott_id))
if (inherits(subtree_result, "phylo")) {
  raw_tree <- subtree_result
} else if (is.list(subtree_result) && inherits(subtree_result$phylo, "phylo")) {
  raw_tree <- subtree_result$phylo
} else {
  stop(sprintf("Unexpected induced-subtree class: %s", paste(class(subtree_result), collapse = ", ")))
}
write.tree(raw_tree, file.path(outdir, "opentree_induced_raw.newick"))

extract_ott <- function(x) {
  hit <- regmatches(x, regexpr("ott[0-9]+", x))
  suppressWarnings(as.numeric(sub("^ott", "", hit)))
}
tip_ott <- vapply(raw_tree$tip.label, extract_ott, numeric(1))
name_by_ott <- setNames(matched$search_string, as.character(as.numeric(matched$ott_id)))
tip_names <- unname(name_by_ott[as.character(tip_ott)])
valid_tip <- is.finite(tip_ott) & !is.na(tip_names) & nzchar(tip_names)
if (any(!valid_tip)) raw_tree <- drop.tip(raw_tree, raw_tree$tip.label[!valid_tip])

tip_ott <- vapply(raw_tree$tip.label, extract_ott, numeric(1))
raw_tree$tip.label <- unname(name_by_ott[as.character(tip_ott)])
if (anyDuplicated(raw_tree$tip.label)) {
  raw_tree <- drop.tip(raw_tree, raw_tree$tip.label[duplicated(raw_tree$tip.label)])
}
if (length(raw_tree$tip.label) < 20) stop("Fewer than 20 uniquely mapped tree tips")

# OpenTree topology is not time calibrated. Grafen branch lengths make this a
# topology-aware sensitivity analysis rather than a dated evolutionary model.
resolved_tree <- multi2di(raw_tree, random = FALSE)
resolved_tree <- compute.brlen(resolved_tree, method = "Grafen", power = 1)
resolved_tree <- ladderize(resolved_tree)
write.tree(resolved_tree, file.path(outdir, "opentree_grafen_resolved.newick"))

model_data <- d[d$canonical_name %in% resolved_tree$tip.label, , drop = FALSE]
model_data <- model_data[match(resolved_tree$tip.label, model_data$canonical_name), , drop = FALSE]
rownames(model_data) <- model_data$canonical_name

zlog <- function(x, label) {
  y <- log1p(pmax(as.numeric(x), 0))
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop(sprintf("Predictor %s is constant or invalid", label))
  as.numeric((y - mean(y, na.rm = TRUE)) / s)
}
for (v in c(components, base_covariates)) {
  model_data[[paste0(v, "_z")]] <- zlog(model_data[[v]], v)
}
write.csv(model_data, file.path(outdir, "phylogenetic_model_dataset.csv"), row.names = FALSE)

results <- list()
coefficients <- list()
comparison <- list()
failures <- list()
for (metric in components) {
  vars <- c("among", paste0(metric, "_z"), paste0(base_covariates, "_z"))
  dd <- model_data[complete.cases(model_data[, vars, drop = FALSE]), , drop = FALSE]
  tr <- keep.tip(resolved_tree, rownames(dd))
  dd <- dd[match(tr$tip.label, rownames(dd)), , drop = FALSE]
  rownames(dd) <- tr$tip.label
  if (nrow(dd) < 20 || length(unique(dd$among)) < 2) {
    failures[[metric]] <- data.frame(metric = metric, reason = "insufficient_complete_binary_data")
    next
  }
  formula_text <- paste(
    "among ~",
    paste(c(paste0(metric, "_z"), paste0(base_covariates, "_z")), collapse = " + ")
  )
  form <- as.formula(formula_text)
  fitted <- tryCatch({
    fit <- phyloglm(form, data = dd, phy = tr, method = "logistic_MPLE", btol = 50)
    sm <- summary(fit)
    cf <- as.data.frame(sm$coefficients)
    cf$term <- rownames(cf)
    rownames(cf) <- NULL
    if (ncol(cf) < 4) stop("Unexpected phyloglm coefficient table")
    names(cf)[1:4] <- c("estimate", "std_error", "z_value", "p_value")
    cf$metric <- metric
    focal <- cf[cf$term == paste0(metric, "_z"), , drop = FALSE]
    if (nrow(focal) != 1 || !all(is.finite(unlist(focal[, c("estimate", "std_error", "p_value")]))) ) {
      stop("Missing or non-finite focal coefficient")
    }
    nonphy <- glm(form, data = dd, family = binomial())
    nonphy_cf <- summary(nonphy)$coefficients[paste0(metric, "_z"), ]
    coefficients[[metric]] <<- cf
    results[[metric]] <<- data.frame(
      metric = metric,
      n_species = nrow(dd),
      n_within = sum(dd$among == 0),
      n_among = sum(dd$among == 1),
      estimate = focal$estimate,
      std_error = focal$std_error,
      odds_ratio = exp(focal$estimate),
      odds_ratio_ci_low = exp(focal$estimate - 1.96 * focal$std_error),
      odds_ratio_ci_high = exp(focal$estimate + 1.96 * focal$std_error),
      p_value = focal$p_value,
      phylogenetic_alpha = fit$alpha,
      logLik = as.numeric(logLik(fit)),
      method = "phyloglm_logistic_MPLE_OpenTree_Grafen",
      stringsAsFactors = FALSE
    )
    comparison[[metric]] <<- data.frame(
      metric = metric,
      n_species = nrow(dd),
      phylogenetic_estimate = focal$estimate,
      phylogenetic_p = focal$p_value,
      nonphylogenetic_same_subset_estimate = nonphy_cf[[1]],
      nonphylogenetic_same_subset_p = nonphy_cf[[4]],
      stringsAsFactors = FALSE
    )
    TRUE
  }, error = function(e) {
    failures[[metric]] <<- data.frame(metric = metric, reason = conditionMessage(e))
    FALSE
  })
  message(sprintf("%s: %s", metric, if (fitted) "complete" else "failed"))
}

if (length(failures)) {
  write.csv(do.call(rbind, failures), file.path(outdir, "phylogenetic_model_failures.csv"), row.names = FALSE)
}
result <- if (length(results)) do.call(rbind, results) else NULL
if (is.null(result) || nrow(result) != length(components)) {
  stop(sprintf("Only %d of %d phylogenetic models were estimable", if (is.null(result)) 0 else nrow(result), length(components)))
}
result$p_holm_eight_components <- p.adjust(result$p_value, method = "holm")
coef_table <- do.call(rbind, coefficients)
comparison_table <- do.call(rbind, comparison)
write.csv(result, file.path(outdir, "phylogenetic_logistic_component_models.csv"), row.names = FALSE)
write.csv(coef_table, file.path(outdir, "phylogenetic_logistic_all_coefficients.csv"), row.names = FALSE)
write.csv(comparison_table, file.path(outdir, "phylogenetic_vs_nonphylogenetic_same_subset.csv"), row.names = FALSE)

manifest <- list(
  status = "complete",
  input_species = nrow(d),
  opentree_matched_unique_species = nrow(matched),
  tree_tip_species = length(resolved_tree$tip.label),
  model_species_min = min(result$n_species),
  model_species_max = max(result$n_species),
  tree_source = "Open Tree of Life induced synthetic subtree",
  branch_length_method = "Grafen power=1; topology-aware sensitivity, not divergence-time calibration",
  polytomy_resolution = "ape::multi2di(random=FALSE)",
  estimator = "phylolm::phyloglm(method='logistic_MPLE')",
  covariates = base_covariates,
  results = result,
  interpretation_guard = paste(
    "Automated spatial labels remain exploratory.",
    "The model asks whether climate-breadth associations persist after topology-based dependence",
    "and the same range/effort covariates; it is not causal or time calibrated."
  )
)
write_json(manifest, file.path(outdir, "phylogenetic_analysis_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
print(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"))
