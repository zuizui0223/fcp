#!/usr/bin/env Rscript

# Topology-aware sensitivity analysis for the automated 107-species dataset.
# Each climate component is fitted independently so one numerical failure does
# not erase successful models or diagnostic artifacts.

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
d$among <- as.integer(d$spatial_scale == "among_population")

# Open Tree's TNRS echoes the query back lower-cased in `search_string`, so every
# join between TNRS output and the model dataset has to go through a normalised
# key rather than through the raw strings.
name_key <- function(x) tolower(trimws(as.character(x)))
d$canonical_name <- trimws(d$canonical_name)
d <- d[!is.na(d$canonical_name) & nzchar(d$canonical_name), , drop = FALSE]
d$name_key <- name_key(d$canonical_name)
d <- d[!duplicated(d$name_key), , drop = FALSE]
if (nrow(d) < 20 || length(unique(d$among)) < 2) stop("Insufficient binary-response data")

message(sprintf("Matching %d species against Open Tree of Life", nrow(d)))
tnrs <- tnrs_match_names(d$canonical_name, context_name = "Land plants")
if (!all(c("search_string", "ott_id") %in% names(tnrs))) {
  stop(sprintf("Unexpected TNRS columns: %s", paste(names(tnrs), collapse = ", ")))
}
tnrs$query_name <- d$canonical_name[match(name_key(tnrs$search_string), d$name_key)]
write.csv(tnrs, file.path(outdir, "opentree_tnrs_audit.csv"), row.names = FALSE)
matched <- tnrs[!is.na(tnrs$ott_id), , drop = FALSE]
matched <- matched[!duplicated(name_key(matched$search_string)), , drop = FALSE]
# Recover the dataset spelling of each match; `search_string` is lower-cased and
# cannot be compared against `canonical_name` directly.
matched$canonical_name <- matched$query_name
unresolved <- matched$search_string[is.na(matched$canonical_name)]
if (length(unresolved)) {
  message(sprintf(
    "Dropping %d TNRS rows that do not map back to a dataset species (e.g. %s)",
    length(unresolved), paste(utils::head(unresolved, 5), collapse = ", ")
  ))
}
matched <- matched[!is.na(matched$canonical_name), , drop = FALSE]
matched <- matched[!duplicated(matched$ott_id), , drop = FALSE]
if (nrow(matched) < 20) stop(sprintf("Only %d species matched OpenTree", nrow(matched)))

message(sprintf("Requesting induced subtree for %d unique OTT ids", length(unique(matched$ott_id))))
subtree <- tol_induced_subtree(ott_ids = unique(matched$ott_id), label_format = "name_and_id")
raw_tree <- if (inherits(subtree, "phylo")) subtree else subtree$phylo
if (is.null(raw_tree) || !inherits(raw_tree, "phylo")) {
  stop(sprintf("tol_induced_subtree returned unsupported class: %s", paste(class(subtree), collapse = ", ")))
}
if (length(raw_tree$tip.label) < 20) stop("OpenTree induced subtree contains fewer than 20 tips")
write.tree(raw_tree, file.path(outdir, "opentree_induced_raw.newick"))

extract_ott <- function(x) {
  out <- sub(".*_ott([0-9]+)$", "\\1", x)
  suppressWarnings(as.numeric(out))
}
map_cols <- intersect(
  c("canonical_name", "search_string", "unique_name", "ott_id", "is_synonym", "flags"),
  names(matched)
)
map <- matched[, map_cols, drop = FALSE]
map$ott_id <- as.numeric(map$ott_id)

original_tip_labels <- raw_tree$tip.label
tip_ott <- vapply(original_tip_labels, extract_ott, numeric(1), USE.NAMES = FALSE)
tip_map <- data.frame(
  original_tip_label = original_tip_labels,
  ott_id = tip_ott,
  canonical_name = trimws(map$canonical_name[match(tip_ott, map$ott_id)]),
  stringsAsFactors = FALSE
)
valid_tip <- is.finite(tip_map$ott_id) & !is.na(tip_map$canonical_name) & nzchar(tip_map$canonical_name)
valid_tip <- valid_tip & !duplicated(tip_map$canonical_name)
tip_map$keep_for_analysis <- valid_tip
write.csv(tip_map, file.path(outdir, "opentree_tip_mapping_audit.csv"), row.names = FALSE)
if (sum(valid_tip) < 20) {
  stop(sprintf(
    "Only %d of %d OpenTree tips mapped back to a dataset species; see opentree_tip_mapping_audit.csv (unmapped examples: %s)",
    sum(valid_tip), nrow(tip_map),
    paste(utils::head(tip_map$original_tip_label[!valid_tip], 5), collapse = ", ")
  ))
}
if (any(!valid_tip)) raw_tree <- drop.tip(raw_tree, tip_map$original_tip_label[!valid_tip])

kept <- match(raw_tree$tip.label, tip_map$original_tip_label)
mapped_names <- tip_map$canonical_name[kept]
if (anyNA(mapped_names) || any(!nzchar(trimws(mapped_names))) || anyDuplicated(mapped_names)) {
  stop("OpenTree tip labels could not be mapped uniquely to canonical names")
}
raw_tree$tip.label <- mapped_names
if (length(raw_tree$tip.label) < 20) stop("Fewer than 20 uniquely mapped tree tips remain")

resolved_tree <- multi2di(raw_tree, random = FALSE)
resolved_tree <- compute.brlen(resolved_tree, method = "Grafen", power = 1)
resolved_tree <- ladderize(resolved_tree)
write.tree(resolved_tree, file.path(outdir, "opentree_grafen_resolved.newick"))

tip_index <- match(name_key(resolved_tree$tip.label), d$name_key)
if (anyNA(tip_index) || anyDuplicated(tip_index)) {
  unmatched <- resolved_tree$tip.label[is.na(tip_index)]
  stop(sprintf(
    paste0(
      "Model dataset does not map one-to-one onto the resolved tree: ",
      "%d of %d tips unmatched, %d duplicated. Unmatched examples: %s"
    ),
    length(unmatched), length(tip_index), sum(duplicated(tip_index[!is.na(tip_index)])),
    if (length(unmatched)) paste(utils::head(unmatched, 5), collapse = ", ") else "none"
  ))
}
model_data <- d[tip_index, , drop = FALSE]
rownames(model_data) <- model_data$canonical_name

z_log <- function(x) {
  y <- log1p(pmax(as.numeric(x), 0))
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(y)))
  as.numeric((y - mean(y, na.rm = TRUE)) / s)
}
for (v in c(components, base_covariates)) model_data[[paste0(v, "_z")]] <- z_log(model_data[[v]])
write.csv(model_data, file.path(outdir, "phylogenetic_model_dataset.csv"), row.names = FALSE)

empty_results <- function() data.frame(
  metric = character(), n_species = integer(), n_within = integer(), n_among = integer(),
  estimate = numeric(), std_error = numeric(), odds_ratio = numeric(),
  odds_ratio_ci_low = numeric(), odds_ratio_ci_high = numeric(), p_value = numeric(),
  phylogenetic_alpha = numeric(), method = character(), stringsAsFactors = FALSE
)
empty_comparison <- function() data.frame(
  metric = character(), n_species = integer(), phylogenetic_estimate = numeric(),
  phylogenetic_p = numeric(), nonphylogenetic_same_subset_estimate = numeric(),
  nonphylogenetic_same_subset_p = numeric(), stringsAsFactors = FALSE
)
empty_coefficients <- function() data.frame(
  estimate = numeric(), std_error = numeric(), z_value = numeric(), p_value = numeric(),
  term = character(), metric = character(), stringsAsFactors = FALSE
)

result <- empty_results()
comparison_table <- empty_comparison()
coef_table <- empty_coefficients()
fit_audit <- data.frame(
  metric = components,
  status = "pending",
  stage = NA_character_,
  n_species = NA_integer_,
  n_within = NA_integer_,
  n_among = NA_integer_,
  formula = NA_character_,
  warning = NA_character_,
  error = NA_character_,
  stringsAsFactors = FALSE
)

safe_write <- function() {
  out_result <- result
  if (nrow(out_result)) out_result$p_holm_eight_components <- p.adjust(out_result$p_value, method = "holm")
  write.csv(out_result, file.path(outdir, "phylogenetic_logistic_component_models.csv"), row.names = FALSE)
  write.csv(coef_table, file.path(outdir, "phylogenetic_logistic_all_coefficients.csv"), row.names = FALSE)
  write.csv(comparison_table, file.path(outdir, "phylogenetic_vs_nonphylogenetic_same_subset.csv"), row.names = FALSE)
  write.csv(fit_audit, file.path(outdir, "phylogenetic_model_fit_audit.csv"), row.names = FALSE)
}
safe_write()

extract_coef_table <- function(fit_summary, metric) {
  x <- as.data.frame(fit_summary$coefficients, stringsAsFactors = FALSE)
  x$term <- rownames(x)
  rownames(x) <- NULL
  # Lower-case first: normalising the other way round strips the capitals out of
  # phylolm's own column names (Estimate -> "stimate", StdErr -> "tdrr").
  norm <- gsub("[^a-z0-9]", "", tolower(names(x)))
  find_col <- function(candidates) {
    hit <- match(candidates, norm, nomatch = 0)
    hit <- hit[hit > 0]
    if (!length(hit)) NA_integer_ else hit[[1]]
  }
  i_est <- find_col(c("estimate", "coef", "coefficient"))
  i_se <- find_col(c("stderr", "stderror", "se"))
  i_z <- find_col(c("zvalue", "z", "tvalue"))
  i_p <- find_col(c("pvalue", "przz", "prz", "prt", "p"))
  if (anyNA(c(i_est, i_se, i_p))) {
    stop(sprintf("Unsupported phyloglm coefficient columns for %s: %s", metric, paste(names(x), collapse = ", ")))
  }
  z <- if (is.na(i_z)) as.numeric(x[[i_est]]) / as.numeric(x[[i_se]]) else as.numeric(x[[i_z]])
  data.frame(
    estimate = as.numeric(x[[i_est]]),
    std_error = as.numeric(x[[i_se]]),
    z_value = z,
    p_value = as.numeric(x[[i_p]]),
    term = x$term,
    metric = metric,
    stringsAsFactors = FALSE
  )
}

for (metric in components) {
  idx <- match(metric, fit_audit$metric)
  warning_messages <- character()
  fit_audit$stage[idx] <- "prepare_data"
  vars <- c("among", paste0(metric, "_z"), paste0(base_covariates, "_z"))
  dd <- model_data[complete.cases(model_data[, vars, drop = FALSE]), , drop = FALSE]
  fit_audit$n_species[idx] <- nrow(dd)
  fit_audit$n_within[idx] <- sum(dd$among == 0)
  fit_audit$n_among[idx] <- sum(dd$among == 1)

  tryCatch({
    if (nrow(dd) < 20 || length(unique(dd$among)) < 2) {
      stop(sprintf("Insufficient complete data: n=%d, classes=%d", nrow(dd), length(unique(dd$among))))
    }
    tr <- keep.tip(resolved_tree, rownames(dd))
    dd <- dd[match(tr$tip.label, rownames(dd)), , drop = FALSE]
    rownames(dd) <- tr$tip.label

    formula_text <- paste("among ~", paste(c(paste0(metric, "_z"), paste0(base_covariates, "_z")), collapse = " + "))
    fit_audit$formula[idx] <- formula_text
    form <- as.formula(formula_text)
    fit_audit$stage[idx] <- "phyloglm"
    message(sprintf("Fitting %s with %d species", metric, nrow(dd)))
    fit <- withCallingHandlers(
      phyloglm(form, data = dd, phy = tr, method = "logistic_MPLE", btol = 50),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    fit_audit$stage[idx] <- "extract_coefficients"
    cf <- extract_coef_table(summary(fit), metric)
    focal <- cf[cf$term == paste0(metric, "_z"), , drop = FALSE]
    if (nrow(focal) != 1 || !is.finite(focal$estimate) || !is.finite(focal$std_error) || !is.finite(focal$p_value)) {
      stop("Focal coefficient is missing or non-finite")
    }

    fit_audit$stage[idx] <- "nonphylogenetic_comparison"
    nonphy <- glm(form, data = dd, family = binomial())
    nonphy_summary <- summary(nonphy)$coefficients
    focal_name <- paste0(metric, "_z")
    if (!(focal_name %in% rownames(nonphy_summary))) stop("Focal coefficient absent from same-subset GLM")
    nonphy_cf <- nonphy_summary[focal_name, ]

    result <- rbind(result, data.frame(
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
      method = "phyloglm_logistic_MPLE_OpenTree_Grafen",
      stringsAsFactors = FALSE
    ))
    coef_table <- rbind(coef_table, cf)
    comparison_table <- rbind(comparison_table, data.frame(
      metric = metric,
      n_species = nrow(dd),
      phylogenetic_estimate = focal$estimate,
      phylogenetic_p = focal$p_value,
      nonphylogenetic_same_subset_estimate = unname(nonphy_cf[[1]]),
      nonphylogenetic_same_subset_p = unname(nonphy_cf[[4]]),
      stringsAsFactors = FALSE
    ))
    fit_audit$status[idx] <- "success"
    fit_audit$stage[idx] <- "complete"
    fit_audit$warning[idx] <- if (length(warning_messages)) paste(unique(warning_messages), collapse = " | ") else NA_character_
  }, error = function(e) {
    fit_audit$status[idx] <<- "failed"
    fit_audit$error[idx] <<- conditionMessage(e)
    fit_audit$warning[idx] <<- if (length(warning_messages)) paste(unique(warning_messages), collapse = " | ") else NA_character_
    message(sprintf("FAILED %s at %s: %s", metric, fit_audit$stage[idx], conditionMessage(e)))
  })
  safe_write()
}

n_success <- sum(fit_audit$status == "success")
n_failed <- sum(fit_audit$status == "failed")
final_status <- if (n_success == length(components)) "complete" else if (n_success > 0) "partial" else "failed"

manifest <- list(
  status = final_status,
  input_species = nrow(d),
  opentree_matched_species = nrow(matched),
  tree_tip_species = length(resolved_tree$tip.label),
  models_requested = length(components),
  models_successful = n_success,
  models_failed = n_failed,
  successful_metrics = fit_audit$metric[fit_audit$status == "success"],
  failed_metrics = fit_audit$metric[fit_audit$status == "failed"],
  tree_source = "Open Tree of Life induced synthetic subtree",
  branch_length_method = "Grafen power=1; topology-aware sensitivity, not divergence-time calibration",
  polytomy_resolution = "ape::multi2di(random=FALSE)",
  estimator = "phylolm::phyloglm(method='logistic_MPLE')",
  covariates = base_covariates,
  interpretation_guard = "Automated spatial labels remain exploratory. Collinearity diagnostics and pre-specified range/effort covariates must be considered when interpreting component coefficients; this is not a causal or time-calibrated evolutionary model."
)
write_json(manifest, file.path(outdir, "phylogenetic_analysis_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
print(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"))

if (n_success == 0) stop("All eight phylogenetic component models failed; inspect phylogenetic_model_fit_audit.csv")
