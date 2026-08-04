#!/usr/bin/env Rscript

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

# OpenTree taxonomic matching is stored in full so exclusions are auditable.
tnrs <- tnrs_match_names(unique(d$canonical_name), context_name = "Land plants")
write.csv(tnrs, file.path(outdir, "opentree_tnrs_audit.csv"), row.names = FALSE)
matched <- tnrs[!is.na(tnrs$ott_id) & tnrs$is_synonym %in% c(TRUE, FALSE), , drop = FALSE]
matched <- matched[!duplicated(matched$search_string), , drop = FALSE]
if (nrow(matched) < 20) stop(sprintf("Only %d species matched OpenTree", nrow(matched)))

subtree <- tol_induced_subtree(ott_ids = matched$ott_id)
raw_tree <- subtree$phylo
write.tree(raw_tree, file.path(outdir, "opentree_induced_raw.newick"))

extract_ott <- function(x) sub(".*_ott([0-9]+)$", "\\1", x)
tip_ott <- suppressWarnings(as.numeric(vapply(raw_tree$tip.label, extract_ott, character(1))))
map <- matched[, c("search_string", "unique_name", "ott_id", "is_synonym", "flags"), drop = FALSE]
map$ott_id <- as.numeric(map$ott_id)
tip_names <- map$search_string[match(tip_ott, map$ott_id)]
valid_tip <- !is.na(tip_names) & nzchar(tip_names)
if (any(!valid_tip)) raw_tree <- drop.tip(raw_tree, raw_tree$tip.label[!valid_tip])
tip_ott <- suppressWarnings(as.numeric(vapply(raw_tree$tip.label, extract_ott, character(1))))
raw_tree$tip.label <- map$search_string[match(tip_ott, map$ott_id)]
raw_tree <- drop.tip(raw_tree, raw_tree$tip.label[duplicated(raw_tree$tip.label)])

# OpenTree supplies topology but not calibrated branch lengths. Grafen lengths are used
# only for a topology-aware sensitivity analysis; this limitation is explicit in outputs.
set.seed(20260804)
resolved_tree <- multi2di(raw_tree, random = FALSE)
resolved_tree <- compute.brlen(resolved_tree, method = "Grafen", power = 1)
resolved_tree <- ladderize(resolved_tree)
write.tree(resolved_tree, file.path(outdir, "opentree_grafen_resolved.newick"))

model_data <- d[d$canonical_name %in% resolved_tree$tip.label, , drop = FALSE]
model_data <- model_data[match(resolved_tree$tip.label, model_data$canonical_name), , drop = FALSE]
rownames(model_data) <- model_data$canonical_name

z <- function(x) {
  y <- log1p(pmax(as.numeric(x), 0))
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(y)))
  as.numeric((y - mean(y, na.rm = TRUE)) / s)
}
for (v in c(components, base_covariates)) model_data[[paste0(v, "_z")]] <- z(model_data[[v]])

results <- list()
coefficients <- list()
comparison <- list()
for (metric in components) {
  vars <- c("among", paste0(metric, "_z"), paste0(base_covariates, "_z"))
  dd <- model_data[complete.cases(model_data[, vars, drop = FALSE]), , drop = FALSE]
  tr <- keep.tip(resolved_tree, rownames(dd))
  dd <- dd[match(tr$tip.label, rownames(dd)), , drop = FALSE]
  rownames(dd) <- tr$tip.label
  if (nrow(dd) < 20 || length(unique(dd$among)) < 2) next
  formula_text <- paste("among ~", paste(c(paste0(metric, "_z"), paste0(base_covariates, "_z")), collapse = " + "))
  form <- as.formula(formula_text)
  fit <- phyloglm(form, data = dd, phy = tr, method = "logistic_MPLE", btol = 50)
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  names(cf)[1:4] <- c("estimate", "std_error", "z_value", "p_value")
  cf$metric <- metric
  coefficients[[metric]] <- cf
  focal <- cf[cf$term == paste0(metric, "_z"), , drop = FALSE]
  nonphy <- glm(form, data = dd, family = binomial())
  nonphy_cf <- summary(nonphy)$coefficients[paste0(metric, "_z"), ]
  results[[metric]] <- data.frame(
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
  comparison[[metric]] <- data.frame(
    metric = metric,
    n_species = nrow(dd),
    phylogenetic_estimate = focal$estimate,
    phylogenetic_p = focal$p_value,
    nonphylogenetic_same_subset_estimate = nonphy_cf[[1]],
    nonphylogenetic_same_subset_p = nonphy_cf[[4]],
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, results)
if (is.null(result) || nrow(result) != length(components)) stop("Not all eight phylogenetic models were estimable")
result$p_holm_eight_components <- p.adjust(result$p_value, method = "holm")
coef_table <- do.call(rbind, coefficients)
comparison_table <- do.call(rbind, comparison)
write.csv(model_data, file.path(outdir, "phylogenetic_model_dataset.csv"), row.names = FALSE)
write.csv(result, file.path(outdir, "phylogenetic_logistic_component_models.csv"), row.names = FALSE)
write.csv(coef_table, file.path(outdir, "phylogenetic_logistic_all_coefficients.csv"), row.names = FALSE)
write.csv(comparison_table, file.path(outdir, "phylogenetic_vs_nonphylogenetic_same_subset.csv"), row.names = FALSE)

manifest <- list(
  status = "complete",
  input_species = nrow(d),
  opentree_matched_species = nrow(matched),
  tree_tip_species = length(resolved_tree$tip.label),
  model_species_min = min(result$n_species),
  model_species_max = max(result$n_species),
  tree_source = "Open Tree of Life induced synthetic subtree",
  branch_length_method = "Grafen power=1; topology-aware sensitivity, not divergence-time calibration",
  polytomy_resolution = "ape::multi2di(random=FALSE)",
  estimator = "phylolm::phyloglm(method='logistic_MPLE')",
  covariates = base_covariates,
  results = result,
  interpretation_guard = "Automated spatial labels remain exploratory. The phylogenetic model tests whether climate-breadth associations persist after topology-based dependence and the same range/effort covariates; it is not a causal or time-calibrated evolutionary model."
)
write_json(manifest, file.path(outdir, "phylogenetic_analysis_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
print(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"))
