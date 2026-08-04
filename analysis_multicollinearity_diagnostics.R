#!/usr/bin/env Rscript

# Multicollinearity diagnostics for the automated 107-species ecological models.
# Separates (1) redundancy among the eight climate-breadth metrics from
# (2) collinearity in each actually fitted one-climate-component model.

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
if (nrow(d) < 20) stop("Insufficient model data")

z_log <- function(x) {
  y <- log1p(pmax(as.numeric(x), 0))
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(y)))
  as.numeric((y - mean(y, na.rm = TRUE)) / s)
}

all_predictors <- c(components, base_covariates)
for (v in all_predictors) d[[paste0(v, "_z")]] <- z_log(d[[v]])
zvars <- paste0(all_predictors, "_z")

pairwise_long <- function(mat, method) {
  cm <- cor(mat, use = "pairwise.complete.obs", method = method)
  out <- expand.grid(variable_1 = rownames(cm), variable_2 = colnames(cm), stringsAsFactors = FALSE)
  out$correlation <- as.vector(cm)
  out$method <- method
  out <- out[out$variable_1 < out$variable_2, , drop = FALSE]
  out$abs_correlation <- abs(out$correlation)
  out[order(-out$abs_correlation), , drop = FALSE]
}

climate_mat <- d[, paste0(components, "_z"), drop = FALSE]
all_mat <- d[, zvars, drop = FALSE]
pearson_climate <- pairwise_long(climate_mat, "pearson")
spearman_climate <- pairwise_long(climate_mat, "spearman")
pearson_all <- pairwise_long(all_mat, "pearson")
spearman_all <- pairwise_long(all_mat, "spearman")
write.csv(pearson_climate, file.path(outdir, "climate_components_pearson_correlations.csv"), row.names = FALSE)
write.csv(spearman_climate, file.path(outdir, "climate_components_spearman_correlations.csv"), row.names = FALSE)
write.csv(pearson_all, file.path(outdir, "all_predictors_pearson_correlations.csv"), row.names = FALSE)
write.csv(spearman_all, file.path(outdir, "all_predictors_spearman_correlations.csv"), row.names = FALSE)

compute_vif <- function(dat, predictors, model_id) {
  cc <- complete.cases(dat[, predictors, drop = FALSE])
  x <- dat[cc, predictors, drop = FALSE]
  if (nrow(x) < length(predictors) + 3) {
    stop(sprintf("Too few complete rows for %s", model_id))
  }
  rows <- lapply(predictors, function(target) {
    others <- setdiff(predictors, target)
    if (!length(others)) return(data.frame())
    fit <- lm(reformulate(others, response = target), data = x)
    r2 <- summary(fit)$r.squared
    vif <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
    data.frame(
      model_id = model_id,
      n_complete = nrow(x),
      predictor = target,
      r_squared_from_other_predictors = r2,
      vif = vif,
      tolerance = if (is.finite(vif) && vif > 0) 1 / vif else 0,
      flag_vif_gt_3 = is.finite(vif) && vif > 3,
      flag_vif_gt_5 = is.finite(vif) && vif > 5,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

condition_diagnostic <- function(dat, predictors, model_id) {
  cc <- complete.cases(dat[, predictors, drop = FALSE])
  x <- as.matrix(dat[cc, predictors, drop = FALSE])
  # Predictors are already standardized. Include intercept for the fitted design.
  design <- cbind(`(Intercept)` = 1, x)
  singular_values <- svd(design, nu = 0, nv = 0)$d
  condition_number <- max(singular_values) / min(singular_values)
  data.frame(
    model_id = model_id,
    n_complete = nrow(x),
    n_predictors = ncol(x),
    condition_number = condition_number,
    flag_condition_gt_10 = is.finite(condition_number) && condition_number > 10,
    flag_condition_gt_30 = is.finite(condition_number) && condition_number > 30,
    stringsAsFactors = FALSE
  )
}

# Diagnostics for the eight models actually fitted: one BIO component at a time.
actual_vif <- list()
actual_condition <- list()
for (metric in components) {
  predictors <- c(paste0(metric, "_z"), paste0(base_covariates, "_z"))
  actual_vif[[metric]] <- compute_vif(d, predictors, paste0("actual_", metric))
  actual_condition[[metric]] <- condition_diagnostic(d, predictors, paste0("actual_", metric))
}
actual_vif <- do.call(rbind, actual_vif)
actual_condition <- do.call(rbind, actual_condition)
write.csv(actual_vif, file.path(outdir, "actual_component_models_vif.csv"), row.names = FALSE)
write.csv(actual_condition, file.path(outdir, "actual_component_models_condition_numbers.csv"), row.names = FALSE)

# Diagnostic only: all eight climate metrics entered together. This is not the
# inferential model and is reported solely to quantify redundancy among metrics.
joint_predictors <- c(paste0(components, "_z"), paste0(base_covariates, "_z"))
joint_vif <- compute_vif(d, joint_predictors, "diagnostic_all_eight_components_joint")
joint_condition <- condition_diagnostic(d, joint_predictors, "diagnostic_all_eight_components_joint")
write.csv(joint_vif, file.path(outdir, "diagnostic_all_eight_components_joint_vif.csv"), row.names = FALSE)
write.csv(joint_condition, file.path(outdir, "diagnostic_all_eight_components_joint_condition_number.csv"), row.names = FALSE)

summary_table <- data.frame(
  diagnostic = c(
    "max_abs_pearson_among_eight_climate_components",
    "max_abs_spearman_among_eight_climate_components",
    "max_vif_across_actual_one_component_models",
    "max_condition_number_across_actual_one_component_models",
    "max_vif_if_all_eight_components_entered_jointly",
    "condition_number_if_all_eight_components_entered_jointly"
  ),
  value = c(
    max(pearson_climate$abs_correlation, na.rm = TRUE),
    max(spearman_climate$abs_correlation, na.rm = TRUE),
    max(actual_vif$vif, na.rm = TRUE),
    max(actual_condition$condition_number, na.rm = TRUE),
    max(joint_vif$vif, na.rm = TRUE),
    joint_condition$condition_number
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(outdir, "multicollinearity_diagnostic_summary.csv"), row.names = FALSE)

cat("Multicollinearity diagnostics complete\n")
print(summary_table)
cat("\nLargest VIFs in actually fitted models:\n")
print(head(actual_vif[order(-actual_vif$vif), ], 12))
