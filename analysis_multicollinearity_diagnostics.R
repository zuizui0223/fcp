#!/usr/bin/env Rscript

# Diagnose collinearity in the exploratory expanded-set range-adjusted models.
# This script never selects or drops predictors; it reports VIF and condition numbers.

suppressPackageStartupMessages(library(jsonlite))

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
base <- c(
  "geographic_hull_area_km2",
  "occupied_geographic_grid_cells",
  "n_occurrence_records",
  "n_supporting_p1_records"
)

d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(components, base)
missing <- setdiff(required, names(d))
if (length(missing)) stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))

zlog <- function(x) {
  y <- log1p(pmax(as.numeric(x), 0))
  s <- sd(y, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(y)))
  as.numeric((y - mean(y, na.rm = TRUE)) / s)
}

# Match the fitted downstream model: every breadth/range/effort predictor is log1p transformed
# and standardized within the analysis dataset.
for (v in c(components, base)) d[[paste0(v, "_z")]] <- zlog(d[[v]])

vif_for <- function(frame, predictors) {
  out <- lapply(predictors, function(term) {
    others <- setdiff(predictors, term)
    if (!length(others)) return(data.frame(term = term, vif = 1))
    fit <- lm(reformulate(others, response = term), data = frame)
    r2 <- summary(fit)$r.squared
    value <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
    data.frame(term = term, vif = value)
  })
  do.call(rbind, out)
}

vif_rows <- list()
condition_rows <- list()
for (metric in components) {
  predictors <- c(
    paste0(metric, "_z"),
    "geographic_hull_area_km2_z",
    "occupied_geographic_grid_cells_z",
    "n_occurrence_records_z",
    "n_supporting_p1_records_z"
  )
  dd <- d[complete.cases(d[, predictors, drop = FALSE]), predictors, drop = FALSE]
  if (nrow(dd) < 20) next
  vv <- vif_for(dd, predictors)
  vv$metric <- metric
  vv$n_species <- nrow(dd)
  vif_rows[[metric]] <- vv[, c("metric", "n_species", "term", "vif")]

  X <- cbind(`(Intercept)` = 1, as.matrix(dd))
  Xs <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))
  kappa_value <- kappa(Xs, exact = TRUE)
  condition_rows[[metric]] <- data.frame(
    metric = metric,
    n_species = nrow(dd),
    n_predictors = length(predictors),
    max_vif = max(vv$vif, na.rm = TRUE),
    condition_number = as.numeric(kappa_value),
    stringsAsFactors = FALSE
  )
}

vif_table <- do.call(rbind, vif_rows)
condition_table <- do.call(rbind, condition_rows)
if (is.null(vif_table) || nrow(condition_table) != length(components)) {
  stop("Not all eight component models had estimable collinearity diagnostics")
}
write.csv(vif_table, file.path(outdir, "multicollinearity_vif.csv"), row.names = FALSE)
write.csv(condition_table, file.path(outdir, "multicollinearity_condition_numbers.csv"), row.names = FALSE)

manifest <- list(
  status = "complete",
  n_models = nrow(condition_table),
  transformation = "log1p followed by z standardization, matching the downstream component models",
  selection_rule = "diagnostic only; no VIF-based variable deletion",
  vif_threshold_flag = 5,
  models_with_max_vif_gt_5 = condition_table$metric[condition_table$max_vif > 5],
  maximum_vif = max(condition_table$max_vif, na.rm = TRUE),
  maximum_condition_number = max(condition_table$condition_number, na.rm = TRUE),
  results = condition_table
)
write_json(manifest, file.path(outdir, "multicollinearity_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
print(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"))
