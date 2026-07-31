#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phylolm)
})

parse_args <- function(args) {
  out <- list(
    primary = NULL,
    paginated = NULL,
    opentree = NULL,
    dated_s1 = NULL,
    dated_s2 = NULL,
    dated_s3 = NULL,
    outdir = NULL,
    replicates = 100L,
    seed = 20260801L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% names(out)) stop("Unknown argument: --", key)
    if (i == length(args)) stop("Missing value for --", key)
    value <- args[[i + 1L]]
    if (key %in% c("replicates", "seed")) value <- as.integer(value)
    out[[key]] <- value
    i <- i + 2L
  }
  required <- c("primary", "paginated", "opentree", "dated_s1", "dated_s2", "dated_s3", "outdir")
  missing <- required[vapply(required, function(x) is.null(out[[x]]), logical(1))]
  if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

normalise_tree <- function(tree, grafen = FALSE, seed = NULL) {
  tr <- tree
  if (grafen) {
    if (!is.null(seed)) set.seed(seed)
    tr <- if (is.binary(tr)) tr else multi2di(tr, random = TRUE)
    tr$edge.length <- NULL
    tr <- compute.brlen(tr, method = "Grafen", power = 1)
  }
  if (is.null(tr$edge.length)) stop("Tree has no branch lengths")
  positive <- tr$edge.length[is.finite(tr$edge.length) & tr$edge.length > 0]
  if (!length(positive)) stop("Tree has no positive branch lengths")
  tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length <= 0] <- min(positive) * 1e-8
  tr
}

prepare_data <- function(path, metric) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("canonical_name", "spatial_scale", "n_climate_cells", metric)
  missing <- setdiff(required, names(d))
  if (length(missing)) stop(path, " missing: ", paste(missing, collapse = ", "))
  d <- d[d$spatial_scale %in% c("within_population", "among_population"), required, drop = FALSE]
  d <- d[!duplicated(d$canonical_name), , drop = FALSE]
  d$among <- as.integer(d$spatial_scale == "among_population")
  d$metric_z <- zscore(d[[metric]])
  d$effort_z <- zscore(log1p(d$n_climate_cells))
  d$tip_label <- gsub(" ", "_", d$canonical_name, fixed = TRUE)
  d <- d[complete.cases(d[, c("among", "metric_z", "effort_z")]), , drop = FALSE]
  rownames(d) <- d$tip_label
  d
}

fit_one <- function(data, tree, occurrence_dataset, metric, model_type, replicate = NA_integer_, scenario = NA_character_) {
  keep <- intersect(tree$tip.label, rownames(data))
  drop <- setdiff(tree$tip.label, keep)
  tr <- if (length(drop)) drop.tip(tree, drop) else tree
  d <- data[tr$tip.label, , drop = FALSE]
  if (nrow(d) < 20L || length(unique(d$among)) < 2L) stop("Not estimable: ", occurrence_dataset, " / ", metric)
  fit <- phyloglm(
    among ~ metric_z + effort_z,
    data = d,
    phy = tr,
    method = "logistic_MPLE",
    boot = 0
  )
  beta <- unname(coef(fit)[["metric_z"]])
  se <- sqrt(diag(vcov(fit)))[["metric_z"]]
  z <- beta / se
  data.frame(
    occurrence_dataset = occurrence_dataset,
    metric = metric,
    model_type = model_type,
    replicate = replicate,
    scenario = scenario,
    n_species = nrow(d),
    n_within = sum(d$among == 0L),
    n_among = sum(d$among == 1L),
    estimate = beta,
    std_error = se,
    odds_ratio = exp(beta),
    odds_ratio_ci_low = exp(beta - 1.96 * se),
    odds_ratio_ci_high = exp(beta + 1.96 * se),
    p_value = 2 * pnorm(-abs(z)),
    alpha = if (!is.null(fit$alpha)) unname(fit$alpha) else NA_real_,
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
metrics <- c("bio12_q95q05", "bio14_q95q05", "bio15_q95q05", "bio17_q95q05")
datasets <- list(primary = args$primary, paginated_quality_filtered = args$paginated)
open_base <- read.tree(args$opentree)
dated <- list(S1 = read.tree(args$dated_s1), S2 = read.tree(args$dated_s2), S3 = read.tree(args$dated_s3))

rows <- list()
k <- 1L
for (dataset_name in names(datasets)) {
  for (metric_index in seq_along(metrics)) {
    metric <- metrics[[metric_index]]
    d <- prepare_data(datasets[[dataset_name]], metric)
    for (r in seq_len(args$replicates)) {
      tr <- normalise_tree(open_base, grafen = TRUE, seed = args$seed + metric_index * 1000L + r)
      rows[[k]] <- fit_one(d, tr, dataset_name, metric, "OpenTree_Grafen", replicate = r)
      k <- k + 1L
    }
    for (scenario in names(dated)) {
      tr <- normalise_tree(dated[[scenario]], grafen = FALSE)
      rows[[k]] <- fit_one(d, tr, dataset_name, metric, "dated_megaphylogeny", scenario = scenario)
      k <- k + 1L
    }
  }
}

fits <- do.call(rbind, rows)
write.csv(fits, file.path(args$outdir, "precipitation_component_phylogenetic_fits.csv"), row.names = FALSE)

summary_rows <- list()
k <- 1L
for (dataset_name in unique(fits$occurrence_dataset)) {
  for (metric in unique(fits$metric)) {
    for (model_type in unique(fits$model_type)) {
      x <- fits[
        fits$occurrence_dataset == dataset_name &
          fits$metric == metric &
          fits$model_type == model_type,
        , drop = FALSE
      ]
      if (!nrow(x)) next
      summary_rows[[k]] <- data.frame(
        occurrence_dataset = dataset_name,
        metric = metric,
        model_type = model_type,
        n_fits = nrow(x),
        n_species_min = min(x$n_species),
        n_species_max = max(x$n_species),
        odds_ratio_min = min(x$odds_ratio),
        odds_ratio_median = median(x$odds_ratio),
        odds_ratio_max = max(x$odds_ratio),
        ci_low_min = min(x$odds_ratio_ci_low),
        ci_high_max = max(x$odds_ratio_ci_high),
        p_value_min = min(x$p_value),
        p_value_max = max(x$p_value),
        negative_fraction = mean(x$estimate < 0),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}
summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(args$outdir, "precipitation_component_phylogenetic_summary.csv"), row.names = FALSE)

manifest <- list(
  status = "complete",
  metrics = metrics,
  occurrence_datasets = names(datasets),
  opentree_replicates = args$replicates,
  seed = args$seed,
  model_formula = "among ~ separately standardised precipitation component + effort_z",
  estimator = "phylolm::phyloglm(method = logistic_MPLE)",
  dated_trees = c("S1", "S2", "S3"),
  results = unname(split(summary, seq_len(nrow(summary)))),
  interpretation_guard = paste(
    "Component-wise phylogenetic models are sensitivity analyses.",
    "They do not establish a causal precipitation mechanism or morph-specific climatic tolerance."
  )
)
writeLines(
  toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"),
  file.path(args$outdir, "precipitation_component_phylogenetic_manifest.json")
)
print(summary)
