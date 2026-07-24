#!/usr/bin/env Rscript

# Dated-megaphylogeny sensitivity analysis for the JBI manuscript.
#
# The script uses the time-scaled GBOTB.extended.LCVP backbone distributed
# with V.PhyloMaker2. Species absent from the backbone are inserted under the
# package's S1, S2 and S3 placement scenarios. Each resulting dated tree is
# analysed with phyloglm(logistic_MPLE). Scenario comparison diagnoses
# uncertainty from placement of taxa not already represented in the backbone.

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phylolm)
  library(V.PhyloMaker2)
})

parse_args <- function(args) {
  out <- list(
    dataset = NULL,
    outdir = NULL,
    analysis_label = NULL,
    min_cells = 20L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    if (!key %in% names(out)) stop("Unknown argument: --", key)
    value <- args[[i + 1L]]
    if (key == "min_cells") value <- as.integer(value)
    out[[key]] <- value
    i <- i + 2L
  }
  if (is.null(out$dataset) || is.null(out$outdir) || is.null(out$analysis_label)) {
    stop("--dataset, --outdir and --analysis_label are required")
  }
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_package_field <- function(package, field) {
  desc <- packageDescription(package)
  value <- desc[[field]]
  if (is.null(value) || length(value) == 0L) NA_character_ else as.character(value)
}

normalise_tree <- function(tree) {
  if (!inherits(tree, "phylo")) stop("V.PhyloMaker2 did not return a phylo object")
  if (is.null(tree$edge.length)) stop("Dated tree has no branch lengths")
  tree$edge.length <- as.numeric(tree$edge.length)
  positive <- tree$edge.length[is.finite(tree$edge.length) & tree$edge.length > 0]
  if (!length(positive)) stop("Dated tree has no positive branch lengths")
  replaced <- sum(!is.finite(tree$edge.length) | tree$edge.length <= 0)
  if (replaced > 0L) {
    tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length <= 0] <-
      min(positive) * 1e-8
  }
  list(tree = tree, branch_lengths_replaced = replaced)
}

fit_tree <- function(tree, data, scenario, analysis_label, placement_counts) {
  norm <- normalise_tree(tree)
  tr <- norm$tree
  keep <- intersect(tr$tip.label, rownames(data))
  if (length(keep) < 20L) {
    stop("Fewer than 20 model species remained in ", scenario, ": ", length(keep))
  }
  drop <- setdiff(tr$tip.label, keep)
  if (length(drop)) tr <- drop.tip(tr, drop)
  d <- data[tr$tip.label, , drop = FALSE]
  if (length(unique(d$among)) < 2L) stop("Only one response class in ", scenario)

  fit <- phyloglm(
    among ~ moisture_z + effort_z,
    data = d,
    phy = tr,
    method = "logistic_MPLE",
    boot = 0
  )
  beta <- unname(coef(fit)[["moisture_z"]])
  covariance <- vcov(fit)
  se <- sqrt(diag(covariance))[["moisture_z"]]
  z <- beta / se
  p <- 2 * pnorm(-abs(z))
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  root_age <- max(node.depth.edgelength(tr)[seq_len(length(tr$tip.label))])

  data.frame(
    analysis_label = analysis_label,
    backbone = "GBOTB.extended.LCVP",
    scenario = scenario,
    estimator = "phylolm::phyloglm logistic_MPLE",
    n_species = nrow(d),
    n_within = sum(d$among == 0L),
    n_among = sum(d$among == 1L),
    species_present_in_backbone = placement_counts[["prune"]],
    species_bound_to_backbone = placement_counts[["bind"]],
    species_failed_to_bind = placement_counts[["fail to bind"]],
    branch_lengths_replaced = norm$branch_lengths_replaced,
    root_age_backbone_units = root_age,
    estimate = beta,
    std_error = se,
    odds_ratio = exp(beta),
    odds_ratio_ci_low = exp(ci_low),
    odds_ratio_ci_high = exp(ci_high),
    p_value = p,
    alpha = if (!is.null(fit$alpha)) unname(fit$alpha) else NA_real_,
    log_likelihood = as.numeric(logLik(fit)),
    status = "complete",
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
data$tip_label <- gsub(" ", "_", data$canonical_name, fixed = TRUE)
data <- data[complete.cases(data[, c("among", "moisture_z", "effort_z")]), , drop = FALSE]
if (nrow(data) < 20L) stop("Fewer than 20 complete species before tree construction")
rownames(data) <- data$tip_label

parts <- strsplit(data$canonical_name, " +")
sp_list <- data.frame(
  species = data$canonical_name,
  genus = vapply(parts, `[[`, character(1), 1L),
  family = data$family,
  species.relative = "",
  genus.relative = "",
  stringsAsFactors = FALSE
)
write.csv(sp_list, file.path(args$outdir, "vphylomaker2_input_species.csv"), row.names = FALSE)

# Explicitly load the LCVP time-scaled backbone and its node information.
data("GBOTB.extended.LCVP", package = "V.PhyloMaker2", envir = environment())
data("nodes.info.1.LCVP", package = "V.PhyloMaker2", envir = environment())
if (!exists("GBOTB.extended.LCVP") || !exists("nodes.info.1.LCVP")) {
  stop("V.PhyloMaker2 LCVP backbone objects were not loaded")
}

made <- phylo.maker(
  sp.list = sp_list,
  tree = GBOTB.extended.LCVP,
  nodes = nodes.info.1.LCVP,
  output.sp.list = TRUE,
  output.tree = FALSE,
  scenarios = c("S1", "S2", "S3")
)
if (is.null(made$species.list)) stop("V.PhyloMaker2 did not return a placement table")
placement <- made$species.list
write.csv(
  placement,
  file.path(args$outdir, "vphylomaker2_species_placement.csv"),
  row.names = FALSE
)
status_values <- as.character(placement$status)
status_levels <- c("prune", "bind", "fail to bind")
placement_counts <- setNames(
  vapply(status_levels, function(x) sum(status_values == x, na.rm = TRUE), integer(1)),
  status_levels
)

scenario_objects <- list(
  S1 = made$scenario.1,
  S2 = made$scenario.2,
  S3 = made$scenario.3
)
summary_rows <- list()
for (scenario in names(scenario_objects)) {
  tree <- scenario_objects[[scenario]]
  if (is.null(tree)) stop("Missing V.PhyloMaker2 tree for ", scenario)
  write.tree(tree, file.path(args$outdir, paste0("vphylomaker2_", tolower(scenario), ".tre")))
  summary_rows[[scenario]] <- fit_tree(
    tree, data, scenario, args$analysis_label, placement_counts
  )
}
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL
write.csv(
  summary,
  file.path(args$outdir, "vphylomaker2_dated_phyloglm_summary.csv"),
  row.names = FALSE
)

manifest <- list(
  analysis_label = args$analysis_label,
  dataset = args$dataset,
  min_cells = args$min_cells,
  n_input_species = nrow(data),
  spatial_scale_counts = as.list(table(data$spatial_scale)),
  backbone = "GBOTB.extended.LCVP",
  node_table = "nodes.info.1.LCVP",
  scenarios = names(scenario_objects),
  placement_status_counts = as.list(placement_counts),
  model_formula = "among ~ moisture_z + effort_z",
  estimator = "phylolm::phyloglm(method = logistic_MPLE)",
  confidence_interval = "Wald 95% interval: coefficient +/- 1.96 standard errors",
  branch_length_guard = "Non-finite or non-positive edges, if present, replaced by 1e-8 times the smallest positive edge and counted in the summary.",
  interpretation_guard = paste(
    "V.PhyloMaker2 trees inherit time scaling from the GBOTB.extended.LCVP backbone,",
    "but placements of taxa absent from the backbone depend on scenarios S1-S3.",
    "The analysis is a dated-megaphylogeny sensitivity analysis and does not establish causation."
  ),
  package_versions = list(
    R = R.version.string,
    ape = as.character(packageVersion("ape")),
    jsonlite = as.character(packageVersion("jsonlite")),
    phylolm = as.character(packageVersion("phylolm")),
    V.PhyloMaker2 = as.character(packageVersion("V.PhyloMaker2")),
    V.PhyloMaker2_RemoteSha = safe_package_field("V.PhyloMaker2", "RemoteSha"),
    V.PhyloMaker2_RemoteRepo = safe_package_field("V.PhyloMaker2", "RemoteRepo"),
    V.PhyloMaker2_RemoteUsername = safe_package_field("V.PhyloMaker2", "RemoteUsername")
  ),
  results = unname(split(summary, seq_len(nrow(summary))))
)
writeLines(
  toJSON(manifest, pretty = TRUE, auto_unbox = TRUE, na = "null"),
  file.path(args$outdir, "vphylomaker2_dated_phyloglm_manifest.json")
)
print(summary)
