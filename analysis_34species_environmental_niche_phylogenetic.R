#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phylolm)
})

METRICS <- c(
  "temperature_breadth",
  "moisture_breadth",
  "climatic_heterogeneity",
  "pca_dispersion",
  "pca_hull_area"
)

parse_args <- function(args) {
  out <- list(dataset=NULL, outdir=NULL, min_cells=20L, replicates=100L, seed=20260724L)
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for ", args[[i]])
    if (!key %in% names(out)) stop("Unknown argument --", key)
    value <- args[[i+1L]]
    if (key %in% c("min_cells","replicates","seed")) value <- as.integer(value)
    out[[key]] <- value
    i <- i + 2L
  }
  if (is.null(out$dataset) || is.null(out$outdir)) stop("--dataset and --outdir are required")
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm=TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm=TRUE)) / s
}

normalise_edges <- function(tree, grafen=FALSE, seed=NULL) {
  tr <- tree
  if (grafen) {
    if (!is.null(seed)) set.seed(seed)
    tr <- if (is.binary(tr)) tr else multi2di(tr, random=TRUE)
    tr$edge.length <- NULL
    tr <- compute.brlen(tr, method="Grafen", power=1)
  }
  if (is.null(tr$edge.length)) stop("Tree has no branch lengths")
  tr$edge.length <- as.numeric(tr$edge.length)
  pos <- tr$edge.length[is.finite(tr$edge.length) & tr$edge.length > 0]
  if (!length(pos)) stop("Tree has no positive branch lengths")
  bad <- !is.finite(tr$edge.length) | tr$edge.length <= 0
  replaced <- sum(bad)
  if (replaced) tr$edge.length[bad] <- min(pos) * 1e-8
  list(tree=tr, replaced=replaced)
}

prepare_data <- function(raw, metric, min_cells) {
  req <- c("canonical_name","family","spatial_scale","classification_source","n_climate_cells",metric)
  miss <- setdiff(req, names(raw))
  if (length(miss)) stop("Missing columns for ", metric, ": ", paste(miss, collapse=", "))
  d <- raw[
    raw$classification_source == "baseline_unambiguous" &
      raw$n_climate_cells >= min_cells &
      raw$spatial_scale %in% c("within_population","among_population"),
    req, drop=FALSE
  ]
  d <- d[!duplicated(d$canonical_name),,drop=FALSE]
  d$among <- as.integer(d$spatial_scale == "among_population")
  d$metric_z <- zscore(d[[metric]])
  d$effort_z <- zscore(log1p(d$n_climate_cells))
  d$tip_label <- gsub(" ", "_", d$canonical_name, fixed=TRUE)
  d <- d[complete.cases(d[,c("among","metric_z","effort_z")]),,drop=FALSE]
  rownames(d) <- d$tip_label
  d
}

fit_tree <- function(tree, d, metric, tree_type, scenario, replicate=NA_integer_, seed=NA_integer_, branch_replaced=0L) {
  keep <- intersect(tree$tip.label, rownames(d))
  tr <- if (length(setdiff(tree$tip.label, keep))) drop.tip(tree, setdiff(tree$tip.label, keep)) else tree
  dd <- d[tr$tip.label,,drop=FALSE]
  if (nrow(dd) < 20L || length(unique(dd$among)) < 2L) {
    return(data.frame(metric=metric, tree_type=tree_type, scenario=scenario, replicate=replicate,
      seed=seed, n_species=nrow(dd), n_within=sum(dd$among==0), n_among=sum(dd$among==1),
      branch_lengths_replaced=branch_replaced, estimate=NA_real_, std_error=NA_real_, odds_ratio=NA_real_,
      ci_low=NA_real_, ci_high=NA_real_, p_value=NA_real_, alpha=NA_real_, status="not_estimable"))
  }
  fit <- try(phyloglm(among ~ metric_z + effort_z, data=dd, phy=tr, method="logistic_MPLE", boot=0), silent=TRUE)
  if (inherits(fit,"try-error")) {
    return(data.frame(metric=metric, tree_type=tree_type, scenario=scenario, replicate=replicate,
      seed=seed, n_species=nrow(dd), n_within=sum(dd$among==0), n_among=sum(dd$among==1),
      branch_lengths_replaced=branch_replaced, estimate=NA_real_, std_error=NA_real_, odds_ratio=NA_real_,
      ci_low=NA_real_, ci_high=NA_real_, p_value=NA_real_, alpha=NA_real_, status=paste0("failed: ",substr(as.character(fit),1,160))))
  }
  beta <- unname(coef(fit)[["metric_z"]])
  se <- sqrt(diag(vcov(fit)))[["metric_z"]]
  z <- beta/se
  p <- if (is.finite(z)) 2*pnorm(-abs(z)) else NA_real_
  data.frame(metric=metric, tree_type=tree_type, scenario=scenario, replicate=replicate,
    seed=seed, n_species=nrow(dd), n_within=sum(dd$among==0), n_among=sum(dd$among==1),
    branch_lengths_replaced=branch_replaced, estimate=beta, std_error=se, odds_ratio=exp(beta),
    ci_low=exp(beta-1.96*se), ci_high=exp(beta+1.96*se), p_value=p,
    alpha=if (!is.null(fit$alpha)) unname(fit$alpha) else NA_real_, status="complete")
}

args <- parse_args(commandArgs(trailingOnly=TRUE))
dir.create(args$outdir, recursive=TRUE, showWarnings=FALSE)
raw <- read.csv(args$dataset, stringsAsFactors=FALSE, check.names=FALSE)

# Collinearity diagnostics for the exact predictor pair used in each model.
coll <- do.call(rbind, lapply(METRICS, function(metric) {
  d <- prepare_data(raw, metric, args$min_cells)
  if (nrow(d) != 34L) stop(metric, ": expected 34 frozen species, found ", nrow(d))
  r <- cor(d$metric_z, d$effort_z)
  vif <- 1/(1-r^2)
  X <- cbind(Intercept=1, metric_z=d$metric_z, effort_z=d$effort_z)
  data.frame(metric=metric, n_species=nrow(d), predictor_correlation=r,
    vif_metric=vif, vif_effort=vif, max_vif=vif,
    condition_number=kappa(X, exact=TRUE), stringsAsFactors=FALSE)
}))
write.csv(coll, file.path(args$outdir,"environmental_niche_collinearity.csv"), row.names=FALSE)

# Open Tree topology: same stored topology as the manuscript, 100 random binary resolutions + Grafen lengths.
open_base <- read.tree("docs/supporting/jbi_opentree_induced_topology.tre")
open_rows <- list()
k <- 1L
for (metric in METRICS) {
  d <- prepare_data(raw, metric, args$min_cells)
  for (i in seq_len(args$replicates)) {
    norm <- normalise_edges(open_base, grafen=TRUE, seed=args$seed+i-1L)
    open_rows[[k]] <- fit_tree(norm$tree, d, metric, "OpenTree_Grafen", "random_polytomy_resolution",
      i, args$seed+i-1L, norm$replaced)
    k <- k+1L
  }
}
open_rep <- do.call(rbind, open_rows)
write.csv(open_rep, file.path(args$outdir,"environmental_niche_opentree_replicates.csv"), row.names=FALSE)
open_complete <- open_rep[open_rep$status=="complete",,drop=FALSE]
open_summary <- do.call(rbind, lapply(METRICS, function(metric) {
  x <- open_complete[open_complete$metric==metric,,drop=FALSE]
  data.frame(metric=metric, tree_type="OpenTree_Grafen", n_species=if(nrow(x)) median(x$n_species) else NA,
    n_replicates_complete=nrow(x), median_odds_ratio=median(x$odds_ratio,na.rm=TRUE),
    min_odds_ratio=min(x$odds_ratio,na.rm=TRUE), max_odds_ratio=max(x$odds_ratio,na.rm=TRUE),
    median_ci_low=median(x$ci_low,na.rm=TRUE), median_ci_high=median(x$ci_high,na.rm=TRUE),
    median_p_value=median(x$p_value,na.rm=TRUE), fraction_negative=mean(x$estimate<0,na.rm=TRUE),
    median_alpha=median(x$alpha,na.rm=TRUE), stringsAsFactors=FALSE)
}))
open_summary$holm_median_p_five_metrics <- p.adjust(open_summary$median_p_value, method="holm")
write.csv(open_summary, file.path(args$outdir,"environmental_niche_opentree_summary.csv"), row.names=FALSE)

# Re-use the archived V.PhyloMaker2 S1-S3 time-scaled trees already generated for these 34 species.
dated_files <- c(
  S1="docs/supporting/jbi_dated_phylogeny_s1.tre",
  S2="docs/supporting/jbi_dated_phylogeny_s2.tre",
  S3="docs/supporting/jbi_dated_phylogeny_s3.tre"
)
dated_rows <- list(); k <- 1L
for (scenario in names(dated_files)) {
  norm <- normalise_edges(read.tree(dated_files[[scenario]]), grafen=FALSE)
  for (metric in METRICS) {
    d <- prepare_data(raw, metric, args$min_cells)
    dated_rows[[k]] <- fit_tree(norm$tree, d, metric, "VPhyloMaker2_dated", scenario,
      NA_integer_, args$seed, norm$replaced)
    k <- k+1L
  }
}
dated <- do.call(rbind, dated_rows)
for (sc in unique(dated$scenario)) {
  idx <- dated$scenario==sc & dated$status=="complete"
  dated$holm_p_five_metrics[idx] <- p.adjust(dated$p_value[idx], method="holm")
}
write.csv(dated, file.path(args$outdir,"environmental_niche_dated_phyloglm.csv"), row.names=FALSE)

manifest <- list(
  status="complete",
  dataset_role="frozen_34_species_baseline_only",
  n_species=34L,
  min_cells=args$min_cells,
  metrics=METRICS,
  model_formula="among ~ metric_z + effort_z",
  collinearity=list(method="VIF from the two-predictor correlation and exact condition number of the standardized design matrix", results=split(coll,seq_len(nrow(coll)))),
  opentree=list(topology="docs/supporting/jbi_opentree_induced_topology.tre", replicates=args$replicates,
    branch_lengths="Grafen", results=split(open_summary,seq_len(nrow(open_summary)))),
  dated_phylogeny=list(source="stored V.PhyloMaker2 GBOTB.extended.LCVP trees from manuscript sensitivity analysis",
    scenarios=c("S1","S2","S3"), results=split(dated,seq_len(nrow(dated)))),
  interpretation_guard="All five climatic niche metrics are analysed symmetrically. Phylogenetic models are sensitivity analyses; occupied climate is not physiological tolerance and associations are not causal.",
  package_versions=list(R=R.version.string, ape=as.character(packageVersion("ape")), phylolm=as.character(packageVersion("phylolm")), jsonlite=as.character(packageVersion("jsonlite")))
)
writeLines(toJSON(manifest, pretty=TRUE, auto_unbox=TRUE, na="null"), file.path(args$outdir,"environmental_niche_phylogenetic_manifest.json"))
print(coll)
print(open_summary)
print(dated)
