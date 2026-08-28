#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(rotl)
  library(phylolm)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing value for", flag))
  args[[i+1]]
}
dataset <- get_arg("--dataset")
outdir <- get_arg("--outdir")
replicates <- as.integer(get_arg("--replicates", "100"))
seed <- as.integer(get_arg("--seed", "20260826"))
if (is.null(dataset) || is.null(outdir)) stop("--dataset and --outdir are required")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

metrics <- c("temperature_breadth","moisture_breadth","climatic_heterogeneity","pca_dispersion","pca_hull_area")
outcomes <- c(C="C_local_coexistence_documented", S="S_spatial_segregation_documented")
d <- read.csv(dataset, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(d)==34)

# Resolve the replacement species independently against Open Tree taxonomy.
matches <- tnrs_match_names(d$canonical_name, context_name="Land plants", do_approximate_matching=TRUE)
write.csv(matches, file.path(outdir,"jbi_cs_opentree_tnrs_matches.csv"), row.names=FALSE)

# Keep one unambiguous non-approximate matched OTT id per input name where possible.
# rotl returns one best row per query in normal use; duplicated OTT ids are resolved below.
if (!all(c("search_string","ott_id") %in% names(matches))) stop("Unexpected rotl TNRS columns")
if (!"is_approximate_match" %in% names(matches)) matches$is_approximate_match <- FALSE
if (!"is_synonym" %in% names(matches)) matches$is_synonym <- FALSE
matches$search_norm <- trimws(as.character(matches$search_string))
matches$ott_id <- suppressWarnings(as.integer(matches$ott_id))
usable <- matches[!is.na(matches$ott_id) & !matches$is_approximate_match, , drop=FALSE]
# Preserve one input name per unique OTT id and flag collisions rather than merging species.
dup_ott <- usable$ott_id[duplicated(usable$ott_id) | duplicated(usable$ott_id, fromLast=TRUE)]
usable <- usable[!(usable$ott_id %in% dup_ott), , drop=FALSE]
usable <- usable[usable$search_norm %in% d$canonical_name, , drop=FALSE]

if (nrow(usable) < 20) stop(sprintf("Too few unambiguous OpenTree matches: %d", nrow(usable)))

# Request labels as OTT ids so mapping back to the input species is deterministic.
tr <- tol_induced_subtree(ott_ids=usable$ott_id, label_format="id")
write.tree(tr, file=file.path(outdir,"jbi_cs_opentree_induced_topology.tre"))

# Map tree tip ottNNNN labels back to canonical names.
extract_ott <- function(x) as.integer(sub("^ott", "", x))
tip_ott <- vapply(tr$tip.label, extract_ott, integer(1))
map <- setNames(usable$search_norm, as.character(usable$ott_id))
tip_names <- unname(map[as.character(tip_ott)])
if (any(is.na(tip_names))) stop("Could not map all induced-tree tips back to canonical names")
tr$tip.label <- gsub(" ", "_", tip_names, fixed=TRUE)

matched_species <- gsub("_", " ", tr$tip.label, fixed=TRUE)
md <- d[d$canonical_name %in% matched_species, , drop=FALSE]
md <- md[match(matched_species, md$canonical_name), , drop=FALSE]
rownames(md) <- tr$tip.label
if (any(is.na(md$canonical_name))) stop("Tree/data matching failed")

z <- function(x) as.numeric(scale(as.numeric(x), center=TRUE, scale=TRUE))
set.seed(seed)
raw <- list(); k <- 1L

for (r in seq_len(replicates)) {
  tree_r <- if (is.binary.tree(tr)) tr else multi2di(tr, random=TRUE)
  tree_r <- compute.brlen(tree_r, method="Grafen", power=1)
  for (short in names(outcomes)) {
    outcome_col <- outcomes[[short]]
    for (metric in metrics) {
      x <- md
      x$outcome <- as.integer(x[[outcome_col]])
      x$metric_z <- z(x[[metric]])
      x$effort_z <- z(log1p(as.numeric(x$n_climate_cells)))
      fit <- tryCatch(
        phyloglm(outcome ~ metric_z + effort_z, phy=tree_r, data=x, method="logistic_MPLE"),
        error=function(e) NULL
      )
      if (is.null(fit)) next
      cf <- summary(fit)$coefficients
      if (!("metric_z" %in% rownames(cf))) next
      b <- as.numeric(cf["metric_z","Estimate"])
      se_col <- grep("Std", colnames(cf), value=TRUE)[1]
      p_col <- grep("Pr", colnames(cf), value=TRUE)[1]
      se <- if (!is.na(se_col)) as.numeric(cf["metric_z",se_col]) else NA_real_
      p <- if (!is.na(p_col)) as.numeric(cf["metric_z",p_col]) else NA_real_
      raw[[k]] <- data.frame(
        replicate=r, outcome_short=short, outcome=outcome_col, metric=metric,
        n_species=nrow(x), n_positive=sum(x$outcome),
        estimate=b, odds_ratio=exp(b), std_error=se,
        ci_low=ifelse(is.finite(se),exp(b-1.96*se),NA_real_),
        ci_high=ifelse(is.finite(se),exp(b+1.96*se),NA_real_),
        p_value=p, alpha=fit$alpha, stringsAsFactors=FALSE
      )
      k <- k + 1L
    }
  }
}
raw_df <- do.call(rbind, raw)
if (is.null(raw_df) || nrow(raw_df)==0) stop("No phylogenetic models completed")
write.csv(raw_df, file.path(outdir,"jbi_cs_opentree_phyloglm_replicates.csv"), row.names=FALSE)

agg_fun <- function(g) {
  data.frame(
    n_replicates=nrow(g),
    n_species=unique(g$n_species)[1],
    median_estimate=median(g$estimate,na.rm=TRUE),
    q025_estimate=quantile(g$estimate,0.025,na.rm=TRUE),
    q975_estimate=quantile(g$estimate,0.975,na.rm=TRUE),
    median_odds_ratio=median(g$odds_ratio,na.rm=TRUE),
    q025_odds_ratio=quantile(g$odds_ratio,0.025,na.rm=TRUE),
    q975_odds_ratio=quantile(g$odds_ratio,0.975,na.rm=TRUE),
    fraction_or_below_1=mean(g$odds_ratio<1,na.rm=TRUE),
    median_p_value=median(g$p_value,na.rm=TRUE),
    fraction_p_below_0_05=mean(g$p_value<0.05,na.rm=TRUE),
    median_alpha=median(g$alpha,na.rm=TRUE)
  )
}
groups <- split(raw_df, interaction(raw_df$outcome_short,raw_df$metric,drop=TRUE))
summary_rows <- lapply(groups, function(g) {
  z <- agg_fun(g)
  cbind(outcome_short=g$outcome_short[1], metric=g$metric[1], z, stringsAsFactors=FALSE)
})
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
write.csv(summary_df, file.path(outdir,"jbi_cs_opentree_phyloglm_summary.csv"), row.names=FALSE)

matched_out <- data.frame(canonical_name=matched_species, tree_tip=tr$tip.label, stringsAsFactors=FALSE)
write.csv(matched_out, file.path(outdir,"jbi_cs_opentree_matched_species.csv"), row.names=FALSE)

manifest <- list(
  status="complete",
  requested_species=nrow(d),
  unambiguous_nonapproximate_tnrs_matches=nrow(usable),
  induced_tree_species=length(tr$tip.label),
  unmatched_or_excluded_species=setdiff(d$canonical_name, matched_species),
  replicates_requested=replicates,
  completed_model_rows=nrow(raw_df),
  expected_models_per_replicate=length(metrics)*length(outcomes),
  branch_lengths="Grafen, power=1, applied after random polytomy resolution",
  phylogenetic_model="phylolm::phyloglm(method='logistic_MPLE')",
  outcomes="C and S positive documented-evidence axes fitted separately",
  seed=seed,
  semantic_guard="OpenTree topology is a phylogenetic sensitivity analysis; unmatched taxa are reported and no historical tree is reused."
)
jsonlite::write_json(manifest, file.path(outdir,"jbi_cs_opentree_phyloglm_manifest.json"), auto_unbox=TRUE, pretty=TRUE)
print(summary_df)
print(manifest)

if (length(tr$tip.label) < 20) stop("Induced tree below minimum matched-species threshold")
if (any(summary_df$n_replicates < max(1, floor(0.8*replicates)))) stop("Too many failed phyloglm replicates")
