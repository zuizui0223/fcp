#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clubSandwich)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[[i+1]]
}
dataset <- get_arg("--dataset")
outdir <- get_arg("--outdir")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

d <- read.csv(dataset, stringsAsFactors=FALSE, check.names=FALSE)
if (nrow(d) != 34) stop(sprintf("Expected 34 species, found %d", nrow(d)))
metrics <- c("temperature_breadth","moisture_breadth","climatic_heterogeneity","pca_dispersion","pca_hull_area")
if (!all(c(metrics,"spatial_scale","n_climate_cells","family") %in% names(d))) stop("Missing required columns")

z <- function(x) as.numeric(scale(as.numeric(x), center=TRUE, scale=TRUE))
rows <- list(); raw_rows <- list()

for (metric in metrics) {
  x <- data.frame(
    among = as.integer(d$spatial_scale == "among_population"),
    metric_z = z(d[[metric]]),
    effort_z = z(log1p(as.numeric(d$n_climate_cells))),
    family = d$family,
    stringsAsFactors=FALSE
  )
  fit <- glm(among ~ metric_z + effort_z, family=binomial(), data=x)
  V <- vcovCR(fit, cluster=x$family, type="CR2")
  ct <- coef_test(fit, vcov=V, test=c("z","naive-t","Satterthwaite"))
  tab <- as.data.frame(ct); tab$term <- rownames(tab); rownames(tab) <- NULL; tab$metric <- metric
  raw_rows[[metric]] <- tab
  target <- tab[tab$term == "metric_z", , drop=FALSE]
  if (nrow(target) != 1) stop(paste("Could not locate metric_z for", metric))
  cn <- names(target)
  pick <- function(patterns) {
    for (pat in patterns) {
      hit <- grep(pat, cn, ignore.case=TRUE, value=TRUE)
      if (length(hit)) return(hit[[1]])
    }
    NA_character_
  }
  est_col <- pick(c("^beta$","estimate","coef")); se_col <- pick(c("^SE$","std","se"))
  df_col <- pick(c("Satt.*df","d.f.*Satt","df_Satt","d.f")); p_col <- pick(c("Satt.*p","p.*Satt"))
  if (any(is.na(c(est_col,se_col,df_col,p_col)))) stop(paste("Unsupported clubSandwich output columns:", paste(cn, collapse=", ")))
  est <- as.numeric(target[[est_col]]); se <- as.numeric(target[[se_col]])
  df_s <- as.numeric(target[[df_col]]); p_s <- as.numeric(target[[p_col]]); crit <- qt(0.975, df=df_s)
  rows[[metric]] <- data.frame(
    metric=metric, n_species=nrow(x), n_families=length(unique(x$family)), estimate=est,
    odds_ratio=exp(est), cr2_se=se, satterthwaite_df=df_s,
    cr2_ci_low=exp(est-crit*se), cr2_ci_high=exp(est+crit*se),
    cr2_satterthwaite_p=p_s, stringsAsFactors=FALSE
  )
}
summary <- do.call(rbind, rows); raw <- do.call(rbind, raw_rows)
write.csv(summary, file.path(outdir,"jbi_34species_cr2_satterthwaite_summary.csv"), row.names=FALSE)
write.csv(raw, file.path(outdir,"jbi_34species_cr2_coef_test_raw.csv"), row.names=FALSE)
json <- paste0(
  '{\n  "status": "complete",\n  "n_species": 34,\n  "n_families": ', length(unique(d$family)),
  ',\n  "estimator": "binomial GLM with clubSandwich CR2 cluster-robust variance",',
  '\n  "test": "Satterthwaite small-sample t test",',
  '\n  "cluster": "plant family",',
  '\n  "interpretation_guard": "CR2/Satterthwaite is a small-cluster sensitivity analysis; permutation, leave-one-family-out and phylogenetic analyses remain complementary."\n}\n'
)
writeLines(json, file.path(outdir,"jbi_34species_cr2_manifest.json"))
print(summary)
