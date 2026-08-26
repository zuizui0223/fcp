#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(clubSandwich))

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
metrics <- c("temperature_breadth","moisture_breadth","climatic_heterogeneity","pca_dispersion","pca_hull_area")
outcomes <- c(C="C_local_coexistence_documented", S="S_spatial_segregation_documented")
required <- c(metrics, outcomes, "n_climate_cells", "family", "canonical_name")
if (!all(required %in% names(d))) stop("Missing required columns")
if (nrow(d) != 34) stop(sprintf("Expected 34 species, found %d", nrow(d)))

z <- function(x) as.numeric(scale(as.numeric(x), center=TRUE, scale=TRUE))
rows <- list(); raw_rows <- list(); k <- 1L

for (short in names(outcomes)) {
  outcome_col <- outcomes[[short]]
  for (metric in metrics) {
    x <- data.frame(
      outcome = as.integer(d[[outcome_col]]),
      metric_z = z(d[[metric]]),
      effort_z = z(log1p(as.numeric(d$n_climate_cells))),
      family = d$family,
      stringsAsFactors=FALSE
    )
    fit <- glm(outcome ~ metric_z + effort_z, family=binomial(), data=x)
    V <- vcovCR(fit, cluster=x$family, type="CR2")
    ct <- coef_test(fit, vcov=V, test=c("z","naive-t","Satterthwaite"))
    tab <- as.data.frame(ct); tab$term <- rownames(tab); rownames(tab) <- NULL
    tab$outcome_short <- short; tab$outcome <- outcome_col; tab$metric <- metric
    raw_rows[[k]] <- tab
    target <- tab[tab$term == "metric_z", , drop=FALSE]
    if (nrow(target) != 1) stop(paste("Could not locate metric_z for", short, metric))
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
    df_s <- as.numeric(target[[df_col]]); p_s <- as.numeric(target[[p_col]])
    crit <- qt(0.975, df=df_s)
    rows[[k]] <- data.frame(
      outcome_short=short, outcome=outcome_col, metric=metric,
      n_species=nrow(x), n_families=length(unique(x$family)),
      n_positive=sum(x$outcome), n_negative=nrow(x)-sum(x$outcome),
      estimate=est, odds_ratio=exp(est), cr2_se=se,
      satterthwaite_df=df_s,
      cr2_ci_low=exp(est-crit*se), cr2_ci_high=exp(est+crit*se),
      cr2_satterthwaite_p=p_s,
      stringsAsFactors=FALSE
    )
    k <- k + 1L
  }
}
summary <- do.call(rbind, rows); raw <- do.call(rbind, raw_rows)
write.csv(summary, file.path(outdir,"jbi_cs_cr2_satterthwaite_summary.csv"), row.names=FALSE)
write.csv(raw, file.path(outdir,"jbi_cs_cr2_coef_test_raw.csv"), row.names=FALSE)
manifest <- paste0(
  '{\n  "status": "complete",',
  '\n  "n_species": ', nrow(d), ',',
  '\n  "n_families": ', length(unique(d$family)), ',',
  '\n  "C_positive": ', sum(d$C_local_coexistence_documented), ',',
  '\n  "S_positive": ', sum(d$S_spatial_segregation_documented), ',',
  '\n  "estimator": "binomial GLM with clubSandwich CR2 cluster-robust variance",',
  '\n  "test": "Satterthwaite small-sample t test",',
  '\n  "cluster": "plant family",',
  '\n  "formulae": ["C ~ metric_z + effort_z", "S ~ metric_z + effort_z"],',
  '\n  "semantic_guard": "C and S are positive documented-evidence outcomes; CR2 quantifies small-cluster uncertainty and does not convert zero into biological absence."\n}\n'
)
writeLines(manifest, file.path(outdir,"jbi_cs_cr2_manifest.json"))
print(summary)

if (nrow(summary) != 10) stop("Expected 10 CR2 result rows")
if (any(!is.finite(summary$cr2_satterthwaite_p))) stop("Non-finite CR2 p-value")
