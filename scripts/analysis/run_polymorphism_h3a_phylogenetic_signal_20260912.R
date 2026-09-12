#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key) {
  hit <- which(args == key)
  if (length(hit) != 1L || hit == length(args)) stop("Missing argument ", key)
  args[[hit + 1L]]
}
TREE_ROOT <- normalizePath(get_arg("--tree-root"), mustWork = TRUE)
COVARIATE_PATH <- normalizePath(get_arg("--covariate"), mustWork = TRUE)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grep("^--file=", args_full)][1])
ROOT <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
OUT <- file.path(ROOT, "results", "polymorphism_h3a_phylogenetic_signal_20260912")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DISC <- file.path(ROOT, "data", "derived", "global_monte_carlo_measured_photos_v1.csv")
RES <- file.path(ROOT, "data", "derived", "rgfca_reserve_replication_measured_photos_v1.csv")
MORPHS <- c("white", "yellow_orange", "red_pink", "blue_purple")
MIN_CLASS <- 40L
NSIM_PHYTOOLS <- 10000L # phytools includes the observed K as sim.K[1], leaving 9,999 randomized maps.
SEED_BASE <- 20260912L

as_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

compute_D <- function(path, cohort) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("species", "morph", "global_classifiable")
  miss <- setdiff(req, names(d))
  if (length(miss)) stop("Missing D columns in ", cohort, ": ", paste(miss, collapse = ", "))
  keep <- as_bool(d$global_classifiable) & d$morph %in% MORPHS & nzchar(d$species)
  x <- d[keep, c("species", "morph"), drop = FALSE]
  tab <- table(x$species, factor(x$morph, levels = MORPHS))
  n <- rowSums(tab)
  tab <- tab[n >= MIN_CLASS, , drop = FALSE]
  n <- rowSums(tab)
  p <- sweep(tab, 1L, n, "/")
  D <- 1 - rowSums(p * p)
  ordp <- t(apply(p, 1L, sort))
  second_fraction <- ordp[, ncol(ordp) - 1L]
  out <- data.frame(
    cohort = cohort,
    species = rownames(tab),
    n_classifiable = as.integer(n),
    D = as.numeric(D),
    D_unbiased = as.numeric(D * n / (n - 1)),
    second_fraction = as.numeric(second_fraction),
    stringsAsFactors = FALSE
  )
  for (j in seq_along(MORPHS)) out[[paste0("p_", MORPHS[[j]])]] <- as.numeric(p[, j])
  out <- out[order(out$species), , drop = FALSE]
  rownames(out) <- NULL
  out
}

rank_center <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  r - mean(r, na.rm = TRUE)
}

add_opportunity_residual <- function(d, cov) {
  z <- merge(d, cov, by = c("cohort", "species"), all.x = TRUE, sort = FALSE)
  if (nrow(z) != nrow(d)) stop("Covariate join changed D row count for ", unique(d$cohort))
  req <- c("n_observers_all_measured", "maximum_span_km_all_measured")
  if (any(!complete.cases(z[, req, drop = FALSE]))) stop("Missing frozen opportunity covariates")
  y <- rank_center(z$D)
  x1 <- rank_center(log1p(z$n_classifiable))
  x2 <- rank_center(log1p(z$n_observers_all_measured))
  x3 <- rank_center(log1p(z$maximum_span_km_all_measured))
  X <- cbind(1, x1, x2, x3)
  if (qr(X)$rank != ncol(X)) stop("Frozen opportunity design matrix is rank deficient")
  fit <- lm.fit(X, y)
  z$D_opportunity_residual <- as.numeric(fit$residuals)
  z
}

prepare_trait <- function(tree, frame, column, expected_n) {
  x <- frame[[column]]
  names(x) <- gsub(" ", "_", frame$species, fixed = TRUE)
  keep <- intersect(tree$tip.label, names(x))
  tr <- if (length(setdiff(tree$tip.label, keep))) drop.tip(tree, setdiff(tree$tip.label, keep)) else tree
  x <- x[tr$tip.label]
  if (length(x) != expected_n || Ntip(tr) != expected_n) {
    stop("Frozen tree/cohort tip fingerprint mismatch for ", column, ": ", length(x), " != ", expected_n)
  }
  if (any(!is.finite(x))) stop("Nonfinite trait values for ", column)
  list(tree = tr, x = x)
}

k_randomization <- function(tree, x, seed) {
  set.seed(seed)
  fit <- phytools::phylosig(tree, x, method = "K", test = TRUE, nsim = NSIM_PHYTOOLS)
  # phytools' sim.K[1] is the observed mapping; the remaining 9,999 are randomized.
  if (length(fit$sim.K) != NSIM_PHYTOOLS) stop("Unexpected phytools sim.K length")
  if (abs(fit$sim.K[[1]] - fit$K) > 1e-10) stop("phytools first sim.K is not observed K")
  p_frozen <- (1 + sum(fit$sim.K[-1] >= fit$K - 1e-15)) / NSIM_PHYTOOLS
  if (abs(as.numeric(fit$P) - p_frozen) > 1e-10) stop("phytools P disagrees with frozen Monte Carlo formula")
  list(K = as.numeric(fit$K), p = as.numeric(p_frozen), null = as.numeric(fit$sim.K[-1]))
}

lambda_test <- function(tree, x) {
  fit <- phytools::phylosig(tree, x, method = "lambda", test = TRUE)
  list(
    lambda = as.numeric(fit$lambda),
    logL = as.numeric(fit$logL),
    logL0 = as.numeric(fit$logL0),
    LR = as.numeric(2 * (fit$logL - fit$logL0)),
    p = as.numeric(fit$P)
  )
}

k_effect_only <- function(tree, x) as.numeric(phytools::phylosig(tree, x, method = "K", test = FALSE))

# Outcome opening starts here, after exact trees/covariates have been frozen externally.
disc <- compute_D(DISC, "discovery")
res <- compute_D(RES, "reserve")
if (nrow(disc) != 369L) stop("Discovery D fingerprint drift: ", nrow(disc))
if (nrow(res) != 363L) stop("Reserve D fingerprint drift: ", nrow(res))
if (length(intersect(disc$species, res$species)) != 0L) stop("Discovery/reserve are not species-disjoint")
if (abs(max(disc$D) - 0.707645) > 1e-5) stop("Discovery D maximum fingerprint drift")
if (abs(mean(disc$second_fraction >= 0.10) - 0.4661) > 5e-4) stop("Discovery second-fraction fingerprint drift")

cov <- read.csv(COVARIATE_PATH, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cohort", "species", "n_images_all_measured", "n_observers_all_measured", "maximum_span_km_all_measured") %in% names(cov))) {
  stop("Frozen covariate panel missing required columns")
}
if (any(cov$n_images_all_measured != 100L)) stop("Frozen all-image design is no longer constant at 100")

disc <- add_opportunity_residual(disc, cov[cov$cohort == "discovery", ])
res <- add_opportunity_residual(res, cov[cov$cohort == "reserve", ])
allD <- rbind(disc, res)
write.csv(allD, file.path(OUT, "h3a_species_D_and_sensitivity.csv"), row.names = FALSE)

scenario_files <- c(S1 = "h3a_s1.tre", S2 = "h3a_s2.tre", S3 = "h3a_s3.tre")
expected_n <- c(discovery = 368L, reserve = 341L)
cohorts <- list(discovery = disc, reserve = res)
rows <- list(); null_rows <- list(); kk <- 1L; nn <- 1L
for (ci in seq_along(cohorts)) {
  cohort <- names(cohorts)[[ci]]
  frame <- cohorts[[ci]]
  for (si in seq_along(scenario_files)) {
    scenario <- names(scenario_files)[[si]]
    tr <- read.tree(file.path(TREE_ROOT, scenario_files[[si]]))
    raw <- prepare_trait(tr, frame, "D", expected_n[[cohort]])
    adj <- prepare_trait(tr, frame, "D_opportunity_residual", expected_n[[cohort]])
    unb <- prepare_trait(tr, frame, "D_unbiased", expected_n[[cohort]])
    rawK <- k_randomization(raw$tree, raw$x, SEED_BASE + ci * 100L + si)
    adjK <- k_randomization(adj$tree, adj$x, SEED_BASE + ci * 100L + 20L + si)
    lam <- lambda_test(raw$tree, raw$x)
    Kub <- k_effect_only(unb$tree, unb$x)
    rows[[kk]] <- data.frame(
      cohort = cohort,
      scenario = scenario,
      n_tips = length(raw$x),
      K_raw = rawK$K,
      p_K_raw = rawK$p,
      lambda_raw = lam$lambda,
      LR_lambda0 = lam$LR,
      p_lambda0 = lam$p,
      K_D_unbiased = Kub,
      K_opportunity_residual = adjK$K,
      p_K_opportunity_residual = adjK$p,
      stringsAsFactors = FALSE
    )
    null_rows[[nn]] <- data.frame(cohort = cohort, scenario = scenario, trait = "raw_D", K_null = rawK$null); nn <- nn + 1L
    null_rows[[nn]] <- data.frame(cohort = cohort, scenario = scenario, trait = "opportunity_residual", K_null = adjK$null); nn <- nn + 1L
    kk <- kk + 1L
  }
}
result_table <- do.call(rbind, rows)
null_table <- do.call(rbind, null_rows)
write.csv(result_table, file.path(OUT, "h3a_signal_by_scenario.csv"), row.names = FALSE)
write.csv(null_table, file.path(OUT, "h3a_K_permutation_nulls.csv.gz"), row.names = FALSE)

reserve_rows <- result_table[result_table$cohort == "reserve", , drop = FALSE]
discovery_rows <- result_table[result_table$cohort == "discovery", , drop = FALSE]
raw_pass_n <- sum(reserve_rows$p_K_raw < 0.05 & reserve_rows$K_raw > 0)
adj_pass_n <- sum(reserve_rows$p_K_opportunity_residual < 0.05 & reserve_rows$K_opportunity_residual > 0)
verdict <- if (raw_pass_n == 3L) {
  "H3A_RESERVE_PHYLOGENETIC_SIGNAL_SUPPORTED"
} else if (raw_pass_n >= 1L) {
  "H3A_PLACEMENT_SENSITIVE_UNRESOLVED"
} else {
  "H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED"
}
opportunity_robust <- identical(raw_pass_n, 3L) && identical(adj_pass_n, 3L)

summary_D <- function(x) list(
  n = nrow(x),
  D_min = min(x$D), D_median = median(x$D), D_mean = mean(x$D), D_max = max(x$D),
  n_classifiable_min = min(x$n_classifiable), n_classifiable_median = median(x$n_classifiable), n_classifiable_max = max(x$n_classifiable)
)

result <- list(
  analysis = "polymorphism_h3a_phylogenetic_signal",
  date_jst = "2026-09-12",
  protocol = "docs/POLYMORPHISM_H3A_PHYLOGENETIC_SIGNAL_PROTOCOL_20260912.md",
  outcome_opened = TRUE,
  D_definition = "1 - sum_k p_k^2 over white/yellow_orange/red_pink/blue_purple among frozen global_classifiable rows; n_classifiable >= 40",
  discovery_role = "calibration_retrospective",
  reserve_role = "fresh_primary_replication",
  species_disjoint = TRUE,
  permutation = list(randomized_tip_maps = 9999L, denominator = 10000L, implementation = "phytools phylosig nsim=10000 includes observed map as sim.K[1]; frozen p=(1 + randomized null >= observed)/10000"),
  D_summary = list(discovery = summary_D(disc), reserve = summary_D(res)),
  scenario_results = split(result_table, seq_len(nrow(result_table))),
  decision = list(
    reserve_raw_scenarios_p_lt_0_05 = raw_pass_n,
    reserve_opportunity_scenarios_p_lt_0_05 = adj_pass_n,
    verdict = verdict,
    opportunity_robust = opportunity_robust,
    label = if (opportunity_robust) paste(verdict, "OPPORTUNITY_ROBUST", sep = "+") else verdict
  ),
  package_versions = list(R = R.version.string, ape = as.character(packageVersion("ape")), phytools = as.character(packageVersion("phytools")), jsonlite = as.character(packageVersion("jsonlite"))),
  hard_nonclaim = "Phylogenetic signal is not a causal mechanism, adaptive conservatism, genetic basis, gain/loss reconstruction, or pollination/life-form mechanism."
)
writeLines(toJSON(result, pretty = TRUE, auto_unbox = TRUE, na = "null"), file.path(OUT, "result.json"))
cat(toJSON(result$decision, pretty = TRUE, auto_unbox = TRUE), "\n")
