#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phylolm)
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
OUT <- file.path(ROOT, "results", "polymorphism_h3b_reserve_span_20260912")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DISC <- file.path(ROOT, "data", "derived", "global_monte_carlo_measured_photos_v1.csv")
RES <- file.path(ROOT, "data", "derived", "rgfca_reserve_replication_measured_photos_v1.csv")
MORPHS <- c("white", "yellow_orange", "red_pink", "blue_purple")
MIN_CLASS <- 40L
N_PERM <- 20000L
SEED <- 20260912L

as_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

compute_D <- function(path, cohort) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- as_bool(d$global_classifiable) & d$morph %in% MORPHS & nzchar(d$species)
  x <- d[keep, c("species", "morph"), drop = FALSE]
  tab <- table(x$species, factor(x$morph, levels = MORPHS))
  n <- rowSums(tab)
  tab <- tab[n >= MIN_CLASS, , drop = FALSE]
  n <- rowSums(tab)
  p <- sweep(tab, 1L, n, "/")
  D <- 1 - rowSums(p * p)
  out <- data.frame(
    cohort = cohort,
    species = rownames(tab),
    n_classifiable = as.integer(n),
    D = as.numeric(D),
    D_unbiased = as.numeric(D * n / (n - 1)),
    stringsAsFactors = FALSE
  )
  out[order(out$species), , drop = FALSE]
}

rank_center <- function(x) {
  r <- rank(x, ties.method = "average")
  r - mean(r)
}

spearman_perm <- function(y, x, seed) {
  ry <- rank_center(y)
  rx <- rank_center(x)
  den <- sqrt(sum(ry^2) * sum(rx^2))
  obs <- sum(ry * rx) / den
  set.seed(seed)
  null <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) null[[i]] <- sum(sample(ry, replace = FALSE) * rx) / den
  p <- (1 + sum(abs(null) >= abs(obs) - 1e-15)) / (N_PERM + 1)
  list(rho = as.numeric(obs), p = as.numeric(p), null = null)
}

partial_rank_perm <- function(y, x, nclass, nobs, seed) {
  ry <- rank_center(y)
  rx <- rank_center(x)
  c1 <- rank_center(log1p(nclass))
  c2 <- rank_center(log1p(nobs))
  X <- cbind(1, c1, c2)
  if (qr(X)$rank != ncol(X)) stop("Partial-rank control matrix is rank deficient")
  ey <- lm.fit(X, ry)$residuals
  ex <- lm.fit(X, rx)$residuals
  obs <- cor(ey, ex)
  set.seed(seed)
  null <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) null[[i]] <- cor(sample(ey, replace = FALSE), ex)
  p <- (1 + sum(abs(null) >= abs(obs) - 1e-15)) / (N_PERM + 1)
  list(rho_partial = as.numeric(obs), p = as.numeric(p), null = null)
}

make_frame <- function(d, cov) {
  z <- merge(d, cov, by = c("cohort", "species"), all.x = TRUE, sort = FALSE)
  if (nrow(z) != nrow(d)) stop("Covariate join changed row count")
  req <- c("n_images_all_measured", "n_observers_all_measured", "maximum_span_km_all_measured")
  if (any(!complete.cases(z[, req, drop = FALSE]))) stop("Missing frozen span/observer covariates")
  if (any(z$n_images_all_measured != 100L)) stop("All-image count no longer fixed at 100")
  z$log1p_span <- log1p(z$maximum_span_km_all_measured)
  z
}

rank_pgls <- function(tree, frame, expected_n) {
  d <- frame
  d$tip_label <- gsub(" ", "_", d$species, fixed = TRUE)
  keep <- intersect(tree$tip.label, d$tip_label)
  tr <- if (length(setdiff(tree$tip.label, keep))) drop.tip(tree, setdiff(tree$tip.label, keep)) else tree
  d <- d[match(tr$tip.label, d$tip_label), , drop = FALSE]
  if (nrow(d) != expected_n || any(is.na(d$species))) stop("PGLS tree/data fingerprint mismatch")
  d$y_rank <- rank_center(d$D)
  d$span_rank <- rank_center(d$log1p_span)
  d$nclass_rank <- rank_center(log1p(d$n_classifiable))
  d$nobs_rank <- rank_center(log1p(d$n_observers_all_measured))
  rownames(d) <- tr$tip.label
  fit <- phylolm::phylolm(
    y_rank ~ span_rank + nclass_rank + nobs_rank,
    data = d,
    phy = tr,
    model = "lambda"
  )
  beta <- unname(coef(fit)[["span_rank"]])
  se <- sqrt(diag(vcov(fit)))[["span_rank"]]
  z <- beta / se
  data.frame(
    n_tips = nrow(d),
    beta_span_rank = beta,
    se_span_rank = se,
    z_span_rank = z,
    p_span_rank = 2 * pnorm(-abs(z)),
    lambda = if (!is.null(fit$optpar)) as.numeric(fit$optpar) else NA_real_,
    stringsAsFactors = FALSE
  )
}

disc <- compute_D(DISC, "discovery")
res <- compute_D(RES, "reserve")
if (nrow(disc) != 369L) stop("Discovery D fingerprint drift: ", nrow(disc))
if (nrow(res) != 363L) stop("Reserve D fingerprint drift: ", nrow(res))
if (length(intersect(disc$species, res$species)) != 0L) stop("Cohorts are not species-disjoint")
if (abs(max(disc$D) - 0.707645) > 1e-5) stop("Discovery D fingerprint drift")

cov <- read.csv(COVARIATE_PATH, stringsAsFactors = FALSE, check.names = FALSE)
disc <- make_frame(disc, cov[cov$cohort == "discovery", , drop = FALSE])
res <- make_frame(res, cov[cov$cohort == "reserve", , drop = FALSE])

analyze <- function(frame, cohort_index) {
  raw <- spearman_perm(frame$D, frame$log1p_span, SEED + cohort_index)
  unb <- spearman_perm(frame$D_unbiased, frame$log1p_span, SEED + 10L + cohort_index)
  partial <- partial_rank_perm(
    frame$D, frame$log1p_span,
    frame$n_classifiable, frame$n_observers_all_measured,
    SEED + 20L + cohort_index
  )
  list(raw = raw, unbiased = unb, partial = partial)
}

a_disc <- analyze(disc, 1L)
a_res <- analyze(res, 2L)

primary_pass <- isTRUE(a_res$raw$rho > 0 && a_res$raw$p < 0.05)
opportunity_robust <- isTRUE(primary_pass && a_res$partial$rho_partial > 0 && a_res$partial$p < 0.05)
verdict <- if (primary_pass) "H3B_SAMPLED_SPAN_REPLICATION_SUPPORTED" else "H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED"

scenario_files <- c(S1 = "h3a_s1.tre", S2 = "h3a_s2.tre", S3 = "h3a_s3.tre")
pgls_rows <- list()
for (i in seq_along(scenario_files)) {
  scenario <- names(scenario_files)[[i]]
  tr <- read.tree(file.path(TREE_ROOT, scenario_files[[i]]))
  x <- rank_pgls(tr, res, 341L)
  x$scenario <- scenario
  pgls_rows[[i]] <- x
}
pgls <- do.call(rbind, pgls_rows)
pgls <- pgls[, c("scenario", setdiff(names(pgls), "scenario")), drop = FALSE]
write.csv(pgls, file.path(OUT, "reserve_rank_pgls_s1_s3.csv"), row.names = FALSE)

summary_rows <- data.frame(
  cohort = c("discovery", "reserve"),
  n_species = c(nrow(disc), nrow(res)),
  rho_D_span = c(a_disc$raw$rho, a_res$raw$rho),
  p_D_span = c(a_disc$raw$p, a_res$raw$p),
  rho_Dunbiased_span = c(a_disc$unbiased$rho, a_res$unbiased$rho),
  p_Dunbiased_span = c(a_disc$unbiased$p, a_res$unbiased$p),
  partial_rho_D_span = c(a_disc$partial$rho_partial, a_res$partial$rho_partial),
  p_partial_D_span = c(a_disc$partial$p, a_res$partial$p),
  stringsAsFactors = FALSE
)
write.csv(summary_rows, file.path(OUT, "h3b_span_summary.csv"), row.names = FALSE)

nulls <- rbind(
  data.frame(cohort = "discovery", test = "raw_D", rho_null = a_disc$raw$null),
  data.frame(cohort = "reserve", test = "raw_D", rho_null = a_res$raw$null),
  data.frame(cohort = "discovery", test = "D_unbiased", rho_null = a_disc$unbiased$null),
  data.frame(cohort = "reserve", test = "D_unbiased", rho_null = a_res$unbiased$null),
  data.frame(cohort = "discovery", test = "partial_rank", rho_null = a_disc$partial$null),
  data.frame(cohort = "reserve", test = "partial_rank", rho_null = a_res$partial$null)
)
write.csv(nulls, file.path(OUT, "h3b_span_permutation_nulls.csv"), row.names = FALSE)

result <- list(
  analysis = "polymorphism_h3b_reserve_span",
  date_jst = "2026-09-12",
  protocol = "docs/POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md",
  hypothesis_status = "discovery_selected_reserve_fresh_replication",
  predictor = "log1p(maximum_span_km_all_measured)",
  predictor_guard = "sampled span from all 100 measured rows; not true range size",
  permutations = N_PERM,
  discovery_calibration = list(
    n = nrow(disc),
    rho = a_disc$raw$rho,
    p_two_sided = a_disc$raw$p,
    rho_D_unbiased = a_disc$unbiased$rho,
    p_D_unbiased = a_disc$unbiased$p,
    partial_rho = a_disc$partial$rho_partial,
    partial_p = a_disc$partial$p
  ),
  reserve_primary = list(
    n = nrow(res),
    rho = a_res$raw$rho,
    p_two_sided = a_res$raw$p,
    rho_D_unbiased = a_res$unbiased$rho,
    p_D_unbiased = a_res$unbiased$p,
    partial_rho = a_res$partial$rho_partial,
    partial_p = a_res$partial$p
  ),
  reserve_rank_pgls = split(pgls, seq_len(nrow(pgls))),
  decision = list(
    primary_positive = a_res$raw$rho > 0,
    primary_p_lt_0_05 = a_res$raw$p < 0.05,
    opportunity_robust = opportunity_robust,
    verdict = verdict,
    label = if (opportunity_robust) paste(verdict, "OPPORTUNITY_ROBUST", sep = "+") else verdict
  ),
  package_versions = list(R = R.version.string, ape = as.character(packageVersion("ape")), phylolm = as.character(packageVersion("phylolm")), jsonlite = as.character(packageVersion("jsonlite"))),
  hard_nonclaim = "Sampled geographic span is observational opportunity/geography, not true range size or a causal environmental mechanism."
)
writeLines(toJSON(result, pretty = TRUE, auto_unbox = TRUE, na = "null"), file.path(OUT, "result.json"))
cat(toJSON(result$decision, pretty = TRUE, auto_unbox = TRUE), "\n")
