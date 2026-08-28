#!/usr/bin/env Rscript

# Conservative cross-source trait coverage audit for FCP display-core-v6.
# Uses taxify's registered harmonization maps but does not use a taxonomic
# backbone: focal names have already been curated and are supplied as
# accepted_name. A consensus value is retained only when all non-missing
# source values agree after taxify harmonization. Conflicts remain explicit.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing argument: ", flag)
  args[[i + 1L]]
}
core_path <- get_arg("--core-species")
outdir <- get_arg("--outdir")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(taxify))

core <- read.csv(core_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(core) == 74L)
required <- c("canonical_name", "family", "organization_state")
stopifnot(all(required %in% names(core)))

base <- core[, required]
base$accepted_name <- base$canonical_name

registry <- taxify::list_traits()
write.csv(registry, file.path(outdir, "cross_source_trait_registry_snapshot.csv"), row.names = FALSE)

targets <- c("dispersal_syndrome", "flowering_start", "flowering_end")
missing_targets <- setdiff(targets, registry$trait)
if (length(missing_targets)) {
  stop("Target traits absent from taxify registry: ", paste(missing_targets, collapse = ", "))
}

informative_states <- c(
  "local_coexistence_only",
  "spatial_segregation_only",
  "coexistence_and_segregation"
)

out <- base
coverage <- list()
source_info <- list()

for (trait in targets) {
  info <- taxify::trait_info(trait)
  # Preserve registry/source provenance in a text-safe form.
  if (is.data.frame(info)) {
    write.csv(info, file.path(outdir, paste0("trait_info_", trait, ".csv")), row.names = FALSE)
  } else {
    capture.output(str(info), file = file.path(outdir, paste0("trait_info_", trait, ".txt")))
  }

  wide <- taxify::add_trait(base, trait, mode = "wide", sources = "all", verbose = TRUE)

  # Trait-source columns are everything named trait_<source>, excluding taxify
  # metadata/provenance columns. The registry returns already-harmonized values.
  pref <- paste0("^", trait, "_")
  cols <- grep(pref, names(wide), value = TRUE)
  meta_suffix <- c("unit", "caution", "sources", "n", "min", "max")
  source_cols <- cols[!sub(pref, "", cols) %in% meta_suffix]
  if (!length(source_cols)) stop("No source columns produced for ", trait)

  # Keep source-specific values for auditability.
  for (cc in source_cols) out[[cc]] <- wide[[cc]]

  consensus <- rep(NA_character_, nrow(wide))
  n_sources <- integer(nrow(wide))
  conflict <- logical(nrow(wide))
  for (i in seq_len(nrow(wide))) {
    vals <- unlist(wide[i, source_cols, drop = FALSE], use.names = FALSE)
    vals <- trimws(as.character(vals))
    vals <- vals[!is.na(vals) & nzchar(vals) & vals != "NA"]
    n_sources[[i]] <- length(vals)
    u <- unique(vals)
    if (length(u) == 1L) consensus[[i]] <- u[[1L]]
    if (length(u) > 1L) conflict[[i]] <- TRUE
  }
  out[[paste0(trait, "_consensus")]] <- consensus
  out[[paste0(trait, "_n_sources")]] <- n_sources
  out[[paste0(trait, "_conflict")]] <- conflict

  inf <- out$organization_state %in% informative_states
  val <- !is.na(consensus)
  coverage[[trait]] <- data.frame(
    trait = trait,
    n_registered_sources = length(source_cols),
    registered_source_columns = paste(source_cols, collapse = ";"),
    matched_core_species = sum(val),
    matched_informative_species = sum(val & inf),
    C_only_n = sum(val & out$organization_state == "local_coexistence_only"),
    S_only_n = sum(val & out$organization_state == "spatial_segregation_only"),
    mixed_n = sum(val & out$organization_state == "coexistence_and_segregation"),
    conflict_core_species = sum(conflict),
    conflict_informative_species = sum(conflict & inf),
    unique_consensus_values = length(unique(consensus[val])),
    stringsAsFactors = FALSE
  )
}

write.csv(out, file.path(outdir, "cross_source_core_trait_coverage.csv"), row.names = FALSE)
coverage_df <- do.call(rbind, coverage)
write.csv(coverage_df, file.path(outdir, "cross_source_core_trait_coverage_summary.csv"), row.names = FALSE)

# Explicit prospective gates. Passing a coverage gate only licenses a later
# model; it does not make a trait result significant.
coverage_df$coverage_gate_basic <- with(
  coverage_df,
  matched_informative_species >= 15 & C_only_n >= 4 & S_only_n >= 4 & mixed_n >= 4 & unique_consensus_values >= 2
)
coverage_df$coverage_gate_interaction <- with(
  coverage_df,
  matched_informative_species >= 18 & C_only_n >= 5 & S_only_n >= 4 & mixed_n >= 5 & unique_consensus_values >= 2
)
write.csv(coverage_df, file.path(outdir, "cross_source_core_trait_coverage_gates.csv"), row.names = FALSE)

qc <- list(
  status = "complete",
  protocol = "cross-source-trait-coverage-v1-concordance-only",
  core_species = nrow(core),
  target_traits = targets,
  harmonization = "taxify registered per-source crosswalks",
  consensus_rule = "retain only when all non-missing harmonized source values agree; conflicting rows remain missing for consensus",
  missing_semantics = "missing is unknown, never biological absence",
  basic_gate = "informative>=15; C-only>=4; S-only>=4; mixed>=4; >=2 consensus values",
  interaction_gate = "informative>=18; C-only>=5; S-only>=4; mixed>=5; >=2 consensus values"
)
jsonlite::write_json(qc, file.path(outdir, "cross_source_core_trait_coverage_qc.json"), pretty = TRUE, auto_unbox = TRUE)
cat(paste(capture.output(print(coverage_df)), collapse = "\n"), "\n")
