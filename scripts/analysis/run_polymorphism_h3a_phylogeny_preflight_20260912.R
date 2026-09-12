#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(V.PhyloMaker2)
})

ROOT <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]), "..", ".."), mustWork = FALSE)
# Rscript exposes --file=<path>; recover the repository root more robustly below.
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grep("^--file=", args_full)][1])
ROOT <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
OUT <- file.path(ROOT, "results", "polymorphism_h3a_phylogeny_preflight_20260912")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DISC <- file.path(ROOT, "data", "derived", "global_monte_carlo_measured_photos_v1.csv")
RES <- file.path(ROOT, "data", "derived", "rgfca_reserve_replication_measured_photos_v1.csv")
DISC_FAMILY <- file.path(ROOT, "results", "polymorphism_species_attributes_step4_preflight_20260910", "covariate_panel_preoutcome.csv")

MORPHS <- c("white", "yellow_orange", "red_pink", "blue_purple")
MIN_CLASS <- 40L
MIN_COVERAGE <- 0.90
MIN_TIPS <- 250L
SEED <- 20260912L

as_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

eligible_species <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                colClasses = c(species = "character", morph = "character"))
  req <- c("species", "morph", "global_classifiable")
  miss <- setdiff(req, names(d))
  if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))
  keep <- as_bool(d$global_classifiable) & d$morph %in% MORPHS & nzchar(d$species)
  tab <- table(d$species[keep])
  sp <- sort(names(tab)[tab >= MIN_CLASS])
  data.frame(species = sp, n_classifiable = as.integer(tab[sp]), stringsAsFactors = FALSE)
}

find_col <- function(df, candidates) {
  nms <- names(df)
  low <- tolower(nms)
  for (cand in tolower(candidates)) {
    hit <- which(low == cand)
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

family_map_from_table <- function(df, source_name) {
  if (is.null(df) || !is.data.frame(df)) return(data.frame())
  gcol <- find_col(df, c("genus", "genus.name", "genus_name"))
  fcol <- find_col(df, c("family", "family.name", "family_name"))
  if (is.na(gcol) || is.na(fcol)) return(data.frame())
  x <- data.frame(
    genus = trimws(as.character(df[[gcol]])),
    family = trimws(as.character(df[[fcol]])),
    source = source_name,
    stringsAsFactors = FALSE
  )
  x <- x[nzchar(x$genus) & nzchar(x$family) & !is.na(x$genus) & !is.na(x$family), , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  # Keep only genera that map uniquely to one family in the source table.
  fams <- split(x$family, x$genus)
  good <- names(fams)[vapply(fams, function(z) length(unique(z)) == 1L, logical(1))]
  out <- data.frame(
    genus = good,
    family = vapply(fams[good], function(z) unique(z)[1], character(1)),
    source = source_name,
    stringsAsFactors = FALSE
  )
  out
}

merge_maps <- function(maps) {
  maps <- maps[vapply(maps, nrow, integer(1)) > 0L]
  if (!length(maps)) return(data.frame(genus = character(), family = character(), source = character()))
  x <- do.call(rbind, maps)
  # Source order is priority order. If later sources disagree, preserve the first
  # mapping but record the conflict separately in the audit table.
  x$priority <- seq_len(nrow(x))
  x <- x[order(x$priority), , drop = FALSE]
  x[!duplicated(x$genus), c("genus", "family", "source"), drop = FALSE]
}

disc <- eligible_species(DISC)
res <- eligible_species(RES)
if (nrow(disc) != 369L) stop("Discovery eligibility fingerprint drift: ", nrow(disc))
if (nrow(res) != 363L) stop("Reserve eligibility fingerprint drift: ", nrow(res))
if (length(intersect(disc$species, res$species)) != 0L) stop("Discovery/reserve are not species-disjoint")

union_sp <- sort(unique(c(disc$species, res$species)))
sp <- data.frame(
  species = union_sp,
  genus = vapply(strsplit(union_sp, " +"), `[[`, character(1), 1L),
  stringsAsFactors = FALSE
)

# Load the exact dated backbone system used previously in the repository.
data("GBOTB.extended.LCVP", package = "V.PhyloMaker2", envir = environment())
data("nodes.info.1.LCVP", package = "V.PhyloMaker2", envir = environment())
if (!exists("GBOTB.extended.LCVP") || !exists("nodes.info.1.LCVP")) stop("LCVP backbone objects unavailable")

maps <- list()
# Highest-priority map: the already-frozen discovery pre-outcome family panel.
if (file.exists(DISC_FAMILY)) {
  fam <- read.csv(DISC_FAMILY, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c("species", "family") %in% names(fam))) {
    fam$genus <- vapply(strsplit(as.character(fam$species), " +"), `[[`, character(1), 1L)
    maps[[length(maps) + 1L]] <- family_map_from_table(fam[, c("genus", "family"), drop = FALSE], "frozen_discovery_panel")
  }
}
maps[[length(maps) + 1L]] <- family_map_from_table(nodes.info.1.LCVP, "nodes.info.1.LCVP")

# V.PhyloMaker2 versions commonly expose a tips.info table. Use it if present,
# without making it a required dependency of the protocol.
try({
  data("tips.info.1.LCVP", package = "V.PhyloMaker2", envir = environment())
  if (exists("tips.info.1.LCVP")) {
    maps[[length(maps) + 1L]] <- family_map_from_table(tips.info.1.LCVP, "tips.info.1.LCVP")
  }
}, silent = TRUE)

fmap <- merge_maps(maps)
sp$family <- fmap$family[match(sp$genus, fmap$genus)]
sp$family_source <- fmap$source[match(sp$genus, fmap$genus)]

# For species already represented in the backbone, a missing family should not
# block a direct prune. Assign a sentinel only for the preflight call; such rows
# are audited separately and are not counted as family-resolved.
backbone_tips <- GBOTB.extended.LCVP$tip.label
sp$tip_label <- gsub(" ", "_", sp$species, fixed = TRUE)
sp$present_in_backbone <- sp$tip_label %in% backbone_tips
sp$family_resolved <- !is.na(sp$family) & nzchar(sp$family)
sp$family_for_maker <- sp$family
sp$family_for_maker[!sp$family_resolved & sp$present_in_backbone] <- "BACKBONE_DIRECT"

usable <- sp$family_resolved | sp$present_in_backbone
sp_use <- sp[usable, , drop = FALSE]
sp_list <- data.frame(
  species = sp_use$species,
  genus = sp_use$genus,
  family = sp_use$family_for_maker,
  species.relative = "",
  genus.relative = "",
  stringsAsFactors = FALSE
)
write.csv(sp, file.path(OUT, "species_family_preflight.csv"), row.names = FALSE)
write.csv(sp_list, file.path(OUT, "vphylomaker2_input_species.csv"), row.names = FALSE)

set.seed(SEED)
made <- phylo.maker(
  sp.list = sp_list,
  tree = GBOTB.extended.LCVP,
  nodes = nodes.info.1.LCVP,
  output.sp.list = TRUE,
  output.tree = FALSE,
  scenarios = c("S1", "S2", "S3")
)
if (is.null(made$species.list)) stop("V.PhyloMaker2 did not return species.list")
write.csv(made$species.list, file.path(OUT, "vphylomaker2_species_placement.csv"), row.names = FALSE)

scenarios <- list(S1 = made$scenario.1, S2 = made$scenario.2, S3 = made$scenario.3)
coverage_rows <- list()
for (nm in names(scenarios)) {
  tr <- scenarios[[nm]]
  if (is.null(tr) || !inherits(tr, "phylo")) stop("Missing phylo tree for ", nm)
  write.tree(tr, file.path(OUT, paste0("h3a_", tolower(nm), ".tre")))
  for (cohort in c("discovery", "reserve")) {
    base <- if (cohort == "discovery") disc else res
    target <- gsub(" ", "_", base$species, fixed = TRUE)
    retained <- sum(target %in% tr$tip.label)
    total <- length(target)
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      scenario = nm,
      cohort = cohort,
      eligible_species = total,
      retained_tips = retained,
      coverage = retained / total,
      gate_coverage_0_90 = retained / total >= MIN_COVERAGE,
      gate_min_250 = retained >= MIN_TIPS,
      gate_pass = retained / total >= MIN_COVERAGE && retained >= MIN_TIPS,
      stringsAsFactors = FALSE
    )
  }
}
coverage <- do.call(rbind, coverage_rows)
write.csv(coverage, file.path(OUT, "tree_coverage.csv"), row.names = FALSE)

placement_status <- if ("status" %in% names(made$species.list)) as.list(table(made$species.list$status, useNA = "ifany")) else list()
all_gate <- all(coverage$gate_pass)

pkg <- packageDescription("V.PhyloMaker2")
result <- list(
  analysis = "polymorphism_h3a_phylogeny_preflight",
  date_jst = "2026-09-12",
  D_association_computed = FALSE,
  phylogenetic_signal_computed = FALSE,
  discovery_eligible_species = nrow(disc),
  reserve_eligible_species = nrow(res),
  species_disjoint = TRUE,
  union_species = length(union_sp),
  family_resolved_species = sum(sp$family_resolved),
  backbone_direct_species = sum(sp$present_in_backbone),
  preflight_usable_species = sum(usable),
  family_unresolved_and_not_backbone = sum(!usable),
  backbone = "GBOTB.extended.LCVP",
  node_table = "nodes.info.1.LCVP",
  seed = SEED,
  placement_status_counts = placement_status,
  coverage = split(coverage, seq_len(nrow(coverage))),
  coverage_gate_all_scenarios_both_cohorts = all_gate,
  package = list(
    version = as.character(packageVersion("V.PhyloMaker2")),
    remote_sha = if (is.null(pkg$RemoteSha)) NA_character_ else as.character(pkg$RemoteSha),
    remote_repo = if (is.null(pkg$RemoteRepo)) NA_character_ else as.character(pkg$RemoteRepo)
  ),
  next_action = if (all_gate) "H3A_SIGNAL_ANALYSIS_ALLOWED" else "STOP_AND_REPAIR_TREE_COVERAGE_BEFORE_OPENING_H3A"
)
writeLines(toJSON(result, pretty = TRUE, auto_unbox = TRUE, na = "null"), file.path(OUT, "result.json"))

cat(toJSON(result, pretty = TRUE, auto_unbox = TRUE, na = "null"), "\n")
if (!all_gate) quit(status = 2L)
