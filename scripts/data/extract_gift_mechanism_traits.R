#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(taxify)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[[i + 1L]]
}
core_path <- get_arg('--core-species')
out_path <- get_arg('--out')
coverage_path <- get_arg('--coverage-out')
if (is.null(core_path) || is.null(out_path) || is.null(coverage_path)) {
  stop('Required: --core-species --out --coverage-out')
}

core <- read_csv(core_path, show_col_types = FALSE) %>%
  transmute(canonical_name = as.character(canonical_name),
            family = as.character(family),
            organization_state = as.character(organization_state),
            C_local_coexistence_documented = as.integer(C_local_coexistence_documented),
            S_spatial_segregation_documented = as.integer(S_spatial_segregation_documented))
stopifnot(nrow(core) == 74L, n_distinct(core$canonical_name) == 74L)

# GIFT's bundled enrichment is species-level and joins by accepted_name. The
# v6 names have already passed the upstream taxon-resolution gate, so the FCP
# canonical name is used as the accepted-name key here; raw values remain
# visible for audit and are never silently imputed.
x <- core %>% mutate(accepted_name = canonical_name)
traits <- c(
  'gift_self_fertilization_1',
  'gift_lifecycle_1',
  'gift_dispersal_syndrome_1',
  'gift_flowering_start',
  'gift_flowering_end',
  'gift_seed_mass_mean'
)

y <- add_gift(x, cols = traits, verbose = TRUE)
missing_cols <- setdiff(traits, names(y))
if (length(missing_cols)) stop('Missing requested GIFT columns: ', paste(missing_cols, collapse=', '))

# Keep only the study keys plus raw GIFT values; no recoding occurs at this
# extraction step. This makes later ecological recoding separately auditable.
out <- y %>% select(all_of(names(core)), all_of(traits))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_csv(out, out_path, na = '')

informative <- out %>% filter(organization_state != 'organization_unresolved')
coverage <- bind_rows(lapply(traits, function(tr) {
  v_all <- out[[tr]]
  v_inf <- informative[[tr]]
  tibble(
    trait = tr,
    n_core_nonmissing = sum(!is.na(v_all) & trimws(as.character(v_all)) != ''),
    n_informative_nonmissing = sum(!is.na(v_inf) & trimws(as.character(v_inf)) != ''),
    n_C_only_nonmissing = sum(informative$organization_state == 'local_coexistence_only' & !is.na(v_inf) & trimws(as.character(v_inf)) != ''),
    n_S_only_nonmissing = sum(informative$organization_state == 'spatial_segregation_only' & !is.na(v_inf) & trimws(as.character(v_inf)) != ''),
    n_mixed_nonmissing = sum(informative$organization_state == 'coexistence_and_segregation' & !is.na(v_inf) & trimws(as.character(v_inf)) != ''),
    raw_levels = paste(sort(unique(as.character(v_inf[!is.na(v_inf) & trimws(as.character(v_inf)) != '']))), collapse=' | ')
  )
}))
write_csv(coverage, coverage_path, na='')

cat('GIFT_MECHANISM_COVERAGE\n')
print(coverage, n = Inf, width = Inf)
