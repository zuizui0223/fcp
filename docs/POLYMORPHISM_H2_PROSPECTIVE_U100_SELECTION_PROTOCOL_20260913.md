# Prospective H2 confirmation — U100 species selection protocol

Frozen: 2026-09-13 JST

## Purpose

The current white-versus-non-white H2 target was isolated after the broad H2 geometry had already been opened in the existing 369-species discovery and 363-species reserve cohorts. The next inferential upgrade is therefore a **new species-level prospective confirmation set selected without using colour outcomes**.

This protocol freezes that selection before any high-depth colour outcome is opened for the new set.

## Parent frame

Start from the already frozen U100 opportunity frame:

- global sampling universe U0 = 42,111 species;
- U100 = species with `after_observer_cap >= 100` in `results/rgfca_42111_breadth_depth_step8b_20260911/species_capacity_census.csv.gz`;
- frozen U100 size = **4,730 species**.

The U100 census was generated from metadata/capacity only. Its generating analysis explicitly records `pixels_opened = false` and `flower_colour_used = false`.

## Exclusion of target-development species

Before prospective selection, exclude the union of species present in the two high-depth datasets that were used to develop, audit, transport, or target the current H1/H2 geometry:

- `data/derived/global_monte_carlo_measured_photos_v1.csv` (discovery; 369 D-eligible species);
- `data/derived/rgfca_reserve_replication_measured_photos_v1.csv` (reserve; 363 D-eligible species).

Only the `species` identity column may be read from these two outcome-bearing files during selection. Morph labels, continuous palette coordinates, D, H2 cluster assignments, geography, or any H3 result must not be read.

Define:

`P100 = U100 \ (discovery_species union reserve_species)`.

The exact size of P100 is unknown until the metadata-only selection preflight is run. Species outside U100 are not substituted if exclusions reduce the pool.

One-anchor global breadth colour measurements do not disqualify a species from P100 because they do not identify within-species high-depth geometry; however, P100 is described specifically as **species-disjoint from the high-depth cohorts used to construct the H2 target**, not as never having any flower image inspected anywhere in the project.

## Prospective confirmation sample P500

To keep the high-depth confirmation tractable while retaining a substantially larger species-level test than either existing H2 cohort, freeze a 500-species sample from P100.

For every P100 species, compute

`selection_hash = SHA256("FCP_H2_PROSPECTIVE_20260913|" + inat_taxon_id + "|" + species)`.

Sort ascending by `selection_hash`, then by `inat_taxon_id` as a deterministic tie-break. The first **500 species** are `P500`.

This ordering uses only frozen taxon identity and the preregistered salt. It does not use colour, family, genus, latitude, range, observer count beyond the U100 capacity gate, or any H1/H2/H3 outcome.

If P100 contains fewer than 500 species, stop as `PROSPECTIVE_SELECTION_UNDERIDENTIFIED`; do not lower the U100 threshold or refill from lower-capacity strata.

## Why 500 species

The number is a resource/design commitment rather than an outcome-adaptive power threshold. At 100 raw eligible photos per species it corresponds to at most 50,000 candidate images, comparable to prior high-depth image workflows while allowing attrition through classifiability, H1, and continuous-mode gates.

No minimum number of eventual H1/H2-positive species is assumed from the legacy cohorts. The actual eligibility funnel must be reported as observed attrition.

## Frozen downstream white-axis target

P500 is not used to discover a new direction. If high-depth measurement becomes available, the prospective H2 confirmation must reuse without alteration:

- four-state H1 admission: `n_classifiable >= 40`, second coarse state fraction >=0.10 primary and >=0.20 strict;
- `mixed_uncertain` excluded as a biological morph;
- normalized nine-colour biological palette;
- Hellinger transform;
- deterministic unlabeled two-means mode construction;
- continuous smaller-mode fraction >=0.10 primary and >=0.20 strict;
- canonical fixed `q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`;
- `W = mean_i (u_i^T q_white)^2`;
- the same coarse-state-preserving structured null used in the frozen targeted H2 test.

The current discovery leading axis must not be refit and substituted for `q_white`.

## Prospective primary decision

The exact Monte Carlo implementation and significance threshold remain those of the frozen targeted white-axis test. The prospective primary test is the P500 primary-gate observed W against the same coarse-state-preserving structured null.

Prospective support requires:

1. at least **100 P500 species** survive H1 plus the primary continuous-mode gate; and
2. structured-null upper-tail **p < 0.05** for the fixed W statistic.

The >=100 H2-species floor is an identifiability gate fixed before P500 outcomes are measured. If fewer than 100 species survive, verdict is `H2_PROSPECTIVE_WHITE_AXIS_UNDERIDENTIFIED`, not biological non-support.

If N >=100 and p <0.05: `H2_PROSPECTIVE_WHITE_AXIS_SUPPORTED`.

If N >=100 and p >=0.05: `H2_PROSPECTIVE_WHITE_AXIS_NOT_SUPPORTED`.

The strict 20% gate is a prespecified sensitivity and cannot rescue primary failure.

## Hard stops

1. Do not alter P500 after any P500 colour outcome is opened.
2. Do not replace failed/unavailable P500 species with the next hashes after outcome opening. Technical acquisition failures are reason-coded attrition.
3. Do not lower U100, H1, continuous-mode, or N>=100 gates post hoc.
4. Do not rotate q_white, discover a new palette direction, or use the P500 leading axis as the primary target.
5. Do not use H3 predictors to stratify, subset, or rescue the prospective H2 test.
6. Do not call an underidentified test negative evidence for the biological axis.
7. Do not estimate global polymorphism prevalence from P500.

## Required pre-outcome receipt

Before any P500 high-depth colour result is opened, save:

- exact U100 count;
- exact legacy exclusion count and overlap with U100;
- exact P100 count;
- ordered P500 taxon IDs/species and selection hashes;
- SHA256 of the P500 CSV;
- confirmation that selection read no colour-outcome columns.