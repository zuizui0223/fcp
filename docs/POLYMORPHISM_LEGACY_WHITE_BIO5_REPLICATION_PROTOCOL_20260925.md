# Legacy high-depth BIO5–white replication protocol — 2026-09-25

Status: **frozen before opening BIO5 outcomes in the original discovery/reserve row tables for this replication test**.

This analysis is separate from the frozen New Phytologist manuscript and from the third-cohort mechanism tests.

## Question

Does the broad within-species association between white floral states and warmer long-term environments (WorldClim BIO5) replicate in the two original high-depth cohorts that are species-disjoint from the later third cohort?

## Frozen source rows

Immutable source commit:
`5142f7951af0dde5364bb047a566d67e8c479e51`

Discovery:
`data/derived/global_monte_carlo_measured_photos_v1.csv`
SHA256: `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`

Reserve:
`data/derived/rgfca_reserve_replication_measured_photos_v1.csv`
SHA256: `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`

Only rows passing `global_classifiable` and assigned to the four frozen biological states are used.

## Environmental predictor

WorldClim 2.1, 10 arc-minute BIO5: maximum temperature of the warmest month.

Prediction fixed from the earlier third-cohort result:

> within species, white records occur at locations with higher BIO5.

No additional environmental predictor is tested.

## Primary species-level contrast

For each cohort and species:

`delta_i = mean(z_within_species(BIO5) | white) - mean(z_within_species(BIO5) | nonwhite)`

Eligibility:
- >=5 white rows;
- >=5 nonwhite rows.

Each cohort is estimable only with >=100 eligible species.

Primary evidence for each cohort:
1. median delta > 0;
2. two-sided Wilcoxon signed-rank p < 0.05.

## Photo-level corroboration

Conditional logistic regression stratified by species:

`white ~ z_within_species(BIO5) | species`

Corroboration requires beta > 0 and p < 0.05.

## Replication decision

`LEGACY_BIO5_WHITE_REPLICATION_SUPPORTED` only if both discovery and reserve:
- are estimable;
- have median delta > 0;
- have Wilcoxon p < 0.05;
- have conditional-logistic beta > 0 and p < 0.05.

Otherwise:
`LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST`.

Discovery and reserve are treated symmetrically here; neither may rescue failure of the other.

## Observer sensitivity

A stricter descriptive sensitivity uses only species×observer strata containing both white and nonwhite rows and fits:

`white ~ z_within_species(BIO5) | species × observer`.

This sensitivity is reported when >=30 species remain. It cannot rescue a failed primary replication and does not invalidate a primary result by itself.

## Hard boundaries

This test:
- does not establish causal heat selection;
- does not distinguish genetic from plastic white states;
- does not remove the known image-exposure coupling of the coarse white classifier;
- uses the same iNaturalist/measurement system as the main FCP programme;
- therefore tests species-disjoint transport of a broad environmental sorting association, not independent-source causation.
