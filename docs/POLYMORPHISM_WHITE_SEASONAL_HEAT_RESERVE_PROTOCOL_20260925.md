# Species-disjoint reserve confirmation of seasonal heat — 2026-09-25

Status: **frozen after the third-cohort seasonal-heat result was opened and before monthly WorldClim tmax was attached to the legacy reserve cohort**.

## Question

Does the third-cohort pattern—white observations occurring in relatively hotter flowering months—replicate in the species-disjoint reserve high-depth cohort?

## Primary cohort

Reserve measured table frozen in repository history:
`data/derived/rgfca_reserve_replication_measured_photos_v1.csv`
at source commit `db2514f622468747a7d5896c2696fcbf13729ef6`.

The workflow must verify zero species overlap with the third-cohort measured table before opening the reserve climate result.

Only globally classifiable rows in the four frozen biological colour states are used.

## Predictor

WorldClim 2.1 10-arc-minute monthly tmax.

For each reserve photograph:
`
seasonal_heat_anomaly = tmax(location, observed calendar month) - mean_12_month_tmax(location)
`

BIO5 is included separately as a long-term spatial heat covariate.

## Technical covariates

The reserve cohort predates the direct near-clipping seal. To reduce image-context confounding without inventing a new classifier:
- use frozen `nuisance_pixel_fraction`;
- compute background white fraction from the frozen background palette counts;
- include both as within-species nuisance covariates in the photo-level model.

No new clipping threshold or image reclassification is allowed.

## Eligibility

Species-level test requires >=5 white and >=5 non-white classifiable reserve rows.
Minimum: 100 species.

## Primary tests

1. Species-level white-minus-nonwhite seasonal anomaly:
   - median delta > 0;
   - two-sided Wilcoxon p < 0.05.

2. Species-stratified conditional logistic regression:
`
white ~ seasonal_heat + BIO5 + background_white + nuisance_pixels | species
`
   - seasonal coefficient > 0 and p < 0.05.

3. Same species × WorldClim 10' cell conditional logistic regression:
`
white ~ seasonal_heat + background_white + nuisance_pixels | species × cell
`
   - estimable with >=50 informative strata and >=30 species;
   - seasonal coefficient > 0 and p < 0.05.

## Confirmation gate

`RESERVE_SEASONAL_HEAT_CONFIRMED` only if all three tests pass.

Otherwise:
- `RESERVE_SEASONAL_HEAT_NOT_CONFIRMED`, or
- `RESERVE_LOCAL_GATE_NOT_ESTIMABLE`.

No threshold, cell size, covariate, or outcome definition is altered after opening.

## Interpretation boundary

A positive result would establish species-disjoint replication of seasonal heat association within the same iNaturalist/measurement programme. It would still not prove temperature-induced plasticity, selection, pigment chemistry, or pigmented->white evolutionary direction.
