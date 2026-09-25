# FCP seasonal-heat mechanism protocol — 2026-09-25

Status: **prospectively frozen after the long-term BIO5 association was observed, but before opening month-specific WorldClim temperature values for the third-cohort photographs**.

This is a new mechanism-discrimination test. It does not alter the frozen New Phytologist manuscript, the H2 white/non-white geometry result, or the already opened long-term white-environment result.

## Motivation

The prospective long-term environment test found that measured white records occur, on average, in warmer BIO5 locations within species. The independent transition-direction test did **not** support a recurrent nonwhite-dominant -> white-dominant macroevolutionary asymmetry.

A remaining hypothesis is therefore not irreversible pigment loss, but a recurrent **seasonal heat association**: white states may be expressed or favored during hotter parts of the flowering season, even if white/non-white transitions are macroevolutionarily reversible.

## Data

Biological source:
- immutable third-cohort prospective measurement run `35177668182`;
- globally classifiable rows in the four frozen coarse states only.

Exposure control:
- response-blind third-cohort highlight technical seal run `35687421896`;
- exact frozen high-clip set excluded;
- continuous near-clip fraction retained as a nuisance covariate.

Climate:
- WorldClim 2.1, 10 arc-minute monthly maximum temperature (`tmax`);
- WorldClim 2.1 BIO5 as the already-known long-term warm-extreme covariate.

WorldClim monthly values are climatological month means, not weather observations for the exact year.

## Frozen predictor

For each photograph j at coordinate c and calendar month m:

`
seasonal_heat_anomaly_j =
    tmax(c,m) - mean_{k=1..12} tmax(c,k)
`

This removes the location-specific annual temperature level and asks whether a record falls in a relatively hotter or cooler month at its own location.

## Primary prediction

`
white -> higher seasonal_heat_anomaly
`

The direction is fixed positive before the monthly tmax values are opened.

## Eligibility

After classifiability, highlight availability and frozen high-clip exclusion, a species is eligible for the species-level test when it retains:
- >=5 white rows;
- >=5 non-white rows;
- valid observation month and monthly tmax values.

Minimum primary support: 100 species.

## Primary species-level test

Within each species, standardize the seasonal heat anomaly and calculate:

`
delta_i = mean(z_anomaly | white) - mean(z_anomaly | nonwhite)
`

Primary evidence:
- median delta > 0;
- two-sided Wilcoxon signed-rank p < 0.05 across species.

## Photo-level corroboration

Fit species-stratified conditional logistic regression:

`
white ~ z_within_species(seasonal_heat_anomaly)
      + z_within_species(BIO5)
      + z_within_species(near_clip_fraction)
      | species
`

The seasonal anomaly coefficient must be positive with p < 0.05.

BIO5 is included because its long-term association has already been opened; this model asks whether seasonal heat contributes beyond that spatial climate effect.

## Same-cell local robustness gate

WorldClim 10-arc-minute raster cell is used as an outcome-blind spatial stratum.

Fit a conditional logistic model within `species × WorldClim cell` groups that contain both white and non-white rows:

`
white ~ seasonal_heat_anomaly + near_clip_fraction
      | species × climate_cell
`

Because long-term climate is constant within a raster cell, this directly asks whether white and non-white observations from the same species and same coarse location differ in flowering-month temperature.

This local model is estimable only when:
- >=50 informative species × cell strata; and
- >=30 species are represented among those strata.

## Decision gate

`SEASONAL_HEAT_ASSOCIATION_SUPPORTED` only if all hold:

1. >=100 primary species;
2. species-level median delta > 0;
3. species-level Wilcoxon p < 0.05;
4. species-stratified seasonal coefficient > 0 and p < 0.05;
5. same-cell local model is estimable;
6. same-cell seasonal coefficient > 0 and p < 0.05.

Otherwise the result is:
- `SEASONAL_HEAT_ASSOCIATION_NOT_SUPPORTED_UNDER_THIS_TEST`, or
- `SEASONAL_HEAT_LOCAL_GATE_NOT_ESTIMABLE`.

No threshold, spatial cell size, month definition, or covariate is changed after opening the result.

## Hard nonclaims

A positive result would support repeated seasonal heat sorting/association of the measured white state. It would not by itself establish:
- temperature-induced plastic colour change;
- selection on pigment chemistry;
- anthocyanin downregulation;
- pigmented -> white evolutionary direction;
- actual weather at image date;
- pollinator causation.

A negative result would leave the already observed long-term BIO5 association as a broad geographic sorting signal rather than evidence for seasonal heat response.
