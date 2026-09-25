# Legacy white–heat hue replication protocol — 2026-09-25

Status: **post-result species-disjoint replication** opened only after the third-cohort hue discriminator in PR #79 was observed. It is not an untouched confirmatory test and cannot rescue the failed overall BIO5 replication.

## Question

Does the third-cohort pattern

`white vs red/pink+blue/purple: warmer`

but not

`white vs yellow/orange: warmer`

transport to the two original high-depth species-disjoint cohorts?

## Frozen colour groups

Copied unchanged from PR #79:

- white: `white`
- anthocyanin-associated hue proxy: `red_pink + blue_purple`
- yellow/orange comparator: `yellow_orange`

These are phenotype proxies, not pigment chemistry.

## Source rows

Immutable source commit:
`5142f7951af0dde5364bb047a566d67e8c479e51`

Discovery:
`data/derived/global_monte_carlo_measured_photos_v1.csv`
SHA256 `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`

Reserve:
`data/derived/rgfca_reserve_replication_measured_photos_v1.csv`
SHA256 `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`

Environment:
WorldClim 2.1, 10 arc-minute BIO5.

## Per-cohort species-level contrasts

Within each species, standardize BIO5.

`delta_A = mean(zBIO5 | white) - mean(zBIO5 | red_pink or blue_purple)`

`delta_Y = mean(zBIO5 | white) - mean(zBIO5 | yellow_orange)`

Eligibility:
- each contrast requires >=5 white and >=5 comparator rows.

Minimum support:
- A: >=100 species;
- Y: >=50 species.

Prediction copied from the opened third-cohort result:
- `median(delta_A) > 0`.

Use one-sided and two-sided Wilcoxon signed-rank tests across species.

## Direct specificity

For species with >=5 rows in white, anthocyanin-proxy and yellow/orange:

`S = delta_A - delta_Y`

Minimum support = 30 species per cohort.
Prediction: median S > 0.
Do not lower this gate.

## Row-level corroboration

For each binary contrast, fit species-stratified conditional logistic regression:

`white ~ within-species-z(BIO5) | species`

No highlight/near-clip technical metric exists for these legacy row tables, so row-level models are explicitly measurement-system-sensitive.

## Replication labels

`ANTHOCYANIN_PROXY_TRANSPORT_SUPPORTED` only if **both discovery and reserve**:
- have >=100 eligible A species;
- median delta_A > 0;
- one-sided Wilcoxon p < 0.05;
- species-stratified BIO5 beta > 0 and p < 0.05.

`HUE_SPECIFICITY_TRANSPORT_SUPPORTED` only if **both cohorts**:
- have >=30 direct-specificity species;
- median S > 0;
- one-sided Wilcoxon p < 0.05.

Failure of the yellow/orange test is reported descriptively and is not itself a required proof of equivalence.

## Hard boundaries

This analysis:
- is post-result replication, not untouched confirmation;
- uses the same iNaturalist/measurement system;
- lacks the third cohort's direct highlight technical control;
- does not establish chemical anthocyanin identity;
- does not establish causal heat selection;
- does not change the failed overall legacy BIO5 replication verdict.
