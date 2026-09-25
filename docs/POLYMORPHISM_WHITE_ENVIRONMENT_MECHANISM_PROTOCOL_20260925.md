# FCP white-environment mechanism protocol — 2026-09-25

Status: **prospectively frozen before opening WorldClim values for the third-cohort photo coordinates**.

This is a separate mechanism-analysis line. It does not alter the frozen New Phytologist manuscript, H2 estimand, H2 verdict, or submission package.

## Biological question

Does the recurrent white-versus-nonwhite flower-colour axis show a common environmental sorting direction within species?

The target is not the evolutionary transition direction. The frozen H2 statistic is sign-invariant and establishes a recurrent achromatic–chromatic axis, not pigmented -> white ancestry. Here, “white-environment association” means that, within a species containing both measured white and non-white photographs, the white photographs occupy systematically different long-term environments.

## Data

Primary biological source:
- third-cohort prospective measurement artifact from workflow run `35177668182`;
- 49,900 frozen rows from 499 species;
- only rows passing `global_classifiable` and assigned to white, yellow/orange, red/pink, or blue/purple are eligible.

Technical exposure source:
- response-blind third-cohort highlight technical seal from workflow run `35687421896`;
- the previously frozen high-clip set is reused exactly;
- no clipping threshold is re-estimated.

Environmental source:
- WorldClim 2.1, 10 arc-minute baseline rasters;
- BIO5: maximum temperature of the warmest month;
- BIO14: precipitation of the driest month;
- mean monthly solar radiation across the 12 WorldClim solar-radiation layers.

WorldClim values are long-term environmental context, not weather at the date of each photograph.

## Pre-specified hypotheses

Three mechanistic predictions are tested as a finite family.

### E1 — heat / pigment suppression

Prediction: within species, white records occur in warmer locations.

Operational prediction:
- white-minus-nonwhite standardized BIO5 contrast > 0.

### E2 — reduced drought stress

If pigmentation is favored under drought stress through pleiotropic protective functions, loss/reduction of pigmentation should be less costly in wetter dry-season environments.

Operational prediction:
- white-minus-nonwhite standardized BIO14 contrast > 0.

### E3 — pigment photoprotection

If pigmentation is favored under chronically high radiation, white records should occur in lower-radiation locations.

Operational prediction:
- white-minus-nonwhite standardized mean solar-radiation contrast < 0.

These are competing or partially compatible ecological filters. No single result is allowed to redefine the other predictions after opening outcomes.

## Primary inferential unit

Species are the primary comparative replicate.

For environmental variable x and species i:

`
delta_i = mean(z_within_species(x) | white) - mean(z_within_species(x) | nonwhite)
`

A species contributes only if, after technical availability and response-blind high-clip exclusion, it retains:
- at least 5 white rows; and
- at least 5 non-white rows.

A minimum of 100 estimable species is required for each mechanism variable.

Primary evidence:
1. median delta has the pre-specified sign;
2. two-sided Wilcoxon signed-rank test across species;
3. Holm correction across E1–E3.

The species-level effect is the main evidence because the biological claim is a recurrent direction across species, not merely a photo-level association driven by large species.

## Row-level corroboration

A conditional logistic regression stratified by species is fit to the same primary panel:

`
white ~ within-species-z(environment) + within-species-z(near_clip_fraction) | species
`

The environmental coefficient must have the same pre-specified sign and p < 0.05 to corroborate a mechanism.

The row-level model cannot rescue a failed species-level test.

## Mechanism gate

A mechanism is labeled `gate_pass` only when all conditions hold:
1. >=100 species are estimable;
2. median species-level delta has the pre-specified sign;
3. Holm-adjusted species-level Wilcoxon p < 0.05;
4. conditional-logistic coefficient has the same sign;
5. conditional-logistic p < 0.05.

Otherwise it is `not_supported_under_this_test` or `not_estimable`.

## Exposure control

Primary inference:
- requires a valid near-clip value;
- removes the exact response-blind high-clip set frozen on 2026-09-22;
- retains near-clip fraction as a continuous within-species covariate in the row-level model.

No alternate clipping threshold is allowed.

## Exploratory result excluded from confirmation

Before this protocol was written, a latitude/date-derived extraterrestrial solar-radiation calculation was inspected. It produced a small positive white association in a photo-level species-stratified model, but the species-level direction test was weak.

Because that exploratory value was already seen, it is not part of the prospective environmental mechanism family and cannot be used to define, reverse, or rescue E1–E3.

## Hard nonclaims

A positive environmental gate does not establish:
- pigmented -> white evolutionary transition direction;
- a mutation or anthocyanin-loss mechanism;
- pollinator causation;
- causal selection rather than environmental sorting;
- weather conditions experienced on the image date;
- an artifact-free white state.

A null result for E1–E3 does not exclude local edaphic, microclimatic, pollinator, demographic, or genetic mechanisms.

## Pollinator line

Pollinator causation is deliberately separate. It should be tested only with independent plant–pollinator interaction or occurrence data after a coverage gate is frozen. Pollinator text in the FCP literature is hypothesis-generating evidence, not independent confirmation.
