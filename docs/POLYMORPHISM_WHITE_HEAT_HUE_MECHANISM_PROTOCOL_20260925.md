# Third-cohort white–heat hue-mechanism discriminator — 2026-09-25

Status: post-confirmatory mechanism-discrimination analysis. It does not alter the frozen New Phytologist manuscript, the H2 axis, or the failed legacy BIO5 replication.

## Question

If thermal suppression of anthocyanin/flavonoid pathway activity contributes to some white states, should the third-cohort BIO5 association be stronger when white is contrasted against anthocyanin-associated chromatic hues than when white is contrasted against yellow/orange?

## Frozen hue groups

Using only the four frozen biological coarse states:

- white = `white`
- anthocyanin-associated hue proxy = `red_pink` + `blue_purple`
- yellow/orange comparator = `yellow_orange`

These groups are fixed before BIO5 values are opened for this analysis.

This is a phenotype proxy, not pigment chemistry. Red/pink/blue/purple flowers are not universally anthocyanic, and yellow/orange flowers are not universally non-anthocyanic.

## Inputs

Biological rows:
- immutable third-cohort prospective measurement artifact from workflow run `35177668182`.

Technical control:
- response-blind highlight seal from workflow run `35687421896`;
- use the exact frozen high-clip set;
- no alternate clipping threshold.

Environment:
- WorldClim 2.1, 10 arc-minute BIO5.

## Primary within-species contrasts

After retaining globally classifiable rows, removing the frozen high-clip set, and requiring a valid near-clip metric:

For each species, standardize BIO5 within species and calculate:

- `delta_A = mean(zBIO5 | white) - mean(zBIO5 | red_pink or blue_purple)`
- `delta_Y = mean(zBIO5 | white) - mean(zBIO5 | yellow_orange)`

Eligibility per contrast:
- >=5 white rows;
- >=5 comparator rows.

Minimum support:
- anthocyanin-proxy contrast: >=100 species;
- yellow/orange comparator: >=50 species.

## Primary mechanistic prediction

Prediction A:

`median(delta_A) > 0`

Test:
- one-sided Wilcoxon signed-rank across eligible species;
- also report two-sided p and fraction of species with delta_A > 0.

## Comparator prediction

Yellow/orange is not predicted to show the same heat-linked effect.

Report `delta_Y` with the same summaries, but a non-significant result is not treated as proof of equivalence.

## Within-species specificity contrast

Among species satisfying all three supports:
- >=5 white;
- >=5 anthocyanin-proxy;
- >=5 yellow/orange;

calculate:

`S_i = delta_A - delta_Y`

Prediction:

`median(S) > 0`.

Minimum support for this direct specificity test: >=30 species.

Use one-sided Wilcoxon. If <30 species, mark `specificity_not_estimable`; do not lower the gate.

## Row-level corroboration

For each binary contrast separately, fit species-stratified conditional logistic regression:

`white ~ within-species-z(BIO5) + within-species-z(near_clip_fraction) | species`

The row-level model is corroborative and cannot rescue a failed species-level result.

## Interpretation

A positive anthocyanin-proxy contrast together with a weak yellow/orange contrast and a positive direct specificity contrast would be compatible with the heat × pigment-network model.

It would still not establish:
- actual pigment chemistry of any photograph;
- causal heat selection;
- a universal BIO5 effect;
- pigmented-to-white evolutionary direction;
- that the same molecular node is responsible in all species.

The legacy discovery/reserve BIO5 replication failure remains unchanged.
