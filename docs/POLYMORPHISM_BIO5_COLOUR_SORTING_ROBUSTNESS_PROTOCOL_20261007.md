# FCP BIO5 colour-sorting robustness protocol — 2026-10-07

## Status and chronology

This protocol is frozen after the selected BIO5 colour-sorting signal transported to the third cohort, and before the robustness statistics below are evaluated.

The analysis is a post hoc technical/measurement robustness audit. It cannot upgrade the result to prospective confirmation or local adaptation.

## Target already established before this audit

The selected non-directional target is:

`flower-colour dissimilarity increases with within-species BIO5 difference after geographic distance is partialled out`.

Coarse four-state matched-vertex support was observed in discovery, validation and the fixed third-cohort transport. The third cohort also supported a continuous nine-colour sensitivity. The earlier directional prediction that white should consistently occupy higher BIO5 did not transport across cohorts.

## R1 — continuous nine-colour flower representation across all three cohorts

Recompute the BIO5-sorting statistic using the nine frozen biological palette coordinates rather than the four coarse colour groups.

- discovery/validation: use `palette_count_{white,yellow,orange,red,pink,magenta,purple,blue,bronze}` and normalize rows;
- third cohort: use the frozen `flower_fraction_*` coordinates directly.

For each species, calculate pairwise Jensen–Shannon flower-colour dissimilarity and partial Spearman(BIO5 difference, colour dissimilarity | geographic distance).

Inference uses 199 matched vertex permutations of complete flower-colour rows among fixed coordinates.

**R1 robustness support:** positive equal-species mean and matched-vertex p < 0.05 in discovery, validation and third cohort.

## R2 — matched flower-minus-background control in discovery and validation

The third frozen measured table contains no matched background palette representation, so this control is evaluable only in the original discovery/validation resources.

For each legacy photo use all 12 frozen palette-count anchors for both flower and same-image background. Calculate pairwise flower JSD and background JSD, then use:

`colour differential = flower JSD - background JSD`.

Test partial Spearman(BIO5 difference, colour differential | geographic distance) with the same matched vertex permutation applied jointly to the flower/background observation row.

**R2 robustness support:** positive equal-species mean and matched-vertex p < 0.05 in both discovery and validation.

Failure does not erase the coarse/continuous signal, but it would prevent claiming that BIO5-associated turnover clearly exceeds same-image background structure.

## R3 — same-observer pair control

To reduce observer/camera confounding, restrict the pairwise statistic to photo pairs taken by the **same observer** within a species. Because the acquisition cap is two retained photos per observer per species, each observer contributes at most one pair.

Use the continuous nine-colour flower representation. A species is evaluable only if it has at least five same-observer pairs, nonconstant geographic distance, nonconstant BIO5 difference and nonconstant flower-colour dissimilarity.

Within each species calculate partial Spearman(BIO5 difference, colour dissimilarity | geographic distance) across same-observer pairs.

For the null, independently swap or retain the two colour rows within each two-photo observer block for each permutation, preserving observer membership, coordinates, BIO5 values and each observer's colour-row multiset.

199 deterministic permutations are used.

**R3 sensitivity support:** positive equal-species mean and matched-null p < 0.05 in each cohort with at least 50 evaluable species.

R3 is a stringent sensitivity and is not allowed to rescue a failed R1.

## Hard nonclaims

Even if R1–R3 support the signal, the allowed interpretation remains environmental sorting/temperature-associated turnover. It does not establish:

- BIO5 as the causal variable;
- colour-dependent fitness;
- local adaptation;
- genetic differentiation rather than phenotypic plasticity;
- a common direction of colour response across species.
