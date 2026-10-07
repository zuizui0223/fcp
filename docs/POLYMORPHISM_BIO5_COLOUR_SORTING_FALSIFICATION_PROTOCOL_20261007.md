# FCP BIO5 colour-sorting falsification protocol — 2026-10-07

## Status

This is a post-outcome falsification audit frozen after the three-cohort continuous-colour BIO5-sorting pattern was observed. It cannot upgrade any result to prospective confirmation and cannot alter frozen H1/H2/D-spatial/BIO5 decisions.

## Motivation

The continuous nine-colour BIO5-distance/flower-colour-distance association is positive under matched vertex nulls in discovery, validation and the third species-disjoint cohort. However, validation flower-minus-background and sparse same-observer-pair sensitivities were not supportive. The following tests ask whether the pattern survives representations and controls closer to biological flower-colour state.

## F1 — discrete morph turnover

For each cohort, replace continuous colour distance by a pairwise mismatch indicator: 1 if the two retained observations have different frozen coarse morph labels (white, yellow/orange, red/pink, blue/purple), otherwise 0.

For each species compute partial Spearman(BIO5 absolute difference, morph mismatch | geographic distance). Use 199 deterministic matched vertex permutations of the complete morph-label vector among fixed coordinates. Aggregate by equal-species mean.

F1 is considered cross-cohort supported only if the equal-species mean is positive with matched-null p < 0.05 in discovery, validation and third cohorts.

## F2 — nonwhite-only morph turnover

Repeat F1 after removing all white observations. A species is evaluable when at least 30 nonwhite classifiable observations remain and the pairwise mismatch vector is nondegenerate.

F2 is a sensitivity asking whether any BIO5 sorting extends beyond white-versus-coloured variation. Failure does not negate F1; success would show that the pattern is not confined to the known exposure-sensitive white class.

## F3 — validation flower signal conditional on background

In the validation cohort, use the frozen 12-colour flower and matched-background palette counts. Compute pairwise flower JSD and pairwise background JSD.

For each species compute a rank-based partial correlation between BIO5 distance and flower JSD after residualizing both on geographic distance and background JSD. Permute complete flower-palette rows among fixed coordinates while leaving geography, BIO5 and background fixed. Aggregate by equal-species mean over 199 deterministic permutations.

This is stricter than subtracting flower and background JSD because it asks whether flower-colour turnover contains BIO5-associated information conditional on same-image background turnover.

## F4 — third-cohort highlight exclusion

Use the permanently archived response-blind high-clip set from the prospective third cohort. Map high-clip measurement IDs to photo IDs through the sealed join key, remove those rows, reapply the >=40 classifiable species gate, and rerun:

- F1 discrete morph turnover;
- continuous nine-colour flower-fraction turnover.

F4 asks whether the third-cohort BIO5 sorting survives removal of the exact exposure subset already used in the H2 validity audit.

## Interpretation ceiling

The strongest allowed wording if F1 and F4 are supported and F3 is also positive is that flower-colour turnover is repeatedly associated with warm-season temperature differences beyond geographic distance and is not readily explained by the coarse classifier, the frozen high-clip subset, or matched background variation in the validation cohort.

Even that result would remain **environment-associated sorting consistent with local adaptation or plasticity**, not evidence of adaptive fitness.

## Hard nonclaims

- no fitness measurement;
- no reciprocal transplant or common-garden test;
- no distinction between genetic local differentiation and phenotypic plasticity;
- no causal attribution to BIO5;
- no universal direction of colour change;
- no universal angiosperm rule.
