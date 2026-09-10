# Polymorphism ambiguity bounds — Step 9 protocol

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-ambiguity-bounds-step9`
Parent: `analysis/polymorphism-measurement-missingness-step8` at `815cec0c551c0e225110e90066eb9406c990c092`

## Why Step 9 exists

Step 8 ruled out sampled geographic span plus clearly technical ROI/flip failure as a sufficient explanation of the replicated genus signal and the flower-specific spatial signal. The remaining dominant measurement uncertainty is `not_evaluable_ambiguous_palette_composition`.

Ambiguous palette rate is strongly associated with the observed classifiable-only four-state diversity D, but the earlier bridge-geometry gate did not justify recoding ambiguous rows as intermediate colours. Therefore Step 9 does not impute a preferred morph and does not numerically correct D.

Instead it asks a narrower sensitivity question:

> if every ambiguous-palette row were in fact an unobserved member of one of the same four discrete morph classes, how far could species-level four-state D move, and do the main genus/spatial conclusions survive the two exact endpoint completion rules?

This is a four-state completion sensitivity analysis, not proof that ambiguous rows are latent four-state morphs.

## Frozen source rows

Discovery:
`data/derived/global_monte_carlo_measured_photos_v1.csv`

Reserve:
`data/derived/rgfca_reserve_replication_measured_photos_v1.csv`

Use all 100 fixed acquired photos per species.

Four admitted morphs:

- `white`
- `yellow_orange`
- `red_pink`
- `blue_purple`

Ambiguous rows eligible for completion:

`measurement_status == not_evaluable_ambiguous_palette_composition`

Rows with `not_evaluable_roi_or_flip_gate` or `not_evaluable_no_biological_palette_mass` remain excluded from the four-state completion denominator. They are not reassigned.

## Exact species-level D endpoints

For a species let admitted four-state counts be `c_1,...,c_4`, admitted count `n=sum(c_k)`, and ambiguous-palette count `A`.

For any integer allocation `x_k >= 0`, `sum(x_k)=A`, define

`D(x) = 1 - sum_k ((c_k + x_k)/(n+A))^2`.

### Lower endpoint `D_min4`

The exact minimum of D is obtained by assigning all A ambiguous rows to a currently largest count. This maximizes concentration:

`D_min4 = 1 - [ (c_max + A)^2 + sum(other c_k^2) ] / (n+A)^2`.

Ties are irrelevant because the resulting D is identical.

### Upper endpoint `D_max4`

The exact maximum of D is obtained by allocating ambiguous rows to minimize the sum of squared final counts. Because the objective is separable convex, repeatedly adding one row to a currently smallest count gives an exact integer optimum. Ties are resolved by fixed morph order only for reproducibility; D is invariant to equivalent tie choices.

The current classifiable-only D is retained as the baseline `D_observed`.

The runner must verify its reconstruction of `D_observed` against the frozen species tables before any endpoint downstream analysis is accepted.

## What the endpoints do and do not identify

For every species under the four-state latent-completion assumption, the true completed four-state diversity lies in `[D_min4, D_max4]`.

However, downstream rank correlations and genus clustering are not guaranteed to attain their global extrema when all species are simultaneously set to D_min4 or D_max4. Therefore endpoint downstream tests are called **uniform endpoint stress tests**, not full worst-case partial identification of the paper-level statistic.

No claim may say that all arbitrary species-specific allocations have been exhausted.

## Analysis A — ambiguity interval description

In discovery and reserve report:

- median and 95th percentile of `D_max4 - D_min4`;
- median `D_min4`, observed D, and `D_max4`;
- fraction of species with A=0;
- Spearman correlation of observed D with D_min4 and D_max4;
- Spearman correlation of ambiguity interval width with ambiguous-palette rate.

## Analysis B — genus uniform-endpoint stress test

Use the exact Step-8 genus statistic and exactly the same genus membership.

For each tranche and each response `D_min4`, `D_max4`:

1. rank-transform the response, sampled log-span, and technical_failure_rate;
2. residualize response rank on intercept + ranked sampled span + ranked technical_failure_rate;
3. run the unchanged equal-genus weighted-pair clustering test;
4. use 20,000 lower-tail permutations.

No `n_classifiable` adjustment is used because Step 8 established that it mixes technical and ambiguity-related measurement outcomes; technical failure is controlled separately.

### Genus endpoint criterion

If the species-disjoint reserve shows positive genus clustering with P<0.05 for both D_min4 and D_max4, call the genus signal **uniform-endpoint robust under four-state ambiguity completion**.

If only one endpoint passes, report endpoint sensitivity and do not use endpoint-robust language.

## Analysis C — spatial uniform-endpoint stress test

Reuse the original 999 within-species spatial null arrays, with no new spatial randomization.

For each response `D_min4`, `D_max4`, compute:

`partial Spearman(D_endpoint, spatial response | sampled log-span, technical_failure_rate)`

for:

1. discovery primary spatial rho;
2. reserve primary spatial rho;
3. reserve matched flower-minus-background differential.

For each of 999 frozen spatial-null realizations, recompute the same partial-rank statistic. Directional P value is `(1 + count(null >= observed))/1000`.

### Spatial endpoint criterion

A spatial relationship is uniform-endpoint robust under four-state ambiguity completion only if the adjusted partial correlation is positive and P<0.05 at both D_min4 and D_max4.

The reserve matched-background differential remains the preferred flower-specific test.

## Analysis D — shared repeated-genus endpoint concordance

Reuse the exact 23 genera fixed by Step 7. For D_min4 and D_max4 separately, compute each genus mean in discovery and reserve and their Spearman correlation with 20,000 two-sided permutations.

This is secondary and cannot rescue a failed reserve genus clustering endpoint test.

## Claim boundary

Even if both endpoints survive, the allowed conclusion is only:

> the reported taxonomic/spatial relationships persist under two exact species-level endpoint completion rules spanning the minimum and maximum possible four-state D for ambiguous-palette rows, while holding clearly technical failures excluded.

This does not prove ambiguous rows are four-state morphs, does not exhaust arbitrary species-specific latent allocations, and does not establish formal phylogenetic signal, genetic determination, adaptation, selection, causal range-size effects, population-genetic differentiation, or a shared global boundary.
