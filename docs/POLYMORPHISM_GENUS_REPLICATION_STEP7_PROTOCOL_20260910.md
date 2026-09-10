# Polymorphism genus clustering — Step 7 independent replication and opportunity control

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-genus-replication-step7`
Parent: `analysis/polymorphism-spatial-span-adjusted-step6` at `316a09858359fc151468056aab874c406bb9f05d`

## Why this test exists

Prospectively frozen Step 4 found two Holm-supported species-attribute results in the 369-species discovery frame:

1. flower-colour polymorphism diversity `D` increases with sampled geographic span;
2. `D` is more similar within genera than expected under unrestricted species-label permutation.

The genus result is taxonomic clustering, not phylogenetic signal. It also remains vulnerable to two alternatives:

- closely related species may have similar sampled geographic opportunity;
- closely related species may have similar measurement/classification yield.

An independent 363-species reserve photo tranche is already available and was species-disjoint from the discovery tranche. Step 7 asks whether the genus-level pattern recurs there, using the **same clustering statistic**, and whether it survives explicit adjustment for sampled span and `n_classifiable`.

This is a post-Step-4 replication/robustness analysis. The reserve result is not allowed to change the statistic, genus definition, minimum group size, or permutation direction.

## Frozen inputs

### Discovery

- `results/polymorphism_species_attributes_step4_association_20260910/species_analysis_table.csv`
- `results/polymorphism_spatial_span_adjusted_step6_20260910/discovery_species_span_adjusted_input.csv`

Expected species: 369.

### Reserve

- `results/polymorphism_spatial_span_adjusted_step6_20260910/reserve_species_span_adjusted_input.csv`

Expected eligible species: 363.

The reserve table already contains outcome-blind sampled span computed from all 100 fixed reserve photo coordinates per species.

## Genus definition

Genus is the first token of the accepted frozen binomial species string. No taxonomic remapping or synonym rescue is performed after outcomes are available.

Only genera represented by at least two species contribute to the clustering statistic.

The test is evaluable only when there are at least 10 repeated genera and at least 30 species inside repeated genera. If the reserve fails this deterministic gate, replication is unavailable rather than rescued by singleton pooling or family-level replacement.

## Exact clustering statistic

Use the exact Step-4 statistic.

For every repeated genus `g` with `m_g >= 2`, enumerate all within-genus species pairs. Give every genus equal total weight and every pair within that genus equal share:

`w_pair = 1 / (number_of_repeated_genera * choose(m_g, 2))`.

Then

`W = sum(w_pair * |z_i - z_j|)`.

Smaller `W` means stronger within-genus similarity.

Primary raw outcome `z = D`.

Permutation null: permute `z` without replacement across all species participating in repeated genera, leave genus membership and pair weights fixed, and recompute `W`.

Use 20,000 permutations. Lower-tail P value:

`(1 + number(W_null <= W_observed)) / 20,001`.

## Tests frozen before reserve genus association is computed

### T1 — discovery exact reproduction

Recompute the Step-4 raw-D genus test and require numerical reproduction of:

- repeated genera = 61;
- species in repeated genera = 169;
- `W_observed = 0.1783075360574797`;
- clustering gain approximately `0.20296304739768833`;
- Step-4 lower-tail P approximately `0.0022998850057497125` under its frozen seed.

The runner may use a separate fixed Step-7 seed, so exact P equality is checked by rerunning the Step-4 seed once and the Step-7 replication seed is reported separately.

### T2 — reserve raw-D replication

Apply the identical statistic to reserve `D`.

Primary replication criterion: positive clustering gain and lower-tail P < 0.05.

### T3 — sampled-opportunity / measurement-yield adjusted clustering

Within discovery and reserve separately:

1. rank-transform `D`, sampled log-span, and `n_classifiable`;
2. regress ranked D on intercept + ranked sampled log-span + ranked `n_classifiable`;
3. use the residual as `z` in the same equal-genus pair statistic;
4. permute residuals among species in repeated genera using the same lower-tail test.

This asks whether genus similarity remains after removing monotone associations with sampled geographic opportunity and classification yield. It is not a causal model.

Predeclared sensitivity: replace `D` by `D_unbiased = D * n/(n-1)` before the same rank-residual procedure.

### T4 — cross-tranche genus-mean concordance

This is secondary. Use only genera represented by at least two species in discovery **and** at least two different species in reserve. For each shared repeated genus, compute mean D in each tranche and Spearman correlation between tranche means.

The test is evaluable only if at least 10 shared repeated genera exist. If evaluable, use a two-sided 20,000-permutation test that permutes reserve genus means across the fixed shared-genus labels.

This test asks whether genus ranking of polymorphism tendency transfers across different species sets. It is not required for the within-reserve clustering replication criterion.

## Descriptive measurement controls

Report, without promoting them to primary hypotheses:

- genus clustering of `n_classifiable`;
- genus clustering of sampled log-span;
- Spearman `rho(D, n_classifiable)` in each tranche;
- Spearman `rho(D, sampled log-span)` in each tranche.

These diagnose whether the taxonomic grouping is also strongly structured in measurement opportunity.

## Decision hierarchy

1. **Raw taxonomic replication** requires T2 clustering gain > 0 and P < 0.05.
2. **Opportunity-robust taxonomic replication** additionally requires adjusted reserve clustering gain > 0 and P < 0.05.
3. Discovery adjusted clustering is reported as robustness of the original Step-4 result, not as a second independent replication.
4. Cross-tranche mean concordance is secondary and cannot rescue a failed T2/T3.

If raw reserve clustering fails, keep the Step-4 discovery result as a single-frame taxonomic pattern and do not headline genus propensity.

If raw reserve clustering passes but adjusted reserve clustering fails, describe the taxonomic pattern as partly opportunity-dependent.

If raw and adjusted reserve clustering both pass, the paper may state that flower-colour polymorphism shows genus-level taxonomic clustering across species-disjoint photo tranches and that the pattern is not accounted for by sampled span and classification yield alone.

## Claim boundary

Even a successful Step 7 does not establish:

- phylogenetic signal in the formal comparative-method sense;
- genetic determination of polymorphism;
- adaptation;
- causal effects of genus membership;
- population-level polymorphism prevalence across all angiosperms.

Allowed language remains **genus-level taxonomic clustering of photo-derived species polymorphism diversity**.
