# Polymorphism spatial organization — Step 6 span-adjusted robustness protocol

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-spatial-span-adjusted-step6`
Parent: `analysis/polymorphism-species-attributes-step4b` at `754fdf074fb0c18294c88fa124ead23da98445bb`

## Why this diagnostic is required

Step 5 / Step 5b support a positive association between species-level flower-colour polymorphism and within-species geographic colour organization. Step 4 subsequently found that polymorphism diversity `D` is also positively associated with the outcome-blind sampled geographic span proxy.

That creates a concrete alternative explanation:

> species sampled across larger geographic spans may have more opportunity both to contain multiple colour states and to show positive distance-colour association.

Step 6 asks whether the `D -> stronger within-species spatial organization` pattern remains after explicitly conditioning on sampled geographic opportunity.

This is a **post-Step-4 robustness diagnostic**, not a new prospective primary claim. It cannot convert the exploratory discovery Step 5 into a preregistered result.

## Frozen inputs

### Discovery frame

- Step 4 species table:
  `results/polymorphism_species_attributes_step4_association_20260910/species_analysis_table.csv`
- Step 5 species spatial metrics:
  `results/polymorphism_spatial_subset_step5_20260909/species_membership_and_observed_rho.csv`

Required species count: 369.

Use the already frozen Step-4 geographic-opportunity proxy:

`log1p_span_primary = log1p(maximum_span_km_after_observer_cap)`.

Primary spatial response:

`spatial_observed_rho`.

### Independent reserve frame

- measured photo ledger:
  `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`
- Step 5b species spatial metrics:
  `results/polymorphism_spatial_reserve_step5b_20260909/reserve_species_membership_and_observed_metrics.csv`

Required eligible species count: 363.

For every reserve species, compute an outcome-blind sampled-span proxy from **all 100 fixed reserve photo coordinates**, before using morph/classifiability fields:

`reserve_log1p_span = log1p(maximum pairwise great-circle distance among the species' fixed 100 photo coordinates)`.

Do not select coordinates using morph, measurement status, `global_classifiable`, D, or spatial outcome. The reserve acquisition already fixed the photo set; this calculation only summarizes its geometry.

Primary reserve spatial response:

`rho_primary`.

Flower-specific reserve diagnostic:

`rho_matched_background_differential`.

## Frozen analyses

Use `N_PERM = 20,000` and seed `20260910` plus fixed offsets.

### A. Continuous span-adjusted gradient

For each dataset, rank-transform:

- `D`;
- spatial response;
- sampled log-span;
- `n_classifiable`.

Residualize ranked `D` and ranked spatial response separately on an intercept + ranked log-span + ranked `n_classifiable` by ordinary least squares.

The statistic is Pearson correlation between the two residual vectors, i.e. partial Spearman correlation:

`rho_partial(D, spatial | span, n_classifiable)`.

Permutation test:

1. retain the observed control design and fitted values;
2. permute the ranked-D residual vector across species;
3. correlate the permuted D residual with the fixed spatial residual;
4. use a two-sided `(1 + exceedances) / (N_PERM + 1)` P value.

This is a robustness diagnostic, not a causal mediation analysis.

Run for:

1. discovery `spatial_observed_rho`;
2. reserve `rho_primary`;
3. reserve `rho_matched_background_differential`.

### B. Span-conditioned discrete-polymorphism contrast

Use the already frozen threshold `second_fraction >= 0.10`.

Within each dataset:

1. divide species into five deterministic quantile bins of sampled log-span using rank-based equal-count assignment;
2. compute the observed difference in mean spatial response between `second_fraction >= 0.10` and its complement;
3. within each span bin, permute the binary membership labels while preserving that bin's exact number of polymorphic species;
4. recompute the global mean difference;
5. use a one-sided upper-tail P value because the Step-5 direction is already fixed as `polymorphic > complement`.

Run for:

1. discovery primary spatial response;
2. reserve primary spatial response;
3. reserve matched-background differential.

No threshold search, alternate bin number, successful-subset selection, or replacement of 0.10 is allowed after results are opened.

## Descriptive checks

Report without treating them as new hypothesis families:

- `rho(D, span)`;
- `rho(spatial, span)`;
- `rho(D, spatial)`;
- span range and median;
- classifiable-count range.

The raw `rho(D, spatial)` values must reproduce Step 5 / Step 5b to numerical tolerance before adjusted results are accepted.

Expected anchors:

- discovery `rho(D, spatial) = 0.08921325988911004`;
- reserve primary `rho(D, spatial) = 0.101601` approximately (exact value is checked against the retained result JSON).

## Decision rule

The span-opportunity alternative is **not sufficient to explain the association** only if:

1. discovery continuous partial correlation remains positive, and
2. reserve primary continuous partial correlation remains positive, and
3. the reserve matched-background continuous partial correlation is not directionally reversed.

P values determine strength of evidence; do not define a hidden all-or-none biological truth from one threshold.

The discrete span-conditioned contrasts are supporting diagnostics. They may strengthen or weaken the threshold-based Step-5 wording but do not override the continuous result by themselves.

If adjusted discovery and reserve effects both collapse toward zero or reverse, Claim 2 must be reframed as **span/opportunity-dependent spatial structure** rather than an intrinsic polymorphism-spatial association.

## Claim boundary

Even if the adjusted association survives, do not claim:

- adaptation;
- local selection;
- population-genetic differentiation;
- causal effects of range size;
- a shared global flower-colour boundary;
- true biological species range size.

The allowed conclusion is narrower:

> within these fixed photo-derived frames, the association between flower-colour polymorphism and within-species geographic organization is not fully accounted for by sampled geographic span and classifiable-photo count.
