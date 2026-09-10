# Polymorphism Step 6 — span-adjusted robustness protocol

Date: 2026-09-10 JST

## Status

This is an explicitly **post-outcome robustness diagnostic** motivated by the completed Step 4 result that flower-colour polymorphism D is positively associated with sampled geographic span. It does not alter or reopen the confirmatory Step 4 six-family analysis.

The diagnostic is frozen before any span-adjusted Step 5/5b result is opened.

## Threat being tested

Step 5/5b found that species with greater discrete flower-colour polymorphism have stronger within-species geographic colour organization. Step 4 then found D to be positively associated with sampled geographic opportunity. A possible nuisance explanation is therefore:

> species sampled across broader geographic extents may make both polymorphism and geographic colour structure easier to detect.

Step 6 asks whether the D–spatial-organization relationship remains after directly conditioning on sampled geographic span and classifiable sample size.

## Frozen source tables

### Discovery

- `results/polymorphism_spatial_subset_step5_20260909/species_membership_and_observed_rho.csv`
  - fixed D, second-morph threshold membership, and species-level observed spatial rho;
- `results/polymorphism_species_attributes_step4_preflight_20260910/covariate_panel_preoutcome.csv`
  - frozen `log1p_span_primary`, `n_classifiable`, genus, and family.

### Reserve

- `results/polymorphism_spatial_reserve_step5b_20260909/reserve_species_membership_and_observed_metrics.csv`
  - reserve D, threshold membership, classifiable n, and reserve species-level primary spatial rho;
- `data/frozen/rgfca_reserve_geometry_audit_v1.csv`
  - reserve-photo-tranche `maximum_span_km`; use `log1p(maximum_span_km)`.

No new image acquisition, taxonomy acquisition, environmental extraction, or trait backfilling is allowed in Step 6.

## Primary Claim-2 robustness family

For discovery and reserve separately:

1. rank-transform D, spatial rho, sampled span, and `n_classifiable`;
2. compute the partial correlation between rank(D) and rank(spatial rho) after linear residualization on rank(span) and rank(n);
3. test the partial correlation with 20,000 Freedman–Lane-type residual permutations of the ranked spatial-rho response under the reduced span+n model;
4. use a two-sided permutation p-value;
5. Holm-adjust the two p-values across discovery and reserve.

The central spatial claim is classified as `SPAN_ROBUST_CONTINUOUS` only if both adjusted partial correlations remain positive and both Holm-adjusted p-values are <0.05.

## Threshold sensitivity family

For the pre-existing second-morph >=10% indicator, discovery and reserve separately:

- regress raw species spatial rho on the fixed indicator plus rank(span) and rank(n);
- report the indicator coefficient as an adjusted spatial-rho difference;
- test with 20,000 Freedman–Lane residual permutations under the reduced span+n model;
- Holm-adjust across the two tranches.

This threshold family is secondary; failure does not overturn a surviving continuous adjusted gradient.

## Step-4 reserve recurrence family

Use reserve D and reserve-photo span to test whether the two Step-4 supported attributes recur in the reserve measurement tranche:

1. `sampled_geographic_span`: Spearman rho(D, log1p reserve maximum span), 20,000 species-label permutations, two-sided;
2. `genus_taxonomic_clustering`: use the same equal-genus-weighted mean absolute within-genus D difference as Step 4, with genera having >=2 species, and 20,000 label permutations.

Holm-adjust these two reserve recurrence p-values. Call a Step-4 attribute `RESERVE_RECURRENT` only when its reserve effect has the same direction and Holm p<0.05.

Because discovery and reserve may share species, reserve results are described as **reserve-photo-tranche recurrence**, not automatically as species-disjoint replication.

## Span-adjusted genus diagnostic

For discovery and reserve separately:

- residualize raw D on rank(span) and rank(n);
- apply the same equal-genus-weighted clustering statistic to D residuals;
- test with 20,000 species-label permutations;
- repeat with finite-sample unbiased Simpson diversity `D_unbiased = D*n/(n-1)` as a sensitivity.

This is a robustness diagnostic and is not folded into the original Step-4 Holm family.

## Additional fixed descriptives

Report, without making them new claim gates:

- discovery/reserve species overlap;
- rho(span, spatial rho) in each tranche;
- rho(n_classifiable, spatial rho) in each tranche;
- raw versus adjusted D–spatial effect sizes.

## Interpretation ceiling

A surviving adjusted result supports only that the observed relationship between flower-colour polymorphism and within-species geographic colour organization is not explained by the measured sampled-span and classifiable-n opportunity proxies.

It does **not** establish adaptation, causal range-size effects, population-genetic differentiation, environmental selection, or a globally shared colour boundary.
