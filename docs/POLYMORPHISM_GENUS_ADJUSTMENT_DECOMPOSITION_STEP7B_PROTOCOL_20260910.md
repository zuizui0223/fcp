# Polymorphism genus clustering — Step 7b adjustment decomposition

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-genus-adjustment-decomposition-step7b`
Parent: `analysis/polymorphism-genus-replication-step7` at `d2abafe54043b1f92540be873738793c262f4c55`

## Problem exposed by Step 7

Step 7 found three facts simultaneously:

1. raw genus clustering of photo-derived polymorphism diversity `D` replicated in the species-disjoint reserve;
2. shared repeated-genus mean D values were concordant across discovery and reserve;
3. reserve genus clustering no longer passed after simultaneous adjustment for sampled span and `n_classifiable`.

The third result cannot yet be interpreted as evidence that genus clustering is merely sampling opportunity. In both 100-photo tranches, `n_classifiable` is not raw acquisition effort: every species entered with a fixed photo denominator and `n_classifiable` is the number of photographs that survive the flower-colour classification gate. In the discovery analysis its correlation with D is exactly the sign-reversed counterpart of the previously documented mixed/uncertain-rate association.

Therefore `n_classifiable` can be both a measurement-quality variable and a downstream consequence of difficult/intermediate/diverse colour observations. Conditioning on it may remove biological or measurement-process variation that belongs to the estimand.

Step 7b decomposes the two controls rather than treating `span + n_classifiable` as one undifferentiated confounder.

## Frozen inputs

- discovery: `results/polymorphism_species_attributes_step4_association_20260910/species_analysis_table.csv`
- reserve: `results/polymorphism_spatial_span_adjusted_step6_20260910/reserve_species_span_adjusted_input.csv`
- Step 7 result: `results/polymorphism_genus_replication_step7_20260910/result.json`
- discovery mixed diagnostic: `results/polymorphism_mixed_step2_20260909/result.json`

Expected species counts: discovery 369; reserve 363.

Reserve `reserve_photo_rows` must equal 100 for every retained species. Discovery Step-2 fingerprint must report raw rows per species min=max=100. If either check fails, the fixed-denominator interpretation is not used.

## Statistic

Reuse the exact equal-genus Step-4 / Step-7 clustering statistic and the same deterministic genus definition. Only genera with >=2 species contribute; every repeated genus has equal total weight across its within-genus pairs.

Use 20,000 permutations per test, lower-tail P values, with fixed Step-7b seed family beginning at `20260918`.

## Frozen analyses

### A. Raw finite-sample-corrected replication

Test genus clustering of

`D_unbiased = n_classifiable/(n_classifiable-1) * D`

in discovery and reserve without conditioning on `n_classifiable`.

This addresses the ordinary finite-sample bias of Gini-Simpson diversity while retaining the fixed photo-derived measurement process.

### B. Sampled-span-only adjustment — primary decomposition

Within each tranche separately:

1. rank-transform D and sampled log-span;
2. regress ranked D on intercept + ranked sampled log-span;
3. use the D residual in the unchanged genus clustering test.

Repeat with `D_unbiased`.

This is the preferred test of whether the Step-4 sampled-geographic-opportunity signal is sufficient to explain genus clustering.

### C. Classification-yield-only adjustment — diagnostic

Repeat the same residual clustering after adjusting D only for ranked `n_classifiable`, and then repeat with `D_unbiased`.

This is explicitly a **measurement-process sensitivity**, not the default biological adjustment, because `n_classifiable` is generated after the fixed 100-photo acquisition stage.

### D. Simultaneous adjustment reconciliation

Do not recompute a new primary conclusion from another model family. Read the already frozen Step-7 simultaneous `span + n_classifiable` results and place them alongside A-C.

### E. Cross-tranche D_unbiased genus-mean concordance

For the exact 23 shared repeated genera already identified by Step 7, compute mean `D_unbiased` separately in discovery and reserve and test their Spearman correlation with 20,000 two-sided permutations.

No genus inclusion rule is changed after outcomes are available.

## Decision hierarchy

The sampled-span alternative is considered insufficient to explain genus clustering if the **reserve span-only D or D_unbiased test** retains positive clustering gain and P<0.05.

The classification-yield diagnostic is interpreted separately:

- if span-only passes but n-only fails, state that the taxonomic signal is not explained by sampled geographic span but is entangled with classification success/missingness;
- if both span-only and n-only pass, taxonomic clustering is robust to both controls separately;
- if span-only fails, do not promote an opportunity-independent genus propensity claim.

A successful raw D_unbiased reserve replication shows that ordinary finite-sample bias in D is not sufficient to explain the raw genus result. It does **not** solve non-random classification missingness.

## Claim boundary

No analysis here turns taxonomic clustering into formal phylogenetic signal or genetic determination. `n_classifiable` adjustment cannot by itself distinguish technical classification failure from genuinely intermediate/diverse floral colour states. The surviving claim, if any, must remain about **replicated taxonomic structure in the photo-derived polymorphism metric, with explicit measurement-process uncertainty**.
