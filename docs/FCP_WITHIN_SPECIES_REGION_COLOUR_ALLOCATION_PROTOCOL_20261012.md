# FCP — multi-region photographic colour-allocation transfer (frozen design)

**Date:** 2026-10-12. **Status:** source-exposed retrospective extension, fixed before examining the new region-centered prediction contrasts. This is not a preregistered study or independent biological replication.

## Question

Do jointly measured environmental covariates predict **where colour variation falls among geographic sectors of the same nominal species**, rather than merely its species-wide or genus-average colour composition? This is a photo-based proxy of among-place variation, **not** a classification of biological local coexistence (C) or spatial segregation (S).

Previous C/S comparison in `docs/supporting/jbi_space_time_hypothesis_screen_v1.csv` already tested regional niche-centroid displacement and Gaussian overlap in only **26 species**, both unsupported. This extension does not repeat or rescue those tests. It asks a distinct, directly held-out predictive question on the 20,786 species in the original photographic all36 source.

## Locked source and opportunity-only gate

- Reuse the untouched 2026-10-11 nonlinear transport ZIP, SHA256 `e51547cf9daea2bac68d0ab1fbffe17200cdecd4b810ebb660040df4c40f05a1`. The original photographs and site-environment archives are identified by SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0` and `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`.
- Read `soil_all_oof_predictions.csv.gz` (31,925 classified all36-complete photos, 20,786 nominal species, 4,853 named genera). `soil_quality_oof_predictions.csv.gz` is the exact pre-existing strict-picture sensitivity. The same four models/three exclusion regimes/five deterministic genus folds remain unchanged; never refit or select models on outcomes.
- Primary species have **at least two original classified photos in each of at least two distinct 90-degree longitude sectors**. The gate uses only `inat_taxon_id, genus, sector, n_photos`, not colour or prediction. Observed opportunity support is 493 species / 348 genera in all36; sensitivity `min_photos_per_sector=3` contains 179 species /140 genera. Strict-picture subset with the primary gate has 52 species /46 genera. These subsets are not a representative sample of flowering plants.
- Primary training regime: `buffer500`, which purges unseen genus, target 90-degree sector, and all original photographed species with an original occurrence within 500 km. `sector_purged` and `genus_only` remain complete descriptive sensitivities.
- Retain **all** qualifying sectors for each qualifying species; discard other species only as defined by the opportunity gate. Do not filter by observed colour variation, predicted skill, environmental distance, or confidence interval. Records with fewer than two sectors are ineligible, not zero variation.

## Paired within-species sector contrasts

For species `i`, let `y_is` be the four-class observed photo share in sector `s`, and `p^m_is` be the source-frozen prediction for model `m`. Define equal-sector means `bar(y)_i` and `bar(p^m)_i`, and deviations `Y_is = y_is - bar(y)_i`, `P^m_is = p^m_is - bar(p^m)_i`. Observed colour labels are used only in scoring, never to fit or calibrate these existing predictors.

Define four-output regional allocation loss `L^m_i = mean_s sum_k (Y_isk - P^m_isk)^2`. Define the oracle-species-mean zero-prediction baseline `L^0_i = mean_s sum_k Y_isk^2`. **This baseline is a diagnostic using the test species' overall colour shares, not an operational out-of-sample predictor**. The primary paired endpoint is genus-equal mean `L^nonlinear_geo_i - L^nonlinear_full_i` for all36/buffer500/opportunity gate 2; positive means environmental inputs improve regional allocation, not local adaptation. The directly necessary checks are `L^0 - L^nonlinear_full`, `L^ridge_geo - L^ridge_full`, `L^ridge_full - L^nonlinear_full`, and both models' absolute losses. Include white-vs-other and remaining-hue orthogonal parts as descriptive decomposition, never a substitution for the primary four-hue loss.

Within each genus, species have equal weights after each species' sectors are averaged; primary point estimate gives each genus one vote. Also give species-equal estimates. Use paired whole-genus bootstrap 1,999 draws, fixed seed 2026101201, common resampled genus IDs across models and contrasts, percentile CI95 and CI97.5. The models/folds are fixed: intervals do not cover model refitting or four independent geographic replicates.

## Exact opportunity- and composition-preserving null

For the **one primary endpoint only**, reconstruct integer counts `n_is * y_isk` (verify rounding <= 1e-8 and per-sector sums) from the original saved photo fractions. For each species, pool the original four-class labels across eligible sectors, shuffle complete labels among the **fixed sector slot sizes**, reconstitute the sector fractions, recenter as above, keep all four previously fitted model prediction vectors unchanged, then recompute the primary paired gain. Run **1,999 independent species-conditional permutations**, fixed seed 2026101202, plus-one upper-tail p. This null conditions on the species' observed full colour inventory and exact geographic photo opportunity; it does **not** preserve within-sector observer/time structure or allow causal inference. Verify total class counts, sector sample sizes, fixed predictions, and whole-genus aggregation for each draw.

A second descriptively separate null evaluates observed `mean_i L^0_i` relative to its matched random-labelling distribution, without selecting taxa by their observed colour diversity. Do not use this as the primary p-value.

## Coverage and interpretation stops

Require >=100 species and >=50 genera for inferential bootstrap in the **primary and min3** cohorts. Under this rule, strict-picture 52-species support is reported descriptively with `HOLD_LOW_SUPPORT`, not as a negative replicated test. Report full available support, all 4 regions, original 31,925 complete photo denominator and 54,590 unclassifiable original photo opportunities from upstream; do not equate absence of source colour with white.

We will **not** claim true local coexistence (C), field population turnover (S), genotype frequency, single-plant response, selective colour advantage, pigmentation investment, photosynthesis, pollination, herbivory, reproductive fitness, or population-level adaptation from these photos. Any positive result is a conditional source-exposed prediction pattern subject to camera/photo sampling, geography, lineage and sample-size bias. Any null does not establish an environmental effect of zero.

## Verification and execution

Require duplicate-free `(species,sector,scheme)` identity, equal sectors/targets across the four models and three schemes, exact four-class probability simplex, no NaNs, original source/outcome hashes, consistency of results with an independent pandas/NumPy aggregation, and explicit source model labels. Test synthetic known 4-class contrasts, exact permutation inventory, deterministic bootstrap, direction of each loss difference, shuffled-target leakage guard, and incomplete source/HOLD behavior. Run the actual artifact-derived analysis locally, save paired species-level ledger, all-conditions JSON, bootstrap/permutation receipts, and an ecological interpretation that includes all negative cases. Retain the completed P170 and NP manuscript decisions unchanged.