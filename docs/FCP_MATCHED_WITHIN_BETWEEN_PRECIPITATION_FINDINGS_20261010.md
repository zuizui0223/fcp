# Same-photograph colour response × precipitation: intra- versus inter-specific predictive structure (2026-10-10)

## Key source-verified result

**Original tested observations are now aligned: 14,776 identical original flower photographs, 3,802 source species, 927 source genera, 123 distinct source geographic cells.** Source 85,337 species×cell photographs, 39,075 classified four-colour observations and 38,968 real-site climate-complete classified observations remain frozen, matched by original taxon, observation and photo IDs to 100,543 source photo-site environmental identities.

Both hierarchical tests use the SAME exact heldout photograph's white/yellow-orange/red-pink/blue-purple source classifier outcome, original absolute latitude, longitude sine/cosine, site elevation, BIO1 and BIO5, then add **precipitation BIO12 and seasonality BIO15**, assessing its multiclass Brier improvement. Whole original geographic cells are excluded from training and test focal taxa are excluded entirely from the *between-species* predictor's training across **five geographic cell folds × five source-taxon folds**. In contrast, the *within-species* predictor learns a source species-colour intercept from other training cells of that same species. Between-source-genus training requires >=3 different source species. These training information sets do not estimate identical slopes or causal effects.

| Same heldout 14,776 original photo outcome | Model source and evaluated numerator | Brier baseline GEO+BIO1/BIO5 | Add BIO12/15 | Rain Brier improvement | Species-cluster 95% CI | Original geography-cell 95% CI |
| --- | --- | ---: | ---: | ---: | --- | --- |
| **Within same source species** | Same named species, photos in OTHER train cells only | 0.395477 | 0.395150 | **+0.00032751** | **[+0.00021639,+0.00044641]** | **[+0.00019838,+0.00046154]** |
| **Between congeneric species** | Other genus members, test target species excluded from ALL train photos | 0.499313 | 0.499317 | **−0.00000437** | **[−0.00045951,+0.00043326]** | **[−0.00037655,+0.00044280]** |

**Important counterexample to simplistic shared-effect claim:** the conditional rain block contributes a small positive predictive improvement within original named species on this exact common photo population; no measurable precipitation increment beyond temperature/geography is observed when attempting to predict *the same photographs* from other species in their genus.

However, the **paired difference of the two precipitation improvements** is only +0.00033187 with species-cluster 95% CI **[−0.00008150,+0.00072693]**, genus CI [−0.00008260,+0.00074304], and original photo geographical-cell CI [−0.00010623,+0.00070626]. All include zero. Therefore the appropriate inference is ***within predictive evidence on the common source cohort and unresolved scale contrast***, **not** “within stronger than between” as a statistically established or evolutionary causal contrast.

## Source population and scale comparability

- The original *separate* between-species global congeneric one-photo study involved 13,307 source species/1,417 genera; its BIO12/15 increment +0.001367 was obtained on different source photos with different taxon eligibility and model geometry. That earlier exploratory signal is not overwritten by this smaller but better-aligned matched-source model.
- The original *separate* within-species 250/1000km true photo-geodesic low-rank spatial field studies involved 20,546 source heldout photos from 5,127 species; their weather gains were +0.000433/+0.000335 for precipitation on the climate cohort. Here shared-photo sample sizes dropped to 14,776 photos/3,802 species/927 genera because both training strategies must be possible.
- Same photo response/independent geographical test folds and same covariate definitions improve comparability, but training inputs differ by design, and this model only includes latitude/longitude/elevation, **not** the exact fine-distance spatial random field or real LCVP dated patristic branch covariance. Species fixed intercept is not an independently identified phylogeny.
- True-source spatially matched 50/100km genus photographic label-shuffle nulls from PR #141, within-species all-variable nearest-pair maxT from PR #147, and directly tipped real phylogenetic covariance results from PR #144 remain negative/HOLD for common robust climate selection. This matched-photo incremental model does not rescue them.
- No individual genotypes, trait plasticity measurements, true pollinator visits, seed fitness, two-blind-human botanical colour labels or independent true-color ground truth were supplied. Source photo colours are the previously measured automated four-state outcomes.

## Reproducibility and status

[Successful original-source GitHub Actions 38048277195](https://github.com/zuizui0223/fcp/actions/runs/38048277195) passed **6 synthetic test guards**, recovered the two original Actions source artifacts, matched exact original taxon/photo/observation IDs and enforced both species- and geographic-fold train/test leakage guards. Full numerical original workflow artifact **11668317965**, ZIP SHA256 `759189bcaf8bcf5c2f4bc980da1a57022782939ed8e25f37d9d4a1a752a0119c`. The source stdout exact output is persisted under [result.json](../results/fcp_matched_within_between_precipitation_20261010/result.json) with run provenance.

**Paper-facing decision:** both scales have credible source-observational *potential* environment predictability in separate samples, but there is **no positive, statistically robust single-source same-photo precipitation commonality** in this matched conditional test. At the same time, the difference between levels is not supported as nonzero. Keep "candidate without robust shared identification" in cross-scale discussions. This is a retrospective draft-PR follow-up, and the original *New Phytologist* main, H1/H2 and unopened future 2,000+730 species are unchanged.
