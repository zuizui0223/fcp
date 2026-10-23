# FCP: continuous nine-hue photographic geography across scales (2026-10-23)

**Status:** implemented as a new, post hoc source-locked analysis branch. Numerical results from the entire original 447,500-photo source must not be claimed until the GitHub Actions authentic-data job succeeds and its generated summary is inspected.

## Why this analysis is needed

The original FCP P2 result is species-equal mismatch depletion for four predefined coarse photographic labels at 50 km: 3,600 species with any local pair (+0.016116, permutation P=0.005); stricter >=40 classifiable photos and >=30 local pairs: 1,761 species (+0.017614902). New colour-specific 50-km tests found positive excesses for all four labels, including month-matched shuffles. None establishes two, three or four discrete biological morphs.

The original three species-disjoint deep cohorts also supported continuous nine-hue Jensen-Shannon divergence (JSD) geographic depletion: discovery +0.015185, validation +0.013664, prospective third +0.013008 (P=0.005 each). These older eligible-subset tests were not a 4,475-taxon complete-source continuous spatial refit.

## Exact new estimand

For a photographed species i and radius r:

    T(i,r) = mean JSD of all spatially-valid conspecific pairs
             - mean JSD of spatially-valid conspecific pairs within r.

The complete image-level nine-hue vector is row-normalized. The distance is base-2 JSD *divergence* (not sqrt JSD), exactly following the earlier continuous method. Radius options: 5, 25, 50, 100 and 250 km. A species must have >=40 original classifiable photos and >=30 valid pairs within 50 km; each radius additionally requires >=30 local pairs. Positional accuracy must be 0-5000 metres. No missing hue or non-estimable species is silently imputed as zero.

For each taxon, the complete set of photo palettes and original positions are preserved:

1. Unrestricted photo-vertex shuffle: shuffle complete nine-channel vectors among within-species fixed observation locations.
2. Calendar-month conditional shuffle: shuffle complete vectors only within the original species-calendar-month cells. Controls photographed monthly hue composition, not true flowering stage or year-month.
3. Colour-blind observer thinning: retain one photographed record per observer and refit both nulls instead of reusing the original null.

There are 199 fixed-source null worlds per condition, plus-one upper-tail P, and exploratory Holm correction over five radii within each null family. A side table gives only **observed photographic colour-label counts** (one through four, and >=5 photographs per colour), not estimated inherited morph counts. The analysis does not force K=2 clustering.

**Source integrity gate:** recompute the original 50-km four-state baseline at strict support, exactly 1,761 species and +0.0176149021685411. If this fails, abort instead of changing criteria.

## Interpretation and failure boundaries

A positive continuous effect would show that photographed geographic similarity survives replacing four hard categories with continuous composition differences. A positive calendar-month-null effect would exclude purely month-composition reallocation as a sufficient explanation under the chosen null. These tests do not distinguish discontinuous geographic patches from smooth clines, identify genetic morph counts, natural selection, true populations or inheritance.

The earlier <=5 km, same-year/month, distinct-observer strict local-photograph null remains unsupported (P=0.300 / 0.275 under its two opportunity gates). A positive 50-km estimate cannot negate that result. The source's 220,393 unclassified photos also require continued missingness caution; illumination, organ identity, horticultural material and taxonomic error may bias the observed morphology.

## Reproducibility routes

- Analysis: scripts/analysis/fcp_all4475_continuous_spatial_20261023.py
- Locked loader: scripts/analysis/fcp_full4475_hypothesis_census_20261021.py
- Locked continuous distance: scripts/analysis/run_polymorphism_distributed_polymorphism_robustness_20261007.py
- Original photo release: fcp-all4730-original-method-actual-colour-results-20261017
- CI: .github/workflows/fcp-all4475-continuous-spatial-20261023.yml
- Output: summary.json, original_per_taxon.csv, one_per_observer_per_taxon.csv (no photo IDs, coordinates, observer IDs)

The synthetic test distinguishes genuine geographic sorting within one month, geography fully confounded with month, and continuous geographical colour differences *within one coarse colour category*. Source-run actual estimates remain pending until authenticated Actions execution and receipt verification.
