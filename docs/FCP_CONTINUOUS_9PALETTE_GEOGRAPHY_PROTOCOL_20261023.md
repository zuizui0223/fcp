# FCP — continuous nine-palette photographic geography (2026-10-23)

**Status: exploratory source-locked analysis protocol.** The frozen four-state P1/P2 results are unchanged. Source-authenticated numerical outcomes must not be claimed until the real-photo workflow finishes successfully.

## Ecological question

The existing 447,500-photo FCP evidence shows modest but repeatable 50-km geographic organisation of the four coarse photographed colour categories, including category-specific calendar-month and observer-sensitive spatial excess. These are **not** confirmed hereditary morphs or independently sampled natural populations.

The new test asks whether complete **continuous nine-hue flower-palette compositions** also show geographic similarity; whether it persists after calendar-month matching; and whether there is **additional geographic organisation within the original coarse labels**. This is a test of photographed phenotype geometry, not spatial DBSCAN clusters or genetic polymorphism.

## Immutable input and estimator

- 4,475 taxa with exactly 100 original source photos each (447,500 total), from the October 17 authenticated release.
- Fixed original compressed photo CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6.
- The 227,107 originally classifiable photo outcomes carry nine palette fractions for white/yellow/orange/red/pink/magenta/purple/blue/bronze. The other 220,393 unclassified photos stay unclassified; no hue or petal reflectance is imputed.
- Source continuous color = the frozen ROI palette composition. **Not** independently calibrated CIELAB, true petal reflectance, pigment concentration, or genotype.
- Pairwise continuous dissimilarity = base-2 Jensen–Shannon divergence of entire nine-vectors (0–1, **not its square root**).
- For each species and radius, observed depletion = species-wide mean photo-pair JSD minus mean photo-pair JSD within radius. Positive implies locally more similar measured colours.
- Radii: 5, 25, 50, 100, 250 km. Coordinates with recorded precision worse than 5,000 m are excluded. The shortest radius is especially vulnerable to location uncertainty.
- Each species has equal aggregation weight; pairs are not naively pooled across taxa.
- Broad tier: >=1 conspecific valid pair within 50 km, >=1 pair at the tested radius. Strict tier: >=40 originally classified photos and >=30 50-km pairs; other tested radii also require >=30 pairs. Observer-thinned strict tier additionally requires >=40 retained classified photos.

The analysis must exactly reproduce BOTH historical 50-km four-state totals before it can report continuous results:
(1) 3,600 taxa, mean difference +0.016115883798023662 (any qualifying pair);
(2) 1,761 taxa, mean +0.0176149021685411 (strict 40 photos/30 pairs).

## Three species-conditioned nulls (199 permutations each)

1. **Full species shuffle:** retain the complete nine-palette multiset and original geographic sampling locations, shuffle whole vectors among photo vertices within a species.
2. **Calendar-month shuffle:** retain the same data and the vector composition within each original calendar month. Contrasting colours photographed in different months alone cannot create a positive month-adjusted geographic result.
3. **Calendar-month × original coarse-label shuffle:** retain the exact original coarse category at each original vertex as well as month. Shuffle continuous vectors only within the same taxon, month and original coarse category. Any excess is within-label photographic hue/shade geography, not evidence of a new class.

For each radius/null, compare the species-mean observed statistic against 199 corresponding species-mean null worlds (plus-one Monte Carlo upper P). Report photo/taxon opportunities, per-genus clustered bootstrap 95% intervals of mean excess (499 draws), Holm across five radii and Holm across the full 15 null × radius family, separately by strict/broad and by observer-thinned/unthinned analysis. Labels representing 1–4 observed categories remain descriptive and cannot be renamed cluster/morph number without mixture-model selection and blind verification.

## Controls, deliverables, and limits

Script: scripts/analysis/fcp_all4475_continuous_9palette_geography_20261023.py.
Synthetic checks: tests/test_fcp_all4475_continuous_9palette_geography_20261023.py.
Continuous result: authenticated workflow artifact with summary.json plus species-only metrics, with **no coordinates, observer IDs or photo IDs**.

The synthetic tests must distinguish (a) within-month true spatial contrasts, (b) contrasts explained entirely by calendar month, (c) contrasts explained entirely by existing four-state labels, and (d) within-one-label hue differences. A passing synthetic test is not an ecological effect estimate.

This does not control true petal identity, horticultural versus wild origin, within-month flower age, regional camera/observer colour calibration, genotype, population boundaries, pollinator fitness or pigment maintenance cost. It cannot overturn the older unsupported **<=5-km same-year/month distinct-observer** test (P=0.300 / P=0.275), which uses different restrictions. Thus even a positive multiscale JSD result is photographic spatial structure, not inherited local polymorphism or adaptation.

The next biological stage after verified all-species continuous results is to test whether color distributions are truly multimodal versus continuous/clinal, using species-specific model selection and photo-level blind calibration before relating any inferred K to geography.
