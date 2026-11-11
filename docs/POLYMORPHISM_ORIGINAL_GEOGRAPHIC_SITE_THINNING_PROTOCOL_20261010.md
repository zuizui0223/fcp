# Source-frozen geographical hard-core photo thinning (2026-10-10)

## Question and source population

The previously supported 50-km FCP photographed colour-local-depletion result survives replacement of repeated 50-km photo pairs with distinct photograph endpoints (PR #152). However, even different original photographs may originate from geographically near-identical positions and could represent the same plant or identical photographic microenvironment. The question here is whether the source-observed colour structure survives **removal of closely colocated photograph observations from each species before calculating the 50-km spatial response and null**.

The source remains exactly the frozen three high-depth species cohorts, byte-checked SHA256. Restrict the analysis first to the same original colour-photo locally evaluable species, original >=40 classifiable rows and >=30 photograph pairs at <=50 km, expected 166/181/204 species with observed original mean depletion +0.020529254583812922/+0.018672971642749295/+0.01468491968437185. All input photographs, biological category labels, metadata and high-depth selection remain immutable.

## Before-source-execution frozen design

1. **Outcome-blind hard-core point thinning:** For each original species and each geographic minimum spacing 1/5/10 km, use deterministic random photo order with a seed depending only on original species/cohort, radius and realization. Keep a source photo only if its actual photographed geodesic coordinates are strictly more than the spacing from **every** previously selected source photo. Never use morph, pigment, original local depletion or colour classification success to choose which classifiable photo survives. The original classifiability gate is fixed upstream.
2. **Twelve sample-thinning realizations per eligible species and spacing.** A retained realization requires >=20 classifiable photo sites remaining and >=10 pair edges between >=separated photographed sites that are <=50 km apart. Species remain analyzable at a radius/policy if >=8 of 12 realizations satisfy this geography-only opportunity gate. Otherwise count as underidentified and never impute zero local depletion.
3. **Primary and photographer-control estimands:** Recompute on the thinned photo sample the standard four-state specieswide mismatch fraction minus 50-km photo-pair mismatch fraction. The main mode includes all retained photo pairs, and stricter mode uses pairs whose original photographer ID is nonmissing and different. Both have one equal vote per original species after averaging the admissible thinning realizations.
4. **Matched photographic null:** For each retained original photo subset, use 199 complete vertex (source photo) label shuffles preserving the retained sample's exact four-state colour composition and fixed geography. Average shuffle statistics within species across the admissible realizations before averaging across species, **not** 12 independent empirical cohorts. Use 1,999 species-level bootstrap replicates and one-sided permutation p with floor .005.
5. **Decision gate and all results:** require >=30 retained distinct source species per cohort at each spacing/policy for positive bounded support; further require a positive mean, p<.05 and 95% species-bootstrap lower endpoint >0. All negative, low-power and high-dropout outcomes must be visible. The primary physical sensitivity is 5 km; 1 km and 10 km bound its opportunity dependence.

## Interpretation limitations

- A minimum pairwise distance between observed public photographic coordinates does **not** prove independent plant individual or site identity, given GPS error, plot geometry, missing individual IDs or transplanted/cultivated status.
- The 50-km ecological group is a photographic neighbourhood rather than an interbreeding population. The less-than-50-km spatial correlation between selected sites, photographer expertise and real plant ecology can remain.
- Conditioned subsamples change retained photos and can select a different cohort of species even though the original source eligible species denominator is fixed. Percent-positive claim is conditional on available photographic opportunity.
- The 12 realizations are resamples of the same photos, not independent image datasets or biological trials. Null permutations are post hoc and unadjusted over 3 spacings x 2 photographer policies.
- Neither source real photos nor expert human colour labels are newly acquired; this cannot measure true flower/petal identity, genetic morphs, local fitness or adaptation.
- Historical New Phytologist H1/H2, primary manuscript numerics and unopened future 2,000+730 source species remain unchanged.

**Pre-result status:** code + synthetic controls, no ecological conclusion until source SHA, original baseline and all terminal workflow checks pass.
