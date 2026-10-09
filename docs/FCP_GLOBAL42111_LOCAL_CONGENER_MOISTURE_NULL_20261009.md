# FCP global photographed flower colour: local congeneric precipitation null

## Why matched source-label randomization is necessary

The source-verified 42,111-species photograph atlas's local congener model (PR #140) found positive rainfall block prediction of photographed colour among different congeneric species whose photo positions all lie within 250km or 500km of each other. However original species/colour frequencies are strongly structured by genus, observation effort, and coarse geographic cell. Spatial and genus-cluster bootstrap intervals alone condition on one fitted prediction and cannot establish that the rainfall increment is unusual under a null preserving the photographic species-composition structure.

The stricter test therefore permutes photographed flower-colour labels ONLY among the different original photographed source species from the **same genus×original 162 geographic cell**. For each group, preserve the exact number of each of the original four coarse photo-colour categories; keep every source photograph taxon ID, actual photographic latitude/longitude, site WorldClim mean, photo ROI classification status, site group, and deterministic genus-cell within-group heldout folds unchanged. Repeat fitting, predicting and Brier scoring, with the same local genus-cell baseline, latitude/longitude/elevation and temperature BIO1/BIO5 covariates for baseline, plus rainfall BIO12/BIO15 covariates for the comparison.

## Immutable original 250 and 500km cohorts

Choose source groups before looking at photo-colour values, using >=3 congeneric source photographs entirely within 250km or within 500km maximal great-circle photographed-site group diameter. Do not choose a favourable distance after observing results; report both overlapping source groups, **which are not independent replication**.

- Original climate-classifiable 42,111 species' fixed sample from prior [distance-filtered CI](https://github.com/zuizui0223/fcp/actions/runs/37893629348): **872 source photographed species in 217 groups** at <=250km; observed rainfall Brier improvement **+0.0085980344** beyond local geography and temperature.
- At <=500km: **1,761 source photographed species in 417 groups**; observed improvement **+0.0051152608**.
- <=100km remains source coverage HOLD (**193 original source photo species**, below frozen 300 threshold) and is explicitly NOT subjected to a post-result threshold relaxation.

Use a fixed 199 permutation schedule per distance cohort, RNG seed from `20261009 + threshold`. Preserve source photo group four-state counts *exactly*, apply the exact original taxon-ID-only 5-fold plan, train the within-group categorical photo-colour means from *other* source species, refit the same ridge 10.0 centred predictor coefficients, then re-predict every heldout original species using only the other training species.

Source-conditional one-sided permutation p = (1+#null improvements >= original observed improvement)/(199+1), resolution 0.005; include null mean, SD and 2.5/50/97.5% levels. Source-photo candidate groups were originally selected on geography alone. No image pixels are downloaded and no unclassifiable photo labels become white or monomorphic.

## Interpretation and inference boundary

If observed rainfall improvement is more extreme than group-preserving source-colour label permutations, the original within-local-congener photographic environmental association is less plausibly an automatic statistical artefact of original fixed photo group size/composition alone. This is nevertheless retrospective, not an independent prospective genus/spatial experiment or proof of adaptively selected flower colours. True flower colours are observed once per species, source classifier misses 23,654 taxa, the 162 geographic cells are broad, actual local sites can differ in multiple environmental and biotic factors, and full phylogenetic context remains unmeasured.

If null is not exceeded, the earlier source-cluster bootstrap intervals are insufficient support for claiming genuine within-local environmental predictability. Do not rescue by changing 199 permutations, thresholds, predictor blocks or regional group definitions.

Never reinterpret original high-depth 1,499 source, FCP main paper H1/H2, or untouched independent 2,000+730 taxon cohorts.
