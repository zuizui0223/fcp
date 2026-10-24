# FCP: do the BIO5 white and BIO12 chromatic directions remain after accounting for correlated climate variables?

**24 October 2026. Source-based, post hoc photo-colour climate comparison, stacked after PR #267.** The primary prior photo-colour results are not revised or deleted by these conditional tests.

## Previously authenticated results (not results of this new model)

From 447,500 actual original measured photos and 4,475 species, 227,107 photographs were eligible for original nine-palette colour assessment. The prior original-photo analysis found a small BIO5/white-palette covariance with broad geographical and calendar-month adjustment, and the cross-species *mean photographic hue-climate vector* was about 93% aligned with a fixed white-versus-coloured axis. The observed mean BIO12/photographic hue direction contained chromatic yellow versus bronze contrasts after projecting away white. These were **separate single-climate observational associations**, not evidence of independent climatic causes: BIO1/BIO4/BIO5/BIO12/BIO14/BIO15/SRAD are correlated across geographic ranges.

Other critical original-source results remain: an independent flower-minus-background BIO5 association was unsupported in earlier validation (rho −0.00328, p=0.82), one source-photo observer-thinned fine geographical sensitivity was unsupported, original coarse photo-classification failure was 220,393/447,500, and no individually matched genetic morphs or fitness were measured. **This analysis cannot resolve those limits.**

## Central question and TWO separate estimands

1. **Domain-level independent contrasts:** within each named species, remove broad 3D site geography and original photo calendar-month trends from seven historical climate normals, then estimate the unique part of BIO5 *after accounting for BIO12 and mean solar radiation*, and the unique part of BIO12 *after accounting for BIO5 and radiation*. Also compute conditional solar radiation. This preserves biologically meaningful thermal/moisture target dimensions while recognizing the other temperature, seasonality, and precipitation indices are *not* held fixed in this domain-level model.
2. **Fully independent seven-climate contrasts:** for each of the seven original WorldClim/solar indices, remove within-species linear correlations with **all other six indices**, as well as broad geographical and photo-month terms. If the unique variation is nearly zero, mark the variable/taxon **not identifiable** rather than infer a climate-independent biological response.

For each predictor and species:
- Standardize original climate using that species' sampled climate SD, then regress out globally linear geographical XYZ source photo-location coordinates and calendar month.
- Obtain the climate axis' residual after projecting onto the other climate variables **in the named model**; require at least 5% of the original within-species standardized climatic variance to remain. A stricter 10% threshold is reported in the species support CSV, not used to select favourable biological direction after results are seen.
- Original-source photograph requirements: >=40 classifiable photo vectors with valid geography/dates and all original climate raster values, >=20 original photo observers and >=5 independent WorldClim 10-arcminute raster cells. No substitute or colour imputation for source photo classification failures. All 4,475 original taxa, including non-estimable ones, are retained in every result manifest.
- Calculate the signed photo climate *covariance* with each of the original nine normalized hue fractions and subtract the **analytically exact source-calendar-month permutation expectation**. Also divide by the actual **unique residual climate variance**, to obtain an *observational partial-slope-like association per original species climate SD*. It is **not** a causal per-degree temperature response.
- Exactly decompose the 9D vector into original coarse photographed-colour **between** and **within** components; between is reported only for taxa with >=2 original photo-colour labels each represented by >=5 photos. These are photograph labelling/measurement constructs, not heritable flower-colour morphs or genotypic plasticity.
- Apply the fixed white-versus-all-eight-nonwhite unit axis q=(8,-1,-1,-1,-1,-1,-1,-1,-1)/sqrt(72). Report q projection and the photo-hue covariance vector's orthogonal chromatic component. The squared q-component fraction of a **species-equal mean source-photo colour-climate vector** is not the fraction of genetic variation explained by temperature.
- Conduct colour-blind one original photo per observer sensitivity. Report BIO5/BIO12 signed white responses in each of five species-disjoint acquisition cohorts, including the fraction of evaluable species with positive observed whiteness association, not just the pooled mean.
- Bootstrap complete genera 999 times, reporting simultaneous **seven climate axis** 95% limits for q-projected whiteness separately for each photographic source panel and model. Model family, source photo class component, and full-9-hue multiple comparisons are still exploratory; no across-all-choices familywise error guarantee.
- **Matched-source climate directional comparison** is made on the exact same eligible taxa for BIO5 and BIO12 within each model. Report the bootstrap uncertainty in mean-vector alignment differences, not just pooled ratios across different sets of source species. If BIO12's mean colour vector is near zero its directional ratio may be unstable and this is reported.

## Inferential decisions / claim ceiling

- A positive independently temperature-conditioned photo-white direction would be interesting **observational evidence** that the BIO5 pattern cannot be wholly reproduced by BIO12 and solar covariation under the specified *linear* projections. It does not remove other correlated temperature/seasonality axes unless the strict seven-climate model also identifies a similar direction.
- A non-estimable seven-climate model due to multicollinearity **does not refute** the original mean photo hue direction or prove no biological effect. The denominator and unique variance distribution are ecological and sampling results in their own right.
- A BIO12 yellow/bronze contrast disappearing after thermal and solar adjustment suggests the earlier **single-climate** pattern was not unique to annual precipitation. A contrast that remains is still not a phenotype-specific or causal dry-versus-wet pigment response.
- Any new model is post hoc, same photo source/assay. Photo ROI may include flowers, bud, fruit, sepals, background, lighting and introduced/cultivated photo effects; prior original flower-minus-background negative remains.
- This work **does not demonstrate** the existence, maintenance or relative frequency of genetic flower colour morphs, anthocyanin biosynthesis costs, local selection, fitness, pollinator change or a molecular mechanism.

## Source immutability and execution

- Original photo release: fcp-all4730-original-method-actual-colour-results-20261017; original source CSV SHA256 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6.
- Original WorldClim BIO archive SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf; monthly SRAD archive SHA256 c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026.
- Code: scripts/analysis/fcp4475_joint_conditional_climate_hue_20261024.py.
- Original source CI: .github/workflows/fcp-4475-joint-seven-climate-hue-20261024.yml.
- Outputs: one JSON with exact source-photo denominators, independent climate axis opportunity and genus bootstrap; four *full 4,475-taxon* support ledgers and four compressed numerical 4,475 × 7 climate × 3 photographic components × 9 original hue photo covariances/conditional slope arrays **without individual photo IDs, source coordinates or observer identity**.
