# FCP global 42,111 original flower photos: different congeneric species WITHIN one geographic cell

## New ecological question

The already-verified 18,413 original climate-complete classified photos demonstrate a small rainfall-associated colour predictor among DIFFERENT congeneric species when their genus category mean is estimated from other species and broad source cells are held out. However, that comparison can still use distant geographic cells, so the result may represent habitat/geographic sorting among congeneric species rather than environmental variation among closely colocated species.

A stricter observational contrast asks:

**When two or more different photographed species are from the SAME genus AND the SAME original 162 equal-area geographic cell, does variation in their actual photographed-site climate still improve prediction of the source four-colour classes?**

This is NOT a within-species contrast, within-population colour polymorphism, phylogenetic evolutionary-rate estimate or selection test. The 162 equal-area cells are deliberately coarse and may encompass substantially different local climates. A photo-position difference inside a cell is an ecological gradient, not proof of a shared community.

## Pre-existing original source and fixed groups

The entire original source remains 42,111 species, 18,457 flower-ROI-classifiable (23,654 unclassifiable). Use exact original photo coordinates and actual WorldClim 2.1 10-minute climate/SoilGrids 5km topsoil already verified in [original environmental Actions 37861582446](https://github.com/zuizui0223/fcp/actions/runs/37861582446).

- First cohort: 18,413 original colour-classifiable climate+elevation-complete species, without soil-completeness conditioning.
- Second cohort: 14,136 original colour-classifiable climate+soil/elevation-complete species (separate observational support).
- Define source group as `(original genus name, geographic equal-area cell from TRUE original photo lat/lon)`. There are exactly 162 possible spatial bins but the number of sampled groups is data dependent and **must be reported**.
- For local group eligibility require at least THREE DIFFERENT original source species. This is an observer-independent, species/photo-ID and source-geography-only sufficiency rule, never selected from flower-colour class or temperature/rainfall responses. If insufficient support, emit a coverage HOLD without rescuing the minimum.
- Expose discarded original-source denominators, eligible genera, eligible cells and local groups. No unclassified original source photo becomes white/monomorphic, and no new imagery/metadata is acquired.

## Fitting and predicted outcome

Original species/photo-ID-only SHA256 ranks within each fixed group assign one of 5 folds independently of photo colour or climate. A test specimen's training genus-cell colour baseline must be calculated from at least two OTHER original species in that same group. Training and test source species never overlap; necessarily other congeners in the SAME cell appear in training by construction (this is **local within-group heldout species prediction**, not out-of-cell spatial generalization).

Use train-only same-genus-cell centring and one shared ridge 10.0 within-group linear response for the original four class probabilities; project them to the probability simplex. Feature families:
- genus×cell categorical baseline alone;
- baseline + actual finer geographic absolute latitude / longitude sine-cosine / elevation;
- plus BIO1/BIO5 temperature block only;
- plus BIO12/BIO15 precipitation amount/seasonality block only;
- plus both temperature and precipitation blocks;
- plus five SoilGrids properties on the fixed soil-complete subpopulation only.

Within each cohort, all feature families must evaluate exactly the SAME original test photo-species and original 5-fold assignments. Outcome is heldout multiclass Brier error; positive incremental improvement means lower error on other photographed congeneric species of the same cell. Report full climate increment, temperature increment conditional on moisture, moisture increment conditional on temperature, and soil increment conditional on climate.

Uncertainty: 999 resamples of the fixed prediction errors clustered separately at (i) original full geographic cell, (ii) source genus and (iii) local genus×cell. These condition on already generated predictions, not refitted models or independent randomization. A robust positive contribution should have a consistent sign and positive intervals across all 3 clusterings, rather than relying on a single high-richness region/lineage.

## Biological limits and possible outcomes

- **Signal persists within same genus×cell:** stronger observational evidence that the original moisture pattern is not solely explained by turnover between very broad geographic cells or genera. Still not a causal adaptation inference; geographical cells are vast, photos are one per species and narrower environment/habitat/relatedness confounding remains.
- **Signal disappears:** earlier weak moisture predictability may largely reside in genus-specific sorting across distant cells, observational/phylogenetic composition or a drop in power/coverage; it does **not prove** a complete absence of local environment selection.
- **Insufficient group support:** an identifiability/representation bottleneck, not a biological zero; do not relax the >=3 local congeneric taxa requirement just to claim a result.

The original same-species 3,182 two-locality photo-pair negative result and source photo-classifiability negative control must remain separate. Original FCP New Phytologist 1,499 deep-photo H1/H2 conclusions and untouched 2,000+730 prospective taxa remain unchanged.

Dedicated source-only CI `.github/workflows/fcp-global42111-local-congeners-20261009.yml` reads ONLY previously verified original 42,111 photo-source climate metadata. Output folder `results/fcp_global42111_local_congeners_20261009`, status and full original species denominator are verified before any explanatory claims.


## Executed original source true-photographic-distance and null results

**Complete [Actions 37893253874](https://github.com/zuizui0223/fcp/actions/runs/37893253874):** the original 6,290 photographed taxon group-level sample had been selected from very large fixed equal-area cells. Requiring all members of each original genus-cell to be within an actual maximum great-circle diameter of 100km left 193 climate-complete and 127 soil-complete species (**HOLD**, below fixed 300). At 250km there were 872 classified climate-complete taxa (217 local groups, 184 genera, 51 original cells), and at 500km 1,761 taxa (417 groups, 306 genera, 69 original cells). The original stricter photo-distance moisture predictor's heldout multiclass Brier improvement was **+0.0085980** at 250km and **+0.0051153** at 500km, with group resampling intervals positive under geographic-cell, genus and local genus-cell cluster policies. These nested subsets are **not independent confirmations**.

**Soil caution:** on the much smaller soil-complete source 620 taxa (250km) and 1,322 taxa (500km), the incremental SoilGrids topsoil gain was +0.000529 and +0.000574, respectively, but the genus/cell/group bootstrap intervals cross zero. Therefore the previously apparent SoilGrids gain in the 4,564-tax on 20-degree-cell congeneric local model is not demonstrably robust after actual site-distance restriction. Do not claim that soil determines petal colour or that its modelled null is universal.

**Independent model-fitting null calibration:** completed [Actions 37894008208](https://github.com/zuizui0223/fcp/actions/runs/37894008208), reading only the identical original source photos. Each of 199 null realizations shuffled four-state photo-colour labels within the SAME source genus×original cell, preserving that group’s exact four-state class counts, species IDs, geographic original photos, WorldClim values, folds and model hyperparameters. Both baseline GEO+BIO1/BIO5 and full GEO+BIO1/BIO5+BIO12/BIO15 models were RE-FITTED on all 199 shuffled realizations. Full [source result JSON](../results/fcp_global42111_local_congener_rain_null_20261009/result.json).

| Source photographed original distance cohort | Observed out-of-species Brier improvement from rain beyond geography+heat | Mean within-group label-shuffle null improvement | Null 2.5–97.5% interval | Nulls >= observed among 199 | One-sided conditional permutation p |
|---|---:|---:|---:|---:|---:|
| 250km, 872 taxa | +0.008598 | +0.001654 | [−0.001334, +0.005938] | 1 | **0.010** |
| 500km, 1,761 taxa | +0.005115 | +0.001602 | [−0.000304, +0.004031] | 2 | **0.015** |

**Source-grounded ecological inference:** original photographed visible flower-colour classes are associated with actual site rainfall predictors among DIFFERENT CONGENERIC species from neighbourhoods limited to 250–500km, beyond coarse genus×cell categorical composition, local longitude/latitude/elevation, and original site heat features. A simple class-frequency-preserving source shuffle rarely recovers the same magnitude of Brier reduction. Still this is retrospective and uses a high-level photo class, large site distance caps, local species turnover and selected climate/classifiable photo support; not an experimental evolutionary adaptation demonstration. No p-value should be sold as independently preregistered or multiplicity-corrected, and the two nested threshold tests are overlapping.

The major unaddressed rival explanations are other spatially correlated ecological gradients (habitat, soil moisture, land use, light, pollinators), species-level phylogenetic inheritance beneath genus, observer/photo-colour classifier bias and selective taxon sampling. This does not contradict the separate same-species 3,182 source pair climate null because that assay targets another biological/inferential scale with one photographed pair per species.
