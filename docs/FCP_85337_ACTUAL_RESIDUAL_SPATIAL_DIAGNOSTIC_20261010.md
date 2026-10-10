# FCP original 85,337 photographed taxon×region records: heldout spatial residuals after all abiotic blocks (2026-10-10)

## Source and reproducibility

Actual verified [GitHub Actions 38022723695](https://github.com/zuizui0223/fcp/actions/runs/38022723695) passed all original 85,337 source records/39,075 classified photo status guards, seven original abiotic blocks, and original source 100,543 real photo site metadata (no image pixels accessed). Output [source numeric receipt](../results/fcp_original85337_trainonly_geodesic_spatial_environment_20261010/result.json).

- **Climate source cohort**: 38,968 original classified original photo species×cell records with complete non-soil WorldClim/solar/wind/vapor/elevation; **20,546 heldout photos from 5,127 different nominal source species** have a training-only species colour mean from other original geographical cells.
- **Climate+soil cohort**: 27,003 original complete classified photo source cells; **11,136 heldout photos from 3,401 species** have a training-only colour mean; soil source variables include five original SoilGrids properties plus CEC, sand, silt, bulk density and coarse fragments.
- Training-only original photo geographical knots fitted independently in each source geographic-cell heldout fold without accessing its PHOTO COLOUR or true test locations in fitting. Two fixed geodesic spatial low-rank exponential decay scales: **250km** and **1000km**. Residual category-probability predictions retain the same originally photographed photo ID train/test fold for all three models: species + nonlinear geographic polynomial; plus real photo-distance low-rank spatial field; plus all abiotic predictors.

### Positive incremental abiotic photographic prediction remains

Positive means lower heldout 4-state source photo-colour Brier:

| Source population / original geodesic spatial kernel | Added all environment beyond original species+nonlinear spatial+photo distance | Original source-region group 95% interval | Original species group 95% interval |
|---|---:|---:|---:|
| Climate; 250km | **+0.001254** | [0.000771,0.001702] | [0.000880,0.001719] |
| Climate; 1000km | **+0.000753** | [0.000265,0.001268] | [0.000360,0.001175] |
| Climate+soil; 250km | **+0.003076** | [0.002259,0.003852] | [0.002332,0.003918] |
| Climate+soil; 1000km | **+0.002333** | [0.001494,0.003295] | [0.001412,0.003142] |

Conditional block-wise Brier increment under both spatial bandwidth settings: precipitation is positive in both the source climate and climate+soil cohorts, and soil is positive in the soil-complete cohort. Temperature is NOT positively distinguished at both bandwidths in both cohorts; solar radiation is NOT independently positive at either scale. All these are source-photo predictive associations conditional on chosen basis, not independently replicated causal effects or local adaptation.

### Newly computed original photo-neighbour residual Moran-like statistic

The next source-specific question is not just whether environment improves Brier. **Are heldout category-probability residuals still spatially similar at local original photographed sites?**

Using actual public original photographed latitudes and longitudes, identify each heldout photo's 10 nearest photographed source neighbours (BallTree haversine geodesic); omit pairs that are two photos of the SAME nominal species; reduce to unique undirected original cross-species neighbor edges. For source PHOTO residual vector `r_i = one-hot(original classified flower-photo colour) - heldout predicted four probabilities`, subtract its class-vector average across all heldout source photos, then calculate a descriptive centered four-class Moran-like index:

`I = n_photos * SUM_(undirected original nearby cross-species pairs) dot(centered_r_i, centered_r_j) / [n_undirected_photo_pairs * SUM_i norm(centered_r_i)^2]`

Report separately for pairs no more than **50km / 100km / 250km** apart, on the EXACT same heldout source photo IDs, for (a) global spherical polynomial alone, (b) real photo-distance spatial field, (c) real photo-distance spatial field + all environment. This is a **bounded 10-nearest-neighbor source graph**, NOT a full all-photos spatial weights matrix, and no permutation-calibrated p-values/uncertainty are computed.

Example, spatial field decay **250km**:

| Original repeated source population | Max original cross-species PHOTO-neighbor distance | n original local neighbor edges | Moran-like I after source spatial field alone | I after spatial field + all abiotic |
|---|---|---:|---:|---:|
| 20,546 climate-eligible heldout source photos | 50km | 86,036 | +0.002039 | **+0.002600** |
| climate | 100km | 107,641 | +0.001929 | **+0.002523** |
| climate | 250km | 127,397 | +0.002308 | **+0.002919** |
| 11,136 soil-complete heldout source photos | 50km | 35,919 | +0.006553 | **+0.007299** |
| soil | 100km | 49,954 | +0.006852 | **+0.007269** |
| soil | 250km | 65,884 | +0.003471 | **+0.004008** |

With the separately fixed 1000km field bandwidth, the after-environment Moran-like source index was 0.003791/0.003527/0.003946 in the climate population across 50/100/250km photo-neighbour radii, and 0.007944/0.007895/0.004518 in the soil-complete population. It therefore never vanishes, and the addition of all abiotic blocks does not uniformly reduce the diagnostic.

**Interpretation caution:** these I values are numerically small and no calibrated spatial null or standard error exists. The minor increase after adding environment cannot be declared statistically significant, nor imply biological reverse causation. But neither can we claim the residual field has zero significant spatial autocorrelation or that space was 'fully corrected'. The finite-rank source spatial kernel only approximates an underlying photo-space covariance, with observational species/photographer and plant community structure still unmeasured. The distinction between Brier predictive improvement and spatial residual independence is vital.

### Status and scientific decision

- **Retrospective conditional photo prediction**: the source same-species repeated photo grain retains small environmental predictive information under a real distance-dependent spatial basis, with precipitation and (if fully observed) modelled soil often contributing unique heldout Brier improvements.
- **Not spatially independent/causal identification**: close photographed species still exhibit small positive centered residual spatial similarity, true phylogenetic covariance in the separate 342/649 direct-tree-tip interspecies cohort did NOT yield robust seven-block abiotic improvements over its combined spatial and true-tree baseline, and the original 872/1761 local congeneric fine spatial permutation null did NOT support a robust rainfall-colour adaptation story.
- **Full genetic/fitness inference HOLD**: neither observation grain establishes an inherited pigment-morph distribution, pollinator fitness tradeoff, developmental canalization or a causal selection response to climatic sunshine/precipitation/soil.

The original higher-depth 1,499 species FCP paper and unopened future 2,000+730 source photo taxa remain untouched. PR #148 is a draft standalone spatial sensitivity to PR #144, not a merge-ready causal ecological finding.
