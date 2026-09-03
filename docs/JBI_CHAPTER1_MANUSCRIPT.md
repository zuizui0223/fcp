# Continuous flower colour is locally organized within species without a universal global transition boundary

**Running title:** Spatial organization of flower colour  
**Article type:** Research Article  
**Target journal:** *Journal of Biogeography*  
**Draft status:** The photograph split, operational colour representation, Stage-A/Stage-B contracts and numerical results are frozen. The prospectively frozen 200-species terminal scale-out is also closed as `not_evaluable` at its measurement-completeness gate. Introduction and Discussion wording remains subject to reference audit and co-author revision.

## Abstract

### Aim

We asked whether continuous flower colour is spatially organized within species after conditioning on each species' sampled range, and whether independent species place their strongest colour transitions in the same global regions.

### Location

Global community-photograph sample.

### Taxon

Angiosperms: six species with documented natural intraspecific flower-colour variation.

### Methods

We acquired 1,200 georeferenced photographs and froze an outcome-blind 480/720 calibration–evaluation split. Species-specific continuous colour vectors were derived from frozen flower regions and standardized using calibration-only parameters. Stage A compared local colour discontinuity on colour-blind within-species nearest-neighbour graphs with 9,999 species-conditioned permutations. After its prospective gate passed, Stage B tested whether ranked within-species transition intensity was geographically concentrated across geometry-defined detectable cells under a second 9,999-permutation null.

### Results

All 720 held-out photographs were processed without localization failure or post-opening rule modification. Primary Stage A discontinuity was lower than its null expectation (1.3911 vs 1.4294; standardized deficit = 2.31; *p* = 0.0113), with the same direction at `k = 3` (*p* = 0.0066) and `k = 8` (*p* = 0.0065). Primary Stage B concentration exceeded its null mean but was not confirmatory (0.00823 vs 0.00568; standardized excess = 1.44; *p* = 0.0906). One of eight non-primary spatial sensitivities was nominally supported.

### Main conclusions

Continuous flower colour was spatially organized within species, but independent species did not show confirmatory concentration of their strongest transitions along one shared global geography. Spatial organization and spatial coincidence are therefore distinct empirical questions; the present evidence supports heterogeneous within-species geography rather than a universal boundary or common mechanism.

**Keywords:** citizen science, colour polymorphism, continuous traits, flower colour, random labelling, spatial boundaries, spatial organization, species-conditioned permutation

## Introduction

Intraspecific phenotypic variation is not only a matter of how many phenotypes occur within a species. It also has a spatial organization. Alternative phenotypes can coexist locally, replace one another among populations, vary continuously across geographic space or combine these patterns at different scales. These configurations matter because mechanisms capable of maintaining local diversity are not identical to those that generate differentiation among populations. Negative frequency-dependent selection, opposing selective agents and microenvironmental heterogeneity can maintain local variation, whereas spatially varying selection, restricted dispersal, demographic history and drift can generate geographic structure (Delph & Kelly, 2014; Narbona et al., 2018).

Flower colour is a tractable system for studying this organization. Floral pigmentation can affect pollinator attraction, but pollinator-mediated selection is neither universal nor the only relevant process. Abiotic stress, herbivory, correlated traits, gene flow, drift and demographic history can all contribute to colour variation (Narbona et al., 2018; Trunschke et al., 2021). Consequently, detecting spatial organization in flower colour does not by itself identify a selective agent, and similar-looking geographic patterns need not share a mechanism.

Broad-scale analyses have related flower-colour properties to environmental and biotic conditions, while individual-species studies have used field surveys or community photographs to describe geographic variation in colour morphs (Dalrymple et al., 2020; Farquhar et al., 2023; Jansen et al., 2025). These approaches establish that colour variation can be geographically structured, but two inferential questions are often conflated. The first is whether colour is non-randomly arranged within a species after accounting for its range and sampling geometry. The second, stronger question is whether independent species place their colour transitions in the same geographic regions.

A pooled global colour map cannot answer either question on its own. Species composition, observation effort and range boundaries are spatially heterogeneous. Unconditional shuffling of colour observations among species would therefore test the geography of species composition rather than the spatial organization of flower colour. Species identity may be visually suppressed in a global display, but it must remain in the null model.

Comparing transition locations among species creates an additional denominator problem. A species can contribute evidence about a geographic cell only where its observation geometry would have allowed a transition to be detected. A cell with no observed transition is not equivalent to a cell in which the species was not sampled on a spatial support capable of revealing one. Shared-boundary analysis must therefore distinguish a detected transition, an evaluable absence and a non-evaluable location.

We developed a held-out, species-conditioned analysis of continuous flower colour from georeferenced community photographs. We asked two ordered questions. **Stage A** tested whether geographically neighbouring observations were more similar in colour than expected after fixing each species' locations and permuting complete continuous colour vectors only within that species. **Stage B**, run only after Stage A passed its prospective gate, tested whether relatively strong transitions from independent species were geographically concentrated in the same globally evaluable cells. We expected Stage A to detect spatial organization if the photographs retained biological geographic signal. We treated Stage B as a distinct, stronger hypothesis: support would indicate repeated geographic coincidence, whereas non-support would imply that spatial organization is heterogeneous among species rather than governed by one universal global boundary.

## Materials and Methods

### Study design and frozen sample

The analysis used six species in the frozen Chapter-1 development set: *Antirrhinum majus*, *Dactylorhiza sambucina*, *Gentiana lutea*, *Ipomoea purpurea*, *Lysimachia arvensis* and *Raphanus sativus*. For each species, 200 georeferenced photographs were acquired, giving 1,200 photographs in total.

Photographs were assigned deterministically and outcome-blind to 80 calibration and 120 held-out evaluation observations per species. The resulting totals were 480 calibration and 720 evaluation photographs. The split, assignment basis and per-species counts were hash-frozen before the evaluation set was opened. Photograph measurements were joined to geographic coordinates by `photo_id`.

The candidate species set and community-photograph acquisition procedure define the inferential scope of this study. The six species are not a random sample of angiosperms, and the analysis does not estimate the global prevalence of spatially organized flower-colour variation.

### Operational continuous colour representation

Calibration feature geometry did not support forcing one common discrete-morph representation across all six species. We therefore froze a species-specific continuous colour-vector representation before inspecting evaluation values.

For each photograph, the operational flower region was localized using the frozen Florence model and prompt path. A predeclared species-specific subset of visible colour features was extracted from that region. Each feature was standardized using the corresponding species' calibration-only mean and population standard deviation. The complete standardized vector was retained as one indivisible observational label.

The vectors are operational descriptions of visible colour structure in community photographs. They are not spectrophotometric measurements, biological morph labels or estimates of pigment concentration. No evaluation-derived feature selection, scaling parameter, region-of-interest rule or vector dimension was introduced after the evaluation set was opened.

### Evaluation integrity

The frozen extraction path was applied unchanged to all 720 held-out photographs in 36 deterministic shards of 20 observations. Evaluation integrity required 120 unique observations per species, unique `photo_id` and blind identifiers, successful feature extraction, absence of calibration observations and no final discrete biological label.

All evaluation records were retained. The analysis did not impute missing colour vectors, discard observations on the basis of their measured values or tune extraction rules after inspecting the evaluation set.

### Stage A: species-conditioned local colour organization

#### Colour-blind spatial graph

Within each species, we constructed a symmetrized spherical nearest-neighbour graph using observation latitude and longitude only. The primary graph used `k = 5`. Graphs with `k = 3` and `k = 8` were frozen as sensitivity analyses. Colour values did not select neighbours, remove edges or choose graph scale.

#### Edge discontinuity and global statistic

For observations `i` and `j` joined by a graph edge, with standardized colour vectors `z_i` and `z_j`, continuous colour discontinuity was

`d_ij = sqrt(mean((z_i - z_j)^2))`.

For species `s`, `Q_s` was the arithmetic mean of its edge discontinuities. The global Stage-A statistic was the arithmetic mean of the six species-specific `Q_s` values. This equal-species aggregation gave each species one vote irrespective of its graph-edge count.

#### Random-labelling null and prospective gate

For each of 9,999 permutations, observation locations and graph edges were held fixed, complete colour vectors were permuted strictly within species, and all species-specific and global statistics were recomputed. No vector component was permuted independently and no vector moved between species.

The confirmatory alternative was lower observed discontinuity than expected under species-conditioned random labelling. A plus-one lower-tail Monte Carlo probability was used. Stage B was run only if the primary Stage-A probability was at most 0.05.

### Stage B: shared continuous transition concentration

#### Geometry-only selection of spatial support

Stage B used the primary `k = 5` graph. Before any observed colour edge score was calculated, we evaluated a fixed ordered set of spatial configurations defined by edge-length caps of 500, 1,000 and 2,000 km and equal-area longitude–sin(latitude) grids of 36×18, 24×12 and 18×9 cells.

For each species and configuration, graph edges longer than the candidate cap were removed and retained edges were assigned to the cell containing their great-circle midpoint. A species was detectable in a cell only when at least two retained geometry edges had midpoints in that cell. Let `D_s(x)` equal one when species `s` was detectable in cell `x` and zero otherwise. The opportunity denominator was

`A(x) = sum_s D_s(x)`.

A configuration passed the frozen geometry criteria only if every species retained at least 30 edges, at least eight cells had `A ≥ 2`, at least two cells had `A ≥ 3`, and at least four species contributed to shared opportunity. The first passing configuration in the predeclared order was selected as primary. All nine configurations passed, so the primary support was the 500-km edge cap and 36×18 grid.

#### Species-cell transition intensity

Within each species and fixed configuration, RMS colour discontinuities were calculated on retained edges and average-rank transformed to `[0,1]`. Ranked edge intensities were then averaged within detectable cells, producing species-cell transition intensity `T_s(x)`. Ranking allowed comparison of the geographic placement of relatively strong and weak transitions without treating raw distances from different species-specific feature spaces as directly commensurate.

For cells with at least two detectable species, shared transition intensity was

`S(x) = sum_s D_s(x) T_s(x) / A(x)`.

Cells with `A < 2` were stored as not evaluable (`NaN`) rather than assigned zero.

#### Concentration statistic and complete null

The primary Stage-B statistic was the `A(x)`-weighted variance of `S(x)` across cells with `A ≥ 2`. Larger values indicated stronger geographic concentration of shared transition intensity.

For each of 9,999 primary permutations, locations, graph geometry, edge filters, grids and detectability masks were held fixed. Complete colour vectors were permuted strictly within species, after which edge discontinuities, within-species ranks, species-cell intensities, the shared surface and the concentration statistic were recomputed. The confirmatory probability was the plus-one upper-tail Monte Carlo probability.

The eight non-primary configurations were evaluated as predeclared descriptive sensitivities with 1,999 permutations each. A sensitivity result could not replace the prospectively selected primary configuration after inspection.

### Reproducibility and claim control

Analysis contracts, frozen inputs, null distributions, result manifests and canonical figure hashes are versioned in the repository. Environmental and historical reference layers were not used to select Stage-A or Stage-B support or to modify the observed colour surface. Because the primary Stage-B gate was not passed, no geographic-reference correspondence analysis was promoted to the confirmatory sequence.

### Prospective terminal scale-up

A separate prospective extension tested whether the image-first measurement system could produce an evaluable colour field at substantially broader taxonomic scale; it did not re-estimate or replace the six-species Stage-A/B results. Eight disjoint cohorts of 25 genus-distinct species were frozen with 300 observations per species, giving 200 species and 60,000 photographs. All 60,000 records, source roles, ROI-v4 estimator, thresholds, cohort requirements, branch order and stopping rules were fixed before terminal outcome inspection, with no replacement, resampling, favourable subset selection or early stopping.

The terminal sequence required complete location-blind measurement before any protected coordinate-colour join. A predeclared `not_evaluable` measurement-completeness outcome stopped the cascade before species-conditioned spatial organization and prohibited all downstream shared-transition, environmental and pollinator analyses. The scale-out therefore tested measurement evaluability, not the biological spatial hypotheses, unless its frozen completeness gate was passed.

## Results

### Held-out colour measurement

All 720 evaluation photographs were processed successfully, with 120 unique records for each species. All records had successful feature status, no localization failure occurred and no final discrete biological label was emitted. The evaluation path was completed without post-opening changes to feature selection, scaling or localization rules.

### Stage A: continuous flower colour was locally organized

The primary equal-species mean edge discontinuity was 1.39114, below the species-conditioned null mean of 1.42943 (Figure C1A). The standardized clustering deficit was 2.311 and the lower-tail Monte Carlo probability was 0.0113. Geographically neighbouring observations were therefore more similar in operational continuous flower colour than expected after preserving each species' observed locations and colour-vector distribution.

The direction was stable across the predeclared graph scales (Figure C1B). At `k = 3`, observed global discontinuity was 1.37661 compared with a null mean of 1.42967 (standardized deficit = 2.552; *p* = 0.0066). At `k = 8`, observed discontinuity was 1.39556 compared with a null mean of 1.42978 (standardized deficit = 2.613; *p* = 0.0065).

Species-specific results were heterogeneous (Figure C2). Primary lower-tail probabilities were 0.2635 for *A. majus*, 0.2560 for *D. sambucina*, 0.7865 for *G. lutea*, 0.2923 for *I. purpurea*, 0.0023 for *L. arvensis* and 0.0080 for *R. sativus*. The global inference was the prospectively specified equal-species test; the species-level results describe heterogeneity rather than define six post-hoc primary hypotheses.

### Stage B: shared transition concentration was not confirmed

Geometry-only selection chose the 500-km edge cap and 36×18 equal-area grid as the primary support. It retained 246, 363, 372, 241, 300 and 329 edges for *A. majus*, *D. sambucina*, *G. lutea*, *I. purpurea*, *L. arvensis* and *R. sativus*, respectively. Detectable-cell counts were 15, 8, 7, 17, 22 and 6. Twenty-four cells had `A ≥ 2`, nine had `A ≥ 3`, four had `A ≥ 4`, and maximum opportunity was all six species (Figure C3; Figure C-S2).

Observed opportunity-weighted concentration was 0.0082315, compared with a null mean of 0.0056757 and null standard deviation of 0.0017763. The standardized concentration excess was 1.4389. The upper-tail Monte Carlo probability was 0.0906 and the descriptive two-sided probability was 0.1372. The primary shared-concentration null was therefore not rejected.

Sensitivity results depended on spatial support (Figure C4). The 500-km/24×12 configuration gave an upper-tail probability of 0.0445, while probabilities for 500-km/18×9, 1,000-km/36×18, 1,000-km/24×12, 1,000-km/18×9, 2,000-km/36×18, 2,000-km/24×12 and 2,000-km/18×9 were 0.3415, 0.2235, 0.4945, 0.1920, 0.4690, 0.5500 and 0.3495, respectively. The isolated nominal result on the coarser 500-km grid was retained as exploratory scale sensitivity and did not replace the frozen primary result.

### Prospective terminal scale-up was not evaluable

The terminal 200-species extension completed all 256 location-blind compute partitions and reassembled exactly 60,000 unique measurement records into 16 semantic shards. The frozen measurement-completeness gate was nevertheless not met. Only 58 of 200 species were measurement-evaluable, with 3, 4, 7, 8, 8, 12, 7 and 9 evaluable species across the eight fixed 25-species cohorts. The gate therefore returned `not_evaluable_scaleout_measurement_completeness`.

The protected coordinate join remained prohibited (`coordinate_join_permitted = false`; `coordinates_opened = false`). Consequently, species-conditioned spatial organization and all downstream shared-transition, environmental-concordance and pollinator-concordance analyses were not run. This is a measurement-evaluability result rather than evidence for or against those biological hypotheses.

## Discussion

### Spatial organization and spatial coincidence are different results

The analysis produced an intentionally asymmetric result. Stage A showed that continuous flower colour was locally organized within species after preserving each species' observed range, sampling geometry and colour-vector distribution. Stage B did not show that independent species concentrated their strongest transitions in the same global cells at the frozen primary support. A finding of non-random within-species organization therefore does not imply one common biogeographic boundary.

This distinction matters for interpretation. Similar spatial autocorrelation can arise from different combinations of environmental gradients, dispersal, demographic history, drift, biotic interactions and image-scale measurement variation. The Stage-A result establishes geographic structure in the operational colour vectors. It does not establish that all six species respond to the same driver or that their transitions should coincide geographically.

### Heterogeneous species contributions

The equal-species global statistic was chosen to test a repeated tendency without allowing species with denser graphs to dominate. The strongest primary species-level signals occurred in *Lysimachia arvensis* and *Raphanus sativus*, whereas the remaining species were not individually resolved. This heterogeneity does not negate the prospective global test, but it limits any claim of universal species-level strength.

Several non-exclusive explanations remain. Species can differ in the scale and shape of their colour gradients, the alignment of sampling locations with those gradients, the environmental or historical processes structuring populations and the reliability with which community photographs represent visible floral colour. Estimating how frequently strong local organization occurs across angiosperms will require a future measurement system whose transfer across taxa is prospectively qualified before a new large-scale biological test is opened.

### Why the Stage-B non-result is informative

The shared-boundary test asked a substantially stronger question than Stage A. It required not only structured colour variation within species, but repeated placement of relatively strong transitions in the same globally evaluable cells. The primary probability of 0.0906, together with inconsistent neighbouring spatial supports, does not justify promoting a universal boundary.

The nominal 500-km/24×12 sensitivity result indicates that some cross-species overlap may be visible at a coarser grid. However, the absence of support at the prospectively selected finer grid and at the remaining seven supports shows that this overlap is scale-dependent. Treating the single nominal sensitivity as the main result would reverse the predeclared selection rule and overstate evidence.

The denominator `A(x)` is central to this interpretation. Strong shared intensity in a cell can reflect agreement among only the species for which transitions were geometrically detectable there. Figure C-S2 makes this opportunity structure explicit and prevents non-evaluable species from being counted as biological zeros or silent agreement.

### Community photographs as comparative trait data

The complete held-out processing demonstrates that a frozen computer-vision path can produce reproducible continuous colour descriptors from large community-photograph samples. Retaining continuous vectors avoided forcing heterogeneous species into a universal discrete-morph scheme unsupported by calibration geometry. Freezing the representation, scaling and analysis before opening the evaluation set also separated method development from confirmatory inference.

These advantages do not convert community photographs into calibrated spectral measurements. Illumination, camera processing, viewing angle, flower condition and imperfect localization can contribute measurement variation. The present test is therefore about spatial organization of operational visible-colour vectors under a fixed pipeline. Future work can compare these vectors with spectrophotometry, standardized field imaging or repeated photographs of the same individuals and populations.

### Limits and next tests

The study has six principal limits. First, the six species are a selected development set rather than a random sample of angiosperms. Second, community photographs provide uneven spatial and temporal sampling even though the null preserves the observed geometry. Third, colour vectors are species-specific operational measurements and cannot be interpreted as directly commensurate pigment or receptor-space coordinates. Fourth, the analysis detects spatial organization but does not estimate morph-specific fitness, gene flow or causal environmental effects. Fifth, Stage-B power is constrained by the number of species and the limited set of cells with shared detectability. Sixth, the prospectively frozen terminal scale-out showed that the current automated ROI-v4 measurement system did not transfer with sufficient completeness across 200 species to authorize a new biological inference.

That terminal extension is informative because its stopping rule was applied before coordinates were joined to measured colour. Exact location-blind measurement was completed for all 60,000 frozen records, but only 58 of 200 species met the fixed measurement-evaluable rule and none of the eight cohorts passed the frozen cohort requirement. We therefore stopped rather than relaxing the gate after observing performance. The scale-out must not be interpreted as a negative test of spatial organization, shared transitions or environmental concordance; those analyses were never opened.

A future confirmatory expansion must therefore begin as a new prospective experiment with broader taxonomic measurement transfer qualified before new outcome pixels are opened. Environmental or historical layers can then be introduced only after a reproducible shared-transition pattern is established, or tested within species under independently frozen hypotheses. In the present six-species analysis, adding geographic overlays after the unsupported Stage-B primary test would be exploratory and cannot rescue the common-boundary hypothesis. The present 60,000 terminal records likewise cannot be re-analysed as a favourable 58-species subset or rescued with substituted post-result thresholds.

## Conclusions

Across 720 held-out community photographs from six angiosperm species, geographically neighbouring observations were more similar in continuous flower colour than expected under species-conditioned random labelling. This result was stable across three nearest-neighbour graph scales. The stronger hypothesis that independent species concentrated their strongest colour transitions in the same global regions was not supported by the frozen primary test. Flower-colour variation can therefore be repeatedly spatially organized within species without being demonstrably governed by one universal global boundary. A separately frozen 200-species extension did not reach biological inference because its predeclared measurement-completeness gate was not met, underscoring that current automated colour measurement does not yet transfer uniformly enough across taxa for this terminal global atlas design.

## Data availability and reproducibility

The frozen split, operational colour representation, evaluation features, Stage-A and Stage-B contracts, null distributions, result manifests and canonical figure products are versioned in this repository. The terminal scale-out design, exact artifact provenance and final `not_evaluable` decision are also versioned separately so that the scale-out cannot be mistaken for a completed biological test. Primary entry points are:

- `data/frozen/jbi_ch1_photo_split_v1.csv`;
- `docs/supporting/jbi_ch1_continuous_colour_representation_v1.json`;
- `data/evaluation/jbi_ch1_florence_evaluation_features_v1.jsonl`;
- `docs/supporting/jbi_ch1_stage_a_continuous_graph_v1.json`;
- `data/evaluation/jbi_ch1_stage_a_primary_null_v1.csv`;
- `docs/supporting/jbi_ch1_stage_b_geometry_audit_v1.json`;
- `docs/supporting/jbi_ch1_stage_b_shared_transition_concentration_v1.json`;
- `data/evaluation/jbi_ch1_stage_b_shared_transition_surface_v1.csv`;
- `data/evaluation/jbi_ch1_stage_b_primary_null_v1.csv`;
- `docs/supporting/jbi_ch1_figure_manifest_v1.json`;
- `docs/JBI_CHAPTER1_TERMINAL_SCALEUP_RESULT.md`;
- `docs/supporting/jbi_atlas_terminal_measurement_v5_receipt.json`.

## Figure legends

### Figure C1

**Species-conditioned spatial organization of continuous flower colour.** (A) The observed equal-species mean edge discontinuity in the primary five-nearest-neighbour graph relative to 9,999 null values generated by permuting complete standardized colour vectors strictly within species while holding observation geometry fixed. Lower values indicate stronger local colour similarity. (B) Standardized clustering deficit across the primary and predeclared graph-degree sensitivities. All three lower-tail tests rejected random labelling at 0.05.

### Figure C2

**Heterogeneous species contributions to the global Stage-A result.** Points show species-specific standardized clustering deficits in the primary five-nearest-neighbour graph; larger positive values indicate lower observed discontinuity than the species-conditioned null. Point size represents retained graph-edge count. The prospectively defined global test gave each species equal weight; individual probabilities describe heterogeneity rather than redefine the global hypothesis.

### Figure C3

**Observed shared-transition intensity under the frozen primary spatial support.** Grey points are the 720 held-out evaluation locations. Coloured cells show mean ranked transition intensity among detectable species for the primary 500-km edge cap and 36×18 equal-area grid. Marker size represents the opportunity denominator `A`, the number of species with sufficient geometry support in the cell. Cells with fewer than two detectable species are not evaluable and are not plotted as zero.

### Figure C4

**Shared-transition concentration was sensitive to spatial support and was not confirmed by the primary test.** Upper-tail Monte Carlo probabilities and standardized concentration excesses are shown for the selected primary configuration and all predeclared sensitivity configurations. The primary 500-km/36×18 analysis did not reject the shared-concentration null. One coarser 500-km sensitivity was nominally below 0.05, whereas the remaining configurations were not.

### Figure C-S1

**Frozen measurement and inferential workflow.** The diagram records the outcome-blind 480/720 split, freezing of the species-specific continuous representation before evaluation, the passed Stage-A gate, the unsupported Stage-B primary gate and the resulting decision not to promote geographic-reference correspondence as a confirmatory analysis.

### Figure C-S2

**Detectability behind the primary shared-transition surface.** The upper panel shows the opportunity denominator for each geographically ordered evaluable cell. The lower matrix shows species-cell ranked transition intensity where the species was detectable; grey entries denote non-evaluable species-cell combinations rather than zero transition intensity.

## References cited in this draft

Dalrymple, R. L. et al. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist*, 228, 1972–1985.

Delph, L. F. & Kelly, J. K. (2014). On the importance of balancing selection in plants. *New Phytologist*, 201, 45–56.

Farquhar, J. E., Pili, A. & Russell, W. (2023). Using crowdsourced photographic records to explore geographical variation in colour polymorphism of an Australian varanid. *Journal of Biogeography*, 50, 1409–1421.

Jansen, N., Pruijn, N. & Mayer, M. (2025). Citizen observations shed new light on geographic variation in colour polymorphism of a widespread reptile. *Journal of Biogeography*, 52, 629–640.

Narbona, E., Wang, H., Ortiz, P. L., Arista, M. & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology*, 20(S1), 8–20.

Trunschke, J., Lunau, K., Pyke, G. H., Ren, Z.-X. & Wang, H. (2021). Flower color evolution and the evidence of pollinator-mediated selection. *Frontiers in Plant Science*, 12, 617851.