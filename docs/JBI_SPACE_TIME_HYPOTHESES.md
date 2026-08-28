# JBI space–time hypothesis registry for flower-colour polymorphism

Status: **exploratory literature-driven registry, v1 (2026-08-27)**.  This document is downstream of the strict v2.2 single-pass C/S evidence audit.  It does not replace the frozen 34-species historical analysis and must not be described as preregistered confirmatory inference.

## 1. Main biogeographic question

The response is the **spatial organization of intraspecific flower-colour variation**, not colour itself:

- **C**: local coexistence documented;
- **S**: spatial segregation / among-place differentiation documented;
- **C+S**: both signals documented;
- unresolved: no positive C/S organization signal under the strict pass, not biological absence.

The working synthesis is that a species-wide climatic niche is too coarse if it is represented by a single breadth number.  The next tests ask how ecological conditions are partitioned through **time and space**.

## 2. Literature synthesis

### 2.1 Within-population maintenance is inherently pluralistic

Sapir, Gallagher & Senden (2021) review balancing mechanisms including multiple pollinators, opposing selection, fluctuating selection, heterozygote advantage and frequency-dependent selection. Narbona et al. (2018) additionally emphasize non-pollinator selection, drift, gene flow, autonomous selfing and clonal reproduction; Mediterranean FCP species commonly form mosaics of monomorphic and polymorphic populations rather than one uniform spatial state.

### 2.2 Temporal variation is a direct candidate mechanism for C

Long-standing examples show that selection on colour can change through time. Frey (2004) found opposing selection via herbivores/pathogens that fluctuated among seasons in *Claytonia virginica*. In *Clarkia xantiana*, pollinator species show different forms of frequency-dependent foraging, and the balance among pollinators can vary through space and time. In 2026, Heinze et al. showed that *Anemone pavonina* morphs differed in flowering peak and pollinator association along elevation; beetle versus bee abundance provides an explicitly spatiotemporal route to maintaining colour variation.

### 2.3 Persistent spatial environmental heterogeneity is a direct candidate mechanism for S

Abiotic selection can covary with colour across landscapes. Vaidya et al. (2018) linked colour to drought/herbivory in *Boechera stricta*. Grossenbacher et al. (2026) found pink *Leptosiphon parviflorus* more frequent under serpentine soil, higher temperature, drought deficit and UV-related conditions. Fetterly et al. (2026) found *Castilleja coccinea* colour associated with habitat conditions across predominantly monomorphic populations. Jaeger et al. (2026) showed pollinator, herbivore and climatic pressures changing across a floral colour transition zone. These studies motivate spatial environmental sorting but also show that a single climatic breadth is not equivalent to spatially structured selection.

### 2.4 Neutral isolation must remain an explicit alternative

Spatial colour differentiation can be generated or reinforced by restricted gene flow and drift. Work on *Iris lutescens* demonstrates that neutral processes can contribute to spatial morph-frequency patterns. Therefore a spatial signal should be decomposed into environmental turnover versus geographic fragmentation/isolation rather than automatically labelled local adaptation.

### 2.5 Hypervolume is useful only if geometry is decomposed

Blonder et al. (2014) provide the n-dimensional hypervolume framework; subsequent debate shows strong sensitivity to sample size, dimensionality, bandwidth and sampling bias (Qiao et al. 2017; Blonder 2017). Hypervolumes are most useful here for **regional centroid separation, overlap and temporal turnover**, not simply species-level total volume. Vilas et al. (2022) illustrate temporal hypervolume comparison using centroid distance, volume change, intersection and Jaccard/Sørensen dissimilarities.

## 3. Hypotheses and predictions

| ID | Mechanism | Prediction for C/S organization | First observable | v1 status |
|---|---|---|---|---|
| H1a | temporal thermal heterogeneity | stronger temporal thermal variability favours local coexistence C over persistent S | mean occupied BIO4 | **candidate only** |
| H1b | temporal precipitation heterogeneity | stronger precipitation seasonality favours C | mean occupied BIO15 | unsupported v1 |
| H2 | geographic fragmentation / restricted connectivity | more fragmented ranges favour S and/or C+S | multiscale disconnected-range index | unsupported v1 |
| H3a | spatial environmental turnover | environmental distance should increase more strongly with geographic distance in S | long-vs-short distance climate turnover | unsupported v1 |
| H3b | regional niche displacement | geographically separated sectors should have more displaced climate centroids in S | normalized regional centroid distance | unsupported v1 |
| H3c | regional niche overlap | geographically separated sectors should overlap less in climate space in S | Gaussian/Bhattacharyya overlap | unsupported v1 |
| H4 | pollinator turnover | temporal turnover of pollinators favours C; persistent spatial turnover favours S | pollinator diversity/turnover by place and season | not yet measured |
| H5 | opposing mutualist–antagonist selection | variable balance among pollinators, herbivores and pathogens favours C/mosaic states | interaction-community traits | not yet measured |
| H6 | breeding system / clonality | selfing/clonality modify the loss or persistence of local morph diversity and effective gene flow | self-compatibility, autonomous selfing, clonality | not yet measured |
| H7 | drift / gene flow | low connectivity with weak environmental correspondence favours S | neutral genetic differentiation / dispersal proxies | not yet measured |
| H8 | pigment-mediated abiotic trade-offs | spatially persistent stress gradients favour S; temporally fluctuating stress can favour C | drought, UV, soil, extremes | partly proxied; needs better data |
| H9 | phenological partitioning | morph-specific flowering time can stabilize local coexistence and alter pollinator exposure | flowering phenology / season length | not yet measured |
| H10 | negative frequency dependence / overdominance | stabilizing rare-morph or heterozygote advantage favours C | reward status, pollinator FDS evidence, genetic architecture | literature-only currently |
| H0 | total niche-size null | total hypervolume size alone should not be treated as the mechanism | rarefied PC1–PC3 hull volume | unsupported as discriminator / retained as negative control |

## 4. v1 tests on display-core-v6

Input state: 74 display-core species; 66 climate-eligible; 32 species with documented C/S organization. Conditional comparisons use the existing IPW L2 multinomial design, geographic-range control and 499 family bootstraps. These are **exploratory tests generated after literature inspection**.

### H1a — temperature seasonality

For a 1-SD increase in mean occupied BIO4:

- S vs C: OR = **0.349**, family-bootstrap 95% CI **0.085–0.909**;
- C+S vs C: OR = **0.401**, CI **0.140–1.002**.

This is the only v1 screen with a bootstrap interval excluding 1 for S vs C. However, the S/C interval crosses 1 after adding mean temperature or absolute latitude:

- + mean BIO1: OR 0.412, CI 0.100–1.227;
- + median absolute latitude: OR 0.425, CI 0.100–1.120;
- + both: OR 0.458, CI 0.112–1.342.

**Decision:** retain H1a as the highest-priority candidate, but do not claim a temperature-seasonality effect. The present BIO4 association may partly encode broad latitudinal/thermal geography. Confirmation requires dynamic, year-specific climate rather than climatological seasonality alone.

### H1b — precipitation seasonality

S vs C OR = 1.50, CI 0.41–5.15. No support for the predicted C increase.

### H2 — geographic fragmentation

The point direction is consistent with the prediction (S vs C OR = 1.80; C+S vs C OR = 1.75), but both intervals are broad and cross 1. Fragmentation is therefore a **mechanistically plausible but unconfirmed** axis.

### H3 — spatial niche structure

Three independent diagnostics are currently null:

- pairwise spatial environmental turnover: S vs C OR 1.37, CI 0.49–7.51;
- two-sector normalized climate-centroid separation: OR 1.62, CI 0.56–3.48;
- two-sector Gaussian climate overlap: OR 1.00, CI 0.48–3.36.

Thus the present WorldClim/GBIF data do not show that S species have clearly stronger large-scale climatic partitioning than C species.

### H0 — total hypervolume size

A PC1–PC3 convex-hull hypervolume was rarefied to 20 climate cells per species over 99 draws. Median hypervolume size does not clearly separate C/S states; S vs C OR = 2.08 (CI 0.46–6.79), while C+S vs C OR = 0.47 (CI 0.26–1.16). This reinforces the original conclusion that **total niche size is not the useful object by itself**.

## 5. Next execution order

1. **Dynamic temporal climate (priority).** Replace WorldClim BIO4 as the temporal proxy with occurrence-year/flowering-season climate anomalies and interannual variability (e.g. ERA5-Land/TerraClimate). Test temporal hypervolume turnover within species.
2. **Availability-corrected regional hypervolumes.** Compare realized regional hypervolumes only after defining accessible environmental background M for each region; quantify centroid displacement, overlap and shape change with rarefaction/sensitivity analyses.
3. **Biotic turnover.** Add pollinator guild/richness/seasonality and, where available, herbivore/florivore interaction evidence. Test temporal turnover -> C versus persistent spatial turnover -> S.
4. **Breeding system and dispersal.** Add self-compatibility/autonomous selfing/clonality and dispersal/gene-flow proxies; test whether they modify fragmentation effects.
5. **Neutral-vs-environmental decomposition.** For systems with population genetics, compare phenotypic spatial differentiation with neutral differentiation rather than inferring selection from geography alone.

## 6. Interpretation boundary

The current v1 result does **not** establish that temporal heterogeneity maintains FCP. It establishes a narrower point: once total niche breadth and total hypervolume fail, temperature seasonality is the most promising existing-data direction, but it is confounded enough that dynamic temporal climate is now the decisive next test.

## Key references

- Narbona E, Wang H, Ortiz PL, Arista M, Imbert E. 2018. Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology* 20(S1):8–20. doi:10.1111/plb.12575.
- Sapir Y, Gallagher MK, Senden E. 2021. What Maintains Flower Colour Variation within Populations? *Trends in Ecology & Evolution* 36:507–519. doi:10.1016/j.tree.2021.01.011.
- Trunschke J, Lunau K, Pyke GH, Ren Z-X, Wang H. 2021. Flower Color Evolution and the Evidence of Pollinator-Mediated Selection. *Frontiers in Plant Science* 12:617851. doi:10.3389/fpls.2021.617851.
- Frey FM. 2004. Opposing natural selection from herbivores and pathogens may maintain floral-color variation in *Claytonia virginica*. *Evolution* 58:2426–2437. doi:10.1111/j.0014-3820.2004.tb00872.x.
- Gigord LDB, Macnair MR, Smithson A. 2001. Negative frequency-dependent selection maintains a dramatic flower color polymorphism in *Dactylorhiza sambucina*. *PNAS* 98:6253–6255. doi:10.1073/pnas.111162598.
- Vaidya P et al. 2018. Ecological causes and consequences of flower color polymorphism in a self-pollinating plant (*Boechera stricta*). *New Phytologist* 218:380–392. doi:10.1111/nph.14998.
- Grossenbacher DL et al. 2026. Soil and climate contribute to maintenance of a flower color polymorphism. *American Journal of Botany* 113:e70018. doi:10.1002/ajb2.70018.
- Fetterly E et al. 2026. Selection maintains floral color polymorphism in scarlet paintbrush, *Castilleja coccinea*, reflecting combined ecological factors. *American Journal of Botany* 113:e70094. doi:10.1002/ajb2.70094.
- Jaeger SL et al. 2026. Pollinator, herbivore, and climatic selective pressures differ across a floral color transition zone. *American Journal of Botany* 113:e70142. doi:10.1002/ajb2.70142.
- Heinze J et al. 2026. Flower color polymorphism in the peacock anemone (*Anemone pavonina*) reflects spatiotemporal variation in pollinator abundance. *American Journal of Botany* 113:e70189. doi:10.1002/ajb2.70189.
- Blonder B et al. 2014. The n-dimensional hypervolume. *Global Ecology and Biogeography* 23:595–609. doi:10.1111/geb.12146.
- Qiao H et al. 2017. A cautionary note on the use of hypervolume kernel density estimators in ecological niche modelling. *Global Ecology and Biogeography* 26:1066–1070. doi:10.1111/geb.12492.
- Blonder B. 2017. Using n-dimensional hypervolumes for species distribution modelling: a response to Qiao et al. *Global Ecology and Biogeography*. doi:10.1111/geb.12611.
- Vilas D et al. 2022. Understanding the temporal dynamics of estimated environmental niche hypervolumes for marine fishes. *Ecology and Evolution*. doi:10.1002/ece3.9604.
