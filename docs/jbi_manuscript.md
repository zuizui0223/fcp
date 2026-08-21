# Climatic niche breadth and spatial organization of flower-colour variation

**Running title:** Flower-colour variation and climate

## Abstract

### Aim

Intraspecific phenotypic variation may be maintained through local coexistence or expressed as geographic differentiation, but comparative studies rarely treat this spatial organization as a distinct trait. We tested whether documented within-population flower-colour polymorphism and geographically structured flower-colour variation differ in species-level occupied climatic niche breadth.

### Location

Global, literature-derived sample.

### Taxon

Angiosperms with documented natural intraspecific flower-colour variation.

### Methods

We combined source-traceable, rule-derived literature classifications with GBIF occurrences and WorldClim 2.1 climate data. The frozen comparison comprised 34 species from 25 families: 20 within-population and 14 geographically structured cases. At a minimum of 20 occupied climate cells, we evaluated temperature breadth, moisture breadth, climatic heterogeneity, PCA dispersion and PCA hull area in separate binomial models controlling for occurrence effort. Inference used family-clustered sandwich uncertainty, 9,999 label permutations, leave-one-family-out refits and Holm multiplicity context. We also evaluated predictor collinearity, Open Tree and time-scaled phylogenetic logistic models, CR2/Satterthwaite finite-cluster inference and design-based power/precision diagnostics.

### Results

All five primary point estimates were below one. Moisture breadth showed the largest contrast (odds ratio = 0.412, 95% confidence interval = 0.180–0.947; family-clustered Wald p = 0.0368; permutation p = 0.0423), followed by PCA hull area (odds ratio = 0.577). The remaining odds ratios were 0.681 for climatic heterogeneity, 0.712 for PCA dispersion and 0.817 for temperature breadth. Every leave-one-family-out estimate remained below one. Predictor VIFs were low (1.095–1.415). After Holm correction, moisture-breadth support weakened (Wald p = 0.184; permutation p = 0.212). Open Tree and time-scaled phylogenetic analyses retained negative point estimates for all five metrics, but all phylogenetic confidence intervals included one. For moisture breadth, CR2/Satterthwaite inference gave an odds ratio of 0.407 and p = 0.063.

### Main conclusions

Across five complementary summaries, point estimates indicated lower odds of geographically structured flower-colour variation as sampled occupied climatic breadth increased, with the strongest observed contrast along the moisture axis. This direction is not what a simple broad-gradient environmental-sorting explanation would predict and is more compatible with an alternative expectation in which broader ecological occupancy accompanies local maintenance of intraspecific variation. Direction was stable to family deletion and two phylogenetic treatments, whereas multiplicity, phylogenetic and finite-sample analyses showed substantial inferential uncertainty. The result is therefore a comparative ecological pattern, not evidence for a uniquely established moisture mechanism. Species-level occupied climatic breadth is not physiological tolerance, and the analysis does not test morph-specific adaptation or climatic causation.

**Keywords:** climatic niche breadth, evidence synthesis, flower-colour variation, GBIF, geographic differentiation, intraspecific polymorphism, macroecology

## Introduction

Intraspecific phenotypic variation can be expressed through local coexistence, geographic differentiation or a combination of both. These configurations are biologically distinct because the processes capable of maintaining diversity within a population are not identical to those that maintain differences among populations. Plant balancing-selection theory distinguishes mechanisms such as negative frequency-dependent selection, which can maintain polymorphism within populations, from spatially varying selection, which can maintain adaptive differences among populations when not overwhelmed by migration (Delph & Kelly, 2014). Local coexistence and geographic differentiation can therefore be treated as alternative spatial organizations of intraspecific variation rather than as interchangeable descriptions of the same phenomenon.

Flower colour provides a tractable system for examining this distinction. Floral pigmentation can affect interactions with pollinators, but colour evolution may also reflect abiotic selection, pleiotropic effects of pigment-pathway genes, genetic drift, gene flow and mating-system processes (Rausher, 2008; Trunschke et al., 2021; Wessinger & Rausher, 2012). Reviews of flower-colour polymorphism emphasize that pollinator-mediated selection is only one component of a broader set of biotic, abiotic and demographic processes, and that the balance among those processes can vary across populations and colour variants (Narbona et al., 2018; Trunschke et al., 2021). Flower colour is therefore both an ecologically meaningful trait and a useful model for asking how intraspecific diversity is arranged geographically.

Terminology is central to this comparison. Flower-colour polymorphism is conventionally defined as the coexistence of at least two discrete flower-colour variants in the same population (Narbona et al., 2018). We use *intraspecific flower-colour variation* as the umbrella term. *Within-population flower-colour polymorphism* denotes documented local coexistence of discrete variants, whereas *geographically structured flower-colour variation* denotes differentiation among populations or regions without retained evidence of local coexistence. Cases supported at both spatial scales are classified as mixed and do not enter the binary comparison.

Previous work has established several pieces of the broader problem without directly comparing these two spatial states across species. Macroecological studies have related floral pigmentation or average flower-colour properties to geography, temperature, precipitation, radiation and biotic context (Dalrymple et al., 2020; Koski & Ashman, 2016). Individual systems show that colour-morph composition or morph-specific realized niches can vary with climate; for example, climatic variables predict geographic morph composition in *Podarcis muralis* and distinguish melanistic and non-melanistic niches in *Vipera berus* (Pérez i de Lanuza et al., 2018; Spaseni et al., 2025). Separately, the niche-variation tradition links ecological breadth with the amount of intraspecific variation, and modern plant syntheses emphasize that intraspecific trait variation can influence performance, distributions and niche occupancy while remaining strongly structured across space (Van Valen, 1965; Westerband et al., 2021). Species-level niche breadth is also positively associated with geographic range size across diverse taxa and niche dimensions (Slatyer et al., 2013). These literatures make climatic breadth a plausible macroecological correlate of the spatial organization of variation, but they do not uniquely predict its direction.

We therefore considered two competing ecological expectations. Under a **broad-gradient environmental-sorting expectation**, species spanning wider macroclimatic conditions provide greater opportunity for spatially varying selection and geographic replacement of colour variants; geographically structured cases could consequently occupy climatic niches at least as broad as, or broader than, within-population cases. Under a **broad-niche coexistence expectation**, inspired by the niche-variation tradition and mechanisms that maintain local polymorphism, broader ecological occupancy may instead accompany the maintenance of multiple phenotypes within populations, predicting broader occupied climatic niches for within-population cases and narrower niches for geographically structured cases (Van Valen, 1965; Delph & Kelly, 2014; Westerband et al., 2021). Neither expectation implies a single causal pathway: geographic differentiation can also arise through fine-scale environmental differences, demographic history, restricted dispersal or drift, whereas local coexistence can reflect frequency dependence, opposing selective agents or microenvironmental heterogeneity.

Our contribution is therefore not another general test of whether flower colour covaries with climate. We ask whether **the spatial organization of intraspecific phenotypic diversity itself** — local coexistence versus geographic differentiation — is associated with species-level occupied climatic niche breadth. We evaluated five climatic summaries symmetrically rather than selecting a single climatic axis a priori. Moisture breadth is consequently interpreted as the strongest observed component of the resulting cross-metric pattern, not as a prospectively preregistered endpoint or an established mechanism. Throughout, we distinguish species-level macroecological association from morph-level mechanism.

## Methods

### Study design and evidence reduction

The study population comprised documented natural cases of intraspecific flower-colour variation identified by the repository evidence workflow. The original discovery path retained 1,075 literature works and mapped them to 664 candidate species from 140 families. Direct-evidence screening produced an initial 72-species review queue; targeted follow-up and evidence aggregation produced a resolved 111-species queue. Mixed, unclear and otherwise non-binary cases were not forced into the final response. After climatic eligibility was imposed, the baseline-unambiguous comparison was frozen at 34 species from 25 families: 20 within-population and 14 geographically structured cases.

The candidate pool was assembled by literature discovery and targeted evidence follow-up rather than sampled randomly from angiosperms. The inferential population is therefore the assembled set of documented cases, not all angiosperms. The analysis does not estimate the global prevalence of flower-colour variation or the prevalence of either spatial configuration.

A later systematic-map search substantially expanded bibliographic search coverage, but it was developed after the original evidence path and its unreviewed expanded species sets are not used in the present downstream analysis. This distinction prevents search-completeness infrastructure from being represented as if it were the deterministic selection path that produced the frozen 34-species sample.

### Literature discovery and classification

Initial automated discovery queried the OpenAlex works index (Priem et al., 2022) with eight English-language expressions: `"flower color polymorphism"`, `"flower colour polymorphism"`, `"floral color polymorphism"`, `"floral colour polymorphism"`, `"flower color variation" pollinator`, `"flower colour variation" pollinator`, `"floral color morph"` and `"floral colour morph"`. Implemented defaults retrieved up to two cursor-paginated pages of 200 records for each query, and work identifiers were deduplicated across queries.

Candidate binomials were matched against the repository angiosperm census. A species mention was retained only when the binomial occurred in the title or near flower-, colour-, pigment-, morph- or polymorphism-related language in available abstract metadata. Strong flower-colour context increased the evidence score, whereas titles indicating contexts unlikely to document natural flower-colour variation were penalized. Source titles, DOI or OpenAlex identifiers and supporting text were retained for traceability. English search phrases and dependence on indexed title and abstract metadata may have produced language and database-coverage bias.

Spatial labels were assigned by rule-based screening of retained source text. A within-population label required a within-population signal without a geographic signal; an among-population label required the converse; simultaneous signals produced `mixed`, and absence of both produced `unclear`. Mixed and unclear records were excluded from binary comparative models. A rule audit identified that the phrase `within populations` had not been recognized by an earlier expression; the rule was corrected and the affected case moved to `mixed` before the final baseline was frozen.

The resolved evidence queue retains `review_status = unreviewed`, and the repository does not document completed independent screening, duplicate review or formal adjudication for every included classification. We therefore describe the retained states as **source-traceable, rule-derived classifications**, not independently human-validated annotations. A blinded review sheet and separate rule key are provided as submission-facing review infrastructure but are not represented as completed review.

### Occurrence data and occupied climate

The occurrence dataset underlying the frozen climatic summaries was produced from GBIF coordinate-bearing present records. The original primary retrieval used a deterministic first-page sample of at most 300 records per taxon. Invalid coordinates, the coordinate origin and duplicate coordinate pairs rounded to 0.001 degrees were removed before climate extraction.

We extracted WorldClim 2.1 bioclimatic data (Fick & Hijmans, 2017). Nine variables were retained: annual mean temperature (BIO1), temperature seasonality (BIO4), maximum temperature of the warmest month (BIO5), minimum temperature of the coldest month (BIO6), annual temperature range (BIO7), annual precipitation (BIO12), precipitation of the driest month (BIO14), precipitation seasonality (BIO15) and precipitation of the driest quarter (BIO17).

Records lacking any retained climate value were removed. Within species, identical nine-variable climate vectors were deduplicated and treated as occupied climate cells. The operational unit is therefore an occupied combination of raster-cell climate values, not an independently estimated biological population. The nine climatic variables were standardized across the occurrence dataset used to construct the summaries, and a principal components analysis was fitted to the standardized values. The first three axes explained 45.41%, 31.50% and 13.89% of total variance, respectively.

### Climatic niche metrics

Five species-level metrics described realized occupied climatic breadth or heterogeneity. Temperature breadth was the mean difference between the 95th and 5th percentiles of BIO1, BIO5, BIO6 and BIO7. Moisture breadth was the corresponding mean percentile range for BIO12, BIO14, BIO15 and BIO17. Climatic heterogeneity was the mean within-species standard deviation across the nine standardized bioclimatic variables. PCA dispersion was the mean Euclidean distance between occupied cells and the species centroid in the first three principal-component dimensions. PCA hull area was the convex-hull area of occupied cells in the first two principal-component dimensions.

The final comparison required at least 20 occupied climate cells and evaluated all five metrics on the same 34 species. These metrics describe realized occupied climate represented in sampled records. They do not estimate fundamental physiological tolerance and were not calculated separately for flower-colour morphs.

### Primary models, permutation and multiplicity

For each climatic metric we fitted a binomial generalized linear model with logit link:

`among ~ metric_z + effort_z`

The response equalled one for geographically structured variation and zero for within-population polymorphism. `metric_z` was the standardized focal climatic metric and `effort_z` was standardized `log1p(n_climate_cells)`. Models were fitted with `statsmodels` 0.14.6 (Seabold & Perktold, 2010). Wald standard errors used family-clustered sandwich covariance, and 95% confidence intervals were calculated on the log-odds scale and exponentiated.

The five climatic metrics were fitted separately, so collinearity among the climatic metrics cannot destabilize coefficients within a fitted model. For each model we nevertheless quantified the correlation between `metric_z` and `effort_z`, variance inflation factors and the condition number of the standardized design matrix.

For each metric we used 9,999 label permutations with fixed seed 20260719. Species were placed in a canonical species-name order before permutation so the finite Monte Carlo sequence is independent of CSV row order. Spatial-category labels were shuffled while climatic metrics, occurrence effort and family labels remained fixed. The two-sided permutation p-value compared the absolute observed coefficient with absolute permuted coefficients.

We refitted each model after omitting each represented plant family in turn. Holm-adjusted p-values across the five climatic metrics are reported for clustered Wald and permutation tests. Because the five metrics are correlated descriptions of climatic space, Holm adjustment is treated as multiplicity context rather than as a substitute for effect sizes, confidence intervals and cross-analysis directional consistency.

### Phylogenetic sensitivity

For topology-based phylogenetic sensitivity, species names were matched to Open Tree Taxonomy without approximate matching. Thirty of 34 species were retained. We induced the Open Tree topology (Hinchliff et al., 2015), randomly resolved polytomies 100 times, assigned Grafen branch lengths and fitted `phyloglm` logistic MPLE models (Ho & Ané, 2014) with the same `among ~ metric_z + effort_z` formula.

We also used V.PhyloMaker2 (Jin & Qian, 2019, 2022) and the time-scaled `GBOTB.extended.LCVP` backbone, which incorporates the broad seed-plant phylogeny of Smith and Brown (2018) and earlier time-calibrated plant-tree information (Zanne et al., 2014). All 34 species were retained. Placement scenarios S1, S2 and S3 were evaluated separately for all five metrics with fixed seed 20260724. These models are sensitivity analyses, not evidence that the residual process necessarily follows either assumed tree model.

### Finite-sample sensitivity

Because the primary comparison contains 34 species clustered in 25 families, we supplemented the family-clustered sandwich analysis with CR2 covariance correction and Satterthwaite degrees of freedom. We also used a design-based simulation that retained the observed 34-species predictor and effort structure to characterize expected sign recovery, nominal Wald coverage, interval width and the probability of p < 0.05 under the observed moisture-breadth effect. These diagnostics quantify finite-sample uncertainty; they are not post-hoc criteria for declaring the sample adequate.

## Results

### Evidence base and frozen comparison

The retained evidence path contained 664 candidate species from 140 families and a resolved 111-species review queue. After the binary spatial-state and climatic-eligibility requirements were applied, the final comparison contained 34 species from 25 families: 20 within-population and 14 geographically structured cases. Because discovery effort was literature-dependent and the candidates were not a random sample of angiosperms, these counts should not be interpreted as prevalence estimates.

The geographic context of the 34 focal species is shown in Figure 1, and species-specific occurrence maps are provided as Supporting Figure S1. Those maps use the broader exact GBIF occurrence subset retained for citation and distribution-context auditing; they are not presented as the exact occurrence sample that generated the checksum-locked climatic metrics, nor are the records labelled by flower-colour morph.

### Five climatic-niche metrics

All five climatic-metric odds ratios were below one (Table 2). Temperature breadth had an odds ratio of 0.817 (95% CI 0.384–1.739; clustered Wald p = 0.6000; permutation p = 0.6131). Moisture breadth showed the largest negative estimate, with an odds ratio of 0.412 (0.180–0.947; clustered p = 0.0368; permutation p = 0.0423). Climatic heterogeneity had an odds ratio of 0.681 (0.294–1.577; p = 0.3700; permutation p = 0.3567), PCA dispersion 0.712 (0.306–1.660; p = 0.4317; permutation p = 0.3859), and PCA hull area 0.577 (0.312–1.067; p = 0.0797; permutation p = 0.2372).

Every leave-one-family-out estimate remained below one for every metric. Odds-ratio ranges were 0.636–0.973 for temperature breadth, 0.306–0.465 for moisture breadth, 0.492–0.782 for climatic heterogeneity, 0.525–0.813 for PCA dispersion and 0.489–0.671 for PCA hull area. Thus, no single represented family generated the shared negative direction.

Figure 4 displays all 125 leave-one-family-out refits across the five climatic metrics. The figure makes the directional stability visible while retaining the important distinction that family deletion diagnoses concentration in individual families and is not a phylogenetic correction.

After Holm adjustment across five metrics, the moisture-breadth clustered Wald p-value was 0.184 and the permutation p-value was 0.212. No metric retained conventional statistical support after this correction. Moisture breadth is therefore the strongest observed association within a broader directional pattern, not a uniquely established climatic driver.

The common direction of the five effect estimates is summarized in Figure 2. Figure 3 shows the 34 species themselves across the five standardized climatic metrics, so the overlap, extreme observations and modest sample size remain visible rather than being represented only by model coefficients.

### Collinearity diagnostics

Correlations between the focal climatic metric and occurrence effort ranged from 0.295 to 0.541. Maximum VIFs ranged from 1.095 to 1.415 and condition numbers from 1.355 to 1.833. These values provide no indication of problematic predictor collinearity with the effort covariate.

### Phylogenetic sensitivity

Open Tree matching retained 30 species and 100 fits completed for each metric. Every point estimate was negative. Median odds ratios were 0.908 for temperature breadth, 0.573 for moisture breadth, 0.726 for climatic heterogeneity, 0.819 for PCA dispersion and 0.654 for PCA hull area. All median 95% confidence intervals included one.

V.PhyloMaker2 retained all 34 species under S1–S3. The five estimates remained below one in every scenario. Across scenarios, temperature-breadth odds ratios were 0.827–0.838, moisture-breadth odds ratios 0.448–0.454, climatic-heterogeneity odds ratios 0.670–0.679, PCA-dispersion odds ratios 0.698–0.707 and PCA-hull-area odds ratios 0.599–0.612. All phylogenetic confidence intervals included one. For moisture breadth, scenario-specific p-values were 0.112–0.115 and Holm-adjusted p-values 0.560–0.574.

### Finite-sample sensitivity

CR2/Satterthwaite estimates retained the negative direction for all five metrics. For moisture breadth, the odds ratio was 0.407 (95% CI 0.156–1.060; p = 0.0632; Satterthwaite df = 11.48). The corresponding CR2 p-values were 0.6245 for temperature breadth, 0.3962 for climatic heterogeneity, 0.4550 for PCA dispersion and 0.1174 for PCA hull area.

Figure 5 compares the primary family-clustered estimates with CR2/Satterthwaite, Open Tree/Grafen and dated-phylogeny treatments for all five metrics. Point estimates remain below one across treatments, whereas confidence intervals broaden under the phylogenetic and finite-cluster analyses.

In 3,000 design-based simulations under the observed moisture-breadth effect, the estimated coefficient was negative in 98.5% of simulations, whereas p < 0.05 occurred in 46.2%. Nominal Wald 95% coverage was 94.6%. This diagnostic indicates that directional recovery under an effect of the observed magnitude can be substantially more stable than conventional significance in a 34-species design.

Supporting Figure S2 expands this design diagnostic across the prespecified odds-ratio grid for all five metrics. It is presented only as a precision diagnostic and not as evidence for the ecological hypothesis.

## Discussion

The principal result is a cross-metric directional pattern rather than a single isolated p-value. In the frozen 34-species comparison, all five point estimates indicated decreasing odds of geographically structured flower-colour variation as sampled occupied climatic breadth increased. Every leave-one-family-out estimate retained the same direction. The largest contrast occurred for moisture breadth, for which the unadjusted family-clustered interval excluded one and the permutation p-value was below 0.05. That support did not survive correction across the five metrics, the CR2 interval included one, and all phylogenetic intervals included one. The data therefore support a directionally coherent comparative signal with moisture as its strongest observed component, not a confirmed moisture-specific mechanism.

The direction of this pattern is ecologically informative because it distinguishes the two expectations posed above. A simple broad-gradient sorting explanation would make geographic differentiation most intuitive when a species spans a broad macroclimatic range across which alternative morphs can be favoured in different places. Morph-level studies demonstrate that such climatic segregation can occur within polymorphic species (Pérez i de Lanuza et al., 2018; Spaseni et al., 2025). In our cross-species comparison, however, geographically structured cases occurred toward the **narrower**, not broader, end of species-level climatic breadth across all five metrics. The result is therefore directionally more compatible with the broad-niche coexistence expectation than with simple broad-gradient climatic sorting. This is a comparison of predictions, not evidence that broad climatic niches themselves maintain polymorphism.

This distinction also changes how geographic flower-colour differentiation should be interpreted. Geographic differentiation need not require large species-level macroclimatic differences. It can arise from environmental contrasts below the resolution of WorldClim, spatially varying biotic interactions, demographic history, founder effects, restricted dispersal, migration–selection balance or drift in subdivided populations (Delph & Kelly, 2014; Narbona et al., 2018). Conversely, local coexistence may be maintained by negative frequency dependence, opposing selective agents, temporal fluctuation or microhabitat heterogeneity. None of those mechanisms was measured here. The present result instead identifies a macroecological contrast that narrows the set of explanations to test with population-level data.

The contrast among climatic metrics adds ecological resolution without isolating a causal axis. Temperature breadth was weakly associated with spatial organization, whereas moisture breadth produced the largest effect size and PCA hull area the second largest. Moisture has unusually direct precedent in flower-colour ecology. Across five polymorphic plant species, pigmented morphs performed relatively better under drought and unpigmented morphs relatively better under well-watered conditions (Warren & Mackenzie, 2001). Long-term work on *Linanthus parryae* found flower-colour selection that changed through time and space and was associated with precipitation, while reciprocal-transplant work linked local morph advantage with environmental differences (Schemske & Bierzychudek, 2001, 2007). In *Boechera stricta*, drought induced floral pigmentation in white-flowered lineages and white flowers had a fecundity advantage under well-watered conditions, alongside effects of herbivory and elevation (Vaidya et al., 2018). These studies make moisture a biologically plausible axis for follow-up, but our GBIF occurrences are not morph-labelled and our moisture metric describes species-level occupied breadth. Moisture is therefore best treated as a **priority axis for mechanistic testing**, not as an established cause of the spatial pattern observed here.

The phylogenetic analyses sharpen rather than overturn this interpretation. Both Open Tree + Grafen and time-scaled V.PhyloMaker2 models retained negative estimates for all five metrics. In the dated analysis, the moisture-breadth odds ratio remained close to the non-phylogenetic estimate (approximately 0.45 versus 0.41), while its confidence interval included one. This persistence argues against describing the observed direction as a simple consequence of one represented lineage, but the models do not establish phylogenetically independent statistical significance.

Finite-sample diagnostics clarify another part of the uncertainty. CR2/Satterthwaite inference moved the moisture result from a conventional p < 0.05 family-clustered Wald result to p = 0.063, and the design-based simulation recovered the negative sign much more often than it crossed a p < 0.05 threshold. This is consistent with the interpretation we adopt throughout: effect direction and magnitude are more informative here than binary significance. The simulation is diagnostic and does not rescue or validate a marginal result.

Collinearity is unlikely to explain the effect ordering. Each climatic metric was analysed in a separate model with occurrence effort, and VIFs were at most 1.415. The cost of separate models is multiplicity, which is why the five estimates are interpreted jointly and Holm-adjusted values are reported transparently.

Classification uncertainty remains an important limitation. The spatial categories are source-traceable and rule-derived, not fully independently adjudicated biological annotations. Literature-derived spatial traits are measurements with error, and source context, extraction rules and review status can affect downstream macroecological inference. The evidence-extraction workflow is useful because it preserves provenance and uncertainty; it should not be treated as creating ground-truth biological labels automatically. Completing the prepared blinded review would strengthen this component before submission.

Occurrence data impose a second boundary. The climatic metrics describe species-level realized occupancy sampled through GBIF, not physiological tolerance, morph-specific niches or local selective environments. The records are not labelled by flower-colour morph. Consequently, the present analysis cannot demonstrate that one morph occupies drier, wetter, warmer or otherwise distinct local conditions. A mechanistic test would require georeferenced morph frequencies or named populations, environmental measurements linked to those states and, ideally, demographic or genomic evidence addressing selection and gene flow.

The study therefore connects evidence synthesis to a general biogeographic question about the **spatial organization of intraspecific diversity**. Previous studies have asked why particular colour morphs occur in particular environments; here, the comparative response is whether colour variants are documented as coexisting locally or differentiated geographically. In the current sample, all five estimates point toward geographically structured variation occurring toward the narrower end of sampled climatic breadth. That direction is more compatible with broad-niche coexistence than with a simple broad-gradient sorting expectation, while remaining agnostic about the causal mechanisms that generate either spatial state. The combination of an explicit spatial trait, reproducible evidence provenance, competing ecological expectations and transparent uncertainty is the central contribution.

## References

Dalrymple, R. L., Kemp, D. J., Flores-Moreno, H., Laffan, S. W., White, T. E., Hemmings, F. A., & Moles, A. T. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist, 228*(6), 1972–1985. https://doi.org/10.1111/nph.16737

Delph, L. F., & Kelly, J. K. (2014). On the importance of balancing selection in plants. *New Phytologist, 201*(1), 45–56. https://doi.org/10.1111/nph.12441

Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology, 37*, 4302–4315. https://doi.org/10.1002/joc.5086

Hinchliff, C. E., Smith, S. A., Allman, J. F., Burleigh, J. G., Chaudhary, R., Coghill, L. M., Crandall, K. A., Deng, J., Drew, B. T., Gazis, R., Gude, K., Hibbett, D. S., Katz, L. A., Laughinghouse, H. D., IV, McTavish, E. J., Midford, P. E., Owen, C. L., Ree, R. H., Rees, J. A., Soltis, D. E., Williams, T., & Cranston, K. A. (2015). Synthesis of phylogeny and taxonomy into a comprehensive tree of life. *Proceedings of the National Academy of Sciences of the United States of America, 112*(41), 12764–12769. https://doi.org/10.1073/pnas.1423041112

Ho, L. S. T., & Ané, C. (2014). A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. *Systematic Biology, 63*(3), 397–408. https://doi.org/10.1093/sysbio/syu005

Jin, Y., & Qian, H. (2019). V.PhyloMaker: An R package that can generate very large phylogenies for vascular plants. *Ecography, 42*(8), 1353–1359. https://doi.org/10.1111/ecog.04434

Jin, Y., & Qian, H. (2022). V.PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. *Plant Diversity, 44*(4), 335–339. https://doi.org/10.1016/j.pld.2022.05.005

Koski, M. H., & Ashman, T.-L. (2016). Macroevolutionary patterns of ultraviolet floral pigmentation explained by geography and associated bioclimatic factors. *New Phytologist, 211*(2), 708–718. https://doi.org/10.1111/nph.13921

Narbona, E., Wang, H., Ortiz, P. L., Arista, M., & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: Occurrence, maintenance and implications for speciation. *Plant Biology, 20*(S1), 8–20. https://doi.org/10.1111/plb.12575

Pérez i de Lanuza, G., Sillero, N., & Carretero, M. Á. (2018). Climate suggests environment-dependent selection on lizard colour morphs. *Journal of Biogeography, 45*, 2791–2802. https://doi.org/10.1111/jbi.13455

Priem, J., Piwowar, H., & Orr, R. (2022). OpenAlex: A fully-open index of scholarly works, authors, venues, institutions, and concepts. *arXiv*. https://doi.org/10.48550/arXiv.2205.01833

Rausher, M. D. (2008). Evolutionary transitions in floral color. *International Journal of Plant Sciences, 169*(1), 7–21. https://doi.org/10.1086/523358

Schemske, D. W., & Bierzychudek, P. (2001). Perspective: Evolution of flower color in the desert annual *Linanthus parryae*: Wright revisited. *Evolution, 55*(7), 1269–1282. https://doi.org/10.1111/j.0014-3820.2001.tb00650.x

Schemske, D. W., & Bierzychudek, P. (2007). Spatial differentiation for flower color in the desert annual *Linanthus parryae*: Was Wright right? *Evolution, 61*, 2528–2543. https://doi.org/10.1111/j.1558-5646.2007.00219.x

Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and statistical modeling with Python. In *Proceedings of the 9th Python in Science Conference* (pp. 92–96).

Slatyer, R. A., Hirst, M., & Sexton, J. P. (2013). Niche breadth predicts geographical range size: A general ecological pattern. *Ecology Letters, 16*, 1104–1114. https://doi.org/10.1111/ele.12140

Smith, S. A., & Brown, J. W. (2018). Constructing a broadly inclusive seed plant phylogeny. *American Journal of Botany, 105*(3), 302–314. https://doi.org/10.1002/ajb2.1019

Spaseni, P., Gherghel, I., Sahlean, T. C., Zamfirescu, Ș. R., Vicol, A. V., & Strugariu, A. (2025). Niche differentiation and spatial segregation of melanistic and non-melanistic colour morphs in a widespread Palearctic reptile. *Journal of Biogeography, 52*, e15150. https://doi.org/10.1111/jbi.15150

Trunschke, J., Lunau, K., Pyke, G. H., Ren, Z.-X., & Wang, H. (2021). Flower color evolution and the evidence of pollinator-mediated selection. *Frontiers in Plant Science, 12*, 617851. https://doi.org/10.3389/fpls.2021.617851

Vaidya, P., McDurmon, A., Mattoon, E., Keefe, M., Carley, L., Lee, C.-R., Bingham, R., & Anderson, J. T. (2018). Ecological causes and consequences of flower color polymorphism in a self-pollinating plant (*Boechera stricta*). *New Phytologist, 218*, 380–392. https://doi.org/10.1111/nph.14998

Van Valen, L. (1965). Morphological variation and width of ecological niche. *The American Naturalist, 99*, 377–390. https://doi.org/10.1086/282379

Warren, J. M., & Mackenzie, S. (2001). Why are all colour combinations not equally represented as flower-colour polymorphisms? *New Phytologist, 151*(1), 237–241. https://doi.org/10.1046/j.1469-8137.2001.00159.x

Wessinger, C. A., & Rausher, M. D. (2012). Lessons from flower colour evolution on targets of selection. *Journal of Experimental Botany, 63*(16), 5741–5749. https://doi.org/10.1093/jxb/ers267

Westerband, A. C., Funk, J. L., & Barton, K. E. (2021). Intraspecific trait variation in plants: A renewed focus on its role in ecological processes. *Annals of Botany, 127*(4), 397–410. https://doi.org/10.1093/aob/mcab011

Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., FitzJohn, R. G., McGlinn, D. J., O'Meara, B. C., Moles, A. T., Reich, P. B., Royer, D. L., Soltis, D. E., Stevens, P. F., Westoby, M., Wright, I. J., Aarssen, L., Bertin, R. I., Calaminus, A., Govaerts, R., . . . Beaulieu, J. M. (2014). Three keys to the radiation of angiosperms into freezing environments. *Nature, 506*(7486), 89–92. https://doi.org/10.1038/nature12872

## Data Accessibility Statement

The canonical downstream dataset is committed as `data/frozen/frozen_34species_five_metric_dataset.csv` (SHA-256 `bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed`) with its freeze metadata in `data/frozen/freeze_manifest.json`. Analysis code, evidence provenance, classification manifests, phylogenetic inputs and submission-facing audits are maintained in the public GitHub repository `zuizui0223/fcp`. The canonical executable entry point is `.github/workflows/34species-paper.yml`.

Generated workflow artifacts are convenience copies of analysis outputs rather than durable citation metadata: workflow-run identifiers change across valid reruns and artifacts expire. Reproduction of the downstream paper analysis instead depends on the committed checksum-locked dataset, versioned code and fixed model seeds recorded in the repository.

**Not verified before submission:** archive the exact final repository release at a permanent DOI and obtain the required citable GBIF occurrence identifier / Derived Dataset registration for the occurrence data underlying the climatic summaries.

## Tables

### Table 1. Operational classification of documented spatial organization

| Category | Operational evidence requirement | Included in binary models |
|---|---|---|
| Within-population flower-colour polymorphism | Explicit retained source evidence that at least two discrete natural flower-colour variants coexist within at least one population | Yes; response = 0 |
| Geographically structured flower-colour variation | Explicit retained evidence of differentiation among populations or regions, without retained evidence of local coexistence | Yes; response = 1 |
| Mixed | Retained evidence supports both within-population coexistence and geographic differentiation | No |
| Unclear | Evidence does not resolve the spatial configuration | No |

### Table 2. Frozen 34-species environmental-niche models

| Metric | Species | Families | Within / among | Odds ratio | Family-clustered 95% CI | Wald p | Permutation p | Holm Wald p | Holm permutation p | Leave-one-family-out OR range |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Temperature breadth | 34 | 25 | 20 / 14 | 0.817 | 0.384–1.739 | 0.6000 | 0.6131 | 1.0000 | 1.0000 | 0.636–0.973 |
| Moisture breadth | 34 | 25 | 20 / 14 | 0.412 | 0.180–0.947 | 0.0368 | 0.0423 | 0.1840 | 0.2115 | 0.306–0.465 |
| Climatic heterogeneity | 34 | 25 | 20 / 14 | 0.681 | 0.294–1.577 | 0.3700 | 0.3567 | 1.0000 | 1.0000 | 0.492–0.782 |
| PCA dispersion | 34 | 25 | 20 / 14 | 0.712 | 0.306–1.660 | 0.4317 | 0.3859 | 1.0000 | 1.0000 | 0.525–0.813 |
| PCA hull area | 34 | 25 | 20 / 14 | 0.577 | 0.312–1.067 | 0.0797 | 0.2372 | 0.3189 | 0.9488 | 0.489–0.671 |

*Note.* Odds ratios are for geographically structured rather than within-population organization per one standard-deviation increase in the focal climatic metric, controlling for standardized `log1p(n_climate_cells)`. The five climatic metrics were fitted separately. Holm adjustments are across the five metrics.

### Table 3. Collinearity diagnostics

| Metric | Correlation with effort | Max VIF | Condition number |
|---|---:|---:|---:|
| Temperature breadth | 0.329 | 1.121 | 1.407 |
| Moisture breadth | 0.298 | 1.097 | 1.360 |
| Climatic heterogeneity | 0.412 | 1.204 | 1.549 |
| PCA dispersion | 0.295 | 1.095 | 1.355 |
| PCA hull area | 0.541 | 1.415 | 1.833 |

### Table 4. Phylogenetic sensitivity

| Metric | Open Tree n | Open Tree median OR (95% CI) | V.PhyloMaker2 S1 OR (95% CI) | S2 OR (95% CI) | S3 OR (95% CI) |
|---|---:|---:|---:|---:|---:|
| Temperature breadth | 30 | 0.908 (0.436–1.893) | 0.827 (0.385–1.777) | 0.838 (0.393–1.789) | 0.827 (0.385–1.777) |
| Moisture breadth | 30 | 0.573 (0.234–1.403) | 0.448 (0.167–1.206) | 0.454 (0.171–1.211) | 0.448 (0.167–1.206) |
| Climatic heterogeneity | 30 | 0.726 (0.322–1.636) | 0.670 (0.290–1.550) | 0.679 (0.295–1.565) | 0.670 (0.290–1.550) |
| PCA dispersion | 30 | 0.819 (0.387–1.731) | 0.698 (0.319–1.527) | 0.707 (0.324–1.542) | 0.698 (0.319–1.527) |
| PCA hull area | 30 | 0.654 (0.261–1.641) | 0.599 (0.232–1.545) | 0.612 (0.242–1.549) | 0.599 (0.232–1.545) |

*Note.* Open Tree values summarize 100 completed fits after polytomy resolution and Grafen branch-length assignment. V.PhyloMaker2 values retain all 34 species under three placement scenarios. All phylogenetic confidence intervals include one.

### Table 5. CR2/Satterthwaite finite-cluster sensitivity

| Metric | OR | 95% CI | p | Satterthwaite df |
|---|---:|---:|---:|---:|
| Temperature breadth | 0.814 | 0.317–2.090 | 0.6245 | 7.29 |
| Moisture breadth | 0.407 | 0.156–1.060 | 0.0632 | 11.48 |
| Climatic heterogeneity | 0.677 | 0.255–1.797 | 0.3962 | 10.46 |
| PCA dispersion | 0.708 | 0.266–1.886 | 0.4550 | 11.09 |
| PCA hull area | 0.572 | 0.276–1.188 | 0.1174 | 8.78 |

*Note.* Exact machine-readable CR2/Satterthwaite estimates are mirrored in `docs/supporting/cr2_satterthwaite_summary.csv` and checked by manuscript consistency CI.

## Figure captions

**Figure 1. Geographic context of the 34 focal species.** The world map shows the broader exact GBIF occurrence subset retained for the same 34 focal taxa, with symbols distinguishing species classified as within-population coexistence or geographic differentiation. The species strip reports the number of occupied climate cells used in the checksum-locked primary climatic analysis. The mapped occurrence archive is used for geographic context and auditability; it is not represented as the exact occurrence sample that created the frozen climatic metrics and is not morph-labelled.

**Figure 2. Five climatic-niche metrics show the same direction of association with spatial organization.** Points are odds ratios from the production binomial models (`among ~ metric_z + effort_z`) and horizontal intervals are 95% Wald confidence intervals from family-clustered sandwich standard errors. Odds ratios below one indicate lower odds of geographically structured rather than within-population flower-colour variation as occupied climatic breadth increases. The figure is selected as the central result because it displays all five prespecified climatic summaries symmetrically rather than privileging the strongest unadjusted p-value.

**Figure 3. Species-level climatic metrics underlying the comparative models.** Each point represents one of the 34 focal species (20 within-population, 14 geographically structured). Climatic metrics are standardized for visualization only so that all five dimensions can be displayed on a common scale; horizontal bars mark category medians. This is a visualization of the frozen observations, not an additional inferential analysis.

**Figure 4. Leave-one-family-out stability across all five climatic-niche metrics.** Grey points are odds ratios from 25 unclustered refits per metric, each omitting one represented plant family; diamonds show the full 34-species point estimates and horizontal grey segments span the deletion range. All 125 deletion estimates remain below one. Family deletion tests whether the shared direction is concentrated in an individual represented family; it does not establish phylogenetic independence.

**Figure 5. Effect direction and uncertainty across inferential treatments.** For each climatic metric, the primary family-clustered estimate is shown alongside CR2/Satterthwaite finite-cluster inference, the median of 100 Open Tree/Grafen polytomy resolutions, and the median dated-phylogeny estimate across V.PhyloMaker2 scenarios S1-S3. Horizontal intervals are the corresponding 95% intervals; for the dated phylogeny the plotted interval is the envelope across the three scenario-specific intervals. All point estimates remain below one, while stricter treatments broaden uncertainty. These treatments are sensitivity analyses of the same data, not independent tests.

**Supporting Figure S1. Geographic occurrence context for each of the 34 focal species.** Species-specific maps show the broader exact GBIF occurrence subset retained for citation and distribution-context quality control, with panel labels indicating the spatial-organization class and occurrence count. These maps provide an auditable view of sampled geographic ranges but do not represent morph-specific ranges or the exact primary occurrence sample used to construct the frozen climatic metrics.

**Supporting Figure S2. Finite-sample design diagnostic across specified effect sizes.** The two panels show, for each climatic metric and simulated true odds ratio, the probability that the fitted coefficient is negative and the probability that the family-clustered Wald test yields p < 0.05 across 3,000 simulations. Simulations retain the observed 34-species predictor, effort and family structure. The figure describes expected precision under specified effects and is not evidence for the ecological hypothesis or a post-hoc adequacy criterion.
