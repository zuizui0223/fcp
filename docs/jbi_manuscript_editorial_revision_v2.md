# Climatic niche breadth and spatial organization of flower-colour variation

**Running title:** Flower-colour variation and climate

## Abstract

### Aim

Intraspecific phenotypic variation may occur through local coexistence or geographic differentiation, but comparative studies rarely treat spatial organization as a distinct property. We tested whether documented within-population flower-colour polymorphism and geographically structured flower-colour variation differ in species-level occupied climatic niche breadth.

### Location

Global, literature-derived sample.

### Taxon

Angiosperms with documented natural intraspecific flower-colour variation.

### Methods

We combined source-traceable, rule-derived literature classifications with GBIF occurrences and WorldClim 2.1 climate data. The frozen baseline comprised 34 species from 25 families (20 within-population and 14 geographically structured cases). At a minimum of 20 occupied climate cells, we evaluated five climatic niche metrics symmetrically: temperature breadth, moisture breadth, climatic heterogeneity, PCA dispersion and PCA hull area. Binomial models controlled for occurrence effort and used family-clustered sandwich uncertainty, 9,999 label permutations and leave-one-family-out refits. We diagnosed predictor collinearity and fitted phylogenetic logistic sensitivity models using both an Open Tree topology with Grafen branch lengths and time-scaled V.PhyloMaker2 trees under three placement scenarios. Holm-adjusted p-values across the five climatic metrics were reported as multiplicity context.

### Results

All five baseline effect estimates were negative. Moisture breadth showed the strongest association (odds ratio = 0.412, 95% confidence interval = 0.180–0.947; family-clustered Wald p = 0.0368; permutation p = 0.0475), followed by PCA hull area (odds ratio = 0.577). The other point estimates were 0.681 for climatic heterogeneity, 0.712 for PCA dispersion and 0.817 for temperature breadth. Every leave-one-family-out estimate remained below one for all five metrics. Predictor VIFs were low (1.095–1.415). After Holm correction across the five metrics, moisture-breadth support weakened (Wald p = 0.184; permutation p = 0.238). Open Tree models retained the negative direction for all five metrics; moisture breadth had a median odds ratio of 0.573 (95% confidence interval = 0.234–1.403). Time-scaled trees retained all 34 species and gave moisture-breadth odds ratios of 0.448–0.454 across placement scenarios; all phylogenetic confidence intervals included one.

### Main conclusions

Geographically structured flower-colour variation was consistently associated with narrower sampled occupied climatic niches across five complementary metrics, with the largest contrast for moisture breadth. The direction was stable to family deletion and two phylogenetic treatments, whereas inferential certainty weakened after multiplicity and phylogenetic correction. The results therefore support a comparative ecological pattern rather than a uniquely established moisture mechanism. Species-level occupied climatic breadth is not physiological tolerance, and the analysis does not test morph-specific adaptation or climatic causation.

**Keywords:** climatic niche breadth, evidence synthesis, flower-colour variation, GBIF, geographic differentiation, intraspecific polymorphism, macroecology

## Introduction

Intraspecific phenotypic variation can be expressed through local coexistence, geographic differentiation or a combination of both. These configurations are biologically distinct. Local coexistence requires multiple forms to persist under at least partly shared demographic and environmental conditions, whereas geographic differentiation may reflect spatial environmental variation, dispersal limitation, demographic history or other forms of regional structure. Combining these configurations into a single category can therefore obscure the spatial scale at which phenotypic variation is documented.

Flower colour provides a tractable system for examining this distinction. Floral pigmentation can affect interactions with pollinators, but colour evolution may also reflect abiotic selection, pleiotropic effects of pigment-pathway genes, genetic drift, gene flow and mating-system processes (Rausher, 2008; Trunschke et al., 2021; Wessinger & Rausher, 2012). These mechanisms need not operate uniformly across a species' range, and their relative contributions can differ among populations and colour variants (Narbona et al., 2018). Flower colour is therefore both an ecologically meaningful trait and a useful system for asking how intraspecific variation is arranged geographically.

Terminology is central to this comparison. Flower-colour polymorphism is conventionally defined as the coexistence of at least two discrete flower-colour variants in the same population (Narbona et al., 2018). We consequently use *intraspecific flower-colour variation* as the umbrella term. *Within-population flower-colour polymorphism* denotes documented local coexistence of discrete variants, whereas *geographically structured flower-colour variation* denotes differentiation among populations or regions without retained evidence of local coexistence. Cases supported at both spatial scales are classified as mixed. This hierarchy avoids calling geographically separated variants polymorphic when coexistence has not been demonstrated.

Previous comparative research has related floral pigmentation or average flower-colour properties to geography, temperature, precipitation, radiation and biotic context (Dalrymple et al., 2020; Koski & Ashman, 2016). Individual systems and regional reviews have also documented geographic variation in colour-morph frequencies (Narbona et al., 2018). These approaches establish that floral colour and its frequencies may vary geographically, but they do not test whether species-level occupied climates differ according to whether colour variants are documented as coexisting locally or differentiated geographically. The contribution of the present study is therefore not another general test of whether flower colour covaries with climate. It is a comparison of the documented spatial organization of intraspecific colour variation across taxa.

We assembled a global, literature-derived sample and asked whether documented within-population polymorphism and geographically structured variation differ in species-level occupied climatic niche breadth. The comparison between the two spatial configurations was theory-led, whereas the environmental characterization was deliberately multivariate: we evaluated temperature breadth, moisture breadth, climatic heterogeneity and two PCA-based measures of occupied climatic extent. Moisture breadth is therefore interpreted as the strongest signal within a five-metric family rather than as a prospectively preregistered endpoint. Our objective was to test whether the *spatial organization* of flower-colour variation covaries with the breadth of climates occupied by species, while keeping the boundary between species-level macroecological association and morph-level mechanism explicit.

## Methods

### Study design and inferential population

The study population comprised species identified by the repository workflow as documented natural cases of intraspecific flower-colour variation. The candidate pool was assembled from literature discovery and targeted evidence follow-up rather than sampled randomly from angiosperms. The inferential population is therefore the assembled set of documented cases. The analyses do not estimate the global prevalence of flower-colour variation or the proportions of all angiosperms showing each spatial configuration.

The resolved evidence base contained 664 candidate species from 140 families. Of these, 111 were retained as validated natural cases and 553 remained deferred. A separate ascertainment analysis related validation status to retained literature effort. Because this analysis concerns validation conditional on candidature, it characterizes the evidence-assembly process rather than biological prevalence.

### Literature discovery and validation

Initial automated discovery queried the OpenAlex works index (Priem et al., 2022) using eight English-language search expressions: `"flower color polymorphism"`, `"flower colour polymorphism"`, `"floral color polymorphism"`, `"floral colour polymorphism"`, `"flower color variation" pollinator`, `"flower colour variation" pollinator`, `"floral color morph"` and `"floral colour morph"`. The implemented defaults retrieved up to two cursor-paginated pages of 200 records for each query. Work identifiers were deduplicated across queries.

Candidate binomials were matched against the repository angiosperm census. A species mention was retained only when the binomial occurred in the title or within a 220-character window around flower-, colour-, pigment-, morph- or polymorphism-related language in an available abstract. Strong flower-colour context increased the evidence score; titles indicating genomes, transcriptomes, checklists, floras or other contexts unlikely to document natural flower-colour variation were penalized. Candidate evidence required a score of at least eight. Source titles, DOI or OpenAlex identifiers and supporting text were retained for traceability. Subsequent targeted scripts searched unresolved or high-priority species and fed the resolved review queue.

Repository history records automated global discovery and targeted follow-up or enrichment activity from 16 to 19 July 2026 (Appendix S1). No explicit language filter was coded, and the English search phrases and dependence on indexed title and abstract metadata may have produced language and database-coverage bias.

### Classification of spatial organization

Spatial labels were generated by rule-based screening of retained source text. The baseline classifier concatenated the retained title, evidence passage and review-reason fields and searched prespecified regular expressions for within-population and geographic-differentiation language. A within-population label was assigned when a within-population signal occurred without a geographic signal; an among-population label was assigned for the converse; simultaneous signals produced `mixed`, and absence of both produced `unclear`. Mixed and unclear rows were flagged as requiring manual review and excluded from binary comparative models. The implementation explicitly describes these automated labels as screening labels rather than final biological annotations.

The resolved evidence queue retained `review_status = unreviewed`, and the repository contains no field documenting completed independent human screening, duplicate review, inter-reviewer agreement or formal adjudication for every included binary classification. We therefore report the within- and among-population states as source-traceable, rule-derived classifications rather than independently human-screened annotations.

Source identifiers, evidence snippets and rule-derived decision notes were propagated into the analysis dataset and checked for traceability. A targeted rule audit detected that the phrase `within populations` had not been recognized by the initial expression. The rule was corrected and the analysis was rerun; a case with signals at both spatial scales was consequently moved to `mixed` and removed from the binary baseline set. The final baseline-unambiguous manifest was frozen without reference to the climatic model results and contains 34 model-eligible species from 25 families: 20 within-population and 14 among-population cases.

A broader evidence set supplemented ambiguous baseline records with automated high-confidence enrichment. This broader set is retained as an evidence-quality sensitivity analysis, whereas the present five-metric comparison is restricted to the frozen 34-species baseline so that all environmental metrics, family-deletion analyses and phylogenetic corrections are evaluated on the same biological sample.

### Occurrence data and sampling sensitivity

The primary GBIF workflow retrieved a deterministic first-page sample of at most 300 records per taxon using `hasCoordinate=true` and `occurrenceStatus=present`. Invalid coordinates, the coordinate origin and duplicated coordinate pairs rounded to 0.001° were removed.

A separate paginated sensitivity workflow for the 34 baseline species resolved strict GBIF backbone matches and queried accepted taxon keys with `hasCoordinate=true`, `hasGeospatialIssue=false` and `occurrenceStatus=present`, using pages of 300 and a cap of 3,000 records per species. Records with invalid coordinates, reported uncertainty above 20 km or duplicated 0.001° coordinate pairs were removed. This paginated analysis remains a sampling sensitivity for the focal moisture result rather than the basis of the symmetric five-metric baseline reported here.

**Not verified:** a citable GBIF derived-dataset DOI is required before submission.

### Climatic variables and occupied-climate summaries

We extracted WorldClim 2.1 bioclimatic data at 10 arc-min resolution (Fick & Hijmans, 2017). The nine retained variables were annual mean temperature (BIO1), temperature seasonality (BIO4), maximum temperature of the warmest month (BIO5), minimum temperature of the coldest month (BIO6), annual temperature range (BIO7), annual precipitation (BIO12), precipitation of the driest month (BIO14), precipitation seasonality (BIO15) and precipitation of the driest quarter (BIO17).

Records lacking any of the nine climatic values were removed. Within a species, records with identical nine-variable climate vectors were deduplicated and treated as occupied climate cells. This operational unit is therefore an occupied combination of raster-cell climate values rather than an independently estimated biological population.

The nine climatic variables were standardised across the occurrence dataset used to construct the climatic summaries. A principal components analysis was fitted to the standardised values. The first three axes explained 45.41%, 31.50% and 13.89% of total variance, respectively.

### Climatic niche metrics

Five species-level metrics described the breadth or heterogeneity of realised occupied climate. Temperature breadth was the mean difference between the 95th and 5th percentiles of BIO1, BIO5, BIO6 and BIO7. Moisture breadth was the corresponding mean percentile range for BIO12, BIO14, BIO15 and BIO17. Climatic heterogeneity was the mean within-species standard deviation across the nine standardised bioclimatic variables. PCA dispersion was the mean Euclidean distance between occupied cells and the species centroid in the first three principal-component dimensions. PCA hull area was the convex-hull area of occupied cells in the first two principal-component dimensions.

The comprehensive baseline comparison used a minimum threshold of 20 occupied climate cells and evaluated all five metrics symmetrically in the same 34 species. Broader specification analyses at minimum thresholds of 10, 20, 30 and 50 cells are retained as additional sensitivity analyses. These metrics describe realised occupied climate represented in sampled records. They do not estimate fundamental physiological tolerance and were not calculated separately for colour morphs.

### Spatial-organization models and collinearity diagnostics

For each metric we fitted a binomial generalized linear model with a logit link:

`among ~ metric_z + effort_z`

The response equalled one for geographically structured variation and zero for within-population polymorphism. `metric_z` was the standardized climatic metric within the model sample, and `effort_z` was standardized `log1p(n_climate_cells)`. Models were fitted in Python 3.12 using `statsmodels` 0.14.6 (Seabold & Perktold, 2010). Wald standard errors were estimated with family-clustered sandwich covariance; 95% confidence intervals were calculated on the log-odds scale and exponentiated for odds-ratio intervals.

The five climatic metrics were not entered simultaneously. Thus, collinearity among temperature, moisture, heterogeneity and PCA-derived metrics cannot destabilize coefficients within a fitted model. For each two-predictor model, we nevertheless quantified the correlation between `metric_z` and `effort_z`, variance inflation factors and the condition number of the standardized design matrix. Maximum VIFs ranged from 1.095 to 1.415 and condition numbers from 1.355 to 1.833, providing no indication of problematic predictor collinearity.

### Permutation, multiplicity and family-deletion analyses

For each of the five metrics in the frozen 34-species baseline, we used 9,999 label permutations with fixed seed 20260719. Spatial-category labels were shuffled among model-eligible species while metrics, occurrence effort and family labels remained fixed. The two-sided permutation p-value compared the absolute observed coefficient with the absolute permuted coefficients.

To assess concentration in individual plant families, we refitted the model after omitting each represented family in turn. We also report Holm-adjusted p-values across the five climatic metrics for both clustered Wald tests and permutation tests. The Holm correction is treated as multiplicity context because the five metrics are correlated descriptions of occupied climatic niche breadth rather than independent biological experiments; effect sizes, confidence intervals and cross-analysis directional consistency remain central to interpretation.

### Phylogenetic sensitivity models

For topology-based phylogenetic sensitivity, names were matched to Open Tree Taxonomy without approximate matching. Thirty of the 34 species were retained. We induced the Open Tree topology (Hinchliff et al., 2015), randomly resolved polytomies 100 times, assigned Grafen branch lengths and fitted `phyloglm` logistic MPLE models (Ho & Ané, 2014) for each of the five climatic metrics using the same formula `among ~ metric_z + effort_z`.

We also used V.PhyloMaker2 (Jin & Qian, 2019, 2022) and the time-scaled `GBOTB.extended.LCVP` backbone, which incorporates the broad seed-plant phylogeny of Smith and Brown (2018) and earlier time-calibrated plant-tree information (Zanne et al., 2014). All 34 species were retained. Using fixed seed 20260724, placement scenarios S1, S2 and S3 were evaluated separately for all five metrics. These models are phylogenetic sensitivity analyses rather than evidence that the residual process itself follows the assumed tree model.

### Candidate-versus-control and coarse spatial sensitivity analyses

We retained the previously implemented candidate-versus-control analysis to test whether documented colour-variable species generally occupied broader climates than taxonomically matched controls. Controls were selected from outside the complete candidate list and were not asserted to be monomorphic. Five metrics were evaluated at four cell thresholds for all controls and same-genus controls.

We also retained coarse range-fragmentation and environmental-turnover sensitivity analyses derived from GBIF point clouds. These components are unsupervised summaries of occurrence geometry, not verified populations, barriers, gene-flow units or colour-morph distributions.

## Results

### Evidence base and ascertainment

The resolved evidence base contained 664 candidate species, of which 111 were validated as natural cases of intraspecific flower-colour variation. Validation probability showed a nonlinear association with retained literature effort, demonstrating strong research-effort dependence in the assembled sample and precluding interpretation as a census of angiosperm prevalence.

### Candidate-versus-control climatic niches

Across the previously implemented conditional-logit specifications, there was no consistent evidence that documented colour-variable species occupied broader climatic niches than taxonomically matched controls. In the same-genus comparison at the 20-cell threshold, odds ratios ranged from 0.830 for temperature breadth to 1.184 for PCA hull area, and all five confidence intervals included one. Thus, the main comparison concerns variation among documented colour-variable species in the spatial scale at which colour variation is expressed, rather than a general climatic distinction between colour-variable species and other congeners.

### Comprehensive five-metric baseline comparison

The frozen baseline contained 34 species from 25 families: 20 within-population and 14 geographically structured cases. All five climatic-metric odds ratios were below one (Table 2). Temperature breadth had an odds ratio of 0.817 (95% CI 0.384–1.739; clustered Wald p = 0.6000; permutation p = 0.6029). Moisture breadth had the largest negative effect, with an odds ratio of 0.412 (0.180–0.947; p = 0.0368; permutation p = 0.0475). Climatic heterogeneity had an odds ratio of 0.681 (0.294–1.577; p = 0.3700; permutation p = 0.3574), PCA dispersion 0.712 (0.306–1.660; p = 0.4317; permutation p = 0.3820), and PCA hull area 0.577 (0.312–1.067; p = 0.0797; permutation p = 0.2382).

Every leave-one-family-out estimate remained below one for every metric. The leave-one-family-out odds-ratio ranges were 0.636–0.973 for temperature breadth, 0.306–0.465 for moisture breadth, 0.492–0.782 for climatic heterogeneity, 0.525–0.813 for PCA dispersion and 0.489–0.671 for PCA hull area. Thus, the shared negative direction was not generated by a single represented family.

After Holm correction across the five metrics, the moisture-breadth clustered Wald p-value was 0.184 and the permutation p-value was 0.238. No other metric retained conventional statistical support after correction. Accordingly, moisture breadth is best described as the strongest association within a broader, directionally consistent climatic-niche pattern rather than as a uniquely established climatic driver.

### Collinearity diagnostics

Predictor correlations between the focal climatic metric and occurrence effort ranged from 0.295 to 0.541. Maximum VIFs were 1.095 for PCA dispersion, 1.097 for moisture breadth, 1.121 for temperature breadth, 1.204 for climatic heterogeneity and 1.415 for PCA hull area. Condition numbers ranged from 1.355 to 1.833. These values provide no indication that the estimated climatic effects were generated by problematic collinearity with the effort covariate.

### Broader evidence and occurrence-sampling sensitivity

The broader evidence set produced weaker moisture-breadth estimates than the frozen baseline, indicating that effect magnitude depends on classification evidence quality. The previously implemented paginated, quality-filtered GBIF sensitivity analysis strengthened the focal moisture association relative to the earlier primary occurrence extraction. We retain this analysis as evidence that the moisture direction is not an obvious artifact of deterministic first-page occurrence sampling, while recognizing that the paginated dataset was itself capped and has not yet been rerun symmetrically for all five metrics.

### Open Tree phylogenetic sensitivity

Open Tree matching retained 30 species, and 100 fits completed for each climatic metric. Every point estimate was negative. Median odds ratios were 0.908 for temperature breadth, 0.573 for moisture breadth, 0.726 for climatic heterogeneity, 0.819 for PCA dispersion and 0.654 for PCA hull area. Median 95% confidence intervals were 0.436–1.893, 0.234–1.403, 0.322–1.636, 0.387–1.731 and 0.261–1.641, respectively. All intervals included one. The negative direction therefore persisted under topology-based phylogenetic correction, while precision declined.

### Time-scaled megaphylogeny sensitivity

V.PhyloMaker2 retained all 34 species under placement scenarios S1–S3. The five climatic effects remained negative in every scenario and were nearly invariant to placement assumptions. Across scenarios, temperature-breadth odds ratios were 0.827–0.838, moisture-breadth odds ratios were 0.448–0.454, climatic-heterogeneity odds ratios were 0.670–0.679, PCA-dispersion odds ratios were 0.698–0.707 and PCA-hull-area odds ratios were 0.599–0.612. All phylogenetic confidence intervals included one. For moisture breadth, scenario-specific p-values were 0.112–0.115 and Holm-adjusted p-values were 0.560–0.574.

Thus, time scaling and retention of all species did not reverse the climatic-niche pattern. Phylogenetic correction attenuated statistical certainty more strongly than it altered effect direction.

### Coarse occurrence-cloud alternatives

Previously implemented fragmentation, connectivity and environmental-turnover summaries did not clearly account for the focal relationship. These analyses cannot exclude environmental sorting among colour morphs because GBIF occurrence records were not labelled by flower-colour state and the derived components are not verified biological populations.

## Discussion

The principal result is a cross-metric pattern rather than a single isolated p-value. In the frozen 34-species comparison, geographically structured flower-colour variation was associated with narrower sampled occupied climatic niches across all five environmental summaries. Temperature, moisture, total climatic heterogeneity and both PCA-based metrics all had odds ratios below one, and every leave-one-family-out estimate retained that direction. The strongest contrast occurred for moisture breadth, where the unadjusted family-clustered confidence interval excluded one and the permutation p-value was below 0.05. However, this support did not survive correction across the five exploratory metrics, and phylogenetic confidence intervals included one. The data therefore support a directionally coherent comparative signal with moisture as its strongest component, not a confirmed moisture-specific mechanism.

This distinction matters biologically. A species in which colour variants coexist locally must maintain multiple phenotypes within at least partly shared environmental and demographic settings. A species in which colour variation is geographically structured can instead exhibit morph differentiation across separated populations or regions. The consistent negative association with occupied climatic breadth is compatible with the idea that these two spatial configurations occur in different macroecological contexts. For example, geographically structured variation could arise more often in species whose sampled climatic occupancy is restricted, whereas local coexistence may be maintained across a broader range of realised conditions. Conversely, the association could arise through dispersal limitation, demographic history, range geometry, sampling processes or correlated life-history traits rather than through direct climatic selection on flower colour.

The contrast among climatic metrics provides additional ecological resolution. Temperature breadth was weakly associated with spatial organization, whereas moisture breadth produced the largest effect size and PCA hull area the second largest. This ordering makes a purely generic “small niche versus large niche” interpretation less complete than one in which hydrological dimensions contribute disproportionately. Nevertheless, the other metrics were also negative and the five metrics are correlated summaries of occupied climatic space. We therefore do not interpret moisture as uniquely causal. Instead, the results identify moisture occupancy as the most prominent axis within a broader climatic-niche pattern that warrants morph-resolved and population-level testing.

The phylogenetic analyses sharpen rather than overturn this interpretation. Both the Open Tree + Grafen analysis and the time-scaled V.PhyloMaker2 analysis retained negative effects for all five metrics. In the time-scaled analysis, the moisture-breadth odds ratio remained close to the non-phylogenetic estimate (approximately 0.45 versus 0.41), while its confidence interval broadened to include one. This is consistent with limited power in a 34-species comparative sample and with shared evolutionary structure contributing uncertainty. The phylogenetic results therefore argue against a claim that the pattern is simply produced by one or a few related lineages, but they do not establish phylogenetically independent statistical significance.

The collinearity diagnostics also reduce a separate modeling concern. Each climatic metric was analysed in a separate model with occurrence effort, rather than entering five correlated niche metrics simultaneously. VIFs were at most 1.415 and condition numbers below 1.84. Thus, coefficient instability caused by multicollinearity with the effort covariate is not a plausible explanation for the observed effect ordering. The price of the separate-model strategy is multiplicity, which is why we report Holm-adjusted values and interpret the five models jointly.

Classification uncertainty remains an important limitation. The spatial categories are source-traceable and rule-derived, not fully independently adjudicated biological annotations. Earlier sensitivity analyses showed that broadening the evidence set attenuated the moisture association. This is scientifically informative: literature-derived spatial traits are themselves measurements with error, and extraction rules, source context and review status affect downstream macroecological inference. The automated evidence-extraction pipeline is therefore a methodological contribution, but its ecological use depends on preserving source provenance and uncertainty rather than treating machine-derived labels as ground truth.

The candidate-versus-control analysis further delimits the inference. Documented colour-variable species did not consistently occupy broader climates than matched controls. The climatic signal is therefore not simply that flower-colour-variable species have unusual niche breadth as a class; it concerns the *spatial organization* of variation among documented colour-variable species. Controls cannot be assumed monomorphic, because exclusion from the candidate list indicates only that they were not identified by the literature pipeline.

Occurrence data impose another boundary. The climatic metrics describe species-level realised occupancy sampled through GBIF, not physiological tolerance, morph-specific niches or local selective environments. The records are not labelled by flower-colour morph, and coarse point-cloud components are not verified populations. Consequently, the present analysis cannot demonstrate that one colour morph occupies drier, wetter, warmer or otherwise distinct local conditions. A direct mechanistic test would require georeferenced morph frequencies, named populations, environmental measurements linked to those states, and ideally demographic or genomic evidence of selection and gene flow.

Despite these constraints, the study establishes a reproducible framework for connecting automated evidence synthesis with comparative ecological analysis. The literature pipeline recovers and classifies documented flower-colour variation at explicit spatial scales; the downstream analysis then asks whether those scales are associated with species-level climatic occupancy. In the current frozen sample, geographically structured variation consistently occurs toward the narrower end of the sampled climatic-niche spectrum, with moisture breadth showing the strongest contrast. The result is robust in direction to family deletion and phylogenetic correction but remains exploratory in inferential strength. That combination—reproducible extraction, transparent uncertainty and an ecologically interpretable comparative pattern—is the central contribution.

## References

Dalrymple, R. L., Kemp, D. J., Flores-Moreno, H., Laffan, S. W., White, T. E., Hemmings, F. A., & Moles, A. T. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist, 228*(6), 1972–1985. https://doi.org/10.1111/nph.16737

Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology, 37*, 4302–4315. https://doi.org/10.1002/joc.5086

Hinchliff, C. E., Smith, S. A., Allman, J. F., Burleigh, J. G., Chaudhary, R., Coghill, L. M., Crandall, K. A., Deng, J., Drew, B. T., Gazis, R., Gude, K., Hibbett, D. S., Katz, L. A., Laughinghouse, H. D., IV, McTavish, E. J., Midford, P. E., Owen, C. L., Ree, R. H., Rees, J. A., Soltis, D. E., Williams, T., & Cranston, K. A. (2015). Synthesis of phylogeny and taxonomy into a comprehensive tree of life. *Proceedings of the National Academy of Sciences of the United States of America, 112*(41), 12764–12769. https://doi.org/10.1073/pnas.1423041112

Ho, L. S. T., & Ané, C. (2014). A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. *Systematic Biology, 63*(3), 397–408. https://doi.org/10.1093/sysbio/syu005

Jin, Y., & Qian, H. (2019). V.PhyloMaker: An R package that can generate very large phylogenies for vascular plants. *Ecography, 42*(8), 1353–1359. https://doi.org/10.1111/ecog.04434

Jin, Y., & Qian, H. (2022). V.PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. *Plant Diversity, 44*(4), 335–339. https://doi.org/10.1016/j.pld.2022.05.005

Koski, M. H., & Ashman, T.-L. (2016). Macroevolutionary patterns of ultraviolet floral pigmentation explained by geography and associated bioclimatic factors. *New Phytologist, 211*(2), 708–718. https://doi.org/10.1111/nph.13921

Narbona, E., Wang, H., Ortiz, P. L., Arista, M., & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: Occurrence, maintenance and implications for speciation. *Plant Biology, 20*(S1), 8–20. https://doi.org/10.1111/plb.12575

Priem, J., Piwowar, H., & Orr, R. (2022). OpenAlex: A fully-open index of scholarly works, authors, venues, institutions, and concepts. *arXiv*. https://doi.org/10.48550/arXiv.2205.01833

Rausher, M. D. (2008). Evolutionary transitions in floral color. *International Journal of Plant Sciences, 169*(1), 7–21. https://doi.org/10.1086/523358

Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and statistical modeling with Python. In *Proceedings of the 9th Python in Science Conference* (pp. 92–96).

Smith, S. A., & Brown, J. W. (2018). Constructing a broadly inclusive seed plant phylogeny. *American Journal of Botany, 105*(3), 302–314. https://doi.org/10.1002/ajb2.1019

Trunschke, J., Lunau, K., Pyke, G. H., Ren, Z.-X., & Wang, H. (2021). Flower color evolution and the evidence of pollinator-mediated selection. *Frontiers in Plant Science, 12*, 617851. https://doi.org/10.3389/fpls.2021.617851

Wessinger, C. A., & Rausher, M. D. (2012). Lessons from flower colour evolution on targets of selection. *Journal of Experimental Botany, 63*(16), 5741–5749. https://doi.org/10.1093/jxb/ers267

Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., FitzJohn, R. G., McGlinn, D. J., O'Meara, B. C., Moles, A. T., Reich, P. B., Royer, D. L., Soltis, D. E., Stevens, P. F., Westoby, M., Wright, I. J., Aarssen, L., Bertin, R. I., Calaminus, A., Govaerts, R., . . . Beaulieu, J. M. (2014). Three keys to the radiation of angiosperms into freezing environments. *Nature, 506*(7486), 89–92. https://doi.org/10.1038/nature12872

## Data Accessibility Statement

Analysis code, source-level evidence fields, frozen classification manifests, correction logs and model outputs are maintained in the public GitHub repository `zuizui0223/fcp`. The comprehensive five-metric 34-species analysis is documented in `docs/supporting/jbi_environmental_niche_comprehensive_34species.md` and is generated by workflow run `31142541223` from commit `b82e7cd71db77b6f184aaa0bd0847d13c14858bd`. The corresponding artifact is `34species-environmental-niche-comprehensive` (artifact ID `8980386463`).

**Not verified:** before submission, archive the exact code and data release at a permanent DOI; obtain the required citable GBIF occurrence identifiers; and ensure that all manuscript tables and figures are generated from the final frozen analysis release.

## Tables

### Table 1. Operational classification of the documented spatial organization of intraspecific flower-colour variation

| Category | Operational evidence requirement | Included in binary models |
|---|---|---|
| Within-population flower-colour polymorphism | Explicit source evidence that at least two discrete natural flower-colour variants coexist within at least one population | Yes; response = 0 |
| Geographically structured flower-colour variation | Explicit source evidence of differentiation among populations or regions, without retained evidence of local coexistence | Yes; response = 1 |
| Mixed | Retained evidence supports both within-population coexistence and geographic differentiation | No |
| Unclear | Evidence does not resolve the spatial configuration | No |

### Table 2. Comprehensive frozen-baseline environmental-niche models at the 20-cell threshold

| Metric | Species | Families | Within / among | Odds ratio | Family-clustered 95% CI | Wald p | Permutation p | Holm Wald p | Holm permutation p | Leave-one-family-out OR range |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Temperature breadth | 34 | 25 | 20 / 14 | 0.817 | 0.384–1.739 | 0.6000 | 0.6029 | 1.0000 | 1.0000 | 0.636–0.973 |
| Moisture breadth | 34 | 25 | 20 / 14 | 0.412 | 0.180–0.947 | 0.0368 | 0.0475 | 0.1840 | 0.2375 | 0.306–0.465 |
| Climatic heterogeneity | 34 | 25 | 20 / 14 | 0.681 | 0.294–1.577 | 0.3700 | 0.3574 | 1.0000 | 1.0000 | 0.492–0.782 |
| PCA dispersion | 34 | 25 | 20 / 14 | 0.712 | 0.306–1.660 | 0.4317 | 0.3820 | 1.0000 | 1.0000 | 0.525–0.813 |
| PCA hull area | 34 | 25 | 20 / 14 | 0.577 | 0.312–1.067 | 0.0797 | 0.2382 | 0.3189 | 0.9528 | 0.489–0.671 |

*Note.* Odds ratios are for geographically structured rather than within-population organization per one standard-deviation increase in the focal climatic metric, controlling for standardized `log1p(n_climate_cells)`. The five climatic metrics were fitted separately. Holm adjustments are across the five metrics within this frozen-baseline family.

### Table 3. Collinearity diagnostics for the comprehensive baseline models

| Metric | Correlation with effort | Max VIF | Condition number |
|---|---:|---:|---:|
| Temperature breadth | 0.329 | 1.121 | 1.407 |
| Moisture breadth | 0.298 | 1.097 | 1.360 |
| Climatic heterogeneity | 0.412 | 1.204 | 1.549 |
| PCA dispersion | 0.295 | 1.095 | 1.355 |
| PCA hull area | 0.541 | 1.415 | 1.833 |

### Table 4. Phylogenetic sensitivity of all five environmental-niche metrics

| Metric | Open Tree n | Open Tree median OR (95% CI) | V.PhyloMaker2 S1 OR (95% CI) | S2 OR (95% CI) | S3 OR (95% CI) |
|---|---:|---:|---:|---:|---:|
| Temperature breadth | 30 | 0.908 (0.436–1.893) | 0.827 (0.385–1.777) | 0.838 (0.393–1.789) | 0.827 (0.385–1.777) |
| Moisture breadth | 30 | 0.573 (0.234–1.403) | 0.448 (0.167–1.206) | 0.454 (0.171–1.211) | 0.448 (0.167–1.206) |
| Climatic heterogeneity | 30 | 0.726 (0.322–1.636) | 0.670 (0.290–1.550) | 0.679 (0.295–1.565) | 0.670 (0.290–1.550) |
| PCA dispersion | 30 | 0.819 (0.387–1.731) | 0.698 (0.319–1.527) | 0.707 (0.324–1.542) | 0.698 (0.319–1.527) |
| PCA hull area | 30 | 0.654 (0.261–1.641) | 0.599 (0.232–1.545) | 0.612 (0.242–1.549) | 0.599 (0.232–1.545) |

*Note.* Open Tree values summarize 100 completed fits after polytomy resolution and Grafen branch-length assignment. V.PhyloMaker2 values retain all 34 species under three placement scenarios. All phylogenetic confidence intervals include one.

## Figure plan

The previous moisture-only figures should not be treated as final figures for this revision. The preferred replacement is a forest plot displaying the five non-phylogenetic effect estimates together with the Open Tree and time-scaled phylogenetic sensitivities, plus a compact panel or supplementary plot showing leave-one-family-out ranges. This avoids visually privileging moisture while preserving its larger effect size.

## Supporting Information

The exact comprehensive results, diagnostics, provenance and interpretation guard are recorded in `docs/supporting/jbi_environmental_niche_comprehensive_34species.md`. Existing supporting tables for the evidence pipeline, candidate-versus-control analyses, occurrence sampling, classification audit and earlier moisture-focused sensitivities remain useful provenance, but manuscript submission tables should be regenerated against the final frozen branch so that no stale 0.426-era values are mixed with the current five-metric analysis.
