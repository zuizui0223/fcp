# Moisture niche breadth is associated with the documented spatial organization of intraspecific flower-colour variation

## Abstract

### Aim

Intraspecific phenotypic variation may be maintained through local coexistence or partitioned geographically among populations. These spatial expressions are often collapsed into a single category, although they imply different ecological and evolutionary contexts. We asked whether documented within-population flower-colour polymorphism and geographically structured flower-colour variation differ in species-level occupied climatic niche breadth.

### Location

Global literature-derived sample.

### Taxon

Angiosperms with documented natural intraspecific flower-colour variation.

### Methods

We assembled a literature-derived evidence base and classified documented flower-colour variation as within-population or among-population. Classifications supported directly by baseline evidence formed a frozen baseline-unambiguous set; less direct literature enrichment formed a broader sensitivity set. We combined these classifications with cleaned GBIF occurrences and WorldClim-derived occupied climatic niche metrics. Logistic models predicted among-population rather than within-population organization from standardized niche metrics while controlling for the number of occupied climate cells. Uncertainty was estimated using family-clustered sandwich standard errors, 9,999 label permutations, and leave-one-family-out analyses. Multiple climatic metrics and occurrence thresholds were retained in the multiplicity record.

### Results

The baseline-unambiguous moisture model included 34 species from 25 families, comprising 19 within-population and 15 among-population cases. Geographically structured variation was associated with narrower species-level realised moisture niche breadth (odds ratio 0.403, family-clustered 95% confidence interval 0.165–0.985; clustered Wald p = 0.0463; two-sided permutation p = 0.0439). The direction was retained after omission of each represented family, with leave-one-family-out odds ratios ranging from 0.284 to 0.456. The estimate attenuated in the broader evidence set (odds ratio 0.566, 95% confidence interval 0.293–1.093; permutation p = 0.0886). Coarse GBIF occurrence-cloud fragmentation and generic climatic turnover among unlabelled occurrence components did not account for the focal association.

### Main conclusions

The documented spatial organization of intraspecific flower-colour variation carries an evidence-sensitive comparative climatic signal. The result is not a test of morph-specific physiological tolerance, local adaptation, or causal climatic effects. Mechanistic inference requires colour-state-labelled localities or named populations linked directly to environmental conditions.

**Keywords:** climatic niche breadth; flower-colour variation; geographic differentiation; intraspecific polymorphism; evidence synthesis; GBIF; macroecology

## Introduction

Intraspecific phenotypic variation can be expressed at more than one spatial scale. Distinct phenotypes may coexist within the same population, be partitioned among populations or regions, or occur in a mixture of local and geographic configurations. These alternatives are not interchangeable. Local coexistence requires that multiple forms persist under at least partly shared demographic and environmental conditions, whereas geographic differentiation can arise under spatially separated environments, restricted dispersal, demographic history, or other forms of regional structure. Treating both situations as a single binary state risks obscuring the spatial scale at which variation is actually maintained or documented.

Flower colour provides an unusually rich system for studying this distinction. Natural flower-colour variation has been associated with pollinator behaviour, abiotic stress, pleiotropic effects, drift, mating systems, and demographic processes [CITATIONS REQUIRED]. Most comparative syntheses, however, focus on whether colour variation is present, which mechanisms maintain polymorphism, or how average floral traits vary with climate. Individual case studies also document geographic colour clines and regional differentiation, but these studies do not establish whether the spatial organization of flower-colour variation shows a broader comparative association with species-level occupied climate.

A central conceptual problem is terminology. We use **intraspecific flower-colour variation** as the umbrella term. **Within-population flower-colour polymorphism** refers to discrete colour variants documented as coexisting within at least one natural population. **Geographically structured flower-colour variation** refers to colour variants documented among populations or regions without demonstrated local coexistence. The latter category should not automatically be called flower-colour polymorphism because coexistence at the population scale has not been established.

Here, we test whether these two documented spatial expressions differ in occupied climatic niche across a global literature-derived sample. The broad comparison between local coexistence and geographic differentiation was theory-led. Several occupied-climate metrics were evaluated, and moisture niche breadth is therefore treated as the focal reported association rather than as a prospectively preregistered endpoint. We explicitly retain the full metric-by-threshold result matrix to avoid presenting the focal association without its multiplicity context.

We also address four alternative explanations and inferential risks. First, the association may depend on classification evidence quality. Second, it may be concentrated in one or a few taxonomic families. Third, occurrence sampling effort may influence niche estimates. Fourth, coarse fragmentation or climatic turnover in unlabelled GBIF occurrence clouds may correlate with the documented spatial category. Our aim is not to infer causation or local adaptation, but to test whether the documented spatial organization of flower-colour variation covaries with species-level realised climatic niche breadth while making evidence quality and ascertainment explicit.

## Methods

### Literature-derived study population

The study population consisted of species identified through the repository's literature acquisition and validation pipeline as documented natural cases of intraspecific flower-colour variation. The candidate pool is not a random sample of angiosperms. Consequently, the inferential population is the assembled set of documented cases rather than all flowering plants, and the analysis does not estimate worldwide prevalence.

The latest resolved evidence base contained 664 candidate species, of which 111 were retained as validated natural cases and 553 remained deferred. Literature effort strongly predicted validation probability in a separate ascertainment analysis. This relationship is reported to establish that the evidence base is shaped by research effort, not as a biological prevalence model.

### Spatial-organization classification

Each retained case was classified using source text as `within_population`, `among_population`, `mixed`, or `unclear`. Within-population cases required explicit evidence that multiple discrete colour variants coexisted within at least one natural population. Among-population cases required evidence of colour differentiation among populations or regions without demonstrated local coexistence. Mixed and unclear cases were excluded from binary spatial-scale models.

The baseline-unambiguous set retained cases whose classification was directly supported by the baseline source record. The pipeline preserves species name, family, assigned spatial category, classification source, source title, DOI or OpenAlex identifier, evidence snippet, and review note. The strict list is sorted deterministically and written with a SHA-256 digest. Any later change must be recorded in a correction log. Classifications were frozen independently of the climatic model output.

A broader evidence set incorporated high-confidence literature enrichment for cases unresolved by baseline evidence. This set was analysed as a generalization sensitivity analysis rather than as a replacement for the baseline-unambiguous set.

### Occurrence data and climatic niche metrics

Species occurrences were obtained from GBIF and cleaned through the repository workflow. Records without valid coordinates were removed, duplicate and unsuitable records were filtered according to the pipeline rules, and occurrences were summarized by occupied climate cells. WorldClim 2.1 bioclimatic layers supplied climatic values [RESOLUTION AND VARIABLES TO BE INSERTED FROM PIPELINE DOCUMENTATION]. Species-level metrics included PCA dispersion, climatic heterogeneity, PCA hull area, temperature breadth, and moisture breadth.

These metrics describe realised occupied climatic space represented by cleaned occurrence records. They do not estimate fundamental physiological tolerance and are not assigned to individual colour morphs.

### Comparative models

For each climatic metric and minimum occupied-cell threshold of 10, 20, 30, or 50 cells, we fitted a binomial generalized linear model with logit link:

`among ~ metric_z + effort_z`

where `among` equalled one for geographically structured variation and zero for within-population polymorphism, `metric_z` was the standardized climatic niche metric, and `effort_z` was the standardized log-transformed number of occupied climate cells. Models were fitted in `statsmodels` using a binomial GLM. Sandwich covariance was clustered by plant family. We report the log-odds coefficient, clustered standard error, odds ratio, Wald 95% confidence interval, clustered Wald p-value, species and family counts, class counts, convergence status, iteration count, and fitted-probability range.

For the focal moisture-breadth comparison, we additionally used 9,999 permutations of spatial-category labels to obtain a two-sided permutation p-value. Leave-one-family-out models quantified whether the direction depended on any represented family. These analyses reduce sensitivity to individual families but do not constitute a phylogenetic comparative analysis.

### Coarse spatial negative controls

We evaluated GBIF point-cloud fragmentation metrics and generic climatic turnover among unsupervised occurrence components. These tests address whether coarse properties of the sampled occurrence cloud account for the focal association. The components are not interpreted as biological populations, barriers, gene-flow units, or morph distributions.

## Results

### Evidence base and ascertainment

The resolved acquisition artifact contained 664 candidates and 111 validated natural cases. Validation probability increased strongly with retained literature effort and showed nonlinear saturation. This confirms that the documented evidence base is research-effort dependent and should not be interpreted as a representative census of angiosperms.

### Candidate-versus-control climatic niche comparison

The matched comparison included 70 focal species and 243 taxonomically matched control species. Across five occupied-climate metrics and multiple occurrence thresholds, there was no clear evidence that documented colour-variable species generally occupied broader climatic niches than matched controls. At the same-genus, 20-cell specification, all confidence intervals included an odds ratio of one.

### Spatial organization and moisture niche breadth

The validated baseline-unambiguous moisture model contained 34 species from 25 families: 19 within-population and 15 among-population cases. The coefficient for standardized moisture breadth was negative, corresponding to an odds ratio of 0.403 for geographically structured rather than within-population variation. The family-clustered 95% confidence interval was 0.165–0.985, and the clustered Wald p-value was 0.0463. The model converged in four iterations.

Permutation inference gave a two-sided p-value of 0.0439 from 9,999 valid permutations. After omitting each represented family in turn, odds ratios ranged from 0.284 to 0.456 and remained below one in every refit. Thus, the negative direction was not produced by any single represented family.

The broader evidence set produced a weaker estimate. Its moisture-breadth odds ratio was 0.566, with a 95% confidence interval of 0.293–1.093 and a permutation p-value of 0.0886. The attenuation indicates that the inferred association depends on classification evidence quality.

### Multiplicity context

All five climatic metrics were evaluated at four minimum occupied-cell thresholds, yielding a complete 20-specification matrix. The moisture-breadth result is presented as the focal reported association because it produced the clearest evidence-sensitive pattern, not because it was the only climatic metric evaluated. The complete matrix is retained in the supplementary outputs.

### Coarse occurrence-cloud alternatives

No examined point-cloud fragmentation metric showed a clear independent association with documented spatial organization after accounting for moisture breadth and sampling effort. Adding fragmentation metrics did not improve model fit relative to the baseline moisture model, and the moisture coefficient remained negative.

Generic climatic turnover among unsupervised occurrence components was also not associated with documented spatial organization. At the primary 100-km, minimum-three-record specification, the turnover odds ratio was 1.007, with a 95% confidence interval of 0.367–2.760 and p = 0.989. These results do not show that fragmentation or environmental differentiation is absent. They show only that the examined coarse, unlabelled occurrence-cloud summaries did not account for the focal association.

## Discussion

Our analysis identifies an evidence-sensitive association between species-level realised moisture niche breadth and the documented spatial organization of intraspecific flower-colour variation. In the baseline-unambiguous set, geographically structured cases were associated with narrower moisture niche breadth than cases with documented within-population coexistence. The direction persisted after removing every represented family, but the estimate weakened when less certain classifications were included.

Several biological interpretations are compatible with this pattern. Geographic differentiation may be more readily documented in species occupying comparatively restricted moisture contexts, or local coexistence may persist across a wider range of realised moisture conditions. Dispersal, demographic structure, range geometry, or research practices may also covary with occupied niche breadth. The present data cannot distinguish among these explanations. In particular, species-level GBIF occurrences lack colour-state labels, so the analysis does not test whether individual morphs occupy different climates or whether local adaptation maintains geographic colour differences.

The evidence-set contrast is itself a major result. Comparative trait syntheses often treat categorical classifications as fixed observations, but literature-derived categories differ in directness and source quality. Here, the focal odds ratio shifted from 0.403 in the baseline-unambiguous set to 0.566 in the broader evidence set. This attenuation demonstrates that conclusions about spatial organization are sensitive to how confidently coexistence or geographic differentiation is documented. Publishing the frozen strict-set audit and the broader sensitivity analysis makes this uncertainty visible rather than absorbing it into a single pooled category.

The family sensitivity results reduce concern that one well-studied lineage generated the association. Nevertheless, family-clustered uncertainty and leave-one-family-out refits are not substitutes for a species-level phylogeny. The result should therefore be interpreted as a cross-species association robust to represented-family deletion, not as a fully phylogenetically independent evolutionary relationship.

The coarse spatial analyses also require restraint. GBIF point-cloud components summarize sampled records rather than biological populations. Their failure to account for the moisture association does not reject morph-specific environmental sorting or geographic differentiation. A decisive mechanistic test would require documented colour states or named populations linked to localities, followed by direct comparison of environments among morph-labelled sites.

The study has additional limitations. The literature-derived sample is non-random and strongly influenced by publication effort. The strict analysis contains only 34 species. Multiple climatic metrics and thresholds were evaluated, and moisture breadth was selected as the focal reported association after comparative inspection. Species-level realised niches are imperfect proxies for the environments experienced by colour morphs. Local coexistence may also be under-documented relative to conspicuous geographic differences. These limitations restrict generality but do not erase the comparative signal observed in the strict evidence set.

The broader contribution is conceptual and methodological. The spatial organization of intraspecific variation should be treated as a comparative property in its own right rather than collapsed into a polymorphic-versus-monomorphic label. Distinguishing local coexistence from geographic differentiation, while preserving source-level evidence and classification uncertainty, provides a more defensible basis for testing how phenotypic variation is distributed across species and environments.

## Data accessibility and reproducibility

All analysis scripts, model specifications, evidence manifests, correction logs, complete metric-by-threshold results, and robustness outputs are maintained in the project repository. The final submission should provide a permanent archived release and DOI.

## Items still required before submission

- insert literature citations throughout the Introduction and Discussion;
- insert exact WorldClim variables and spatial resolution;
- verify author list, affiliations, acknowledgements, and funding;
- add figure captions and supplementary table references;
- archive a release and replace repository-only references with a permanent DOI;
- perform final language and journal-format editing.
