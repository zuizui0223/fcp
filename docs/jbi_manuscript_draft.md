# Evidence-sensitive association between moisture niche breadth and the documented spatial organization of intraspecific flower-colour variation

## Abstract

### Aim

Intraspecific phenotypic variation may be maintained through local coexistence or partitioned geographically among populations. These spatial expressions are often collapsed into one category despite representing different ecological contexts. We asked whether documented within-population flower-colour polymorphism and geographically structured flower-colour variation differ in species-level occupied climatic niche breadth.

### Location

Global literature-derived sample.

### Taxon

Angiosperms with documented natural intraspecific flower-colour variation.

### Methods

We classified documented flower-colour variation as within-population, among-population, mixed or unclear using retained source text. A frozen baseline-unambiguous set was audited with source titles, DOI or OpenAlex identifiers and evidence snippets. We combined binary classifications with cleaned GBIF occurrences and WorldClim-derived occupied-climate metrics. Binomial generalized linear models predicted among-population organization from standardized niche metrics while controlling for occupied-cell count. We used family-clustered sandwich uncertainty, 9,999 label permutations and leave-one-family-out refits. Five climatic metrics and four occurrence thresholds were retained in the multiplicity record.

### Results

The audited baseline-unambiguous moisture model included 34 species from 25 families: 20 within-population and 14 among-population cases. Geographically structured variation showed a negative association with realised moisture niche breadth (odds ratio 0.426, family-clustered 95% confidence interval 0.184–0.985; clustered Wald p = 0.0460). Permutation support was borderline (two-sided p = 0.0556). The direction remained negative after omission of each represented family, with leave-one-family-out odds ratios of 0.317–0.481. The estimate attenuated in the broader evidence set (odds ratio 0.563, 95% confidence interval 0.292–1.085; permutation p = 0.0944).

### Main conclusions

The documented spatial organization of flower-colour variation carries an evidence-sensitive climatic signal, but support depends on uncertainty method and classification quality. The result is exploratory and does not test morph-specific tolerance, local adaptation or causal climatic effects.

**Keywords:** climatic niche breadth; flower-colour variation; geographic differentiation; intraspecific polymorphism; evidence synthesis; GBIF; macroecology

## Introduction

Intraspecific phenotypic variation can coexist locally, be partitioned among populations or occur in mixed configurations. These alternatives are not interchangeable. Local coexistence requires multiple forms to persist under partly shared conditions, whereas geographic differentiation can arise under spatially separated environments, dispersal limitation, demographic history or regional structure.

Flower colour is a well-documented system in which variation has been linked to pollinators, abiotic stress, pleiotropy, drift, mating systems and demography [CITATIONS REQUIRED]. Comparative syntheses usually ask whether colour variation is present or how mean floral traits correlate with climate. They rarely treat the spatial organization of variation as a comparative property.

We use **intraspecific flower-colour variation** as the umbrella term. **Within-population flower-colour polymorphism** requires documented coexistence of discrete variants within at least one natural population. **Geographically structured flower-colour variation** refers to variants differentiated among populations or regions without demonstrated local coexistence. Cases with evidence for both are classified as mixed.

We tested whether these documented spatial expressions differ in occupied climatic niche across a literature-derived sample. Several climatic metrics were evaluated, so moisture niche breadth is treated as the focal reported association rather than a prospectively preregistered endpoint. We also assessed classification sensitivity, family concentration, occurrence effort and coarse occurrence-cloud alternatives.

## Methods

### Literature-derived study population

The evidence base consisted of species identified by the repository pipeline as documented natural cases of intraspecific flower-colour variation. The candidate pool is not a random sample of angiosperms, and inference is limited to the assembled documented cases.

### Spatial-organization classification

Each case was classified as `within_population`, `among_population`, `mixed` or `unclear`. Within-population classification required explicit local coexistence. Among-population classification required geographic differentiation without retained evidence of coexistence. Text indicating both within- and among-population structure produced a mixed classification and exclusion from binary models.

The classifier was audited after source identifiers and evidence snippets were propagated into the analysis dataset. This audit revealed that plural `within populations` had not been recognized. The rule was corrected before the final rerun. The frozen manifest preserves species, family, category, source title, identifier, evidence text and decision note, together with a SHA-256 digest and correction log.

### Occurrence data and climatic niche metrics

Species occurrences were obtained from GBIF and cleaned through the repository workflow. Records were summarized by occupied climate cells and linked to WorldClim 2.1 bioclimatic data [RESOLUTION AND VARIABLES TO BE INSERTED]. Metrics included PCA dispersion, climatic heterogeneity, PCA hull area, temperature breadth and moisture breadth. These represent realised occupied climate, not physiological tolerance or morph-specific climate.

### Comparative models

For each metric and minimum occupied-cell threshold of 10, 20, 30 or 50 cells, we fitted:

`among ~ metric_z + effort_z`

where `among` equalled one for geographically structured variation, `metric_z` was the standardized niche metric and `effort_z` was standardized `log1p(n_climate_cells)`. Models used `statsmodels` binomial GLMs with logit link and family-clustered sandwich covariance.

For focal moisture comparisons, we additionally used 9,999 label permutations and leave-one-family-out refits. These analyses reduce dependence on individual families but are not phylogenetic comparative models.

### Coarse spatial negative controls

GBIF point-cloud fragmentation and generic climatic turnover among unsupervised occurrence components were examined as coarse alternative explanations. Components were not interpreted as populations, barriers or morph distributions.

## Results

### Evidence base and ascertainment

The resolved dataset contained 664 candidates and 111 validated natural cases. Validation probability increased strongly with retained literature effort, confirming that the evidence base is research-effort dependent.

### Candidate-versus-control comparison

The matched comparison included 70 focal species and 243 taxonomically matched controls. Across five metrics and four thresholds, there was no clear evidence that documented colour-variable species generally occupied broader climatic niches.

### Audited spatial organization and moisture breadth

The final baseline-unambiguous model contained 34 species from 25 families: 20 within-population and 14 among-population cases. The standardized moisture-breadth coefficient was -0.854, corresponding to an odds ratio of 0.426. The family-clustered 95% confidence interval was 0.184–0.985 and the clustered Wald p-value was 0.0460. The model converged in four iterations.

The two-sided permutation p-value was 0.0556 from 9,999 valid permutations. Leave-one-family-out odds ratios ranged from 0.317 to 0.481 and remained below one in every refit. Thus, the negative direction was not attributable to one represented family, although permutation support was borderline.

The broader evidence set contained 51 species from 29 families and produced a weaker estimate: odds ratio 0.563, 95% confidence interval 0.292–1.085 and permutation p = 0.0944. This attenuation demonstrates sensitivity to classification evidence.

### Multiplicity context

All five metrics were evaluated at four thresholds, yielding 20 specifications. Moisture breadth is presented as the focal association because it produced the clearest evidence-sensitive pattern, not because it was the sole metric examined.

### Coarse occurrence-cloud alternatives

No examined fragmentation or generic turnover metric clearly accounted for the moisture association. These analyses do not demonstrate absence of morph-specific environmental sorting because occurrence records lack colour-state labels.

## Discussion

The audited analysis identifies a negative, evidence-sensitive association between realised moisture niche breadth and the documented spatial organization of flower-colour variation. The clustered interval narrowly excluded one, family-deletion direction was stable and permutation support was borderline. This combination warrants interpretation as suggestive comparative evidence rather than binary confirmation.

Several explanations remain compatible with the pattern. Geographic differentiation may be more frequently documented in species occupying restricted moisture contexts, local coexistence may persist across broader realised conditions, or dispersal, demography, range geometry and research practice may covary with both variables. The present data cannot distinguish these mechanisms.

The classification audit materially changed the analysis. Recognizing plural within-population language moved a species with evidence at both scales out of the binary strict set. This illustrates why literature-derived categories should not be treated as error-free observations and why source-level manifests and correction logs are central analytical outputs.

The broader-set attenuation reinforces this point. The strict odds ratio was 0.426, whereas the broader estimate was 0.563 and its interval included one. Classification quality is therefore part of the substantive inference, not merely a data-cleaning concern.

Family-clustered uncertainty and family-deletion refits reduce sensitivity to single lineages but do not establish phylogenetic independence. Similarly, GBIF components are sampled point-cloud summaries rather than populations. Mechanistic tests require morph-labelled localities or named populations linked to environmental data.

Limitations include the non-random literature sample, publication-effort dependence, modest strict-set size, multiple evaluated metrics, post-analysis focal selection, incomplete phylogenetic control and lack of morph-labelled occurrences. These limitations restrict generality but do not remove the value of treating spatial organization as a comparative trait.

## Data accessibility and reproducibility

Analysis scripts, source-level manifests, correction logs, model tables, complete metric-by-threshold outputs and robustness analyses are maintained in the repository. The final submission should archive a permanent release with a DOI.

## Items still required before submission

- insert literature citations;
- insert exact WorldClim variables and resolution;
- verify authorship, affiliations, acknowledgements and funding;
- finalize figure captions and supplementary cross-references;
- archive a release and add a permanent DOI;
- complete journal-format and language editing.