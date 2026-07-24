# Moisture niche breadth and spatial organization of flower-colour variation

**Running title:** Flower-colour variation and climate

## Abstract

### Aim

Intraspecific phenotypic variation may occur through local coexistence or geographic differentiation, but comparative studies rarely treat spatial organization as a distinct property. We tested whether documented within-population flower-colour polymorphism and geographically structured flower-colour variation differ in occupied climatic niche breadth.

### Location

Global, literature-derived sample.

### Taxon

Angiosperms with documented natural intraspecific flower-colour variation.

### Methods

We classified documented flower-colour variation as within-population, among-population, mixed or unclear from retained source text. The audited baseline-unambiguous set preserved source identifiers, evidence snippets and classification decisions. We combined binary classifications with cleaned GBIF occurrences and WorldClim 2.1 data at 10 arc-min resolution. Binomial generalized linear models related among-population organization to standardized niche metrics while controlling for occupied-cell count. We used family-clustered sandwich covariance, 9,999 label permutations and leave-one-family-out refits. Five metrics were evaluated at four occurrence thresholds.

### Results

The baseline-unambiguous moisture model included 34 species from 25 families: 20 within-population and 14 among-population cases. Geographically structured variation was negatively associated with realized moisture niche breadth (odds ratio = 0.426, family-clustered 95% confidence interval = 0.184–0.985; clustered Wald p = 0.0460). Permutation support was borderline (two-sided p = 0.0556). The association remained negative after each represented family was omitted, with leave-one-family-out odds ratios of 0.317–0.481. The estimate was weaker in the broader evidence set (odds ratio = 0.563, 95% confidence interval = 0.292–1.085; permutation p = 0.0944).

### Main conclusions

The spatial organization of flower-colour variation showed an evidence-sensitive association with occupied moisture niche breadth. Because support differed between uncertainty procedures and weakened with less certain classifications, the result is exploratory. The analysis does not test morph-specific tolerance, local adaptation or climatic causation.

**Keywords:** climatic niche breadth, evidence synthesis, flower-colour variation, GBIF, geographic differentiation, intraspecific polymorphism, macroecology

## Introduction

Intraspecific phenotypic variation can be expressed through local coexistence, geographic differentiation or a combination of both. These spatial configurations are biologically distinct. Local coexistence requires multiple forms to persist under at least partly shared demographic and environmental conditions, whereas geographic differentiation may reflect spatial environmental variation, dispersal limitation, demographic history or other forms of regional structure. Combining these configurations into a single category can therefore obscure the spatial scale at which phenotypic variation is documented.

Flower colour provides a tractable system for examining this distinction because natural colour variation has been studied in relation to pollinators, abiotic conditions, pleiotropy, genetic drift, mating systems and demographic processes. **Not verified:** supporting citations have not been supplied in the latest manuscript and must be inserted from the primary literature before submission. Existing comparative work commonly asks whether colour variation occurs, which mechanisms maintain it or whether floral colour traits covary with climate. The spatial organization of intraspecific flower-colour variation has received less explicit treatment as a comparative species-level property. **Not verified:** this literature-gap statement requires a documented literature search and appropriate citations.

Here, *intraspecific flower-colour variation* is used as the umbrella term. *Within-population flower-colour polymorphism* denotes the documented coexistence of discrete colour variants within at least one natural population. *Geographically structured flower-colour variation* denotes differentiation among populations or regions without retained evidence of local coexistence. Cases supported at both spatial scales were classified as mixed and excluded from binary models. This distinction avoids applying “polymorphism” to cases in which local coexistence has not been demonstrated.

We tested whether these documented spatial configurations differ in species-level occupied climatic niche breadth across a literature-derived sample. The comparison between local coexistence and geographic differentiation was theory-led, but five climatic metrics were evaluated at four occurrence thresholds. Moisture niche breadth is therefore reported as the focal association identified within a broader model set rather than as a prospectively preregistered endpoint. We evaluated sensitivity to classification quality, taxonomic family, occurrence sampling effort and coarse properties of the sampled occurrence cloud. Our objective was to characterize a comparative association, not to infer climatic causation, morph-specific tolerance or local adaptation.

## Methods

### Literature-derived study population

The study population comprised species identified by the repository pipeline as documented natural cases of intraspecific flower-colour variation. The candidate pool was literature-derived rather than randomly sampled from angiosperms. Inference is therefore restricted to the assembled set of documented cases and does not estimate the prevalence or spatial organization of flower-colour variation across flowering plants as a whole.

The resolved dataset contained 664 candidate species, of which 111 were retained as validated natural cases. Validation probability increased with retained literature effort. This relationship characterized ascertainment in the evidence base rather than biological prevalence.

**Not verified:** the latest manuscript does not provide a complete narrative description of candidate discovery, search dates, bibliographic databases, search strings, screening personnel, duplicate handling, language restrictions or exclusion criteria.

### Spatial-organization classification

Each retained case was classified as `within_population`, `among_population`, `mixed` or `unclear` from retained source text. Within-population classification required explicit evidence that multiple discrete colour variants coexisted within at least one natural population. Among-population classification required documented geographic differentiation without retained evidence of local coexistence. Evidence for both configurations produced a mixed classification. Mixed and unclear cases were excluded from binary models.

The classifier was audited after source identifiers and evidence snippets were propagated into the analysis dataset. The audit identified an unrecognized plural phrase, `within populations`; the classification rule was corrected before the final analysis was rerun. The frozen manifest preserves species, family, spatial category, source title, DOI or OpenAlex identifier, evidence text and decision note, together with a SHA-256 digest and correction log.

**Not verified:** whether all classifications were made blind to climatic results is not stated explicitly. The submission should state the classification freeze date or commit and describe manual adjudication.

### Occurrence data and climatic variables

Species occurrences were obtained from GBIF through the repository workflow. The workflow requested up to 300 records per taxon. Records were summarized by unique combinations of species identity and extracted climatic values, which were treated as occupied climate cells.

WorldClim 2.1 bioclimatic layers at 10 arc-min resolution supplied annual mean temperature (BIO1), temperature seasonality (BIO4), maximum temperature of the warmest month (BIO5), minimum temperature of the coldest month (BIO6), annual temperature range (BIO7), annual precipitation (BIO12), precipitation of the driest month (BIO14), precipitation seasonality (BIO15) and precipitation of the driest quarter (BIO17).

**Not verified:** the latest manuscript does not fully document GBIF download dates, query parameters, taxonomic backbone, coordinate filters, basis-of-record filters, spatial deduplication, country-coordinate checks, cultivated-record handling or geospatial outlier exclusion.

### Climatic niche metrics

Climatic values were standardized across the combined occurrence dataset, and a three-axis principal components analysis was fitted to the nine variables. Temperature breadth was the mean 95th–5th percentile range of BIO1, BIO5, BIO6 and BIO7. Moisture breadth was the mean 95th–5th percentile range of BIO12, BIO14, BIO15 and BIO17. Climatic heterogeneity was the mean within-species standard deviation across the nine standardized variables. PCA dispersion was the mean Euclidean distance from the species centroid in the first three principal-component dimensions. PCA hull area was the convex-hull area in the first two dimensions. The five metrics were evaluated at minimum occupied-cell thresholds of 10, 20, 30 and 50, yielding 20 specifications.

These metrics describe realized occupied climatic space represented by the cleaned GBIF records. They do not estimate fundamental physiological tolerance and were not calculated separately for colour morphs.

### Comparative models

For each climatic metric and threshold, we fitted a binomial generalized linear model with a logit link:

`among ~ metric_z + effort_z`

The response equalled one for geographically structured variation and zero for within-population polymorphism. `metric_z` was the standardized climatic metric and `effort_z` was standardized `log1p(n_climate_cells)`. Models were fitted in `statsmodels`, with family-clustered sandwich covariance for Wald uncertainty.

For focal moisture comparisons, we additionally used 9,999 permutations of spatial-category labels and leave-one-family-out refits. These procedures assess sensitivity to category assignment and individual families but do not constitute a phylogenetic comparative analysis.

**Not verified:** final software versions, the exact permutation statistic and random seed, treatment of failed permutations, missing-value rules, confidence-interval construction and model-separation diagnostics are not reported in the latest manuscript.

### Coarse occurrence-cloud alternatives

GBIF point-cloud fragmentation and generic climatic turnover among unsupervised occurrence components were examined as coarse alternative explanations. Components were treated as summaries of the sampled occurrence cloud and not as populations, barriers, gene-flow units or morph distributions.

**Not verified:** exact formulas, thresholds and model estimates for these analyses require direct transcription from final artifacts.

## Results

### Evidence base and ascertainment

The resolved evidence base contained 664 candidate species and 111 validated natural cases. Validation probability increased with retained literature effort, indicating that the evidence base was research-effort dependent. The dataset is therefore a literature-derived sample rather than a representative census of angiosperms.

**Not verified:** exact effect estimates and uncertainty for the ascertainment analysis are absent from the latest manuscript.

### Candidate-versus-control comparison

The matched comparison included 70 focal species and 243 taxonomically matched controls. Across five metrics and four thresholds, there was no clear evidence that documented colour-variable species generally occupied broader climatic niches than matched controls.

**Not verified:** the complete estimates, confidence intervals, matching procedure and specification-level results must be provided in a table or Supporting Information.

### Audited spatial organization and moisture niche breadth

The baseline-unambiguous model included 34 species from 25 families: 20 within-population and 14 among-population cases. The standardized moisture-breadth coefficient was −0.854, corresponding to an odds ratio of 0.426. The family-clustered 95% confidence interval was 0.184–0.985, and the clustered Wald p-value was 0.0460. The model converged in four iterations.

The two-sided permutation p-value was 0.0556 from 9,999 valid permutations. Leave-one-family-out odds ratios ranged from 0.317 to 0.481 and remained below one after each represented family was omitted. The estimated direction was not attributable to a single represented family, although permutation support was borderline.

The broader evidence set included 51 species from 29 families and produced a weaker estimate (odds ratio = 0.563, 95% confidence interval = 0.292–1.085; permutation p = 0.0944). This attenuation indicates sensitivity to classification evidence quality.

### Multiplicity context

Five metrics were evaluated at four thresholds, producing 20 specifications. Moisture breadth is reported as the focal association because it showed the clearest evidence-sensitive pattern, not because it was the sole metric examined or a prospectively preregistered endpoint.

**Not verified:** the complete 20-specification results, including coefficients, confidence intervals, p-values, sample sizes and diagnostics, are absent from the latest manuscript.

### Coarse occurrence-cloud alternatives

The examined fragmentation and generic environmental-turnover metrics did not clearly account for the focal moisture-breadth association. These analyses do not establish the absence of spatial environmental sorting because GBIF records were not labelled by flower-colour state.

**Not verified:** exact estimates and uncertainty are absent from the latest manuscript.

## Discussion

Across the audited baseline-unambiguous evidence set, geographically structured flower-colour variation was associated with narrower realized moisture niche breadth than documented within-population polymorphism. The estimated direction remained negative in every leave-one-family-out analysis, but the family-clustered confidence interval only narrowly excluded an odds ratio of one and the permutation p-value was 0.0556. The association also weakened when less certain classifications were added. These results provide suggestive, evidence-sensitive comparative support rather than confirmatory evidence for a general climatic relationship.

Several explanations remain compatible with the pattern. Geographic differentiation may be documented more often in species occupying restricted moisture contexts; local coexistence may persist across broader realized conditions; or dispersal, demographic history, range geometry and research practice may covary with both variables. Because occurrences were not labelled by colour morph and the analysis used species-level occupied climates, the data cannot distinguish among these explanations or test local adaptation.

The classification audit was consequential for inference. Correcting an unrecognized plural phrase moved a case supported at both scales out of the binary strict set. Literature-derived trait categories should therefore be treated as observations with measurement uncertainty. Source-level evidence, transparent rules, frozen manifests and correction logs are part of the analytical evidence rather than merely data-cleaning documentation.

The broader-set attenuation reinforces this conclusion. The baseline-unambiguous odds ratio was 0.426, whereas the broader estimate was 0.563 and its confidence interval included one. This difference does not establish that one dataset is correct and the other incorrect; it shows that the estimate depends on how confidently spatial organization can be assigned from the literature.

The focal association arose within 20 evaluated specifications. The full metric-by-threshold matrix is therefore essential. The clustered Wald result and permutation result fall on opposite sides of the conventional 0.05 threshold, making a binary significant/non-significant description unhelpful. The defensible interpretation is that the estimated effect was negative and stable to family deletion, but support was borderline and sensitive to the inferential procedure and classification set.

Family-clustered covariance and family-deletion refits reduce sensitivity to single families but do not establish phylogenetic independence. Likewise, coarse GBIF component analyses evaluate sampled occurrence-cloud properties, not populations or morph distributions. Stronger tests require a phylogenetic comparative model and georeferenced records linked to verified colour states or named populations.

Additional limitations include the non-random, research-effort-dependent sample; modest baseline-unambiguous sample size; post-analysis selection of the focal metric; multiple evaluated specifications; incomplete phylogenetic control; and species-level realized climatic niches. Local coexistence may also be under-documented when studies report regional differentiation without exhaustive population-level sampling. These limitations restrict generality and preclude mechanistic inference.

Despite these constraints, the study identifies a useful comparative distinction. Treating local coexistence and geographic differentiation as separate spatial configurations exposes a dimension of intraspecific phenotypic variation lost when all cases are grouped as flower-colour polymorphism. The evidence suggests that this distinction may covary with occupied moisture niche breadth, while also showing that inference depends on source quality. Future tests using phylogenetic models and morph-labelled environmental data can determine whether the pattern generalizes and which processes generate it.

## Acknowledgements

**Not verified.** Contributors, institutional support, material support and funding sources are not identified in the latest manuscript.

## Author contributions

**Not verified.** A CRediT statement requires the verified author list and agreed contribution roles.

## Conflict of interest

**Not verified.** All authors must confirm the final statement.

## Data accessibility

Analysis scripts, source-level manifests, correction logs, model tables, complete metric-by-threshold outputs and robustness analyses are maintained in the GitHub repository.

**Not verified:** no permanent archive DOI or reviewer-private repository link is provided. Archive the exact analysis release before submission and replace this note with a complete Data Accessibility Statement.

## References

**Not verified.** No reference list is present. Every literature claim and in-text citation must be checked against primary sources.

## Tables

**Not verified.** No final tables, titles, footnotes or in-text citations are present.

## Figure legends

**Not verified.** No final legends or verified figure-number citations are present.

## Supporting Information

**Not verified.** Final Appendix, Table and Figure numbering and cross-references have not been supplied. The complete model matrix, classification audit, correction log, diagnostics, alternative-explanation analyses and species manifests require verified Supporting Information identifiers.