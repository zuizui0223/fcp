# Journal of Biogeography manuscript outline

## Working title

**Climatic niche breadth differs with the spatial organization of intraspecific flower-colour variation**

Alternative conservative title:

**Evidence-sensitive associations between moisture niche breadth and the spatial organization of flower-colour variation**

## Abstract logic

1. **Background** — Flower-colour variation is widely studied, but comparative work rarely distinguishes local coexistence from geographic differentiation.
2. **Question** — Do species differ in occupied climatic niche according to how intraspecific colour variation is spatially organized?
3. **Methods** — Build a global literature-derived evidence set, classify documented variation as within-population or among-population, calculate species-level occupied climatic niches from cleaned GBIF occurrences and climate, and test evidence and family sensitivity.
4. **Results** — The strict evidence subset shows a negative moisture-breadth association for among-population organization; the direction is stable to family omission but attenuates in the broader evidence-enriched set. Coarse fragmentation and generic environmental turnover are null.
5. **Conclusion** — Spatial organization of flower-colour variation carries a comparative climatic signal, but its strength depends on classification evidence. Mechanistic inference requires morph-labelled locality data.

## Introduction

### Paragraph 1 — General problem

Intraspecific phenotypic variation can be maintained locally, partitioned geographically or expressed as a mixture of both. These alternatives are biologically distinct because they imply different combinations of environmental heterogeneity, dispersal, gene flow and local demographic processes.

### Paragraph 2 — Why flower colour

Flower colour is unusually well documented and can respond to pollinators, abiotic stress, pleiotropy, drift and mating-system effects. However, most syntheses focus on why polymorphism exists, which mechanisms maintain it or how individual morphs vary geographically.

### Paragraph 3 — Gap

Existing macroecological studies relate flower colour to climate, and individual systems document geographic morph clines. What remains unresolved is whether the **spatial organization** of intraspecific flower-colour variation differs systematically with species-level occupied climate across taxa.

### Paragraph 4 — Conceptual distinction

Define:

- within-population flower-colour polymorphism: discrete variants coexist locally;
- geographically structured flower-colour variation: variants are documented among populations without demonstrated local coexistence.

Explain that treating both as interchangeable FCP obscures the biological scale of coexistence.

### Paragraph 5 — Predictions

Primary prediction: geographically structured cases occupy narrower realised moisture niches than within-population cases.

Alternative/null explanations:

- the association is an artefact of classification evidence;
- it is concentrated in one or a few families;
- GBIF sampling or coarse range fragmentation accounts for the result.

### Final paragraph — Study contribution

State the global literature-derived comparison, evidence audit, occurrence-derived niche calculation, family sensitivity and interpretation boundary.

## Methods

### 1. Literature-derived study set

Describe candidate acquisition, validation and the non-random nature of the documented sample. Define the estimand as applying to documented cases rather than all angiosperms.

### 2. Spatial organization classification

Operational definitions for `within_population`, `among_population`, `mixed` and `unclear`. Explain baseline-unambiguous and evidence-enriched evidence sets and manual source audit.

### 3. Occurrence data and climate

Describe GBIF record cleaning, coordinate filtering, deduplication, climate-cell occupancy, minimum occupied-cell thresholds and climate source.

### 4. Niche metrics

Make moisture niche breadth the primary response/predictor of interest. List multivariate and temperature metrics as secondary.

### 5. Comparative models

- logistic response: among-population versus within-population;
- standardized niche metric predictor;
- occurrence/research effort controls where applicable;
- family-clustered uncertainty;
- permutation inference;
- leave-one-family-out sensitivity;
- minimum sample-size rule and non-estimable handling.

### 6. Secondary negative controls

Describe point-cloud fragmentation and unsupervised environmental turnover as tests of coarse sampling-derived alternatives, not population processes.

## Results

### 1. Evidence base and classification

Report the number of documented cases in each class, number with source-level matching evidence and number requiring manual review. Keep acquisition effort results short unless needed to establish ascertainment bias.

### 2. Candidate-versus-control niche breadth

Report the null result briefly. This rules out the broad claim that documented FCP species generally occupy broader niches than matched controls.

### 3. Primary spatial-organization result

Lead with baseline-unambiguous moisture breadth:

- sample size and class counts;
- odds ratio and confidence interval;
- permutation p-value;
- leave-one-family-out direction and range.

Avoid describing a single p-value as the discovery criterion.

### 4. Evidence-quality sensitivity

Report attenuation in the evidence-enriched set. Treat this as a major result about classification uncertainty, not an inconvenience.

### 5. Alternative explanations

Report null or weak effects for coarse fragmentation and generic environmental turnover. State explicitly what these analyses cannot test.

## Discussion

### Paragraph 1 — Main answer

The spatial organization of documented flower-colour variation is associated with realised moisture niche breadth in the strict evidence set, with geographically structured cases tending toward narrower breadth.

### Paragraph 2 — Biological interpretations

Discuss compatible explanations:

- environmental sorting across regions;
- local coexistence persisting across broader realised climatic contexts;
- differences in dispersal or demographic structure correlated with climate breadth.

Do not select among these without morph-labelled data.

### Paragraph 3 — Relation to FCP literature

Explain why the result extends rather than repeats existing work: prior studies focus on single systems, average colour traits, or presence of polymorphism, whereas this paper compares the spatial scale at which colour variation is documented.

### Paragraph 4 — Evidence sensitivity

Classification attenuation demonstrates that comparative trait synthesis is sensitive to source quality. Emphasize the value of source-level audit and strict/broad evidence sets.

### Paragraph 5 — Null results

The absence of generic fragmentation and turnover effects argues against simple explanations based on sampled occurrence-cloud structure. It does not reject morph-specific environmental differentiation.

### Paragraph 6 — Limitations

- literature-derived non-random sample;
- modest strict-evidence sample;
- incomplete phylogenetic correction;
- species-level realised niches;
- possible under-documentation of local coexistence;
- no colour-state labels in GBIF.

### Final paragraph — Contribution

End with the general message that spatial organization should be treated as a core comparative property of intraspecific phenotypic variation, not collapsed into a binary polymorphic/monomorphic label.

## Main figures

1. **Conceptual figure** — local coexistence versus geographic differentiation and the species-level moisture niche hypothesis.
2. **Evidence flow** — candidate set, validated variation, baseline-unambiguous classification and evidence-enriched classification.
3. **Primary effect plot** — moisture breadth effect for strict and enriched evidence sets with uncertainty.
4. **Robustness plot** — leave-one-family-out estimates.

## Main tables

1. Classification definitions and evidence requirements.
2. Primary and sensitivity model results.
3. Sample sizes and class composition across analysis sets.

## Supplement

- complete source-level audit;
- candidate-versus-control analyses;
- all alternative niche metrics;
- fragmentation and turnover analyses;
- workflow manifests;
- species lists and family-composition checks.

## Claims checklist

Before submission, verify that the manuscript never claims:

- a global census of flower-colour polymorphism;
- climate causation;
- demonstrated local adaptation;
- morph-specific climatic niches;
- full phylogenetic control;
- that geographically separated variants necessarily coexist locally.
