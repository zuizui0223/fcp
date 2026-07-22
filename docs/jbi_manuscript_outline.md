# Journal of Biogeography manuscript outline

## Working title

**Moisture niche breadth is associated with the documented spatial organization of intraspecific flower-colour variation**

More conservative alternative:

**Evidence-sensitive associations between moisture niche breadth and the spatial organization of flower-colour variation**

The title must not imply causation, preregistered confirmation or a universal climatic rule.

## Abstract logic

1. **Background** — Flower-colour variation is widely studied, but comparative work rarely distinguishes local coexistence from geographic differentiation.
2. **Question** — Do species differ in occupied climatic niche according to how documented intraspecific colour variation is spatially organized?
3. **Methods** — Assemble a global literature-derived evidence set, classify documented variation as within-population or among-population, calculate species-level occupied climatic niches from cleaned GBIF occurrences and climate, and evaluate evidence-set and family sensitivity.
4. **Results** — In the baseline-unambiguous evidence set, geographically structured cases were associated with narrower realised moisture niche breadth. The estimate retained its direction after omission of each represented family but attenuated in the broader evidence-enriched set. The examined coarse occurrence-cloud metrics did not account for the association.
5. **Conclusion** — The spatial organization of documented flower-colour variation carries an evidence-sensitive comparative climatic signal. Mechanistic inference requires morph-labelled locality data.

The Abstract must report the strict-set sample size, odds ratio, 95% confidence interval, permutation support and attenuation in the broader evidence set. It must not describe non-significant secondary models as proof of no effect.

## Introduction

### Paragraph 1 — General problem

Intraspecific phenotypic variation can be maintained locally, partitioned geographically or expressed as a mixture of both. These alternatives are biologically distinct because they may arise under different combinations of environmental heterogeneity, dispersal, gene flow and local demographic processes.

### Paragraph 2 — Why flower colour

Flower colour is unusually well documented and can respond to pollinators, abiotic stress, pleiotropy, drift and mating-system effects. Most syntheses, however, focus on why polymorphism exists, which mechanisms maintain it or how individual morphs vary geographically.

### Paragraph 3 — Gap

Existing macroecological studies relate flower colour to climate, and individual systems document geographic morph clines. What remains unresolved is whether the **documented spatial organization** of intraspecific flower-colour variation differs systematically with species-level occupied climate across taxa.

### Paragraph 4 — Conceptual distinction

Define:

- within-population flower-colour polymorphism: discrete variants coexist locally;
- geographically structured flower-colour variation: variants are documented among populations without demonstrated local coexistence.

Treating both as interchangeable FCP obscures the biological scale of coexistence.

### Paragraph 5 — Study question and analysis status

The broad comparison between local coexistence and geographic differentiation is theory-led. Multiple climatic niche metrics were evaluated. Moisture niche breadth is therefore presented as the **focal reported association**, not as a preregistered directional hypothesis unless a dated analysis record establishes prospective specification.

Alternative explanations evaluated include:

- sensitivity to classification evidence;
- concentration in one or a few families;
- GBIF sampling effort;
- coarse sampled-distribution fragmentation; and
- generic climatic turnover among unlabelled occurrence components.

### Final paragraph — Contribution

State the literature-derived comparison, source-level evidence audit, occurrence-derived niche calculation, family sensitivity and strict interpretation boundary.

## Methods

### 1. Literature-derived study set

Describe candidate acquisition, validation and the non-random nature of the documented sample. Define the estimand as applying to documented cases rather than all angiosperms.

### 2. Spatial-organization classification

Provide operational definitions for `within_population`, `among_population`, `mixed` and `unclear`.

For the baseline-unambiguous set, report:

- explicit source-level inclusion criteria;
- whether classification was performed blind to climatic results;
- the date or commit at which the list was frozen;
- the retained species manifest; and
- the correction log for any later changes.

Explain the evidence-enriched set as a generalization sensitivity analysis rather than a competing primary dataset.

### 3. Occurrence data and climate

Describe GBIF record cleaning, coordinate filtering, deduplication, climate-cell occupancy, minimum occupied-cell thresholds and climate source.

### 4. Niche metrics and multiplicity

List every evaluated climatic niche metric and threshold. Explain that moisture niche breadth became the focal reported metric after comparative evaluation unless prospective documentation establishes otherwise. Provide the complete metric-by-threshold result matrix in the supplement.

### 5. Comparative models

Report the exact model formula and distinguish clearly among:

- model covariates;
- family-clustered standard errors;
- leave-one-family-out sensitivity; and
- any future phylogenetic model.

For every reported model provide:

- number of species and families;
- within/among counts;
- coefficient, odds ratio and 95% CI;
- confidence-interval method;
- package and version;
- permutation procedure and valid permutation count;
- convergence and separation diagnostics; and
- missing-data rules.

### 6. Secondary negative controls

Describe point-cloud fragmentation and unsupervised environmental turnover as tests of coarse sampling-derived alternatives, not population processes. State that non-significance does not establish absence.

## Results

### 1. Evidence base and classification

Report the number of documented cases in each class, number with source-level matching evidence, number requiring manual review and the final frozen strict-set composition.

### 2. Candidate-versus-control niche breadth

Report the null result briefly. This argues against the broad claim that documented colour-variable species generally occupy broader climatic niches than matched controls.

### 3. Focal spatial-organization association

Lead with the baseline-unambiguous moisture-breadth estimate and report together:

- 35 species: 20 within-population and 15 geographically structured;
- coefficient on the log-odds scale;
- odds ratio;
- 95% confidence interval;
- two-sided permutation p-value based on 9,999 valid permutations;
- number of represented families; and
- leave-one-family-out odds-ratio range.

The 95% CI is a submission blocker and must be added from the model output before manuscript completion.

Avoid presenting the permutation p-value as the sole discovery criterion.

### 4. Evidence-quality sensitivity

Report the broader evidence-enriched estimate beside the strict-set estimate. Treat attenuation as a central finding about classification uncertainty.

### 5. Multiplicity context

Report the complete set of evaluated niche metrics and make clear which did and did not show comparable associations.

### 6. Alternative explanations

Use the wording:

> The examined coarse occurrence-cloud fragmentation and unsupervised environmental-turnover metrics did not account for the focal association.

Do not state that fragmentation or turnover was absent.

## Discussion

### Paragraph 1 — Main answer

In the baseline-unambiguous evidence set, geographically structured flower-colour variation was associated with narrower species-level realised moisture niche breadth than documented within-population polymorphism. The estimate retained its direction after omission of each represented family but weakened when less certain classifications were added.

### Paragraph 2 — Interpretation boundary

Discuss compatible explanations, including environmental sorting, dispersal or demographic structure correlated with climate breadth, and persistence of local coexistence across broader realised contexts. Do not select among them without morph-labelled data.

### Paragraph 3 — Exploratory status and multiplicity

State transparently that several climatic metrics were evaluated and that moisture breadth is the focal reported association. Emphasize effect size, uncertainty, permutation support and evidence sensitivity rather than a binary significance claim.

### Paragraph 4 — Relation to flower-colour literature

Explain that the contribution is the comparison of the documented spatial scale of variation, rather than another test of whether climate correlates with colour.

### Paragraph 5 — Evidence sensitivity

Classification attenuation shows that comparative trait synthesis depends on source quality. Present the strict/broad contrast as a substantive methodological result.

### Paragraph 6 — Negative controls

The examined coarse GBIF metrics did not explain the focal association. They cannot reject morph-specific environmental differentiation because records lack colour-state labels.

### Paragraph 7 — Limitations

- literature-derived non-random sample;
- modest strict-evidence sample;
- exploratory focal metric unless prospective documentation is found;
- multiple evaluated niche metrics;
- incomplete phylogenetic correction;
- species-level realised niches;
- possible under-documentation of local coexistence;
- no colour-state labels in GBIF.

### Final paragraph — Contribution

End with the message that spatial organization should be treated as a core comparative property of intraspecific phenotypic variation, while avoiding causal or universal claims.

## Main figures

1. Conceptual figure: local coexistence versus geographic differentiation.
2. Evidence flow: candidate set, validated variation, frozen baseline-unambiguous set and evidence-enriched set.
3. Effect plot: strict and enriched moisture-breadth estimates with 95% CIs.
4. Leave-one-family-out estimate plot.
5. Optional multiplicity panel showing all evaluated niche metrics without selectively displaying the focal result alone.

## Main tables

1. Classification definitions and evidence requirements.
2. Primary and sensitivity model results, including exact formula and estimator.
3. Sample sizes, family counts and class composition across analysis sets.

## Supplement

- frozen source-level classification audit;
- correction log;
- candidate-versus-control analyses;
- all niche metrics and thresholds;
- model diagnostics;
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
- prospective confirmation of moisture breadth without dated evidence;
- proof of no fragmentation or turnover effect; or
- that geographically separated variants necessarily coexist locally.
