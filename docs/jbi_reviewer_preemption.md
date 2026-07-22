# Journal of Biogeography reviewer-preemption checklist

This document records the issues most likely to trigger major revision and the repository evidence required before submission.

## 1. Prospective versus exploratory claims

### Reviewer risk

The manuscript could imply that moisture niche breadth and its predicted direction were specified before the results were inspected.

### Required response

- Describe the spatial-organization question as theory-led.
- Describe moisture niche breadth as the focal reported association unless a dated preregistration or analysis plan proves prospective specification.
- Report all evaluated niche metrics in the supplement.
- Do not write that the primary hypothesis was confirmed solely because the strict-set permutation p-value is below 0.05.

## 2. Strict evidence subset

### Reviewer risk

The baseline-unambiguous subset may appear to have been selected because it produced the strongest result.

### Required response

- Publish the exact inclusion criteria.
- Preserve a frozen species list with source-level support and decision notes.
- State that classification was performed without reference to climatic results.
- Log every post-freeze correction.
- Present the evidence-enriched estimate beside the strict-set estimate.

## 3. Model transparency

### Reviewer risk

Terms such as `family-adjusted` may conceal whether family was treated through fixed effects, random effects or clustered uncertainty.

### Required response

For each reported model provide:

- formula;
- package and version;
- n species and n families;
- within/among counts;
- coefficient, odds ratio and 95% CI;
- confidence-interval method;
- number of valid permutations;
- convergence and separation checks;
- missing-data rules.

The current implementation should be described as family-clustered uncertainty plus leave-one-family-out sensitivity unless the code is changed.

## 4. Small-sample inference

### Reviewer risk

The strict analysis contains only 35 species split 20/15, so asymptotic uncertainty may be fragile.

### Required response

- Lead with effect size and interval, not the p-value.
- Retain permutation inference.
- Show the leave-one-family-out distribution.
- Report whether any model shows separation or influential observations.
- Avoid numerical claims of broad generality.

## 5. Multiplicity

### Reviewer risk

Moisture breadth may be perceived as one positive result selected from several climatic metrics.

### Required response

- State how many metrics and thresholds were evaluated.
- Include the complete result matrix in the supplement.
- Label secondary analyses consistently.
- Do not use a single 0.05 threshold as the discovery rule.

## 6. Negative controls

### Reviewer risk

Non-significant fragmentation and environmental-turnover models could be overstated as proving no role for spatial structure.

### Required response

Use: `the examined coarse occurrence-cloud metrics did not account for the association`.

Do not use: `fragmentation and turnover were absent` or `the null hypothesis was proved`.

Explicitly state that unlabelled GBIF components are not populations, barriers or morph distributions.

## 7. Phylogenetic interpretation

### Reviewer risk

Family-clustered uncertainty is not a phylogenetic comparative method.

### Required response

- Describe the result as robust to represented-family deletion.
- Do not call it phylogenetically independent.
- Add a species-level phylogenetic analysis only if it can be implemented transparently without shrinking or distorting the evidence set.
- Keep incomplete phylogenetic control in the main limitations.

## 8. Scope of inference

### Reviewer risk

A literature-derived sample may be described as globally representative.

### Required response

The inferential population is documented cases of intraspecific flower-colour variation in the assembled evidence base, not all angiosperms and not worldwide prevalence.

## 9. Submission blockers

Do not submit until the following are present:

- complete manuscript text rather than an outline;
- strict-set 95% CI in the results table;
- explicit model formula and estimator;
- frozen classification manifest;
- complete alternative-metric table;
- figure showing strict and enriched estimates together;
- leave-one-family-out figure;
- source-level classification audit in the supplement.

## Recommended result sentence

> In the baseline-unambiguous evidence set, geographically structured flower-colour variation was associated with narrower species-level realised moisture niche breadth than documented within-population polymorphism. The estimate retained its direction after omission of each represented family but attenuated when less certain classifications were added. Because multiple climatic metrics were evaluated and morph-labelled localities were unavailable, we interpret this as an evidence-sensitive comparative association rather than a confirmatory mechanistic test.
