# Current-data research results

## Journal-facing scope

The study concerns **documented intraspecific flower-colour variation**. This umbrella deliberately includes two distinct spatial expressions:

1. **within-population flower-colour polymorphism**, where discrete colour morphs coexist within at least one population; and
2. **geographically structured flower-colour variation**, where colour states are differentiated primarily among populations or regions.

The second category should not be described as flower-colour polymorphism without qualification. The comparative question is whether species differ in occupied climatic niche according to which spatial expression is documented.

The dataset is a **global literature-derived comparative sample of documented cases**, not a random sample of angiosperms and not an estimate of worldwide prevalence.

## Analysis snapshot

The latest resolved acquisition artifact contains:

- 664 candidate species;
- 111 validated cases of natural intraspecific flower-colour variation; and
- 553 deferred candidates.

For the candidate-validation analysis, the estimand is **validation probability conditional on entering the candidate set**. For the spatial-scale analysis, the estimand is the association between documented spatial expression and species-level occupied climatic niche within the literature-derived sample.

## Ascertainment and publication effort

A binomial GLM with family-clustered standard errors showed a strong relationship between retained literature effort and validation:

- linear log-effort coefficient: 1.737;
- odds ratio per one-unit increase in `log(1 + n_works)`: 5.68;
- 95% CI: 3.66–8.81; and
- p < 1e-13.

A quadratic effort model fitted materially better than the linear model:

- linear AIC: 523.45;
- quadratic AIC: 512.98;
- ΔAIC: 10.46;
- quadratic term: −0.912; and
- p = 0.00028.

This pattern indicates saturation in confirmation probability with increasing retained literature. The decline predicted at the extreme should not be interpreted biologically because it likely reflects a small mixture of highly studied but ambiguous cases.

Using only P0/P1 cases gave 97 strict cases, and the nonlinear effort pattern remained:

- linear effort term p = 1.0e-6; and
- quadratic term p = 0.0032.

Thus, ascertainment is a major property of the evidence base and is not driven only by weaker P2/P3 cases.

## Effort-adjusted family signals

For families represented by at least five candidate species, the largest positive standardized residuals after equalising literature effort were:

| Family | Candidates | Observed cases | Expected cases | Standardized residual |
|---|---:|---:|---:|---:|
| Cactaceae | 7 | 5 | 0.93 | 4.61 |
| Liliaceae | 12 | 4 | 1.46 | 2.34 |
| Polemoniaceae | 10 | 5 | 3.04 | 1.51 |
| Ranunculaceae | 22 | 6 | 3.92 | 1.28 |
| Papaveraceae | 6 | 3 | 1.72 | 1.24 |
| Brassicaceae | 25 | 9 | 6.83 | 1.09 |

Largest negative residuals were:

| Family | Candidates | Observed cases | Expected cases | Standardized residual |
|---|---:|---:|---:|---:|
| Rosaceae | 20 | 1 | 3.56 | −1.63 |
| Lamiaceae | 13 | 0 | 1.77 | −1.49 |
| Apiaceae | 8 | 0 | 1.43 | −1.43 |
| Theaceae | 7 | 0 | 1.33 | −1.37 |
| Araceae | 5 | 0 | 0.96 | −1.26 |

These are hypothesis-generating deviations after effort standardisation, not phylogenetically corrected biological effects.

## Occupied climatic niche breadth

### Candidate-versus-control comparison

The matched analysis included 70 focal species with documented intraspecific flower-colour variation and 243 taxonomically matched control species. Across five occupied-niche metrics and multiple occurrence thresholds, there was no clear evidence that documented colour-variable species generally occupy broader climatic niches than matched species outside the candidate set. At the primary same-genus, 20-cell specification, all confidence intervals included an odds ratio of one.

This null result argues against treating broad climatic tolerance as a general property of documented flower-colour variation.

### Spatial expression of colour variation

The biologically more informative comparison concerns the spatial expression of variation. The evidence-enriched classification contained 61 species with climatic metrics: 39 documented as within-population polymorphism and 22 as geographically structured variation. In the primary 20-cell models, geographically structured variation tended to be associated with narrower occupied climatic space:

| Metric | Species | Odds ratio for geographically structured variation | 95% CI | p |
|---|---:|---:|---:|---:|
| PCA hull area | 55 | 0.698 | 0.450–1.081 | 0.107 |
| Moisture breadth | 55 | 0.599 | 0.347–1.032 | 0.0648 |

These full-classification estimates are suggestive rather than conventionally significant.

The pre-specified strict sensitivity subset retaining only baseline-unambiguous classifications produced a stronger moisture-breadth association:

- 35 species: 20 within-population and 15 geographically structured;
- moisture-breadth odds ratio: 0.368;
- two-sided permutation p = 0.0226 using 9,999 valid permutations;
- leave-one-family-out odds-ratio range: 0.258–0.416; and
- the direction was retained after omitting every represented family.

Within this strict evidence subset, species with geographically structured colour variation occupied narrower **species-level realised moisture niches** than species with documented within-population polymorphism.

This result does not establish that moisture caused colour differentiation, that individual colour morphs differ physiologically, or that geographically structured systems are locally adapted. Those mechanisms require colour-state-labelled locality data.

## Evidence-quality sensitivity

The association weakened when less certain enriched classifications were added. This sensitivity is a result rather than a nuisance to conceal: the inferred relationship depends on how confidently the spatial expression of colour variation is documented.

Accordingly, the result hierarchy for the manuscript is:

1. baseline-unambiguous classification as the primary inferential subset;
2. all available classifications as the generalisation sensitivity analysis; and
3. high-confidence enrichment as an explicitly separate subset, reported as non-estimable when sample size is inadequate.

The manuscript should describe the moisture association as a **robust but preliminary comparative signal concentrated in the strict evidence subset**.

## Observed spatial fragmentation

The spatial-fragmentation extension retained 59 classified species with at least 20 coordinate-bearing GBIF records: 37 within-population and 22 geographically structured cases from 33 families.

No sampled-distribution fragmentation metric showed a clear independent association with geographically structured colour variation after controlling for moisture breadth and sampling effort:

| Added GBIF point-cloud metric | Odds ratio | 95% CI | p | Model AIC |
|---|---:|---:|---:|---:|
| Median nearest-neighbour distance | 0.936 | 0.444–1.969 | 0.861 | 80.67 |
| 95% spatial extent | 0.886 | 0.343–2.289 | 0.802 | 80.61 |
| Number of 50 km components | 1.444 | 0.509–4.096 | 0.490 | 80.29 |
| Number of 100 km components | 1.145 | 0.422–3.104 | 0.791 | 80.61 |
| Largest-component fraction at 100 km | 0.642 | 0.293–1.407 | 0.269 | 79.43 |
| Occupied 1-degree grid cells | 1.375 | 0.393–4.806 | 0.618 | 80.48 |

The baseline moisture model had AIC 78.70, lower than every fragmentation-augmented model. The moisture-breadth coefficient remained negative after adding each fragmentation metric.

These analyses do not support the interpretation that geographically structured colour variation is explained simply by more fragmented GBIF occurrence clouds. The connected components are sampling-derived summaries rather than biological populations, dispersal barriers, or gene-flow estimates.

## Environmental turnover among occurrence components

The primary 100 km, minimum-three-record specification retained 51 species with at least two GBIF occurrence components: 34 within-population and 17 geographically structured cases.

The primary environmental-turnover association was null:

- environmental-turnover odds ratio: 1.007;
- 95% CI: 0.367–2.760; and
- p = 0.989.

Across six pre-specified threshold and component-size combinations, turnover odds ratios ranged from 0.900 to 1.273, none had p < 0.05, and estimates occurred on both sides of one.

This is a robust null for generic climatic separation among unsupervised occurrence clusters. It is not a test of environmental differentiation among documented colour states because the GBIF records are not labelled by morph or population colour state.

## Taxonomic non-independence

Family-clustered standard errors and leave-one-family-out analyses reduce sensitivity to single-family dominance, but they do not constitute a full phylogenetic comparative analysis. The current evidence therefore supports a cross-species association that is robust to represented-family deletion, not a fully phylogenetically independent evolutionary relationship.

This limitation must be stated explicitly unless a species-level phylogeny is added.

## Journal-facing interpretation

The present evidence supports four conclusions:

1. Documentation of intraspecific flower-colour variation is strongly shaped by nonlinear publication effort.
2. Documented colour-variable species do not generally occupy broader climatic niches than matched controls.
3. In the strict evidence subset, geographically structured variation is associated with narrower species-level realised moisture niches than within-population polymorphism.
4. Coarse GBIF fragmentation and generic climatic turnover among unlabelled occurrence clusters do not explain that association.

The central contribution is therefore not a universal climatic rule for flower-colour polymorphism. It is the first comparative test in this dataset of whether the **spatial expression of documented intraspecific flower-colour variation** covaries with species-level climatic niche breadth, while explicitly tracking evidence quality and ascertainment.

The next mechanistic step would require documented colour states or named populations linked to localities. Further clustering of unlabelled GBIF occurrences is not an adequate substitute.