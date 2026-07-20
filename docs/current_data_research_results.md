# Current-data research results

## Analysis snapshot

This analysis uses the latest resolved acquisition artifact available during PR development:

- 664 candidate species
- 111 validated natural flower-colour polymorphism cases
- 553 deferred candidates

The estimand is **validation probability conditional on entering the candidate set**. It is not global prevalence among angiosperms.

## Main ascertainment result

A binomial GLM with family-clustered standard errors showed a strong relationship between retained literature effort and validation:

- linear log-effort coefficient: 1.737
- odds ratio per one-unit increase in `log(1 + n_works)`: 5.68
- 95% CI: 3.66–8.81
- p < 1e-13

A quadratic effort model fitted materially better than the linear model:

- linear AIC: 523.45
- quadratic AIC: 512.98
- ΔAIC: 10.46
- quadratic term: −0.912
- p = 0.00028

This indicates saturation: additional papers strongly increase confirmation probability at low effort, but marginal gain falls at high effort.

Predicted validation probabilities from the quadratic model were:

| Retained works | Predicted probability |
|---:|---:|
| 1 | 0.064 |
| 2 | 0.184 |
| 3 | 0.306 |
| 5 | 0.466 |
| 10 | 0.581 |
| 20 | 0.521 |

The decline at the extreme should not be interpreted biologically; it likely reflects a mixture of difficult, highly studied ambiguous cases and the small number of very high-effort species.

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

Cactaceae is the only large positive deviation in this snapshot. It is a hypothesis-generating signal, not yet a phylogenetically corrected biological conclusion.

## Sensitivity

Using only P0/P1 cases gave 97 strict cases. The nonlinear effort pattern remained:

- linear effort term p = 1.0e-6
- quadratic term p = 0.0032

Therefore the ascertainment result is not driven only by weaker P2/P3 cases.

## Occupied climatic niche breadth

The matched candidate-versus-control analysis included 70 focal FCP species and 243 taxonomically matched control species. Across five occupied-niche metrics and multiple occurrence thresholds, there was no clear evidence that verified FCP species generally occupy broader climatic niches than matched species outside the candidate set. At the primary same-genus, 20-cell specification, all confidence intervals included an odds ratio of one.

The more informative result concerned the **spatial scale of colour variation among verified FCP species**. The evidence-enriched classification contained 61 species with climatic metrics: 39 classified as within-population variation and 22 as among-population variation. In the primary 20-cell models, among-population variation tended to be associated with narrower occupied climatic space:

| Metric | Species | Odds ratio for among-population variation | 95% CI | p |
|---|---:|---:|---:|---:|
| PCA hull area | 55 | 0.698 | 0.450–1.081 | 0.107 |
| Moisture breadth | 55 | 0.599 | 0.347–1.032 | 0.0648 |

These full-classification results are suggestive rather than conventionally significant. However, the pre-specified strict sensitivity subset retaining only baseline-unambiguous classifications produced a stronger moisture-breadth result:

- 35 species: 20 within-population and 15 among-population
- moisture-breadth odds ratio: 0.368
- two-sided permutation p = 0.0226 using 9,999 valid permutations
- leave-one-family-out odds-ratio range: 0.258–0.416
- direction retained after omitting every represented family

Thus, within the strict evidence subset, species whose colour variation is partitioned among populations occupied narrower moisture niches than species maintaining colour polymorphism within populations. This result is compatible with stronger geographic environmental sorting or local differentiation in among-population systems, but it does not establish causation, physiological tolerance, or niche expansion.

## Observed spatial fragmentation

The spatial-fragmentation extension retained 59 classified species with at least 20 coordinate-bearing GBIF records: 37 within-population and 22 among-population cases from 33 families. All eight pre-specified models completed successfully.

No sampled-distribution fragmentation metric showed a clear independent association with among-population variation after controlling for moisture breadth and sampling effort:

| Added GBIF point-cloud metric | Odds ratio | 95% CI | p | Model AIC |
|---|---:|---:|---:|---:|
| Median nearest-neighbour distance | 0.936 | 0.444–1.969 | 0.861 | 80.67 |
| 95% spatial extent | 0.886 | 0.343–2.289 | 0.802 | 80.61 |
| Number of 50 km components | 1.444 | 0.509–4.096 | 0.490 | 80.29 |
| Number of 100 km components | 1.145 | 0.422–3.104 | 0.791 | 80.61 |
| Largest-component fraction at 100 km | 0.642 | 0.293–1.407 | 0.269 | 79.43 |
| Occupied 1-degree grid cells | 1.375 | 0.393–4.806 | 0.618 | 80.48 |

The baseline moisture model had AIC 78.70, lower than every fragmentation-augmented model. The integrated 100 km connectivity model had AIC 81.08. Therefore the available GBIF point-cloud fragmentation summaries did not improve model fit.

The moisture-breadth association remained negative after adding every fragmentation metric. It was strongest in the model containing largest-component fraction at 100 km (odds ratio 0.490, 95% CI 0.261–0.920, p = 0.0264) and in the integrated 100 km model (odds ratio 0.503, 95% CI 0.266–0.951, p = 0.0346). These are model-dependent estimates rather than multiplicity-adjusted confirmatory tests, but they show that the moisture result is not removed by simple observed-distribution connectivity controls.

This analysis does **not** support the claim that among-population colour variation is explained merely by more fragmented GBIF occurrence clouds. Connected components here are sampling-derived summaries, not biological populations, dispersal barriers, gene flow estimates, or causal mechanisms.

## Environmental turnover among occurrence components

The environmental-turnover analysis measured WorldClim-standardized climatic separation among distance-threshold components of each species' GBIF occurrence cloud. The primary 100 km, minimum-three-record specification retained 51 species with at least two components: 34 within-population and 17 among-population cases.

The primary turnover association was null:

- environmental-turnover odds ratio: 1.007
- 95% CI: 0.367–2.760
- p = 0.989

A pre-specified sensitivity analysis varied both the distance threshold and the minimum component size:

| Threshold | Minimum records per component | Species | Turnover OR | 95% CI | p | Moisture OR | Moisture p |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 km | 3 | 53 | 1.273 | 0.504–3.215 | 0.610 | 0.471 | 0.0759 |
| 50 km | 5 | 51 | 1.135 | 0.502–2.568 | 0.761 | 0.511 | 0.0799 |
| 100 km | 3 | 51 | 1.011 | 0.369–2.769 | 0.983 | 0.552 | 0.153 |
| 100 km | 5 | 47 | 0.900 | 0.379–2.142 | 0.812 | 0.589 | 0.129 |
| 200 km | 3 | 46 | 1.030 | 0.425–2.499 | 0.948 | 0.555 | 0.111 |
| 200 km | 5 | 38 | 0.966 | 0.371–2.514 | 0.943 | 0.589 | 0.149 |

All six specifications completed, none had p < 0.05, and the turnover odds-ratio range was 0.900–1.273. Four estimates were above one and two below one. The conclusion is therefore not sensitive to the chosen 50–200 km threshold or to requiring three versus five records per component.

This is a robust null for **generic climate separation among unsupervised GBIF occurrence clusters**. It does not test whether documented colour states occupy different environments, because the occurrence records are not labelled by colour morph or population state.

## Scientific interpretation

The first substantive result is that discovery effort is a major, nonlinear component of the observed evidence pattern. Any ecological model that ignores publication effort will be biased. Family differences remain after effort standardisation and identify candidate clades for the macroecological and phylogenetic model.

The climatic analysis does not support a simple claim that FCP species universally occupy broad climatic niches. Instead, it points to a more specific and biologically useful pattern: among-population colour variation is associated with narrower occupied moisture breadth in the strict evidence subset. Coarse GBIF fragmentation and generic climate turnover among unsupervised occurrence clusters do not explain that pattern.

Further threshold expansion of the same unlabelled occurrence-cluster analysis is not justified. The next inferential step must connect **documented colour states or named populations to localities**, then test whether those state-labelled localities differ in moisture conditions. Without morph-labelled locality evidence, additional clustering of generic GBIF records cannot adjudicate environmental sorting of flower-colour states. Island status, latitude, life history, pollination system, breeding system, range size, phylogeny identifiers and literature effort remain required covariates for the broader macroecological model.
