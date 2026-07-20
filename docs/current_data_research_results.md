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

## Scientific interpretation

The first substantive result is that discovery effort is a major, nonlinear component of the observed evidence pattern. Any ecological model that ignores publication effort will be biased. Family differences remain after effort standardisation and identify candidate clades for the macroecological and phylogenetic model.

The climatic analysis does not support a simple claim that FCP species universally occupy broad climatic niches. Instead, it points to a more specific and biologically useful hypothesis: the spatial organization of colour variation may differ with moisture-niche breadth. This distinction should guide the next analysis toward explicit geographic structure, population separation, and environmental turnover rather than another undifferentiated FCP-versus-control comparison.

The next inferential step is therefore to model spatially explicit environmental turnover and population separation for the evidence-classified species, while retaining literature effort and classification uncertainty. Island status, latitude, life history, pollination system, breeding system, range size, GBIF occurrence structure and phylogenetic identifiers remain required covariates for the broader macroecological model.
