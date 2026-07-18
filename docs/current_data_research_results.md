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

## Scientific interpretation

The first substantive result is that discovery effort is a major, nonlinear component of the observed evidence pattern. Any ecological model that ignores publication effort will be biased. Family differences remain after effort standardisation and identify candidate clades for the macroecological and phylogenetic model.

The next inferential step is to populate the species covariate table with island status, latitude, life history, pollination system, breeding system, range size, GBIF occurrences and phylogeny identifiers. The implemented model will then estimate ecological associations while retaining the effort term.
