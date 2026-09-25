# White-environment mechanism status — 2026-09-25

## Prospective result

The frozen three-variable WorldClim mechanism test completed successfully in GitHub Actions run 36101574640 (artifact 10848504388).

The estimable panel contained 281 species after the frozen white/non-white support gate, response-blind high-clip exclusion and near-clip availability filter.

- BIO5 (maximum temperature of warmest month): gate PASS. Median within-species white-minus-nonwhite contrast = +0.069 SD; 57.3% of species were positive; Wilcoxon p = 0.0118; Holm p = 0.0354. Conditional logistic corroboration with continuous near-clip adjustment gave OR = 1.073 per within-species SD (95% CI 1.029–1.119), p = 0.000919.
- BIO14 (precipitation of driest month): not supported. Holm p = 0.692; OR = 0.993, p = 0.749.
- Mean WorldClim solar radiation: not supported. Holm p = 0.692; OR = 0.987, p = 0.544.

The primary signal is therefore specifically associated with warmest-month heat among the three pre-specified environmental filters.

## Post-result robustness

These diagnostics were performed only after the prospective result was opened and do not upgrade or redefine the original gate.

### Single-species influence

Leave-one-species-out Wilcoxon p-values for BIO5 ranged from 0.00866 to 0.01521. No individual species is required for the primary unadjusted species-level association.

Aggressive removal of species with the largest absolute BIO5 contrasts weakens the Wilcoxon result, indicating that effect magnitude is heterogeneous across species rather than uniformly shifted.

### Observer-paired sensitivity

Restricting inference to species-observer strata containing both white and non-white outcomes yielded 144 paired observer strata across 106 species. The BIO5 effect was not supported: conditional OR = 0.787 (95% CI 0.544–1.138), p = 0.202; species-level Wilcoxon p = 0.485.

An observer-balanced broader sensitivity requiring at least three observers per colour retained the predicted median direction across 352 species but was weaker (median delta +0.054 SD; Wilcoxon p = 0.0750; sign p = 0.0484).

### Local geographic matching

Greedy white/non-white matching within species also weakened the signal:
- <=25 km: 901 pairs, 246 species, species-level Wilcoxon p = 0.188; paired conditional p = 0.241.
- <=50 km: 1,372 pairs, 321 species, Wilcoxon p = 0.215; paired conditional p = 0.0765.
- <=100 km: 1,964 pairs, 405 species, Wilcoxon p = 0.611; paired conditional p = 0.135.
- <=250 km: 2,729 pairs, 450 species, Wilcoxon p = 0.642; paired conditional p = 0.0787.

## Current interpretation

The data support a repeatable broad-scale association between measured white states and warmer warmest-month climates across the third cohort. They do **not** yet support a universal local heat-selection mechanism. The signal weakens when comparisons are forced within the same observer or within short geographic distances, so observer geography, regional population structure, phenology, or other spatially correlated factors remain plausible explanations.

Accordingly the current mechanistic status is:

**broad-scale heat sorting candidate; local causal heat mechanism unresolved.**

This line remains separate from the frozen New Phytologist manuscript and does not alter its claim ceiling.
