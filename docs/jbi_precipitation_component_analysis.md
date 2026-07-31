# Precipitation-component decomposition of the legacy moisture-breadth result

## Question

The manuscript defined `moisture_breadth` as the arithmetic mean of the 95th–5th percentile ranges of BIO12, BIO14, BIO15 and BIO17. This diagnostic asks whether the reported association is shared across those variables or carried by one component.

BIO12, BIO14 and BIO17 are precipitation amounts in millimetres. BIO15 is precipitation seasonality expressed as a coefficient of variation. The unstandardised arithmetic mean is therefore not unit-homogeneous.

## Methods

The four component breadths were analysed separately for the same 34 baseline-unambiguous species, using:

- `among ~ component_z + effort_z`;
- standardised `log1p(n_climate_cells)` as the effort covariate;
- family-clustered sandwich uncertainty;
- 9,999 common label permutations within each occurrence dataset;
- leave-one-family-out refits;
- Holm correction across the four component tests;
- primary and paginated quality-filtered occurrence datasets;
- Open Tree topology and fixed dated-megaphylogeny sensitivity models.

The legacy arithmetic-mean moisture breadth was retained only as a comparator.

## Main result

The legacy composite is almost a re-expression of BIO12 annual-precipitation breadth:

- primary: Pearson r = 0.990, R² = 0.980;
- paginated: Pearson r = 0.988, R² = 0.977.

This follows from raw-scale dominance: the primary-sample SD of BIO12 breadth was 482.4 mm, compared with much smaller raw SDs for BIO14, BIO15 and BIO17.

## Non-phylogenetic models

| Dataset | Component | OR | 95% CI | Wald p | Holm p | Permutation p | Holm permutation p | LOO OR range |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| primary | bio12_q95q05 | 0.325 | 0.119–0.891 | 0.0289 | 0.1156 | 0.0200 | 0.0800 | 0.228–0.380 |
| primary | bio14_q95q05 | 0.796 | 0.387–1.638 | 0.5362 | 1.0000 | 0.5824 | 1.0000 | 0.671–1.022 |
| primary | bio15_q95q05 | 1.279 | 0.573–2.851 | 0.5482 | 1.0000 | 0.5338 | 1.0000 | 1.068–1.666 |
| primary | bio17_q95q05 | 0.766 | 0.378–1.552 | 0.4594 | 1.0000 | 0.5146 | 1.0000 | 0.643–0.960 |
| paginated_quality_filtered | bio12_q95q05 | 0.218 | 0.080–0.593 | 0.0029 | 0.0116 | 0.0055 | 0.0220 | 0.174–0.269 |
| paginated_quality_filtered | bio14_q95q05 | 0.695 | 0.341–1.419 | 0.3184 | 0.8085 | 0.3639 | 0.9942 | 0.544–0.861 |
| paginated_quality_filtered | bio15_q95q05 | 0.763 | 0.395–1.475 | 0.4210 | 0.8085 | 0.4864 | 0.9942 | 0.586–0.870 |
| paginated_quality_filtered | bio17_q95q05 | 0.676 | 0.338–1.354 | 0.2695 | 0.8085 | 0.3314 | 0.9942 | 0.531–0.833 |

BIO12 was the only component with clear non-phylogenetic support. In the primary sample its nominal tests were supported but did not remain below 0.05 after four-component Holm correction. In the paginated sample BIO12 remained supported after Holm correction:

- OR = 0.218;
- 95% CI = 0.080–0.593;
- Wald p = 0.0029, Holm p = 0.0116;
- permutation p = 0.0055, Holm permutation p = 0.0220;
- all leave-one-family-out estimates were negative, OR = 0.174–0.269.

BIO14, BIO15 and BIO17 showed no comparable support. BIO15 even changed direction between the primary and paginated non-phylogenetic analyses, indicating that the result is not a general precipitation-seasonality pattern.

## Phylogenetic sensitivity

BIO12 retained a negative direction under both phylogenetic treatments, but intervals included one:

- primary Open Tree: OR 0.543, 95% CI 0.210–1.409, p = 0.209;
- primary dated megaphylogeny: OR 0.412–0.417, CI envelope 0.143–1.189, p = 0.101–0.102;
- paginated Open Tree: OR 0.420, 95% CI 0.137–1.284, p = 0.128;
- paginated dated megaphylogeny: OR 0.290–0.291, CI envelope 0.083–1.019, p = 0.0520–0.0535.

The other components were also phylogenetically unresolved.

## Interpretation

The defensible interpretation is not that a broad composite “moisture niche” predicts spatial organization. The analysis instead identifies an exploratory association with **species-level occupied annual-precipitation breadth (BIO12)**. Species documented as geographically differentiated had lower odds of that category as annual-precipitation breadth increased.

This remains observational and cannot identify morph-specific precipitation tolerance, local adaptation or causation. It is stable to stronger occurrence sampling and family deletion, but statistically unresolved after phylogenetic correction.

## Manuscript implication

Recommended revision:

1. Remove the arithmetic-mean `moisture_breadth` as the focal biological variable.
2. Report BIO12, BIO14, BIO15 and BIO17 separately.
3. Treat BIO12 annual-precipitation breadth as the exploratory focal component, with all four-component tests and Holm adjustment shown.
4. Replace “moisture niche breadth” with “occupied annual-precipitation breadth” in the title, abstract and main interpretation if this component-focused framing is adopted.
5. Retain the legacy composite only as a diagnostic illustrating why the earlier result was BIO12-dominated.
