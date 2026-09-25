# Moricandia anthocyanin-branch specificity analysis — 2026-09-25

## Status

This is a post-publication mechanistic follow-up built on the already opened Moricandia temperature-dose result. It is **not** an untouched confirmatory test and cannot rescue the failed cross-cohort BIO5 replication.

The question is narrower:

> Under the same summer photoperiod, does increasing temperature depress the anthocyanin branch more strongly than the UV-absorbing flavonol branch?

This distinguishes an anthocyanin-specific / branch-allocation response from a nonspecific collapse of floral flavonoid metabolism.

## Public source

Gómez JM et al. 2020, Nature Communications 11:4019, DOI 10.1038/s41467-020-17875-1.

Source Data sheet: `Figure 3F-J`.

The experimental sequence already established in the temperature-dose reanalysis is:
- spring: weighted temperature 14.1667 C;
- mild summer: weighted temperature 23.75 C, 16/8 h photoperiod;
- hot summer: weighted temperature 28.75 C, 16/8 h photoperiod.

Primary panel:
- the same 15 individuals with the sequence spring -> mild summer -> hot summer.

The individual is the inferential replicate; replicate flowers are averaged within individual × period.

## Traits

Primary anthocyanin trait:
- `Cyanidin` concentration (cyanidin-3-glucoside equivalents per fresh weight).

Comparator branch:
- the unique source-data column whose name contains `flavonol` (case-insensitive), corresponding to UV-absorbing flavonols / kaempferol-3-glucoside equivalents per fresh weight.

The script must fail if no unique flavonol column is present.

## Primary branch-specificity statistic

For each individual and each trait:

`
R_trait = hot / mild
`

All mild and hot means must be strictly positive; no pseudocount is allowed.

Transform:

`
L_trait = log(R_trait)
`

Primary within-individual contrast:

`
B_i = L_cyanidin - L_flavonol
`

Prediction:

`
B_i < 0
`

meaning anthocyanin declines proportionally more strongly than flavonol.

Report:
- median and mean hot/mild ratio for each trait;
- median `B_i`;
- fraction of individuals with `B_i < 0`;
- paired Wilcoxon signed-rank test on `B_i`, one-sided alternative < 0;
- two-sided Wilcoxon companion;
- directional sign-test p-value.

## Trait-level descriptive tests

For cyanidin and flavonol separately, report paired mild vs hot:
- mean and median concentration;
- mean and median hot-minus-mild delta;
- hot/mild ratio of cohort means;
- paired two-sided Wilcoxon p-value.

The cyanidin result is a reproduction of the already opened temperature-dose result and is not a new discovery.

## Decision label

`ANTHOCYANIN_BRANCH_SPECIFIC_HEAT_RESPONSE_SUPPORTED` only if:
1. all 15 primary individuals have positive mild and hot values for both traits;
2. median `B_i < 0`;
3. at least 12/15 individuals have `B_i < 0`;
4. one-sided Wilcoxon p < 0.05.

Otherwise:
`ANTHOCYANIN_BRANCH_SPECIFICITY_NOT_SUPPORTED_UNDER_THIS_TEST`.

No equivalence test for flavonols and no post hoc effect-size threshold is allowed.

## Interpretation ceiling

A positive result supports:

> hotter summer conditions shift floral flavonoid allocation away from anthocyanin more strongly than from the UV-absorbing flavonol branch in Moricandia.

It does not prove:
- flux through a particular enzymatic branch by isotope tracing;
- that every FCP white system uses this allocation mechanism;
- that temperature alone explains spring-vs-summer transcriptomic differences;
- that global BIO5 is a replicated driver;
- evolutionary nonwhite -> white transition asymmetry.
