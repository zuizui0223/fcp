# Moricandia floral anthocyanin temperature-dose reanalysis — 2026-09-25

## Status

This is a **post-publication, post-preflight reanalysis** of published Source Data from Gómez et al. (2020), DOI 10.1038/s41467-020-17875-1. The raw anthocyanin values were already visible in the public source-data preflight, so this analysis is not an untouched confirmatory test.

It cannot rescue the failed cross-cohort FCP BIO5 replication. Its purpose is mechanistic:

> when photoperiod is held at the summer setting, does a hotter temperature regime further reduce floral anthocyanin within the same Moricandia individuals?

## Source-data structure

Published Source Data sheet `Figure 3F-J` contains:
- individual ID;
- period;
- season label;
- cyanidin-equivalent anthocyanin concentration;
- weighted experimental temperature.

The observed complete sequence classes are:
1. 15 individuals: 14.1667 -> 23.75 -> 28.75 C;
2. 14 individuals: 14.1667 -> 14.1667 -> 28.75 C;
3. 14 individuals: 14.1667 -> 23.75 -> 14.1667 C.

Two flower measurements per individual-period are first averaged, so the individual, not the flower, is the inferential replicate.

The published Methods identify the summer regimes as:
- mild summer: 30/20 C day/night, mean 23.8 C, 16/8 h day/night;
- hot summer: 35/25 C day/night, mean 28.8 C, 16/8 h day/night.

Thus the 23.75 -> 28.75 contrast holds the stated summer photoperiod constant while increasing temperature severity.

## Primary temperature-dose contrast

Use only individuals with the source-data sequence:

`spring 14.1667 -> summer 23.75 -> extrasummer 28.75 C`.

For each individual:
- average replicate flowers within each period;
- calculate `delta_hot_minus_mild = cyanidin_hot - cyanidin_mild`.

Report:
- n individuals;
- mean and median cyanidin at mild and hot temperature;
- mean and median delta;
- fraction with delta < 0;
- paired Wilcoxon signed-rank test, one-sided alternative `hot < mild`;
- two-sided paired Wilcoxon as a descriptive companion;
- exact sign-test probability for the number of negative individual deltas.

No outlier removal or transformation is allowed.

## Order/reversibility checks

These do not rescue the primary contrast.

### Mild -> spring reversal

For individuals with:

`14.1667 -> 23.75 -> 14.1667 C`

test whether cyanidin increases from period 2 to period 3.

### Spring -> hot

For individuals with:

`14.1667 -> 14.1667 -> 28.75 C`

test whether cyanidin decreases from period 2 to period 3.

### Change-direction contrast

Compare period-2-to-period-3 cyanidin changes between:
- mild -> hot;
- mild -> spring reversal;

using a two-sided Mann-Whitney test.

This checks whether a generic third-period/time effect can explain the decline.

## Interpretation

A strong hot < mild result supports:

> increasing temperature severity under the same summer photoperiod can further depress floral anthocyanin in this Moricandia system.

Together with the public transcriptome result (MYB90, CHS, U75C1/U78D2 and TTG1 all lower in summer-white flowers), this provides a within-species bridge from thermal severity to pigment-network output.

## Boundaries

This analysis does not establish:
- a universal cross-species temperature effect;
- that BIO5 is a replicated FCP driver;
- that temperature rather than all other chamber differences explains every spring-versus-summer transcriptomic change;
- that all natural white/pigmented FCP systems are thermally plastic;
- evolutionary nonwhite -> white transition asymmetry.

The clean general hypothesis remains a genotype/morph x temperature interaction, not a universal temperature main effect.
