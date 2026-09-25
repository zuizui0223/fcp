# Moricandia temperature-dose reanalysis — 2026-09-25

## Status

Post-publication public-data reanalysis. The published article already states that floral plasticity increases with summer-condition severity, and the source-data schema has been inspected. This is therefore **not an independent confirmatory test**.

The analysis asks a narrower mechanistic question:

> with summer photoperiod held constant, does a hotter treatment further suppress floral anthocyanin in the same individuals?

## Public source

Gómez JM et al. 2020, Nature Communications 11:4019.
DOI: 10.1038/s41467-020-17875-1.

Source Data workbook:
`41467_2020_17875_MOESM11_ESM.xlsx`, sheet `Figure 3F-J`.

The sheet contains:
- individual ID;
- period;
- season;
- floral cyanidin-equivalent anthocyanin concentration;
- weighted temperature.

Two flowers are measured per individual-period. The individual, not the flower, is the inferential unit.

## Fixed sequence groups after schema preflight

The source-data schema contains three repeated-measures sequences:

1. `mild_then_hot`: spring (14.17 C weighted mean) -> mild summer (23.75 C) -> hot/extra-summer (28.75 C);
2. `spring_control_then_hot`: spring -> spring control -> hot/extra-summer;
3. `mild_then_recovery`: spring -> mild summer -> spring recovery.

The exact labels and temperatures are read from the workbook rather than manually assigned by individual ID.

## Primary dose-response contrast

Primary mechanistic contrast:

`mild_then_hot: Period Third - Period Second`

This compares the same individuals at:
- mild summer: 23.75 C weighted mean;
- hot summer: 28.75 C weighted mean.

Both summer conditions use the same 16/8 h day/night photoperiod in the published experimental design.

For each individual:
1. average the two flower-level cyanidin measurements within period;
2. calculate `delta_hot_minus_mild`;
3. report mean, median, number decreasing/increasing;
4. use paired Wilcoxon signed-rank test, two-sided.

A negative delta is the predicted direction.

No flower-level pseudoreplication is allowed.

## Triangulation contrasts

These are validation/interpretive contrasts, not separate confirmatory hypotheses:

- `spring_control_then_hot`: Period Third - Period Second;
- `mild_then_recovery`: Period Third - Period Second;
- spring -> mild change in both groups that receive mild summer.

Expected qualitative pattern:
- spring -> mild: anthocyanin down;
- mild -> hot: further down;
- spring-control -> hot: down;
- mild -> spring recovery: up.

## Interpretation ceiling

A strong paired mild->hot decline supports a **temperature dose effect on floral anthocyanin under a fixed summer photoperiod in Moricandia arvensis**.

It does not establish:
- a universal cross-species BIO5 effect;
- that temperature alone explains all spring-to-summer transcriptional change;
- that the RNA-seq pathway changes were measured in the hot treatment (the published RNA-seq used spring vs mild summer);
- evolutionary nonwhite->white transition asymmetry.

## Link to FCP

The dose-response result is mechanistically relevant because the published spring->mild-summer RNA-seq shows significant downregulation of four FCP natural-white node classes: MYB, CHS, UGT and TTG1/WD40.

The combined inference is:
- spring -> mild summer: pigment-network repression plus anthocyanin loss;
- mild -> hotter summer: further anthocyanin loss under the same photoperiod.

This is a same-species mechanistic bridge, not a global ecological main effect.
