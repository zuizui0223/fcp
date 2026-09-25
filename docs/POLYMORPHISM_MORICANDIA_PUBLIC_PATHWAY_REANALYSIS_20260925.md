# Moricandia public-data pigment-pathway reanalysis — 2026-09-25

## Status

This is a **post-publication reproducibility / mechanistic-triangulation analysis**. It is not outcome-blind: the published paper already reports that summer white flowers show reduced expression of several anthocyanin-pathway genes. Therefore this workflow cannot be used as an independent confirmatory test or to rescue the failed FCP BIO5 replication.

Its purpose is narrower:

> quantify, from the authors' public supplementary/source-data files, how the experimentally induced spring-lilac -> summer-white transition maps onto the same normalized pigment-network nodes identified in natural FCP white/pigmented systems.

## Public source

Gómez JM et al. 2020. *Within-individual phenotypic plasticity in flowers fosters pollination niche shift*. Nature Communications 11:4019. DOI: 10.1038/s41467-020-17875-1.

Experimental design relevant here:
- the same five plants (MAR-70, MAR-81, MAR-83, MAR-98, MAR-120) were sampled under experimental spring conditions and then mild-summer conditions;
- spring conditions: 20/10 C day/night, 10/14 h day/night;
- mild-summer conditions: 30/20 C day/night, 16/8 h day/night;
- the paper reports paired floral RNA-seq for these 10 samples;
- raw reads are BioProject PRJNA604514.

Public files used in this audit:
- Supplementary Dataset 2:
  `41467_2020_17875_MOESM5_ESM.xlsx`
- Source Data:
  `41467_2020_17875_MOESM11_ESM.xlsx`

## Frozen FCP node panel

The normalized FCP natural-white node classes are fixed before the workbook contents are parsed:

1. MYB
2. CHS
3. F3H
4. ANS
5. UFGT / anthocyanin-relevant UDP-glycosyltransferase class
6. FLS
7. bHLH
8. WD40 / TTG1-class regulator

No node may be added to the FCP panel after workbook inspection.

Context-only pathway genes may also be reported:
- PAL
- 4CL
- CHI
- DFR

They do not count toward FCP-node overlap.

## Extraction rules

Every workbook sheet is read without manual row selection.

For each FCP node:
- search all textual cells using a frozen alias regex;
- retain every matching row and its original sheet/row identity;
- detect any numerical column whose header contains `logFC`, `log2FC`, `fold change`, or equivalent;
- detect FDR / adjusted-p / p-value columns when present;
- never choose among duplicate transcripts based on direction or significance.

The published paper defines the transcriptomic comparison as summer relative to spring. Therefore:
- negative logFC = lower expression in summer-white flowers;
- positive logFC = higher expression in summer-white flowers.

If the workbook does not expose a machine-readable fold-change value, the node is recorded as `direction_not_machine_readable` rather than inferred from prose.

## Node-level summary

For a node with >=1 machine-readable matched transcript:
- report all transcript-level effects;
- node direction = sign of the median available summer-vs-spring logFC;
- node is `summer_down` only when median logFC < 0;
- node is `summer_up` only when median logFC > 0.

No significance threshold is required for the descriptive node direction.

A separate `significant_summer_down` flag is allowed only when the row itself contains an adjusted p/FDR < 0.05 and logFC < 0.

## Descriptive convergence quantities

Report:
- number of the 8 FCP nodes found in the workbook;
- number with machine-readable direction;
- number with median summer-down direction;
- number with at least one significant summer-down transcript;
- exact matched transcript rows.

These are descriptive reproducibility quantities. No enrichment p-value is allowed.

## Interpretation ceiling

This workflow can support:

> the experimentally induced white-flower state in Moricandia suppresses overlapping parts of the same pigment-regulatory network implicated by naturally segregating white/pigmented systems.

It cannot establish:
- temperature alone as the cue, because temperature and photoperiod changed together in the RNA-seq comparison;
- a universal positive BIO5 effect across species;
- evolutionary nonwhite -> white transition asymmetry;
- that the exact same ortholog or causal nucleotide is involved across taxa;
- that all white flowers use anthocyanins rather than betalains/carotenoids or structural coloration.

## Decisive future test

The clean causal experiment remains a factorial design separating temperature from photoperiod and crossing environment with natural colour genotype/morph:

`genotype/morph x temperature x photoperiod`

with pigment chemistry and targeted pathway expression measured on the same flowers.
