# FCP white–bee mechanism coverage protocol — 2026-09-25

Status: **outcome-blind coverage audit**. No flower-colour outcome may be opened in this stage.

## Purpose

Assess whether the 2026 curated global bee–plant interaction dataset has enough exact plant-species coverage in the frozen FCP third cohort to justify a separate pollinator-mechanism test.

## Sources

FCP species identities:
- third-cohort prospective biological artifact, workflow run `35177668182`;
- only the unique `species` column is read;
- `morph`, palette fractions, D, H2 vectors, coordinates and other colour outcomes are not read.

Independent bee–plant source:
- Noori et al. 2026 curated GloBI bee–plant dataset;
- Zenodo record `18303036`, version 3.1;
- use `GloBI_bee_plant_Interactions_Summary.csv` for the first coverage audit.

## Name matching

The audit selects the summary column with the largest exact normalized overlap with the frozen FCP species-name set. Selection uses species identity only.

Normalization:
- trim whitespace;
- collapse repeated whitespace;
- lowercase for matching;
- no genus-only substitution;
- no fuzzy matching;
- no synonym rescue in the primary coverage count.

## Gate

A species-level bee mechanism test may proceed only if:
- >=100 of the 499 frozen third-cohort species have exact species-name coverage in the independent bee summary.

If the exact gate fails:
- report `NOT_ESTIMABLE_CURRENT_BEE_DATASET`;
- do not inspect white/non-white outcomes against bee metrics;
- synonym harmonization may be developed later under a separately frozen outcome-blind protocol.

If the gate passes:
- freeze the bee metric and inferential model before opening flower-colour outcomes.

## Hard boundary

This coverage audit cannot support or reject pollinator causation. It only decides whether the independent bee dataset is sufficiently represented for a downstream test.
