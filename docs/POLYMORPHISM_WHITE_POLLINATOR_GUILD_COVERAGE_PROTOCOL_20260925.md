# FCP pollinator-guild coverage protocol — 2026-09-25

Status: **outcome-blind coverage audit**. No flower-colour outcome may be opened.

## Purpose

Determine whether a versioned, independent GloBI flower-visitation corpus has enough coverage in the frozen FCP third cohort to support a downstream test of pollinator-guild sorting of white versus non-white floral variation.

This follows a separate bee-breadth test. The bee-breadth null does not answer whether different pollinator guilds are associated with white flowers.

## Data sources

FCP:
- third-cohort prospective biological artifact from workflow run `35177668182`;
- read only the unique `species` column.

Independent interaction source:
- Global Biotic Interactions Review Dataset Corpus v3;
- Zenodo record `18064921`, published 2025-12-26;
- file `interactions.tsv.gz`;
- only flower visitation/pollination claims are used.

Accepted interaction types:
- `visitsFlowersOf`
- `flowersVisitedBy`
- `pollinates`
- `pollinatedBy`

Direction is normalized so the FCP taxon is always treated as the plant and the partner as the visitor/pollinator.

## Guild classification

Using the visitor/pollinator taxon path:

- `Sphingidae` -> hawkmoth;
- other `Lepidoptera` -> other_lepidoptera;
- `Diptera` -> diptera;
- `Hymenoptera` -> hymenoptera;
- `Aves` -> bird;
- `Chiroptera` -> bat;
- otherwise -> other_or_unresolved.

This is a coverage taxonomy, not yet an ecological outcome model.

## Exact plant-name matching

Primary matching is exact after:
- trimming whitespace;
- collapsing repeated whitespace;
- case normalization.

No fuzzy matching, genus substitution, or synonym rescue is allowed in this gate.

## Coverage outputs

For every exact-matched FCP plant species report:
- total flower-visitation/pollination records;
- number of distinct visitor taxa;
- record counts by guild;
- number of represented major guilds.

Report cohort-wide:
- species with >=1 flower interaction record;
- species with >=5 records;
- species with >=10 records;
- species with >=20 records;
- species with at least one Lepidoptera record;
- species with at least one Sphingidae record;
- species with both Hymenoptera and at least one non-Hymenoptera record.

## Gate

A general pollinator-guild composition test may be designed only if >=100 frozen FCP species have >=10 independent flower-visitation/pollination records.

A Lepidoptera-specific test may be designed only if >=50 FCP species have >=10 total flower-interaction records and at least one Lepidoptera record.

A hawkmoth-specific test may be designed only if >=30 FCP species have >=10 total flower-interaction records and at least one Sphingidae record.

Failure of a gate means `NOT_ESTIMABLE`, not evidence against that pollinator mechanism.

## Hard boundary

No morph, white fraction, D, H2 vector, colour palette, coordinate, climate variable, or other FCP biological outcome is read in this coverage stage.
