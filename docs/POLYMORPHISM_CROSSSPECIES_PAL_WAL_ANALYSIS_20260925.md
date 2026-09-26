# Cross-species PAL versus WAL frequency analysis rules — 2026-09-25

## Status
Post-publication comparative synthesis using Supplementary Table S1 from Del Valle et al. 2019. The table values were visible during the schema preflight; therefore this is not an untouched confirmatory analysis.

## Coverage
The source contains 13 PAL rows and 13 WAL rows. Both phenotype classes span more than five source-listed taxonomic families/groups, so the frozen coverage gate is passed.

## Parsing rules

### PAL
- exact numeric x -> interval [x, x];
- range a-b -> interval [a, b];
- parenthetical mean values are recorded separately but are not used to collapse the interval.

### WAL
- exact numeric x -> [x, x];
- '< x' -> censored interval [0, x], with x retained as an upper bound;
- qualitative 'rare' or 'extremely rare' -> qualitative only, no numeric imputation;
- rows explicitly labeled greenhouse are excluded from the natural numeric panel but retained in the source registry.

No midpoint imputation is used.

## Descriptive comparison

Report:
1. PAL lower-bound median and upper-bound median;
2. WAL numeric upper-bound median and maximum;
3. number of PAL instances whose lower bound exceeds the maximum quantified WAL upper bound;
4. minimum PAL upper endpoint and maximum WAL upper endpoint;
5. counts of qualitative WAL descriptions ('rare', 'extremely rare');
6. full source-level registry.

No cross-species p-value is primary because the literature table is strongly ascertained and contains interval-censored and qualitative values.

## Claim ceiling

The analysis may support:
> across the literature systems assembled in Table S1, petal-restricted anthocyanin-loss polymorphisms commonly attain frequencies one to two orders of magnitude higher than whole-plant anthocyanin-loss variants.

It cannot establish that tissue restriction causally causes the frequency difference; mutation rate, demographic history, study ascertainment, and other selection can contribute.
