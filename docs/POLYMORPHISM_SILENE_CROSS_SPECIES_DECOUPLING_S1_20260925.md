# Cross-species PAL versus WAL frequency synthesis — 2026-09-25

## Role
Post-publication reproducibility of Supplementary Table S1 in Del Valle et al. 2019 (DOI 10.1186/s12870-019-2082-6). This is not independent confirmation and does not alter the frozen New Phytologist manuscript.

## Source
Supplementary Table S1:
`12870_2019_2082_MOESM1_ESM.docx`

The table explicitly separates:
- species with petal anthocyanin loss (PAL) phenotypes;
- species with whole-plant anthocyanin loss (WAL) spontaneous white mutants.

## Frozen parsing
Every row in the original DOCX table is parsed.

Frequency strings are handled conservatively:
- exact number: documented maximum = number;
- numeric range: documented maximum = upper endpoint;
- '< x': numeric upper bound = x, retained as censored rather than replaced by x/2;
- qualitative 'rare' / 'extremely rare': no numeric value is imputed;
- parenthetical text such as 'greenhouse' is retained as a flag.

## Primary descriptive contrast
No new inferential p-value is calculated because reporting formats and sampling designs differ among source studies.

Report:
1. number of PAL and WAL species;
2. number with numeric frequency information;
3. minimum, median and maximum documented PAL maxima;
4. minimum, median and maximum numeric WAL upper bounds;
5. number of qualitative WAL rare/extremely-rare entries;
6. whether the minimum documented PAL maximum exceeds the largest numeric WAL upper bound;
7. conservative separation factor = minimum PAL maximum / largest numeric WAL upper bound.

This statistic asks whether the frequency scales overlap even under conservative use of the reported ranges.

## Tissue-scope metadata
The table footnotes are preserved:
- PAL species generally retain anthocyanins in vegetative tissues, with stated exceptions where information is unavailable;
- WAL spontaneous mutants generally lack anthocyanins in vegetative tissues, with stated biological exceptions.

These footnotes define the literature synthesis used by the source paper; they are not independently re-adjudicated in this reproduction.

## Claim ceiling
A separated frequency scale supports the descriptive generalization that petal-restricted anthocyanin loss is reported at substantially higher natural frequencies than whole-plant anthocyanin-loss mutants across the systems assembled by Del Valle et al.

It does not establish:
- unbiased species sampling;
- identical survey effort;
- a causal effect of tissue restriction;
- absence of publication/ascertainment bias;
- an independently preregistered comparative meta-analysis.
