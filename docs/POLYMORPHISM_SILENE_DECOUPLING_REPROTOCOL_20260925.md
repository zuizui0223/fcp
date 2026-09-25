# Silene littorea natural decoupling reproducibility protocol — 2026-09-25

## Role
Post-publication reproducibility and mechanistic triangulation of Del Valle et al. 2019, BMC Plant Biology 19:496, DOI 10.1186/s12870-019-2082-6. This is not an untouched confirmatory study and does not alter the frozen New Phytologist manuscript.

## Biological contrast
The paper distinguishes:
- fully pigmented plants;
- PAL (petal anthocyanin loss): white petals while anthocyanins remain in photosynthetic tissues;
- WAL (whole-plant anthocyanin loss): white petals plus broad anthocyanin loss from photosynthetic tissues.

The mechanistic contrast is PAL versus WAL within the same species.

## Public source
Supplementary Table S2, file 12870_2019_2082_MOESM2_ESM.docx.

## Reproducibility rules
1. Download the original DOCX from Springer Nature.
2. Parse every DOCX table using python-docx; do not manually transcribe selected rows.
3. Preserve original table/row order and cell text.
4. Identify population, year, population-size and PAL/WAL frequency columns only from the parsed schema.
5. Convert frequency cells only when unambiguous.
6. Never infer individual counts from rounded percentages unless exact-count fields exist.

## Descriptive outputs
For PAL and WAL separately report:
- surveyed populations and population-year observations;
- populations ever positive and population-years positive;
- positive-frequency minimum, median and maximum;
- populations positive in at least 2 and at least 3 surveyed years;
- per-population trajectories;
- the longest repeated PAL series and its observed range.

## Inference boundary
No new p-value is required to reproduce the published persistence pattern. A paired PAL-versus-WAL test may be added only if both quantities are unambiguous on the same population-year rows and its statistic is specified before calculation. It cannot be described as independent confirmation.

## Chemistry/tissue evidence
Main-paper/source tables may separately verify that PAL retains anthocyanins in photosynthetic tissues, WAL lacks them broadly, and flavone production is not simply absent in WAL.

## Claim ceiling
This analysis can support a within-species natural analogue of floral-specific pigment decoupling: restriction of anthocyanin loss to petals is associated with substantially greater persistence/frequency than whole-plant anthocyanin loss within S. littorea.

It does not establish tissue restriction as the sole causal reason for persistence.
