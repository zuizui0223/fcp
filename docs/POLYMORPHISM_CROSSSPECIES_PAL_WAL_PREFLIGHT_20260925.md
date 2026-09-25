# Cross-species PAL versus WAL natural-frequency preflight — 2026-09-25

## Purpose
Post-publication comparative reproducibility using Supplementary Table S1 from Del Valle et al. 2019 (DOI 10.1186/s12870-019-2082-6).

The biological question is whether naturally occurring flower-color variants with anthocyanin loss restricted to petals (PAL-like) tend to reach higher natural frequencies than white variants with anthocyanin loss extending through the whole plant (WAL-like).

This analysis is separate from the frozen New Phytologist manuscript.

## Source
Supplementary Table S1:
`12870_2019_2082_MOESM1_ESM.docx`

The published table summarizes frequencies of:
- petal anthocyanin loss (PAL) in natural flower-colour polymorphisms;
- whole-plant anthocyanin loss (WAL) spontaneous white mutants.

## First stage: schema-only preflight
Before calculating any cross-species frequency statistic:
1. download the original DOCX;
2. parse all tables and all rows with python-docx;
3. preserve the source text exactly in a row-level CSV;
4. report only table dimensions and headers;
5. do not calculate PAL/WAL medians, ratios or p-values.

## Eligibility for downstream comparison
A species/instance may enter only if the source table provides an unambiguous natural frequency or frequency range and a clear PAL-versus-WAL class.

Horticultural-only or experimentally induced frequencies are excluded if identifiable from the source table.

Ranges are not collapsed until a rule is frozen after schema inspection.

## Downstream comparison gate
Proceed only if at least:
- 10 PAL frequency instances;
- 10 WAL frequency instances;
- 5 plant families in each class
are extractable without interpretation beyond the source table.

If coverage fails, do not lower the gate.

## Hard boundaries
This is literature-derived comparative evidence, not a random sample of angiosperms. Publication bias and discovery bias are expected.

The analysis can address natural frequency/persistence, not directly fitness or causal pleiotropic mechanism.
