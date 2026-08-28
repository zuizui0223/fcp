# JBI Wave 0 independent-review handoff

## Purpose

This handoff is the only non-computational gate currently blocking the canonical v2.2 record-screening pipeline.

Two independent reviewers must complete the fixed 384-record Wave 0 calibration before B01–B13 may start. The same person or AI agent cannot count as both reviewers.

## Reviewer package

Use the GitHub Actions artifact produced by workflow:

`JBI v2.2 reviewer-only Wave 0 package`

Artifact name:

`jbi-v22-wave0-reviewer-only`

The artifact must contain exactly:

- `FCP_Wave0_Reviewer1.csv`
- `FCP_Wave0_Reviewer2.csv`
- `JBI_V22_RECORD_SCREENING_CODEBOOK.md`
- `README_REVIEWERS.md`

It must **not** contain any coordinator key, historical labels, query membership, hidden calibration stratum, automated taxon hints, climatic variables, or the other reviewer's decision columns.

## Assignment

- Reviewer 1 receives only `FCP_Wave0_Reviewer1.csv` plus the codebook/README.
- Reviewer 2 receives only `FCP_Wave0_Reviewer2.csv` plus the codebook/README.
- Reviewers work independently and do not exchange decisions before both files are locked.
- Reviewer 1 fills only `reviewer_1_*` fields.
- Reviewer 2 fills only `reviewer_2_*` fields.
- Neither reviewer edits `record_review_id` or bibliographic fields.

## Conservative review rule

When title/abstract evidence is insufficient, use `uncertain` and/or `full_text_required=yes`; do not force an exclusion from missing information.

## Prespecified pass gate

The gate is fixed before completed reviewer sheets are inspected:

- 384/384 records independently double-coded on all four gated fields;
- raw agreement >= 0.90 for record relevance;
- raw agreement >= 0.90 for natural intraspecific variation;
- raw agreement >= 0.90 for floral display colour;
- raw agreement >= 0.85 for full-text requirement;
- Cohen's kappa >= 0.60 for each gated field when estimable;
- undefined kappa caused by single-category marginals is reported, with the raw-agreement gate still enforced.

`focal_taxon_text` agreement is diagnostic rather than a formal gate.

## Return contract

Return exactly the two completed CSVs without renaming their reviewer columns:

- `FCP_Wave0_Reviewer1.csv`
- `FCP_Wave0_Reviewer2.csv`

The coordinator then runs the already implemented return pipeline:

1. strict merge by `record_review_id`;
2. immutable metadata-drift validation;
3. reviewer-column contamination validation;
4. raw agreement and Cohen's kappa;
5. prespecified PASS / FAIL / NOT_READY decision;
6. disagreement-only adjudication queue;
7. consensus auto-fill for independently identical fields;
8. explicit adjudication for remaining disagreements;
9. retained / full-text / excluded record queues;
10. B01 release only after Wave 0 PASS and completed required adjudication.

## If Wave 0 fails

Do not relax thresholds post hoc. Inspect disagreement patterns, version any material codebook change, retain the failed Wave 0 for provenance, and generate a new blinded calibration wave before full screening.

## Scientific boundary

Wave 0 evaluates reproducibility of **record screening only**. It does not decide within / among / mixed spatial state. Spatial evidence is coded later at source level after record adjudication and taxon resolution.
