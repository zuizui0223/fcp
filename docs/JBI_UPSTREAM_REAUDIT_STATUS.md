# JBI upstream re-audit status

## Current stage

The historical 34-species freeze remains unchanged. Independent duplicate reviewer screening was intentionally waived on **2026-08-26** and is no longer a blocking gate.

The canonical replacement path is now:

**v2.2 retrieval → all 12,064 records retained → strict positive-evidence single pass → conservative GBIF taxon resolution → source-level local/geographic positive axes → species-level mixed-preserving aggregation → versioned replacement freeze → climatic analysis.**

The replacement analysis must be described as a **single-pass conservative documented-evidence audit**, not as an independently double-reviewed systematic review.

Detailed current protocol: `docs/JBI_SINGLE_PASS_STRICT_PROTOCOL.md`.

## Completed retrieval boundary

OpenAlex title/abstract OQL v2.2 remains fixed:

- 15 query blocks;
- 13,911 query memberships;
- 12,064 deduplicated works;
- 0 truncated v2.2 blocks;
- historical exact source recovery 34/34;
- benchmark recovery 34/34 by exact/version/direct-citation diagnostics.

All 12,064 records remain visible in the audit universe. Automatic taxon extraction cannot silently remove a record.

## Duplicate-review infrastructure — archived, not required

The repository still contains:

- 384-record Wave 0 calibration materials;
- independent Reviewer 1/Reviewer 2 packages;
- agreement/kappa gate code and boundary tests;
- full B01–B13 duplicate-review packages;
- merge/adjudication safeguards.

These are retained for provenance or optional sensitivity work only. Issue #18 was closed as `not_planned` after the project decision to proceed without independent duplicate review.

## Strict positive-evidence rule

A spatial axis can be positive only from explicit source wording.

### Local coexistence

`local_coexistence_documented = 1` requires discrete floral-display colour variation plus explicit same-population/site coexistence or equivalent direct wording.

Generic `color/colour morph`, `polymorphic`, multiple colour forms, or multiple sampling locations cannot by themselves create a local positive.

### Geographic structure

`geographic_structure_documented = 1` requires explicit among-population/site/region/island differentiation or explicit geographic/spatial/clinal/latitudinal/altitudinal/elevational colour structure.

A broad range or many sampled localities without reported colour differentiation is insufficient.

### Absence is never inferred

Failure to detect a positive phrase means **not documented by the strict pass**, not biological absence. Ambiguous cases remain `unresolved`.

## Eligibility boundary

The strict pass stores explicit conflicts rather than force-classifying them. Cultivar/horticultural-line, induced mutation/transgenic/editing/tissue-culture, and ontogenetic colour-change wording creates `conflict_unresolved` unless independent natural evidence is available elsewhere.

Only sources with floral-display colour context, natural/field/population/geographic context, and relevant variation evidence enter `eligible_high_confidence`.

## Taxon-resolution boundary

Automatic binomial hints are used only for potentially informative records and are resolved against GBIF.

A source contributes to a species state only when it maps unambiguously to one accepted plant species.

The historical manifest is allowed only as an exact **source → taxon rescue map** for the 34 known benchmark papers. Historical spatial labels, climate-cell counts, model results, and effect directions are not supplied to the strict evidence builder.

## Mixed retained

Species-level documented states are generated from independent positive axes:

- local=1, geographic=0 → `within_evidence_only`;
- local=0, geographic=1 → `among_evidence_only`;
- local=1, geographic=1 → `mixed_evidence`;
- insufficient evidence → `unresolved`.

The two mixed axes may come from different papers and different parts of the range.

## Running now

Workflow: **`JBI v2.2 strict single-pass evidence audit`**.

Input is the fixed 12,064-record canonical artifact. The workflow is currently generating:

- `v22_single_pass_source_evidence.csv`;
- `v22_single_pass_species_states.csv`;
- `v22_single_pass_gbif_cache.json`;
- `v22_single_pass_summary.json`.

The first run is being treated as a diagnostic freeze candidate, not yet as a manuscript result. After it completes, the immediate checks are:

1. number of high-confidence positive sources;
2. within-only / among-only / mixed / unresolved species counts;
3. historical-34 source taxon rescue = 34/34;
4. conflict/unresolved burden;
5. whether strict wording is so conservative that one state becomes too sparse for stable modelling.

## Legacy diagnostics remain sensitivity only

Previous P1/P2/P3 outputs remain useful for explaining the original discrepancy but do not define the replacement state:

- P1: 209 GBIF-valid species; automated 109 among / 41 mixed / 26 within / 33 unresolved;
- P2/P3: 522 GBIF-valid species, all unresolved under the old directional flag system;
- P1–P3 union: 656 candidate species.

These numbers are not the new biological result.

## Manuscript boundary

No historical manuscript conclusion, figure, DOCX, or frozen climatic dataset is overwritten yet.

After the strict single-pass freeze is built and QC passes, the next steps are:

1. freeze the new source/species evidence tables with checksums;
2. compute literature-effort/observation-process variables;
3. recompute climatic metrics for the resolved replacement species set;
4. fit two parallel positive-evidence models and, if estimable, a non-ordinal within/among/mixed model;
5. rerun robustness checks, figures, Main/SI, references/data-source appendix, and DOCX;
6. retain the historical 34-species binary model as a provenance-preserving sensitivity analysis.
