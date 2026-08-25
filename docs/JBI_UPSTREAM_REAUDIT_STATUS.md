# JBI upstream re-audit status

## Current stage

The historical 34-species freeze remains unchanged. The re-audit has now moved **upstream of automatic taxon extraction**.

The canonical path is:

**v2.2 retrieval → all-record duplicate blind screening → record adjudication/full-text escalation → focal taxon validation → source-level natural/spatial evidence review → species-level two-axis aggregation → new freeze → climatic analysis.**

No new climatic model has been fit and no historical manuscript conclusion has been overwritten.

## Completed: replacement search surface

OpenAlex title/abstract OQL v2.2 is now the completeness-defining search layer.

- query blocks: **15**;
- raw query memberships: **13,911**;
- deduplicated works: **12,064**;
- duplicate memberships removed: **1,847**;
- truncated v2.2 query blocks: **0**;
- historical exact source recovery: **34/34**;
- historical exact-or-version recovery: **34/34**;
- direct/citation benchmark recovery: **34/34**;
- seven prespecified review seeds resolved: **7/7**.

The legacy 19 truncated shards (15 Crossref + 4 OpenAlex) remain archived for provenance but are no longer an active blocker.

## Completed: diagnosis of the automatic-taxon bottleneck

A diagnostic priority-independent taxon/source pass was attempted directly on the 12,064 v2.2 works.

Its build step produced approximately:

- **1,684** GBIF-valid candidate species;
- **3,358** species×source rows;
- historical benchmark species represented: **31/34**.

The missing benchmark taxa were *Convolvulus arvensis*, *Disa porrecta*, and *Lupinus perennis* even though the v2.2 search contains all 34 exact source papers.

Therefore automatic title/abstract binomial extraction itself can create **silent false-negative source-to-taxon attribution**. This diagnostic universe is not canonical.

## Completed: canonical all-record blind screening universe

Workflow: `JBI v2.2 canonical record screening universe`.

Successful canonical record-screen output:

- input v2.2 records: **12,064**;
- reviewer-facing blind rows: **12,064**;
- coordinator-key rows: **12,064**;
- records with a hidden automatically detected binomial hint: **9,071**;
- records with **no** detected binomial but retained for human review: **2,993**;
- historical exact source records retained: **34/34**;
- review batches: **13** (target 1,000 records; final batch 64);
- reviewer-facing query-membership columns: **0**;
- reviewer-facing historical-status columns: **0**;
- reviewer-facing automatic-taxon-hint columns: **0**.

Additional audit of the blind artifact found **1,259 records without abstracts**. These records are retained; missing abstract is a reason for conservative/full-text handling, not automatic exclusion.

This closes the major upstream failure mode in which absent/abbreviated taxon names could silently remove relevant source papers before review.

## Completed: Wave 0 calibration infrastructure

Workflow: `JBI v2.2 record-screening Wave 0 calibration`.

Successful run: `32835886607`.

The calibration workflow builds **384 unique blind records** with hidden mutually exclusive strata:

- historical benchmark sources: **34**;
- no automatically detected binomial: **100**;
- detected binomial: **100**;
- non-English: **100**;
- missing abstract: **50**.

Review order is deterministic but independent of stratum. Reviewer-facing material contains no stratum, query membership, historical status, or automatic taxon hints.

`score_v22_record_calibration.py` was smoke-tested on the blank 384-row sheet and correctly reports zero double-coded rows with null agreement coefficients. Once two independent reviewers complete Wave 0, it will calculate raw agreement and Cohen's kappa for record relevance, natural intraspecific variation, floral-display colour, and full-text requirement, plus normalized exact agreement for focal taxon text.

The coding rules are fixed in `docs/JBI_V22_RECORD_SCREENING_CODEBOOK.md`.

## Legacy priority diagnostics — sensitivity only

These results are retained to quantify how strongly the old screening path shaped the candidate universe. They no longer define canonical inclusion.

### P1 diagnostic

- P1 records: **543**;
- primary taxon strings checked against GBIF: **874**;
- GBIF-valid species: **209**;
- automated navigation states: 109 among-only / **41 mixed** / 26 within-only / 33 unresolved;
- blind species×source rows: **383**;
- historical species recovered through the P1 candidate path: **22/34**.

Direct historical-source rescue showed that all **34/34** source papers existed in the archived corpus and **15/34** were outside P1. Thus P1 cannot be an inclusion criterion.

### P2/P3 diagnostic

The 1,996 legacy P2/P3 records were processed through the corrected source-attribution and GBIF path.

- primary taxon names checked: **2,978**;
- GBIF-valid species: **522**;
- blind source-review rows: **1,196**;
- navigation state: all **522 unresolved**.

The unresolved-only state is expected: legacy P2/P3 records did not carry the directional within/among signals used by the old navigation classifier. It is not evidence that these species lack spatial structure.

### P1–P3 union diagnostic

- GBIF-valid species: **656**;
- species added beyond P1: **447**;
- automated navigation states: **124 among-only / 41 mixed / 30 within-only / 461 unresolved**;
- blind source-review rows: **1,625**;
- new mixed-navigation species from P2/P3: **0**;
- P1 species promoted to mixed by P2/P3: **0**.

The unchanged mixed count is an artifact of legacy directional-signal semantics, not a biological conclusion. The large +447 species expansion is further evidence that legacy priority strata cannot define the replacement dataset.

## Why the old result and the re-audit diverged

The divergence is now traceable to multiple upstream mechanisms rather than one bad species label:

1. the historical 34-species freeze and later 79,242-record map followed different discovery chains;
2. the old within regex accepted generic morph terminology as local coexistence evidence;
3. geographic evidence vocabulary was comparatively narrow;
4. mixed candidates were removed before exploratory ecological modelling;
5. P1 priority missed exact papers already known to be relevant;
6. automatic binomial extraction can miss relevant source papers even after search retrieval is complete.

The replacement pipeline therefore removes **all of those automated steps as irreversible biological filters**.

## Spatial representation to be adjudicated

For every eligible taxonomically resolved species:

1. `local_coexistence_documented` — direct evidence of discrete natural floral-display colour variants within at least one same population;
2. `geographic_structure_documented` — direct evidence of spatial differentiation among populations/regions/sites/islands or other geographic units.

Species aggregation:

- local only → within-only evidence;
- geographic only → among-only evidence;
- both → mixed/multiscale evidence;
- unresolved → outside the primary state analysis until resolved.

A zero means **not documented in reviewed evidence**, not proof of biological absence.

## Current blocker and next scientific action

The computational/search infrastructure is no longer the limiting step. The next irreducible scientific gate is **genuinely independent duplicate review**.

1. Two independent reviewers code the 384-record Wave 0 using the fixed codebook.
2. Agreement and disagreements are scored before any rule change.
3. Material ambiguities trigger codebook revision plus a new calibration wave.
4. After calibration passes, screen all 12,064 records in the 13 preassigned blind batches.
5. Adjudicate record relevance and full-text needs.
6. Only retained/uncertain records proceed to taxon validation and source-level spatial evidence coding.
7. Only after source adjudication is complete is a new species freeze created and climatic inference rerun.

The same person or AI agent is **not** counted as two independent reviewers. Automated triage may assist coordinators but cannot populate independent reviewer fields or adjudicated biological states.

## Manuscript boundary

Until the new record/source review and adjudication are complete:

- the historical 34-species results remain reproducible **historical/sensitivity** results;
- the current Main/SI statistical conclusions are not promoted to final canonical inference;
- no new figures/DOCX should be presented as final results from the replacement pipeline;
- the new mixed-preserving models remain preregistered/planned inference, not completed biological results.
