# JBI upstream re-audit status

## Current stage

The historical 34-species freeze remains unchanged. The re-audit has moved **upstream of automatic taxon extraction** and all computationally resolvable upstream gates are now in place.

Canonical path:

**v2.2 retrieval → all-record duplicate blind screening → record adjudication/full-text escalation → focal taxon validation → source-level natural/spatial evidence review → species-level two-axis aggregation → new freeze → climatic analysis.**

No replacement climatic model has been fit and no historical manuscript conclusion has been overwritten.

## Completed: replacement search surface

OpenAlex title/abstract OQL v2.2 is the completeness-defining search layer.

- query blocks: **15**;
- raw query memberships: **13,911**;
- deduplicated works: **12,064**;
- duplicate memberships removed: **1,847**;
- truncated v2.2 query blocks: **0**;
- historical exact source recovery: **34/34**;
- historical exact-or-version recovery: **34/34**;
- direct/citation benchmark recovery: **34/34**;
- seven prespecified review seeds resolved: **7/7**.

The legacy 19 truncated shards remain archived for provenance but are no longer an active blocker.

## Completed: automatic-taxon failure diagnosis

A diagnostic automatic source-to-taxon pass produced about **1,684** GBIF-valid candidate species and **3,358** species×source rows but represented only **31/34** historical benchmark species, even though all 34 source papers are in v2.2.

Therefore automatic title/abstract binomial extraction is not a canonical inclusion step.

## Completed: canonical all-record blind screening universe

- reviewer-facing records: **12,064/12,064**;
- coordinator-key rows: **12,064**;
- hidden detected-binomial hints: **9,071** records;
- no detected binomial but retained: **2,993** records;
- historical exact source records retained: **34/34**;
- records without abstracts retained: **1,259**;
- preassigned blind batches: **13**;
- reviewer-facing query membership, historical status, and automatic taxon hints: **0 columns**.

Thus absent/abbreviated taxon names can no longer silently remove a source before review.

## Completed: Wave 0 calibration infrastructure

Successful workflow: `JBI v2.2 record-screening Wave 0 calibration`, run **`32835886607`**.

Wave 0 contains **384 unique blind records** with hidden, mutually exclusive strata:

- historical benchmark sources: **34**;
- no detected binomial: **100**;
- detected binomial: **100**;
- non-English: **100**;
- missing abstract: **50**.

All CI checks passed, including stratum/blinding validation and the agreement-scorer smoke test.

Reviewer 1 and Reviewer 2 receive separate copies of the same blind sheet. Reviewer 1 fills only `reviewer_1_*`; Reviewer 2 fills only `reviewer_2_*`. Neither receives the coordinator key or the other's completed sheet. Adjudication remains blank until both independent copies are locked.

## Legacy priority diagnostics — sensitivity only

### P1

- **543** records;
- **209** GBIF-valid species;
- automated navigation states: 109 among-only / **41 mixed** / 26 within-only / 33 unresolved;
- **383** blind source rows;
- historical species recovered through P1 candidate path: **22/34**;
- exact historical source papers outside P1: **15/34**.

### P2/P3

- **1,996** records;
- **2,978** primary names checked;
- **522** GBIF-valid species;
- **1,196** blind source rows;
- all **522 unresolved** under the old navigation signal system.

### P1–P3 union

- **656** GBIF-valid candidate species;
- **+447** species beyond P1;
- automated states: **124 among-only / 41 mixed / 30 within-only / 461 unresolved**;
- **1,625** blind source rows.

The old mixed count does not increase because P2/P3 did not carry directional within/among flags. These are screening-path diagnostics, not biological results.

## Why the old result diverged

The mismatch is a stack of upstream filters:

1. different discovery chains;
2. generic morph language interpreted as local coexistence;
3. narrower geographic-evidence vocabulary;
4. mixed candidates removed before modelling;
5. priority-dependent entry into candidate analysis;
6. automatic taxon extraction missing relevant source papers.

The replacement path removes all six as irreversible biological filters.

## Biological state to be adjudicated later

For every eligible resolved species:

- `local_coexistence_documented`: discrete natural floral-display colour variants documented within at least one same population;
- `geographic_structure_documented`: colour or morph-frequency differentiation documented among geographic units.

Aggregation:

- local only → within-only evidence;
- geographic only → among-only evidence;
- both → mixed/multiscale evidence;
- unresolved → not entered into primary state analysis until resolved.

A zero means **not documented in reviewed evidence**, not biological absence.

## Current blocker — exactly one external input

The computational/search side has reached its honest stopping point. The next required input is:

**two independently completed copies of the 384-record Wave 0 blind sheet.**

After both are returned, the pipeline can resume automatically:

1. merge reviewer columns by `record_review_id`;
2. calculate raw agreement and Cohen's kappa;
3. identify disagreements;
4. inspect hidden strata only after coding is locked;
5. revise the codebook only if systematic ambiguity is demonstrated;
6. run a second calibration wave if a material rule changes, otherwise start the 13-batch full screen;
7. after record adjudication, proceed to taxon resolution and source-level spatial evidence review.

The same person or AI agent is **not** counted as two independent reviewers.

## Manuscript boundary

Until the new record/source review and adjudication are complete:

- historical 34-species results remain reproducible **historical/sensitivity** results;
- current Main/SI statistical conclusions are not replacement canonical inference;
- replacement figures/DOCX are not final;
- mixed-preserving climatic models remain planned inference.
