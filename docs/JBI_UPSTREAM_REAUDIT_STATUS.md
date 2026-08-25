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

The legacy 19 truncated shards (15 Crossref + 4 OpenAlex) remain archived for provenance but are no longer an active blocker.

## Completed: diagnosis of the automatic-taxon bottleneck

A diagnostic priority-independent taxon/source pass directly on the 12,064 v2.2 works produced approximately:

- **1,684** GBIF-valid candidate species;
- **3,358** species×source rows;
- historical benchmark species represented: **31/34**.

It missed *Convolvulus arvensis*, *Disa porrecta*, and *Lupinus perennis* despite the exact source papers being present in v2.2. Therefore automatic title/abstract binomial extraction can create **silent false-negative source-to-taxon attribution** and is not the canonical inclusion step.

## Completed: canonical all-record blind screening universe

Successful workflow: `JBI v2.2 canonical record screening universe`.

- input v2.2 records: **12,064**;
- reviewer-facing blind rows: **12,064**;
- coordinator-key rows: **12,064**;
- records with hidden automatically detected binomial hint: **9,071**;
- records without detected binomial but retained for review: **2,993**;
- historical exact source records retained: **34/34**;
- blind review batches: **13** (12 × 1,000 + final 64);
- reviewer-facing query membership: **0 columns**;
- reviewer-facing historical status: **0 columns**;
- reviewer-facing automatic taxon hints: **0 columns**.

Additional blind-artifact audit: **1,259 records have no abstract**. Missing abstract triggers conservative/full-text handling, never automatic exclusion.

## Completed: Wave 0 calibration gate infrastructure

Successful workflow: `JBI v2.2 record-screening Wave 0 calibration`, run **`32835886607`**.

All CI steps passed: canonical input verification, 384-record construction, blinding/stratum validation, blank agreement-scorer smoke test, and artifact upload.

Wave 0 contains **384 unique records** with hidden mutually exclusive strata:

- historical benchmark sources: **34**;
- no detected binomial: **100**;
- detected binomial: **100**;
- non-English: **100**;
- missing abstract: **50**.

Reviewer order is independent of stratum. Reviewer-facing files contain no stratum, query membership, historical status, old spatial label, or automatic taxon hint.

`docs/JBI_V22_RECORD_SCREENING_CODEBOOK.md` now explicitly requires separate reviewer copies: Reviewer 1 fills only `reviewer_1_*`, Reviewer 2 fills only `reviewer_2_*`, and adjudication fields stay blank until both independent copies are locked.

`score_v22_record_calibration.py` computes raw agreement and Cohen's kappa for record relevance, natural intraspecific variation, floral-display colour, and full-text requirement, plus normalized exact agreement for focal taxon text.

## Legacy-priority diagnostics — sensitivity only

These diagnostics quantify screening-path bias but do not define canonical inclusion.

### P1

- records: **543**;
- GBIF-valid species: **209**;
- automated states: 109 among-only / **41 mixed** / 26 within-only / 33 unresolved;
- blind source rows: **383**;
- historical species recovered through P1 candidate path: **22/34**;
- exact historical source papers outside P1: **15/34**.

### P2/P3

- records: **1,996**;
- primary taxon names checked: **2,978**;
- GBIF-valid species: **522**;
- blind source rows: **1,196**;
- navigation states: **522 unresolved**.

The unresolved-only slice is expected because the legacy P2/P3 screen did not carry directional within/among signals.

### P1–P3 union

- GBIF-valid species: **656**;
- additional species beyond P1: **447**;
- navigation states: **124 among-only / 41 mixed / 30 within-only / 461 unresolved**;
- blind source rows: **1,625**;
- mixed count added by P2/P3: **0**.

The unchanged mixed count is a property of the old signal definitions, not a biological result.

## Why the old result diverged

The mismatch is not one mislabeled species. It is a stack of upstream filters:

1. different discovery chains for the historical freeze and later systematic map;
2. generic morph terminology accepted as local coexistence;
3. narrower geographic-evidence vocabulary;
4. explicit removal of mixed candidates;
5. priority-dependent candidate admission;
6. automatic taxon extraction capable of missing relevant source papers.

The replacement pipeline removes all six as irreversible biological filters.

## Biological state to be adjudicated later

For each eligible resolved species:

- `local_coexistence_documented`: discrete natural floral-display colour variants documented within at least one same population;
- `geographic_structure_documented`: colour/morph-frequency differentiation documented among geographic units.

Aggregation:

- local only → within-only evidence;
- geographic only → among-only evidence;
- both → mixed/multiscale evidence;
- unresolved → not entered into the primary state analysis until resolved.

A zero means **not documented in the reviewed evidence**, not biological absence.

## Current blocker: external independent duplicate review

The next main step cannot be completed honestly by repeating the same reviewer/AI twice. The required input is now exactly one thing:

**two independent completed copies of the 384-record Wave 0 blind sheet.**

After both are returned:

1. merge only the appropriate reviewer columns by `record_review_id`;
2. calculate raw agreement and Cohen's kappa;
3. inspect disagreements by field and hidden stratum only after independent coding is locked;
4. revise the codebook only if a systematic ambiguity is demonstrated;
5. if rules change materially, run a new blind calibration wave;
6. otherwise begin the 13-batch full 12,064-record screen;
7. only after record adjudication proceed to focal taxon resolution and spatial source review.

The same person or AI agent is **not** counted as two independent reviewers. Coordinator automation may organize records and score agreement but cannot populate independent reviewer/adjudication decisions.

## Manuscript boundary

Until record/source review and adjudication are complete:

- historical 34-species results remain reproducible **historical/sensitivity** results;
- current Main/SI statistical conclusions are not promoted to replacement canonical inference;
- no replacement figures/DOCX should be labeled final;
- mixed-preserving climatic models remain planned inference, not completed biological results.
