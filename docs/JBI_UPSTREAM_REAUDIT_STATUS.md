# JBI upstream re-audit status

## Current stage

The historical 34-species freeze remains unchanged. A separate upstream re-audit is running from the archived systematic-search corpus. Stage-1 P1 taxon validation, mixed-preserving evidence aggregation, truly blinded source-material generation, historical-source rescue, and the replacement search-surface audit have all completed successfully.

Latest successful upstream review workflow: `32825713465` (`JBI upstream spatial evidence re-audit`).

Latest successful replacement-search workflow: `32825713331` (`JBI search v2 scope and retrieval`).

## Established facts

### Legacy archived search corpus

- archived systematic search: **105,249 raw records** and **79,242 deduplicated records**;
- Stage 1 contains **543 P1 records** (`P1_high_natural_itv` + `P1_high_population_itv`);
- Stage 2 contains **1,996 P2/P3 records** that remain explicitly pending rather than silently excluded;
- the legacy search log contains **19 truncated query/database shards**: **15 Crossref** and **4 OpenAlex**;
- those 19 legacy truncations are retained as provenance but are no longer the active completeness boundary because the search surface has been replaced by v2.2 below;
- final natural eligibility and spatial classification were never human-adjudicated in the archived corpus.

### Replacement search v2.2

The broad legacy search was not repaired by increasing caps. The search surface itself was corrected:

1. OpenAlex is the completeness-defining retrieval source;
2. retrieval is restricted to **title/abstract OQL**, rather than ordinary search that can include full text;
3. Crossref is retained for bibliographic metadata / DOI resolution, not as a completeness-defining high-recall source;
4. non-English exact phrases were replaced by bounded language-specific organ × colour × variation/population blocks after scope probes;
5. generic lexical omissions found during benchmark inspection were repaired without adding any focal species name: singular/plural `color/colour morph(s)` and `geography/geographic structure` variants were added.

Successful v2.2 retrieval:

- query blocks: **15**;
- raw query memberships: **13,911**;
- deduplicated works: **12,064**;
- duplicate memberships removed: **1,847**;
- truncated v2.2 query blocks: **0**;
- historical 34 exact classification sources recovered directly by the v2.2 queries: **34/34**;
- historical 34 recovered by exact identifier or conservative same-title/version matching: **34/34**;
- prespecified seven review seeds resolved: **7/7**;
- historical benchmark sources recovered by direct v2.2 query or one-generation prespecified citation chasing: **34/34**.

The old 19 truncated shards are therefore treated as a **resolved legacy-search defect**, not a blocker for the replacement search. Search-query membership remains retrieval evidence only and never determines biological inclusion or spatial state.

### Stage-1 taxon validation

The P1 corpus contained many permissive two-word candidate strings. Source attribution is therefore separated from taxon validation before any species-level evidence aggregation.

- contextual candidate strings in P1 records: **4,613**;
- primary source candidates (name in the title or first 2,000 abstract characters): **874**;
- contextual strings deferred as non-primary rather than treated as taxa: **3,739**;
- GBIF queries sent: **874**;
- strict GBIF-valid candidate species: **209**;
- GBIF resolution errors: **0**.

These 209 rows are **taxonomically valid review candidates**, not adjudicated biological inclusions.

### Automated navigation states after taxon validation

The 209 Stage-1 candidate species currently have the following automated evidence states:

- `among_evidence_only`: **109**;
- `mixed_evidence`: **41**;
- `within_evidence_only`: **26**;
- `unresolved`: **33**.

These are navigation states only. `mixed_evidence` means that the archived record set contains at least one automated within signal and at least one automated among signal after taxon validation and source aggregation. It does **not** mean that human review has already established a true multiscale biological state.

### Source-level review workload and blinding

- Stage-1 species-by-source review rows: **383**;
- species with any automated within signal: **67**;
- species with any automated among signal: **150**;
- species with both signals in the same source: **30**;
- historical 34-species cases recovered in the Stage-1 P1 candidate queue: **22/34**.

The source-review materials are split into two files:

1. `systematic_source_review_blind.csv` — reviewer-facing bibliographic information, evidence excerpt, and blank reviewer/adjudication fields only;
2. `systematic_source_review_coordinator_key.csv` — automated within/among/natural/cultivated/induced signals, screen priority, and source-link metadata.

The workflow asserts that the reviewer-facing sheet contains **zero automated-signal columns**. The coordinator key must not be supplied to independent reviewers before coding is complete.

### Historical 34-species direct source rescue

The historical sample is no longer dependent on P1 recovery. Every frozen source identifier is looked up directly against the full archived systematic corpus, independently of search priority.

- historical species retained for rescue review: **34/34**;
- exact historical classification sources recovered in the archived systematic corpus: **34/34**;
- exact historical sources that fall outside Stage-1 P1: **15**;
- historical sources missing from the archived systematic corpus: **0**.

This establishes that the earlier 22/34 P1 recovery was a **priority/screening-path issue, not literature absence**.

Historical rescue materials are also split:

1. `historical_34_source_review_blind.csv` — reviewer-facing source material with no old spatial label;
2. `historical_34_source_review_coordinator_key.csv` — old frozen label, source recovery state, archived screen priority, and automated signals.

The workflow asserts that the reviewer-facing historical sheet contains **zero historical-label columns**.

## Why the old and new states diverge

The historical 34-species analysis and the later systematic map were built through different evidence paths. The old paper freeze came from the earlier 1,075-work discovery chain, whereas the systematic map was acquired later and contains 79,242 deduplicated records. In addition, the old automated classification logic could infer a within-population state from generic morph terminology and the exploratory systematic-analysis code explicitly removed mixed within/among candidates before ecological modelling.

The direct source rescue adds a stronger diagnosis: **15/34 historical source papers are outside P1 despite being the exact papers used for the frozen labels**. Therefore P1 recall cannot be used as a biological inclusion criterion. Search priority is now only review triage.

The v2.2 benchmark additionally shows that several apparent search misses were lexical/indexing artifacts rather than biological evidence gaps. Generic vocabulary repair restored **34/34 exact historical source recovery** without species-specific query terms.

The new audit does not treat disagreement with the historical label as evidence that either label is automatically correct. It returns to source-level evidence before climatic data are inspected.

## Spatial representation retained for adjudication

The re-audit stores two evidence axes separately:

1. `local_coexistence_documented` — direct source evidence that discrete natural floral-colour variants coexist in at least one population;
2. `geographic_structure_documented` — direct source evidence for geographic/among-population differentiation in colour or morph frequency.

After source adjudication:

- local = 1, geographic = 0 → within-only evidence;
- local = 0, geographic = 1 → among-only evidence;
- local = 1, geographic = 1 → mixed/multiscale evidence;
- unresolved evidence remains outside the primary inferential freeze.

A zero is interpreted as **not documented in the reviewed evidence**, not as proof of biological absence.

## Review priority now fixed

### Priority A — mixed and historical cases

- all **41 automated mixed candidates**;
- all **34 historical species** via direct historical-source rescue;
- overlap between these sets is retained rather than collapsed at the species level because the review unit is species × source;
- coordinator batch plan currently contains **195 source-review units**.

The automated mixed list and historical labels are used only by the coordinator to assemble the batch; reviewers receive the blind source sheets without the priority reason.

### Priority B — directional candidates with multiple independent sources

Within-only or among-only species supported by more than one source.

### Priority C — eligibility conflicts and unresolved cases

Cultivated/induced/ontogenetic conflicts, non-display floral traits, or unresolved spatial evidence.

### Stage 2 — P2/P3 expansion

**1,996 records**, processed after Stage-1 review infrastructure is stable. P2/P3 records remain pending, not excluded.

Priority is a workload triage variable only. It cannot determine inclusion or spatial state.

## Remaining gates before a replacement climatic analysis

1. complete duplicate blind source-level review and adjudication for Priority A, then B/C;
2. process Stage-2 P2/P3 records through the corrected taxon/source-review path and document what they add;
3. freeze the adjudicated natural/display-colour dataset with separate local and geographic evidence axes;
4. only then calculate the climatic variables and fit the new two-axis / mixed-preserving models.

The legacy 19-shard truncation problem is no longer a remaining gate: it has been replaced by the v2.2 title/abstract search with zero truncation and 34/34 benchmark recovery.

Until the review/adjudication gates pass, the historical 34-species analysis remains a reproducible historical/sensitivity result, not the final canonical inference.
