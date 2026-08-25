# JBI upstream re-audit status

## Current stage

The historical 34-species freeze remains unchanged. A separate upstream re-audit is now running from the archived systematic-search corpus, and the Stage-1 P1 taxon/source audit has completed successfully.

Successful workflow run: `32822067112` (`JBI upstream spatial evidence re-audit`).

## Established facts

### Search corpus

- archived systematic search: **105,249 raw records** and **79,242 deduplicated records**;
- Stage 1 contains **543 P1 records** (`P1_high_natural_itv` + `P1_high_population_itv`);
- Stage 2 contains **1,996 P2/P3 records** that remain explicitly pending rather than silently excluded;
- the search log contains **19 truncated query/database shards**: **15 Crossref** and **4 OpenAlex**;
- the old manifest did not reliably expose those truncations, so the new gate reads `itv_fcp_search_log.csv` directly;
- final natural eligibility and spatial classification were never human-adjudicated in the archived corpus.

### Stage-1 taxon validation

The P1 corpus contained many permissive two-word candidate strings. Source attribution was therefore separated from taxon validation before any species-level evidence aggregation.

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

### Source-level review workload

- Stage-1 species-by-source review rows: **383**;
- species with any automated within signal: **67**;
- species with any automated among signal: **150**;
- species with both signals in the same source: **30**;
- historical 34-species cases recovered in the Stage-1 P1 candidate queue: **22/34**.

The 12 historical species not recovered in P1 are not considered absent from the literature. They are now routed to a direct historical-source rescue lane using the classification source identifier in the frozen manifest, independent of automated search priority.

## Why the old and new states diverge

The historical 34-species analysis and the later systematic map were built through different evidence paths. The old paper freeze came from the earlier 1,075-work discovery chain, whereas the systematic map was acquired later and contains 79,242 deduplicated records. In addition, the old automated classification logic could infer a within-population state from generic morph terminology and the exploratory systematic-analysis code explicitly removed mixed within/among candidates before ecological modelling.

The new audit therefore does not treat disagreement with the historical label as evidence that either label is automatically correct. It returns to source-level evidence before climatic data are inspected.

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

1. **Priority A — mixed and historical cases**: all 41 automated mixed candidates plus every historical-34 source record;
2. **Priority B — directional candidates with multiple independent sources**: within-only or among-only species supported by more than one source;
3. **Priority C — eligibility conflicts and unresolved cases**: cultivated/induced/ontogenetic conflicts, non-display floral traits, or unresolved spatial evidence;
4. **Stage 2 — P2/P3 expansion**: 1,996 records, processed only after Stage-1 review infrastructure is stable.

Priority is a workload triage variable only. It cannot determine inclusion or spatial state.

## Remaining gates before a replacement climatic analysis

1. complete direct rescue of all 34 historical classification sources regardless of P1/P2/P3 priority;
2. complete duplicate species-by-source review for Priority A, then B/C, with adjudication;
3. repair or explicitly bound the 19 truncated search shards using a narrower systematic-search v2 rather than simply increasing the old broad caps;
4. process Stage-2 P2/P3 records and document what they add;
5. freeze the adjudicated natural/display-colour dataset with separate local and geographic evidence axes;
6. only then calculate the climatic variables and fit the new two-axis / mixed-preserving models.

Until those gates pass, the historical 34-species analysis remains a reproducible historical/sensitivity result, not the final canonical inference.
