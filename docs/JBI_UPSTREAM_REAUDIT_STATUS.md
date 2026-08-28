# JBI upstream re-audit status

## Current stage

The evidence-side rebuild is complete through a **versioned coexistence/segregation freeze**. The historical 34-species binary freeze remains unchanged and is retained as a provenance/sensitivity dataset.

Independent duplicate review was intentionally waived on **2026-08-26**. The replacement evidence synthesis must therefore be described as a **single-pass conservative documented-evidence audit**, not as an independently double-reviewed systematic review.

Current canonical path:

**OpenAlex v2.2 retrieval → all 12,064 records retained → strict positive-evidence audit → GBIF taxon resolution → C/S species aggregation → versioned C/S evidence freeze → fresh GBIF occurrence reconstruction → WorldClim 2.1 climate metrics → C/S models.**

The current active task is the fresh occurrence/climate reconstruction.

## Retrieval boundary — complete

OpenAlex title/abstract OQL v2.2:

- query blocks: **15**;
- query memberships: **13,911**;
- deduplicated works: **12,064**;
- truncated v2.2 blocks: **0**;
- historical exact-source recovery: **34/34**;
- benchmark recovery by exact/version/direct-citation diagnostics: **34/34**.

All 12,064 records remain visible in the upstream audit. Automatic taxon extraction cannot silently remove a bibliographic record from the retrieval universe.

## Biological representation — frozen

The old `within_population` / `among_population` binary response has been replaced for the new analysis by two independent positive-evidence axes.

### C — local coexistence

`C_local_coexistence_documented = 1` means that discrete natural floral-display colour variants are explicitly documented within the **same population/site**.

Accepted positive evidence includes direct coexist/co-occur statements, within-population flower-colour polymorphism, mixed-colour/polymorphic populations, or population-level morph frequencies establishing more than one morph with positive frequency.

Generic `colour morph`, species-level `polymorphic`, multiple localities, or silence cannot create C=1.

### S — spatial segregation

`S_spatial_segregation_documented = 1` means that floral-colour morph occurrence/frequency is explicitly spatially structured among populations, sites, regions, islands, elevation bands, or other geographic units.

Accepted positive evidence includes explicit geographic/spatial differentiation, morph-frequency contrasts among populations, restriction/replacement of morphs, monomorphic-versus-polymorphic population contrasts, and explicit altitudinal/geographic changes in polymorphism.

Multiple sampling localities or a broad range without a reported colour contrast cannot create S=1.

### Zero semantics

C=0 and S=0 mean **not documented by the strict pass**, never biological absence. S is not defined as the negation of C.

## C/S evidence freeze — complete

Frozen input:

`data/frozen/frozen_34species_coexistence_segregation_v22.csv`

Provenance manifest:

`data/frozen/frozen_34species_coexistence_segregation_v22_manifest.json`

Positive-source index:

`docs/supporting/jbi_cs_positive_source_index_v22.csv`

Final v5 evidence workflow:

- workflow run: **32926071541**;
- artifact: **9591546337**;
- artifact SHA-256: `0a2b7c7f20e6093daedb6b4b317fd8e7b8a89b818c2df042866ddb6bcf660951`;
- source record-universe run: **32835034885**;
- source record-universe artifact SHA-256: `d2c2ccd196043ddff41c7a4794453ac5a311db49d34ca8e94a2b3f334e150dbb`.

Frozen species states:

- **local coexistence only (C=1,S=0): 11**;
- **spatial segregation only (C=0,S=1): 8**;
- **coexistence + segregation (C=1,S=1): 15**;
- strict primary evidence freeze: **34 species**.

Evidence-source diagnostics:

- strict positive bibliographic records: **50**;
- positive source records contributing to the freeze after taxon resolution: **49**;
- C-positive source records: **37**;
- S-positive source records: **29**;
- positive source records left taxonomically unresolved: **1**.

The final pass explicitly excludes colour polymorphism restricted to non-display sexual organs such as anther/pollen/androecium colour. `Epimedium pubescens` therefore does not enter the floral-display freeze.

## Historical comparison

Only **10 species** overlap between the old binary 34-species freeze and the new C/S 34-species freeze; **24 species differ**. Therefore old climatic metrics cannot be treated as the new dataset and old odds ratios cannot be interpreted as estimates for the C/S question.

The old frozen dataset and its checksums remain untouched.

## Climate reconstruction — running

The old climate-method history confirms that occupied climate cells were defined at **WorldClim 2.1, 10 arc-minute resolution**. Nine BIO variables are retained:

`BIO1, BIO4, BIO5, BIO6, BIO7, BIO12, BIO14, BIO15, BIO17`.

Occurrence reconstruction deliberately matches the previous paginated-GBIF sensitivity settings:

- maximum retrieved records per species: **3,000**;
- page size: **300**;
- `hasCoordinate=true`;
- `hasGeospatialIssue=false`;
- `occurrenceStatus=present`;
- accepted basis of record: preserved specimen, material sample, human observation, machine observation, observation;
- maximum coordinate uncertainty: **20,000 m**;
- coordinate deduplication: latitude/longitude rounded to **3 decimals**.

Active workflow:

`JBI C/S occurrence and climate rebuild`

Current run:

**32926741272**

The workflow will not force the evidence freeze to remain n=34. Species with fewer than **20 occupied 10-arc-minute climate cells** remain in the evidence freeze but are excluded only from the climatic model matrix and reported explicitly.

Five climatic summaries are recomputed from scratch:

1. temperature breadth;
2. moisture breadth;
3. climatic heterogeneity;
4. PCA dispersion;
5. PCA hull area.

These are realised occupied-climate summaries, not physiological tolerances or morph-specific niches.

## Planned primary inference — implemented, awaiting climate artifact

New model entry point:

`scripts/run_cs_models.py`

For each of the five climate metrics, two symmetric primary models are fitted:

- `C ~ metric_z + occurrence_effort_z`;
- `S ~ metric_z + occurrence_effort_z`.

A C+S species contributes a positive case to **both** outcomes.

Primary uncertainty/robustness already implemented:

- family-clustered sandwich covariance;
- 9,999 outcome-label permutations;
- leave-one-family-out refits;
- Holm correction across the five metrics separately within C and S;
- sensitivity adding `z(log1p(n_resolved_sources))` as a documentation/source-effort covariate.

The three-state `local-only / spatial-only / coexistence+segregation` comparison is secondary and will be attempted only after the climate-eligible state counts are known.

## Manuscript boundary

No replacement manuscript conclusion, figure, or DOCX is final yet. The next valid inferential checkpoint is:

1. climate workflow completes;
2. >=20-cell analysis sample and state counts are fixed;
3. C and S five-metric models run;
4. effect directions, finite-sample stability, and documentation-effort sensitivity are inspected;
5. only then are figures/Main/SI rewritten.

The historical binary result remains available as a sensitivity/provenance analysis and is not silently rewritten as the new result.
