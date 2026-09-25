# New Phytologist submission readiness audit — 2026-09-25

Status: **scientific, narrative and provenance package ready; administrative metadata, permanent external archive DOI, and final formatted submission file remain.**

Active submission files:

- `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`
- `docs/POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md`
- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`
- `docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`
- `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`

Official guideline source rechecked on 2026-09-25:

- New Phytologist Author Guidelines:
  `https://nph.onlinelibrary.wiley.com/hub/journal/14698137/about/author-guidelines`

## 1. Scientific claim state

**PASS**

The decisive biological result remains:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

The evidential sequence is now explicit:

1. **Measurement validity** — species-level four-state diversity D is reproducible under the first-frozen observer-disjoint rule and transports strongly across a later fresh image set within the same measurement system.
2. **Target discovery/localization** — the original discovery/reserve cohorts localize recurrent continuous colour geometry to a white-versus-nonwhite axis, but that named axis is retrospective in those cohorts.
3. **Untouched prospective confirmation** — a separately selected species-disjoint H2 cohort tests the already frozen q_white/W target with unchanged gates and structured null.
4. **Post-H2 validity/environment follow-up** — only after prospective H2 terminalization are the same physical third-cohort rows reused for highlight validity and the pre-specified BIO5/BIO14/solar analysis.
5. **Spatial organization** — greater D is associated with stronger within-species geographic colour organization; current RGFCA identifiability does not decide whether cross-species spatial realization is shared, partly shared or species-specific.
6. **Alternative-explanation filters** — broad reserve phylogenetic signal is unsupported and the discovery sampled-span association fails species-disjoint replication.

### Required claim boundaries

The manuscript explicitly preserves:

- species-disjoint prospective H2 confirmation, but not independent-source replication;
- no global polymorphism-prevalence claim;
- no evolutionary transition-direction claim;
- no universal pigment mechanism;
- no pollinator causation;
- no universal climate causation;
- no claim that cross-species geographic maps are species-specific;
- measured coarse white remains exposure-coupled rather than artifact-cleared.

### Secondary BIO5 result

The secondary environmental result is reported rather than hidden:

- third-cohort eligible species = **281**;
- median white-minus-nonwhite BIO5 contrast = **+0.0690 SD**;
- Holm-adjusted p = **0.0354**;
- conditional OR = **1.073**, p = **0.000919**.

But the later fixed species-disjoint transport rule fails:

- discovery species-level p = **0.743**;
- reserve species-level p = **0.0541**;
- verdict = `LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST`.

Temperature is therefore discussed as a context-dependent mechanistic clue, not a universal driver.

## 2. Narrative clarity

**PASS**

The manuscript now distinguishes physical datasets from inferential uses.

### Legacy resource

Discovery and reserve are species-disjoint halves of the original high-depth RGFCA measurement resource and support H1, legacy H2, spatial organization and H3 under different frozen roles.

### Prospective H2 cohort

One physical dataset:

- 499 species;
- 49,900 terminal rows;
- 377 measurement-evaluable species.

Two chronological uses:

1. untouched prospective H2 confirmation;
2. only later, post-H2 secondary highlight/environment analyses.

The manuscript Table 1 and Supporting Information display these as separate inferential executions, preventing BIO5 from being mistaken for part of the untouched prospective H2 design.

## 3. Full Paper format audit

The current New Phytologist guidance describes Full Papers as usually approximately 6,500–7,500 words with 6–8 display items and a 200-word bulleted Summary. Initial submissions require the standard research-paper sections, 1.5-line spacing, page and continuous line numbering, title/author/correspondence metadata, section word counts, and 5–8 keywords. Cover letters answer three editorial questions in no more than 50 words each.

### Current manuscript measurements

| Requirement | Current state | Decision |
|---|---:|---|
| Title approximately <=130 characters | 114 characters | PASS |
| Summary <=200 words | 180 words | PASS |
| Summary structure | 4 bullets | PASS |
| Keywords | 6, alphabetical | PASS |
| Introduction | 541 words | RECORDED |
| Materials and Methods | 2,858 words | RECORDED |
| Results | 1,807 words | RECORDED |
| Discussion | 2,087 words | RECORDED |
| Main text, Introduction–Discussion | **7,293 words** | PASS |
| Discussion share | **28.6%** | PASS (<30%) |
| Main figures | 5 | PASS |
| Main tables | 1 | PASS |
| Total display items | 6 | PASS (usual 6–8) |
| Required sections | present | PASS |
| Figure legends 1–5 | present | PASS |
| Supporting legends S1–S9 | present | PASS |

The previous 7,712-word version was trimmed by moving acquisition/firewall implementation detail to Supporting Information and the data-lineage map. Statistical decision rules and numerical results were not removed.

## 4. Cover-letter audit

**PASS**

Current answer lengths:

1. What hypotheses or questions does this work address? — **40 words**
2. How does this work advance our current understanding of plant science? — **49 words**
3. Why is this work important and timely? — **48 words**

All satisfy the <=50-word rule.

## 5. Data lineage and reproducibility

**PASS, with one external-byte archival action remaining**

Reader-facing routing map:

- `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`

Machine-readable audit:

- `results/polymorphism_data_lineage_audit_20260925/result.json`

Permanent Git archive:

- `archive/fcp_submission_20260925/`

### Current traceability grades

| Component | Grade | Reason |
|---|---|---|
| H1 | A | protocol/result plus immutable legacy input hashes |
| Legacy H2 | A | protocol/result plus immutable legacy inputs |
| Prospective H2 | A | selection, protocol, result, immutable measured table/hash and terminal artifact |
| Highlight validity | A | exact former Actions-only technical inputs now copied into Git |
| H3a | A | exact trees, pre-outcome covariate panel, signal outputs and permutation nulls copied into Git |
| H3b | A | exact pre-outcome inputs and full permutation/PGLS outputs copied into Git |
| Spatial organization | A- | reporting receipt points to immutable source commit/blobs |
| Secondary BIO5 | A- | biological/technical inputs fixed; WorldClim bytes checksum-pinned but not Git-mirrored |
| BIO5 transport | A- | legacy inputs fixed; WorldClim bytes checksum-pinned but not Git-mirrored |

### Expiring Actions artifacts

**Resolved.**

The exact files formerly dependent on Actions retention are now committed byte-for-byte under:

- `archive/fcp_submission_20260925/highlight/`
- `archive/fcp_submission_20260925/h3a_tree/`
- `archive/fcp_submission_20260925/h3a_covariate/`
- `archive/fcp_submission_20260925/h3a_signal/`
- `archive/fcp_submission_20260925/h3b/`

Source run/artifact IDs and per-file SHA256 values are recorded in:

- `archive/fcp_submission_20260925/ARCHIVE_MANIFEST.md`

### WorldClim byte identity

A dedicated checksum run froze the current exact WorldClim 2.1 10-arc-minute archives and analysis rasters:

- BIO zip SHA256: `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`
- SRAD zip SHA256: `c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026`

BIO5, BIO14 and all 12 SRAD TIFF hashes are stored in:

- `archive/fcp_submission_20260925/worldclim_checksums.txt`

The environmental workflows now hard-fail if downloaded bytes do not match the frozen checksums.

### Remaining provenance action

For provider-independent bit-for-bit replay, mirror the two checksum-pinned WorldClim zip archives themselves in the final Zenodo/institutional release. This is the only material external-byte gap remaining in the headline chain.

## 6. Automated guards

**PASS on the archival-hardening revision merged to main**

The latest verified revision passed:

- New Phytologist submission guard;
- polymorphism manuscript claim guard;
- WorldClim checksum freeze;
- checksum-enforced third-cohort white-environment rerun;
- checksum-enforced legacy BIO5 transport rerun.

The reruns reproduce the same frozen BIO5 positive third-cohort result and the same legacy transport failure.

## 7. Figure package

**PASS**

Canonical figure directory:

`docs/figures/polymorphism_20260918/`

Main display sequence:

1. measurement frame;
2. H1 observer-disjoint reproducibility;
3. legacy H2 target localization;
4. prospective H2 confirmation;
5. spatial organization plus bounded alternative-explanation tests.

Main figures = 5; main tables = 1; total display items = 6.

No new BIO5 main figure is required. BIO5 remains a concise secondary result in text/SI so that the manuscript's visual spine remains measurement -> discovery -> prospective confirmation -> spatial ecology.

## 8. Literature and mechanistic interpretation

**PASS for submission drafting**

The literature layer remains interpretation-only and cannot alter frozen empirical results.

The Discussion may state that the recurrent achromatic–chromatic axis is consistent with many-to-one accessibility of floral pigment networks and that temperature can be one context-dependent input into such networks.

It must continue to state that:

- H2 is sign-invariant and does not identify pigmented -> white evolutionary direction;
- the third-cohort BIO5 association failed cross-cohort transport;
- no universal thermal whitening rule is established;
- the measured white state remains exposure-coupled.

## 9. Remaining formal-submission blockers

These require author/external input rather than new biological analysis.

### A. Authorship metadata

Current placeholders remain for:

- final author order;
- full institutional affiliations;
- corresponding author name/email;
- ORCIDs.

### B. Acknowledgements and funding

Complete:

- funding sources/grant numbers;
- institutional/technical support;
- data/provider acknowledgements;
- non-author contributions.

### C. Competing interests

A final declaration is required.

### D. Author contributions

Complete after final authorship is frozen.

### E. Permanent external archive DOI

Create the final Zenodo/institutional release and include:

- the repository submission state;
- `archive/fcp_submission_20260925/`;
- the two checksum-pinned WorldClim zip archives;
- the final archive manifest/checksums.

Add the DOI/version to Data availability.

### F. Final formatted submission file

Export the frozen manuscript to a review-ready DOCX/PDF with:

- 1.5-line spacing;
- page numbers;
- continuous line numbering;
- consistent font;
- completed title-page metadata.

## 10. Submission decision

**Do not run new biological analyses to improve this submission.**

The scientific and provenance packages are sufficient for submission. The remaining path is:

`author metadata -> acknowledgements/conflicts/contributions -> mirror WorldClim bytes + permanent DOI -> final formatted DOCX/PDF -> submit`

No H2 target, threshold, null, cohort definition, BIO5 predictor family or H3 predictor should be reopened during packaging.
