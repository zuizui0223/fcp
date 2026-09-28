# New Phytologist submission readiness audit — 2026-09-25

Status: **scientific, narrative and reproducibility/provenance package ready; administrative author metadata and the final formatted submission file remain. A DOI-bearing mirror may be added for citation but is no longer required to make the numerical evidence chain recoverable.**

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

The current manuscript has one reader-facing inferential sequence:

1. **Outcome-blind opportunity frame** — metadata-only iNaturalist discovery defines 42,111 candidate species before flower colour is examined; 4,730 can support >=100 retained photographs.
2. **Discovery + species-disjoint validation** — 500 + 500 species, 100 photographs per species, provide the high-depth resource for D reliability, target discovery/localization and replicated spatial organization.
3. **Prospective confirmation** — because the white-versus-nonwhite target was identified only after the original geometry was opened, the target and inference machinery were frozen before 499 previously unused species and 49,900 new photographs were measured.
4. **Post-confirmatory annotations** — highlight metrics, WorldClim climate and phylogenetic placements test measurement coupling and simple alternative explanations without upgrading the prospective H2 status.

The decisive biological result remains `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`. The ecological mainline is **measurement validity → recurrent phenotype-space geometry → prospective confirmation → geographic organization**, with phylogeny, sampled span and climate treated as bounded explanatory filters.

### Required claim boundaries

The manuscript explicitly preserves:

- species- and photo-disjoint prospective H2 confirmation, but not independent-source replication;
- no global polymorphism-prevalence claim;
- no evolutionary transition-direction claim;
- no universal pigment, pollinator or climate mechanism;
- no claim that cross-species geographic maps are species-specific;
- measured coarse white remains exposure-coupled rather than artifact-cleared.

The secondary BIO5 association is reported but does not transport under its frozen rule: discovery p = **0.743**, validation p = **0.0541**. The prespecified cross-cohort transport criterion was therefore not met. Temperature is a context-dependent clue, not a universal driver.

## 2. Narrative and data-lineage clarity

**PASS**

The manuscript now explains every data resource with the same four questions: **where it came from, how it was selected, why the step was necessary, and what it can establish**.

The physical flower-colour resources are deliberately simple:

- **100,000 photographs / 1,000 species** for discovery and species-disjoint validation;
- **49,900 new photographs / 499 species** for prospective confirmation.

Fresh-image D transport remeasures new photo IDs under the same system. Highlight, climate and phylogenetic data are annotations or reacquisitions of these resources, not additional independent biological cohorts. Development-only execution history is retained in repository provenance rather than presented as part of the manuscript cohort architecture.

Table 1 is now titled **“Data sources, lineage and inferential necessity”**, and Figure 1 shows the same sequence visually: outcome-blind frame → discovery/validation → target freeze → prospective confirmation → dashed post-H2 annotations.

## 3. Full Paper format audit

The current New Phytologist guidance describes Full Papers as usually approximately 6,500–7,500 words with 6–8 display items and a 200-word bulleted Summary.

### Current manuscript measurements

| Requirement | Current state | Decision |
|---|---:|---|
| Title approximately <=130 characters | 105 characters | PASS |
| Summary <=200 words | **179 words** | PASS |
| Summary structure | 4 bullets | PASS |
| Keywords | 6, alphabetical | PASS |
| Introduction | **725 words** | RECORDED |
| Materials and Methods | **1,894 words** | RECORDED |
| Results | **1,809 words** | RECORDED |
| Discussion | **1,678 words** | RECORDED |
| Main text, Introduction–Discussion | **6,106 words** | RECORDED — concise relative to the usual range |
| Discussion share | **27.5%** | PASS (<30%) |
| Main figures | 5 | PASS |
| Main tables | 1 | PASS |
| Total display items | 6 | PASS (usual 6–8) |
| Required sections | present | PASS |
| Figure legends 1–5 | present | PASS |
| Supporting legends S1–S8 | present | PASS |

The data-provenance rewrite initially expanded the manuscript. Duplicate implementation detail was then removed while preserving, for every main stage, its source, selection logic, inferential necessity and claim boundary. The resulting main text is shorter than the journal's stated usual range, but retains the complete Full Paper structure and all frozen numerical results and decision rules.

## 4. Cover-letter audit

**PASS**

Current answer lengths:

1. What hypotheses or questions does this work address? — **40 words**
2. How does this work advance our current understanding of plant science? — **49 words**
3. Why is this work important and timely? — **48 words**

All satisfy the <=50-word rule.

## 5. Data lineage and reproducibility

**PASS — all headline analysis inputs are permanently recoverable through Git history or checksum-verified release assets**

Reader-facing routing map:

- `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`

Machine-readable audit:

- `results/polymorphism_data_lineage_audit_20260925/result.json`

Permanent Git archive:

- `archive/fcp_submission_20260925/`

### Current traceability grades

| Component | Grade | Reason |
|---|---|---|
| H1 | A | protocol/result plus immutable discovery/validation input hashes |
| Original-cohort H2 target localization | A | protocol/result plus immutable discovery/validation inputs |
| Prospective H2 | A | selection, protocol, result, immutable measured table/hash and terminal artifact |
| Highlight validity | A | exact former Actions-only technical inputs now copied into Git |
| H3a | A | exact trees, pre-outcome covariate panel, signal outputs and permutation nulls copied into Git |
| H3b | A | exact pre-outcome inputs and full permutation/PGLS outputs copied into Git |
| Spatial organization | A- | reporting receipt points to immutable source commit/blobs |
| Secondary BIO5 | A | biological/technical inputs fixed; exact WorldClim BIO/SRAD zip bytes mirrored as checksum-verified release assets |
| BIO5 transport | A | discovery/validation inputs fixed; exact WorldClim BIO zip mirrored as checksum-verified release asset |

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

### WorldClim permanent byte archive

**Resolved.**

The checksum-pinned WorldClim 2.1 10-arc-minute archives are mirrored under GitHub Release tag:

`fcp-worldclim-2.1-10m-20260925`

- BIO asset ID: `588322903`, 49,869,449 bytes, SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`
- SRAD asset ID: `588322905`, 16,233,364 bytes, SHA256 `c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026`

Canonical receipt:

- `archive/fcp_submission_20260925/WORLDCLIM_RELEASE_RECEIPT.md`

The committed SHA256 manifest defines canonical byte identity and the analysis workflows enforce those hashes before fitting.

### Self-contained manuscript provenance snapshot

**Resolved.**

A single end-to-end verification package is maintained under GitHub Release tag:

`fcp-np-provenance-20260926`

- asset: `fcp-np-provenance-20260926.tar.gz`
- receipt: `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`

The Git-tracked receipt records the source commit, asset byte count, SHA256 and packaged-file count for the currently published package. The package contains the active manuscript/SI/figures, frozen protocols/results/code/tests, exact discovery/validation measured tables, the exact prospective-H2 measured table, frozen spatial receipts, fresh-D provenance, permanent highlight/H3 inputs and the checksum-pinned WorldClim BIO/SRAD archives. A per-file SHA256 manifest is included.

The package is refreshed when the active manuscript/evidence surface changes; a refresh repackages frozen evidence and does not alter biological results. A later Zenodo/institutional DOI mirror would provide a citation identifier, not missing analytical evidence.


## 6. Automated guards

**PASS on the archival-hardening revision merged to main**

The latest verified revision passed:

- New Phytologist submission guard;
- polymorphism manuscript claim guard;
- WorldClim checksum freeze;
- checksum-enforced prospective-cohort white-environment rerun;
- checksum-enforced BIO5 transport rerun.

The reruns reproduce the same frozen BIO5 positive prospective-cohort result and the same transport failure.

## 7. Figure package

**PASS**

Canonical figure directory:

`docs/figures/polymorphism_20260918/`

Main display sequence:

1. measurement frame;
2. H1 observer-disjoint reproducibility;
3. original-cohort H2 target localization;
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
- the prospective-cohort BIO5 association failed cross-cohort transport;
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

### E. DOI-bearing mirror for citation

The complete verification bytes are already frozen in GitHub releases and the self-contained provenance snapshot above. A Zenodo/institutional mirror may be registered to provide a DOI for long-term citation. If registered before submission, add that DOI/version to Data availability; it does not require rerunning or changing any analysis.

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

`author metadata -> acknowledgements/conflicts/contributions -> optional DOI-bearing mirror -> final formatted DOCX/PDF -> submit`

No H2 target, threshold, null, cohort definition, BIO5 predictor family or H3 predictor should be reopened during packaging.
