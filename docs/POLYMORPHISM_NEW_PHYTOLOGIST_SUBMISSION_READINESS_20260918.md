# New Phytologist submission readiness audit — 2026-09-18

Status: **scientific package ready; administrative metadata still blocks formal submission.**

This audit is for the active flower-colour polymorphism paper:

- `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`
- `docs/POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md`
- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`

Current official author-guideline source checked on 2026-09-18:

- New Phytologist Author Guidelines:
  `https://nph.onlinelibrary.wiley.com/hub/journal/14698137/about/author-guidelines`

## 1. Scientific claim state

**PASS**

The paper's decisive biological claim is frozen as:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

The current manuscript preserves the required boundary:

- species-disjoint prospective confirmation;
- same iNaturalist opportunity/source universe;
- same measurement system;
- not independent-source replication;
- replicated D–spatial organization presented as a structural correlate, not a causal mechanism;
- no global prevalence claim;
- no pigment, transition-direction, pollinator or climate mechanism claim.

Authoritative/result-reporting files:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`
- `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`
- `results/polymorphism_spatial_organization_clue_20260918/result.json` (reporting-only synthesis of previously frozen PR #32 results)

## 2. New Phytologist Full Paper format audit

The official guidelines currently describe Full Papers as usually approximately 6,500–7,500 words with 6–8 display items and a 200-word bulleted Summary. Initial submissions use free format but should include the required research-paper sections, title/authors/addresses/correspondence information, section word counts, and counts of figures/tables/supporting information. Cover letters must answer three editorial questions in no more than 50 words each.

### Current manuscript measurements

| Requirement | Current state | Decision |
|---|---:|---|
| Title approximately <=130 characters | 111 characters | PASS |
| Summary <=200 words | 171 words | PASS |
| Summary structure | 4 bullets | PASS |
| Keywords | 6, alphabetical | PASS |
| Introduction | 541 words | RECORDED |
| Materials and Methods | 2,251 words | RECORDED |
| Results | 1,025 words | RECORDED |
| Discussion | 1,167 words | RECORDED |
| Main text, Introduction–Discussion | 4,984 words | RECORDED |
| Discussion share of main text | 23.4% | PASS (<30%) |
| Main figures | 5 | PASS |
| Main tables | 1 | PASS |
| Total display items | 6 | PASS (guideline range 6–8) |
| Required article sections | present | PASS |
| Figure legends 1–5 | present | PASS |
| Supporting legends S1–S9 | present | PASS |

The main text is shorter than the journal's usual 6,500–7,500-word region. This is not itself a formal failure: the guideline describes a usual range rather than a minimum. Do not add filler solely to approach the range. Any expansion should improve biological framing, methodological reproducibility or interpretation.

## 3. Cover-letter audit

The three required editor questions are present.

Measured answer lengths:

1. What hypotheses or questions does this work address? — **40 words**
2. How does this work advance our current understanding of plant science? — **39 words**
3. Why is this work important and timely? — **40 words**

All satisfy the <=50-word rule.

The cover letter also preserves the same-source/species-disjoint boundary and does not overstate H2 as independent-source replication.

## 4. Figure package

**PASS**

Canonical figure directory:

`docs/figures/polymorphism_20260918/`

Present in both PNG and PDF:

1. Figure 1 — measurement frame;
2. Figure 2 — H1 observer-disjoint reproducibility;
3. Figure 3 — legacy H2 target localization;
4. Figure 4 — prospective H2 confirmation;
5. Figure 5 — replicated spatial clue plus bounded alternative-explanation tests.

Canonical manifest:

`docs/figures/polymorphism_20260918/polymorphism_figure_manifest_20260918.json`

Manifest state:

- `status = generated_from_frozen_results`
- `scientific_claims_changed = false`

Figure 4 contains the decisive prospective result:

- primary: 158 species, W = 0.5172457461, p = 0.001;
- strict: 86 species, W = 0.5329282123, p = 0.001;
- support-evaluable species = 377;
- verdict = `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

## 5. Supporting Information and evidence chain

**PASS**

Supporting evidence map:

`docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`

It records:

- cohort roles;
- H1 chronology;
- legacy H2 target localization;
- third-cohort selection and pre-opening qualification;
- prospective biological run and immutable result;
- replicated D–spatial organization reporting receipt and robustness checks;
- P500 terminal non-verdict;
- H3a/H3b alternative-explanation boundaries;
- canonical figures;
- manuscript claim guard;
- literature-positioning audit.

The most recent post-architecture manuscript claim-guard run recorded in the Supporting Information is:

- run `35299258800`;
- conclusion: success.

## 5b. Current automated guards

**PASS**

Latest manuscript claim guard:

- workflow run: `35314652195`;
- job: `105503609508`;
- result: **8 / 8 tests passed**.

Latest New Phytologist submission-format guard:

- workflow run: `35314606151`;
- job: `105503472346`;
- result: **9 / 9 tests passed**.

The submission-format guard now checks live section word counts, Discussion <=30% of main text, 6–8 total display items, alphabetical keywords, four-bullet Summary <=200 words, cover-letter question lengths, required sections and frozen H2 claim boundaries.

## 6. Literature layer

**PASS for submission drafting**

Canonical audit:

`docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`

The literature layer is restricted to background and interpretation. It does not modify the frozen empirical claim ceiling.

The journal-specific manuscript currently cites direct precedents for:

- intraspecific trait variation;
- flower-colour polymorphism;
- image-derived colour measurement;
- iNaturalist observer structure;
- validation design;
- Monte Carlo p-value resolution.

## 7. Remaining formal-submission blockers

These items require author input and must not be inferred from repository metadata.

### Blocker A — authorship metadata

Current placeholders:

- `[AUTHOR LIST TO CONFIRM]`
- `[AFFILIATIONS TO INSERT]`
- `[NAME / EMAIL TO INSERT]`

Required before submission:

- final author order;
- full institutional affiliations;
- corresponding author name and email;
- corresponding-author ORCID.

### Blocker B — acknowledgements/funding

Current placeholder:

`[TO COMPLETE BEFORE SUBMISSION: funding, institutional support, data-provider acknowledgements, and individual contributions that do not meet authorship criteria.]`

Required:

- funding sources and grant numbers where applicable;
- institutional or technical support;
- data/provider acknowledgements;
- non-author contributions.

### Blocker C — competing interests

Current placeholder:

`[TO CONFIRM BEFORE SUBMISSION.]`

A final declaration is required even if the answer is no competing interests.

### Blocker D — author contributions

Current placeholder:

`[TO COMPLETE AFTER FINAL AUTHOR LIST.]`

This must be completed only after final authorship is known.

### Blocker E — permanent archive

Current Data availability text states that a permanent archival DOI/version should be added before submission.

A Zenodo, institutional repository or equivalent permanent version should be frozen and cited in the final Data availability statement.

### Blocker F — final formatted submission file

The repository source is Markdown. The formal initial-submission manuscript should be exported to a review-ready word-processing/PDF format with:

- 1.5-line spacing;
- page numbers;
- continuous line numbering;
- a consistent font;
- the completed title-page author/correspondence metadata.

This is a packaging requirement and does not change the frozen scientific content.

## 8. Current submission decision

**Do not run new biological analyses to improve this submission.**

The scientific package already contains the highest-value prospective H2 upgrade identified in the earlier architecture.

The remaining route to a formal New Phytologist submission is:

`author metadata -> acknowledgements/conflicts/contributions -> permanent archive DOI -> final formatted submission file (1.5 spacing + page/continuous line numbering)`

No H2 target, threshold, null, cohort or H3 predictor should be reopened during this packaging stage.
