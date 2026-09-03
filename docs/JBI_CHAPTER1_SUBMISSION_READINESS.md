# JBI Chapter 1 — submission readiness audit

Audit date: 2026-09-03  
Target: *Journal of Biogeography*, Research Article  
Scope: submission packaging only. Frozen scientific decisions, thresholds, metrics and denominators are not reopened by this audit.

## Current scientific package

Primary biological result:

- six species / 1,200 photographs;
- 480 calibration + 720 held-out evaluation;
- Stage A primary `p = 0.0113` — supported local within-species organization;
- Stage B primary `p = 0.0906` — universal shared-transition concentration not confirmed.

Separate terminal scale-out:

- 200 species / 60,000 frozen photographs;
- exact 256-partition location-blind measurement/reassembly completed;
- 58 / 200 species measurement-evaluable;
- all eight fixed cohorts `not_evaluable`;
- coordinate join prohibited and never opened;
- no species-conditioned spatial, shared-transition, environmental or pollinator result exists for the terminal scale-out.

The terminal scale-out is reported only as a prospectively stopped measurement-transfer/evaluability extension. It cannot replace or re-estimate the six-species result.

## Journal-format audit

Current Journal of Biogeography author guidance for a Research Article calls for a title of at most 115 characters, a running title shorter than 40 characters, a structured abstract of at most 300 words using Aim / Location / Taxon / Methods / Results / Main conclusions, usually no more than 6,000 words of main text, and roughly up to six main tables/figures. The journal uses double-anonymous review and requires data/code access for reviewers at submission. A cover letter is currently optional.

| Requirement | Current state | Decision / action |
|---|---|---|
| Research Article format | PASS | Manuscript uses Introduction, Materials and Methods, Results and Discussion. |
| Title <=115 characters | PASS | Current title is 107 characters including spaces. |
| Running title <40 characters | PASS | Revised to `Spatial organization of flower colour` (37 characters). |
| Structured abstract <=300 words | PASS | Revised abstract is about 250 words including heading text. |
| Required abstract headings | PASS | Aim, Location, Taxon, Methods, Results, Main conclusions are present. |
| Main biological claim matches frozen results | PASS | Stage A support and Stage B primary non-support are unchanged. |
| Main illustrative material | PASS | Four main result figures (C1–C4); C-S1 and C-S2 remain Supporting Information. |
| Main text <=6,000 words | VERIFY AT EXPORT | Draft is written to Research Article scale; run the final word count on the anonymized submission file after title-page separation and reference cleanup. |
| Main-file order | PARTIAL | Final anonymous export should place References before the Data Accessibility Statement, then tables/figure legends/embedded figures as required by current guidance. |
| Double-anonymous main text | PARTIAL | Author names/affiliations are absent, but the reviewer data/code route must not expose identifying repository ownership. |
| Separate title page | TODO | Add author names, one corresponding author, affiliations, ORCIDs and present addresses outside the anonymized main text. |
| Acknowledgements | TODO AUTHOR INPUT | Add funding/material support and contributor acknowledgements in the submission package while preserving anonymous review. |
| Conflict of Interest | TODO AUTHOR INPUT | Add final author-approved declaration in the required submission metadata/file. |
| Author Contributions | TODO AUTHOR INPUT | Add final CRediT-style or journal-compatible contribution statement. |
| Data availability for peer review | BLOCKING PACKAGING ITEM | Create a reviewer-accessible anonymized/private archive. JBI currently offers private-for-peer-review Dryad deposition at submission. Do not rely on an owner-identifying public GitHub URL in the anonymous manuscript. |
| Permanent public archive at publication | TODO | Archive data and code with a persistent DOI; add the DOI/link to the final Data Availability Statement. |
| Supporting Information separation | TODO PACKAGING | Export C-S1/C-S2 and detailed terminal scale-out receipts/contracts to SI rather than inflating the main paper. |
| Cover letter | OPTIONAL | Current JBI guidance makes the cover letter optional. If supplied, keep it concise and focused on the biogeographic conceptual advance. |
| Taxon image | TODO | Supply a legally reusable image of a focal study taxon as required by current JBI submission guidance. |
| Species authorities on first main-text use | TODO NOMENCLATURE | Add taxonomic authorities after verifying the frozen taxon names against a stable nomenclatural source. This is a formatting/nomenclature step only; it must not change the frozen sample identity. |

## Nomenclature check for final formatting

Current Plants of the World Online records support the following author citations for five of the names directly and give the standard authorship for the sixth sampled name:

- *Antirrhinum majus* L.;
- *Dactylorhiza sambucina* (L.) Soó;
- *Gentiana lutea* L.;
- *Ipomoea purpurea* (L.) Roth;
- *Lysimachia arvensis* (L.) U.Manns & Anderb.;
- *Raphanus sativus* L.

POWO currently treats *Raphanus sativus* L. as a synonym of *Raphanus raphanistrum* subsp. *sativus*. Because the frozen analysis used the sampled taxon identity *Raphanus sativus*, do not silently rename or remap that taxon in the manuscript or data. Resolve nomenclatural presentation explicitly at final copy-edit while preserving the frozen analysis identity.

## Anonymous-review data strategy

The current repository is excellent as the durable development and audit record but is not, by itself, an ideal double-anonymous reviewer endpoint because repository ownership can identify the authors. The submission package should therefore separate:

1. **review archive:** anonymized/private reviewer-accessible data, code, frozen contracts, primary null distributions and figure-generation inputs;
2. **development provenance:** this GitHub repository and PR history, retained unchanged;
3. **publication archive:** DOI-backed public data/code release after acceptance/publication according to journal policy.

The review archive must preserve the same frozen files and hashes. Anonymization is packaging, not a scientific rewrite.

## Main-paper versus Supporting-Information boundary

Keep in the main paper:

- the six-species held-out design;
- operational continuous-colour representation;
- Stage A test and result;
- Stage B geometry/opportunity logic and primary result;
- concise sensitivity statement;
- one short terminal-scale-out Methods paragraph;
- one short terminal-scale-out Results paragraph;
- one limitation paragraph explaining measurement NE.

Move or retain in Supporting Information:

- detailed scale-out provenance and exact artifact hashes;
- cohort-by-cohort terminal gate records;
- all frozen preimage/source-role/qualification contracts;
- complete Stage-B sensitivity tables/null diagnostics;
- reproducibility manifests;
- C-S1 and C-S2.

The terminal scale-out must not acquire a main biological figure because no biological coordinate-colour test was opened.

## Submission claim ceiling

Allowed central claim:

> Continuous flower colour was locally organized within species in the held-out six-species analysis, whereas the stronger hypothesis of one shared global geography of the strongest transitions was not confirmed by the frozen primary test.

Allowed terminal-scale-out claim:

> A prospectively frozen 200-species extension completed exact 60,000-record location-blind measurement but was not evaluable under its predeclared measurement-completeness gate, so no downstream biological spatial test was opened.

Do not claim:

- that the terminal scale-out found no spatial organization;
- that the terminal scale-out refuted shared transitions, climate or pollinators;
- that 58 evaluable species constitute a valid post-result biological subset;
- that the one nominal Stage-B sensitivity replaces the primary `p = 0.0906` result;
- that the six species establish the prevalence of spatial organization across angiosperms;
- that community-photo colour vectors are spectrophotometric or causal pigment measurements.

## Remaining blocking items before a real submission upload

1. reference audit and completion of the bibliography;
2. species-authority/nomenclature formatting without changing frozen taxon identity;
3. anonymized reviewer data/code archive;
4. final Supporting Information assembly;
5. separate title page, acknowledgements, conflict-of-interest and author-contribution statements;
6. taxon image with appropriate reuse rights;
7. final main-file ordering, exported word count and anonymous-file inspection.

An optional cover letter can be prepared after these blockers are cleared. No additional confirmatory atlas inference is a submission blocker. The analysis is scientifically closed; the remaining work is editorial, nomenclatural and packaging work.
