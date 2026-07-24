# Journal of Biogeography editorial submission check

## Verified source state

- Source manuscript reviewed: `docs/jbi_manuscript_draft.md`, blob `91d061710613790091b0489e4035246bbdffd6ee`.
- Primary analysis: workflow run `29972327794`, artifact digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`.
- GBIF and Open Tree sensitivity: workflow run `30067762848`, artifact digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`.
- Fixed-seed dated-megaphylogeny sensitivity: workflow run `30076757379`, artifact `8590190840`, digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`.
- Automated literature chronology: repository commits document global discovery and recorded follow-up/enrichment activity from 16 to 19 July 2026; the reconstruction is preserved in `docs/jbi_literature_search_provenance.md`.
- Submission-package validation: workflow run `30079894645` passed; artifact `8591349630`, digest `sha256:3158e1f58dc8eaec9aa623fffbcb0574b37411eba6f8c814f8e4f57b85b96335`.
- Taxon-image candidate: the *Ipomoea purpurea* three-colour photograph is documented as CC0 on Wikimedia Commons; source, author, caption and interpretation boundary are recorded in `docs/jbi_taxon_image_candidate.md`.
- No author detail, funding statement, archive DOI, GBIF DOI, human-screener identity, declaration or molecular-phylogeny result has been invented.

## Journal-format check

| Item | Status | Verified finding |
|---|---|---|
| Title | ✓ Pass | 74 characters and non-causal. |
| Running title | ✓ Pass | 35 characters. |
| Structured abstract | ✓ Pass | Uses Aim, Location, Taxon, Methods, Results and Main conclusions; validated length is 244 words. |
| Keywords | ✓ Pass | Seven alphabetized keywords. |
| Main-text length | ✓ Pass | Introduction–Discussion contains 4,424 words in the validated branch state. |
| Introduction | ✓ Pass | Defines the spatial-organization question and separates polymorphism from geographic differentiation. |
| Automated literature chronology | ✓ Pass | Initial global output, deferred follow-up, evidence aggregation and ambiguous-case enrichment are dated from 16–19 July 2026 by repository commits. |
| Manual screening documentation | △ Needs author confirmation | The repository does not identify human screener names or number, independent duplicate screening, agreement statistics or a formal disagreement-resolution procedure. |
| Classification | ✓ Pass | Four classes, source audit, freeze rule and strict manifest are documented. |
| GBIF methods | ✓ Pass analytically | Primary and paginated sensitivity workflows are described exactly. A citable GBIF derived-dataset DOI is missing. |
| Climate and metrics | ✓ Pass | WorldClim variables, resolution, PCA and five metrics are reproducible. |
| Statistical methods | ✓ Pass | Formula, covariance, confidence intervals, seeds, permutations, family deletion and diagnostics are reported. |
| Open Tree sensitivity | ✓ Pass with limitation | Two 100-replicate topology-based analyses completed; the model retained 30 species and used Grafen branch lengths. |
| Dated-megaphylogeny sensitivity | ✓ Pass with limitation | Six fixed-seed V.PhyloMaker2 models completed on time-scaled `GBOTB.extended.LCVP` trees; all retained 34 species, but six species were inserted. |
| Results | ✓ Pass | Primary, broader, matched-control, occurrence-sampling and both phylogenetic treatments are reported with uncertainty. |
| Discussion | ✓ Pass | Causal overstatement is removed; agreement in direction and uncertainty in inference are both explicit. |
| References | ✓ Pass for cited text | Fifteen references support the current text. The GBIF data citation awaits a DOI. |
| Tables and figures | ✓ Pass | Four main tables, two figure legends and Tables S1–S17 are cross-referenced. |
| Supporting provenance | ✓ Pass | Analysis tables, literature chronology, placement audit, Open Tree topology, three dated trees and manifests are indexed or directly referenced. |
| Automated package QA | ✓ Pass | File presence, figure links, required estimates, forbidden claims, S1–S17 row counts and all indexed SHA-256 values passed CI. |
| Data Accessibility | △ Needs revision | Permanent code/data DOI and GBIF DOI remain `Not verified`. |
| Title page and declarations | ✗ Missing author confirmation | Authors, affiliations, ORCIDs, funding, CRediT, conflicts, acknowledgements and biosketch remain `Not verified`. |
| Taxon image | △ Candidate verified | A high-resolution CC0 image and caption are identified; final author approval, download and submission upload remain. |

## Current statistical interpretation

### Primary audited model

- 34 species from 25 families; 20 within-population and 14 among-population cases.
- OR = 0.426; family-clustered 95% CI = 0.184–0.985; Wald p = 0.0460.
- Permutation p = 0.0556; leave-one-family-out OR = 0.317–0.481.

### Stronger occurrence sampling

- All 34 species matched the GBIF backbone exactly and retained at least 20 records.
- 58,455 coordinates were retained; 20,859 climate-linked records remained after climate-vector deduplication.
- OR = 0.300; 95% CI = 0.133–0.675; Wald p = 0.00361; permutation p = 0.0164.
- Leave-one-family-out OR = 0.247–0.359.

### Open Tree topology-based sensitivity

- 30 species were retained in the synthetic topology.
- All 100 polytomy-resolution fits completed for each occurrence dataset.
- Primary data: OR = 0.592; 95% CI = 0.244–1.434; p = 0.246.
- Paginated data: OR = 0.472; 95% CI = 0.175–1.272; p = 0.138.

### Time-scaled megaphylogeny sensitivity

- V.PhyloMaker2 retained all 34 species: 28 were already present in the backbone and six were inserted; no species failed to bind.
- Fixed seed 20260724 produced identical S1–S3 trees for the primary and paginated datasets.
- Primary data across S1–S3: OR = 0.464–0.470; 95% CI envelope = 0.176–1.231; p = 0.121–0.124.
- Paginated data across S1–S3: OR = 0.366–0.369; 95% CI envelope = 0.124–1.081; p = 0.0677–0.0691.
- All six estimates were negative and all confidence intervals included one.

The defensible synthesis is:

> The association was stable to stronger occurrence sampling and remained negative under both topology-based and time-scaled phylogenetic treatments, but all phylogenetic confidence intervals included one.

Do not describe the result as confirmed, phylogenetically robust, independent of ancestry or causal.

## Round-six change history

| Before | After | Reason |
|---|---|---|
| Automated search execution dates marked `Not verified` | Reconstructed global discovery and follow-up/enrichment activity from repository commits dated 16–19 July 2026 | Replaced a broad uncertainty statement with auditable provenance. |
| Search and human review treated as one unresolved block | Separated verified automated chronology from unrecorded human screening arrangements | Prevented commit timestamps from being misrepresented as reviewer metadata. |
| No durable search-history supplement | Added Appendix S1 with commit SHA, timestamp, stage and epistemic boundary | Made the evidence-assembly chronology reviewable. |
| No package-level automated audit | Added CI validation for journal structure, numerical guardrails, reference presence, causal wording, files, row counts and SHA values | Prevents silent drift during final author edits. |
| Taxon image entirely unresolved | Identified a high-resolution *I. purpurea* image released under CC0 and drafted the credit, caption and interpretation boundary | Reduced the image task to final author approval and upload. |
| Title-page word counts stale | Replaced with validated values: Abstract 244 words and Introduction–Discussion 4,424 words | Synchronized submission metadata with the current manuscript. |

## Editor Check

### Provisional decision: **Major revision before submission**

The manuscript now addresses occurrence-sampling sensitivity and phylogenetic non-independence using two distinct phylogenetic constructions, documents the automated evidence-search chronology and passes an automated package audit. Directional consistency is strong: the moisture-breadth estimate remained negative in every family-deletion refit, every Open Tree replicate and every dated-megaphylogeny scenario. Inferential certainty remains limited because all phylogenetic confidence intervals include one, moisture breadth was selected from a 20-specification matrix and the literature-derived response can contain classification error.

The remaining `Major revision` label no longer reflects missing search dates, a missing phylogenetic analysis, an unlocated image or an internally inconsistent submission package. It reflects exploratory focal selection, observational scale mismatch, non-random evidence assembly, incomplete documentation of human review, unresolved author-controlled declarations and missing permanent data citations.

### Likely reviewer concerns, in priority order

1. **Post-analysis focal selection.** Moisture breadth was selected from five metrics and four thresholds.
2. **Classification error.** Among-population status may reflect under-documentation of local coexistence.
3. **Scale mismatch.** Species-level realised niches cannot identify morph-specific climatic sorting.
4. **Phylogenetic residual uncertainty.** Both phylogenetic treatments preserve the negative direction, but every interval includes one; six taxa were inserted into the dated backbone.
5. **Non-random literature sampling.** Research effort, English queries and metadata availability shape inclusion.
6. **Human review documentation.** The repository does not establish whether screening was single-reviewer, duplicated or formally adjudicated.
7. **Occurrence-data limits.** The paginated sample remains capped and lacks a citable GBIF download DOI.
8. **Small comparative sample.** The strict comparison contains 34 species.
9. **Mechanistic non-identifiability.** Occurrences are not labelled by flower-colour morph.

## Remaining submission blockers

- Confirm the people who performed human screening/classification, the number of screeners and how disagreements were handled; state explicitly when duplicate screening or agreement assessment was not performed.
- Archive the exact submission release and add its permanent DOI.
- Create and cite a GBIF derived-dataset DOI.
- Verify authors, affiliations, corresponding author, ORCIDs, funding and CRediT roles.
- Obtain conflict-of-interest confirmations.
- Complete Acknowledgements and the biosketch.
- Approve the CC0 *I. purpurea* image candidate, download the original-resolution file, retain the licence record and upload the final image.
- Confirm whether the current time-scaled megaphylogeny sensitivity is adequate for submission or whether a bespoke species-level molecular tree is feasible; this is an editorial-strengthening choice rather than an unaddressed analytical omission.
