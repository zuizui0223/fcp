# Journal of Biogeography editorial submission check

## Verified source state

- Source manuscript reviewed: `docs/jbi_manuscript_draft.md`, blob `91d061710613790091b0489e4035246bbdffd6ee`.
- Primary analysis: workflow run `29972327794`, artifact digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`.
- GBIF and Open Tree sensitivity: workflow run `30067762848`, artifact digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`.
- Fixed-seed dated-megaphylogeny sensitivity: workflow run `30076757379`, artifact `8590190840`, digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`.
- No author detail, funding statement, archive DOI, GBIF DOI, search date, declaration or molecular-phylogeny result has been invented.

## Journal-format check

| Item | Status | Verified finding |
|---|---|---|
| Title | ✓ Pass | 74 characters and non-causal. |
| Running title | ✓ Pass | 35 characters. |
| Structured abstract | ✓ Pass | Uses Aim, Location, Taxon, Methods, Results and Main conclusions and remains below 300 words. |
| Keywords | ✓ Pass | Seven alphabetized keywords. |
| Introduction | ✓ Pass | Defines the spatial-organization question and separates polymorphism from geographic differentiation. |
| Literature discovery | △ Needs revision | Automated queries and filters are documented; original dates, manual screeners and adjudication remain `Not verified`. |
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
| Supporting provenance | ✓ Pass | Analysis tables, placement audit, Open Tree topology, three dated trees and manifests are indexed with SHA-256 values. |
| Data Accessibility | △ Needs revision | Permanent code/data DOI and GBIF DOI remain `Not verified`. |
| Title page and declarations | ✗ Missing author confirmation | Authors, affiliations, ORCIDs, funding, CRediT, conflicts, acknowledgements and biosketch remain `Not verified`. |
| Taxon image | ✗ Missing | Image, rights and caption remain `Not verified`. |

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

## Round-four change history

| Before | After | Reason |
|---|---|---|
| Open Tree topology and Grafen branch lengths only | Added fixed-seed V.PhyloMaker2 S1–S3 models using a time-scaled vascular-plant backbone | Tested whether the result depended on the absence of dated branch lengths. |
| Open Tree model retained 30 species | Dated-megaphylogeny models retained all 34 species | Removed sample loss as the sole explanation for broad phylogenetic intervals. |
| One phylogenetic construction | Compared topology-based and time-scaled treatments | Separated directional consistency from inferential certainty. |
| “A dated species phylogeny is needed” | Replaced with a limitation concerning six inserted taxa and the value of a species-specific molecular phylogeny | Removed an obsolete blocker without overstating the megaphylogeny. |
| Tables S1–S15 | Added Table S16, Table S17, three fixed-seed dated trees and a provenance manifest | Made the dated analysis auditable and reproducible. |

## Editor Check

### Provisional decision: **Major revision before submission**

The manuscript now addresses occurrence-sampling sensitivity and phylogenetic non-independence using two distinct phylogenetic constructions. Directional consistency is stronger than before: the estimated moisture-breadth coefficient remained negative in every family-deletion refit, every Open Tree replicate and every dated-megaphylogeny scenario. Inferential certainty remains limited because all phylogenetic confidence intervals include one, moisture breadth was selected from a 20-specification matrix and the literature-derived response can contain classification error.

The remaining `Major revision` label no longer reflects a missing phylogenetic analysis. It reflects a combination of exploratory focal selection, observational scale mismatch, non-random evidence assembly, unresolved author-controlled declarations and missing permanent data citations.

### Likely reviewer concerns, in priority order

1. **Post-analysis focal selection.** Moisture breadth was selected from five metrics and four thresholds.
2. **Classification error.** Among-population status may reflect under-documentation of local coexistence.
3. **Scale mismatch.** Species-level realised niches cannot identify morph-specific climatic sorting.
4. **Phylogenetic residual uncertainty.** Both phylogenetic treatments preserve the negative direction, but every interval includes one; six taxa were inserted into the dated backbone.
5. **Non-random literature sampling.** Research effort, English queries and metadata availability shape inclusion.
6. **Occurrence-data limits.** The paginated sample remains capped and lacks a citable GBIF download DOI.
7. **Small comparative sample.** The strict comparison contains 34 species.
8. **Manual review documentation.** Search dates, screeners and disagreement resolution remain missing.
9. **Mechanistic non-identifiability.** Occurrences are not labelled by flower-colour morph.

## Remaining submission blockers

- Record the original search dates, manual screeners and adjudication procedure.
- Archive the exact submission release and add its permanent DOI.
- Create and cite a GBIF derived-dataset DOI.
- Verify authors, affiliations, corresponding author, ORCIDs, funding and CRediT roles.
- Obtain conflict-of-interest confirmations.
- Complete Acknowledgements and the biosketch.
- Select and clear the required taxon image.
- Confirm whether the current time-scaled megaphylogeny sensitivity is adequate for submission or whether a bespoke species-level molecular tree is feasible; this is now an editorial-strengthening choice rather than an unaddressed analytical omission.
