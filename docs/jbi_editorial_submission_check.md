# Journal of Biogeography editorial submission check

## Verified source state

- Source manuscript reviewed: `docs/jbi_manuscript_draft.md`, blob `91d061710613790091b0489e4035246bbdffd6ee`.
- Primary analysis: workflow run `29972327794`, artifact digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`.
- GBIF and phylogenetic sensitivity: workflow run `30067762848`, artifact digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`.
- No author detail, funding statement, archive DOI, GBIF DOI, search date or declaration has been invented.

## Journal-format check

| Item | Status | Verified finding |
|---|---|---|
| Title | ✓ Pass | 74 characters and non-causal. |
| Running title | ✓ Pass | 35 characters. |
| Structured abstract | ✓ Pass | Required six headings and fewer than 300 words. |
| Keywords | ✓ Pass | Seven alphabetized keywords. |
| Introduction | ✓ Pass | Defines the spatial-organization question and separates polymorphism from geographic differentiation. |
| Literature discovery | △ Needs revision | Automated queries and filters are documented; original dates, screeners and adjudication remain `Not verified`. |
| Classification | ✓ Pass | Four classes, source audit, freeze rule and strict manifest are documented. |
| GBIF methods | ✓ Pass analytically | Primary and paginated sensitivity workflows are described exactly. A citable GBIF derived-dataset DOI is missing. |
| Climate and metrics | ✓ Pass | WorldClim variables, resolution, PCA and five metrics are reproducible. |
| Statistical methods | ✓ Pass | Formula, covariance, CI, seed, permutations, family deletion and diagnostics are reported. |
| Phylogenetic sensitivity | ✓ Pass with limitation | Two 100-replicate Open Tree topology analyses completed; this is not a dated phylogeny. |
| Results | ✓ Pass | Primary, broader, matched-control, occurrence-sampling and phylogenetic estimates are reported with uncertainty. |
| Discussion | ✓ Pass | Causal overstatement is removed and the phylogenetic uncertainty is explicit. |
| References | ✓ Pass for cited text | Eleven references support the current text. GBIF citation awaits DOI. |
| Tables and figures | ✓ Pass | Four main tables, two figure legends and Tables S1–S15 are cross-referenced. |
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

### Topology-based phylogenetic sensitivity

- 30 species were retained in the Open Tree synthetic topology.
- All 100 polytomy-resolution fits completed for each occurrence dataset.
- Primary occurrence data: OR = 0.592; 95% CI = 0.244–1.434; p = 0.246.
- Paginated data: OR = 0.472; 95% CI = 0.175–1.272; p = 0.138.

The defensible synthesis is:

> The negative association was stable to stronger occurrence sampling but remained statistically unresolved after topology-based phylogenetic correction.

Do not describe the result as confirmed, phylogenetically robust or causal.

## Round-three change history

| Before | After | Reason |
|---|---|---|
| Primary GBIF first-page sample only | Added a paginated, quality-filtered 3,000-record sensitivity analysis | Tested whether the 300-record cap generated the focal direction. |
| Family clustering and family deletion only | Added two topology-based phylogenetic logistic analyses | Directly addressed shared-ancestry sensitivity. |
| “Phylogenetic correction has not been performed” | Replaced with exact Open Tree methods, results and limitations | Removed an obsolete statement after verified analysis completion. |
| GBIF limitation stated generically | Distinguished primary sampling from capped paginated sensitivity | Preserved transparency without ignoring the stronger analysis. |
| Tables S1–S7 | Added Tables S8–S15, QC manifests and induced topology | Made the new analyses fully auditable. |

## Editor Check

### Provisional decision: **Major revision**

The manuscript now addresses the two largest analytical reviewer concerns: occurrence-sampling sensitivity and phylogenetic non-independence. The negative estimate strengthened with improved occurrence sampling, but topology-based phylogenetic confidence intervals included one. The manuscript is therefore substantially stronger, yet the central comparative signal remains phylogenetically unresolved rather than confirmatory.

### Likely reviewer concerns, in priority order

1. **Phylogenetic uncertainty.** The Open Tree analysis is topology-based, uses Grafen branch lengths and retains 30 species; a dated species phylogeny could alter uncertainty.
2. **Post-analysis focal selection.** Moisture breadth was selected from five metrics and four thresholds.
3. **Classification error.** Among-population status may reflect under-documentation of local coexistence.
4. **Scale mismatch.** Species-level realised niches cannot identify morph-specific sorting.
5. **Non-random literature sampling.** Research effort, English queries and metadata availability shape inclusion.
6. **Occurrence-data limits.** The paginated sample remains capped and lacks a GBIF download DOI.
7. **Small comparative sample.** The strict model has 34 species and the phylogenetic model 30.
8. **Manual review documentation.** Search dates, screeners and disagreement resolution remain missing.

## Remaining submission blockers

- Record the original search dates, manual screeners and adjudication procedure.
- Archive the exact submission release and add its permanent DOI.
- Create and cite a GBIF derived-dataset DOI.
- Verify authors, affiliations, corresponding author, ORCIDs, funding and CRediT roles.
- Obtain conflict-of-interest confirmations.
- Complete Acknowledgements and the biosketch.
- Select and clear the required taxon image.
- Decide whether the topology-based model is sufficient or a dated species phylogeny is required before submission.
