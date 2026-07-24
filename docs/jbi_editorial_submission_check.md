# Journal of Biogeography editorial submission check

## Verified source state

- Repository: `zuizui0223/fcp`
- Default branch reviewed: `main`
- Latest source manuscript reviewed: `docs/jbi_manuscript_draft.md`
- Reviewed manuscript blob SHA: `91d061710613790091b0489e4035246bbdffd6ee`
- Editorial manuscript: `docs/jbi_manuscript_editorial_revision.md`
- Final analysis source: successful `PR climatic niche breadth` workflow run 57, run ID `29972327794`
- Analysis commit: `5c7bd0f22aef54fee9b690e8697e9d8f580f82b2`
- Analysis artifact digest: `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`

No author detail, funding statement, archive DOI, GBIF derived-dataset DOI, search date, reviewer identity or phylogenetic result has been invented. Unresolved items remain marked **Not verified**.

## Journal-format check

Journal standards were checked against the Wiley Journal of Biogeography Author Guidelines accessed on 24 July 2026. The journal uses double-anonymous review; the anonymized main file must contain the title, running title, abstract and keywords, Introduction, Methods, Results, Discussion, References, Data Accessibility Statement, tables, figure legends and embedded figures. The title limit is 115 characters, the running-title limit is fewer than 40 characters, and the structured abstract limit is 300 words.

| Item | Status | Verified finding |
|---|---|---|
| Title | ✓ Pass | 74 characters; non-causal and contains the principal terms. |
| Running title | ✓ Pass | 35 characters. |
| Abstract structure | ✓ Pass | Uses Aim, Location, Taxon, Methods, Results and Main conclusions. |
| Abstract length | ✓ Pass | 280 words excluding headings and keywords. |
| Keywords | ✓ Pass | Seven alphabetized keywords. |
| Introduction | ✓ Pass | Reorganized around the spatial-organization question and supported by nine verified references. |
| Literature-discovery Methods | △ Needs revision | OpenAlex queries, matching rules and automated evidence filters are documented. Original execution dates, manual screeners and disagreement resolution remain Not verified. |
| Classification Methods | ✓ Pass | Four classes, audit rule, freeze rule, strict-set composition and SHA-256 manifest are specified. |
| GBIF Methods | △ Needs revision | Query parameters, first-page cap, coordinate rules, absence of basis filter and 0.001° deduplication are explicit. A citable derived-dataset DOI is missing. |
| Climate and metric Methods | ✓ Pass | WorldClim version, 10 arc-min resolution, BIO variables, deduplication, PCA and all five metric formulas are documented. |
| Statistical Methods | ✓ Pass | Formula, covariance, CI method, minimum sample, missing-row rule, seed, permutation statistic, valid permutation count, software version and diagnostics are reported. |
| Results | ✓ Pass | Focal, broader, matched-control, multiplicity and alternative-control results now include exact estimates and uncertainty. |
| Discussion | ✓ Pass | Causal and confirmatory overstatement removed; classification sensitivity and multiplicity integrated. |
| References | ✓ Pass for cited text | Nine cited references are listed in APA form. A GBIF data citation remains missing because its DOI is Not verified. |
| Tables | ✓ Pass | Three self-contained main tables added. |
| Figure legends | ✓ Pass | Two legends correspond to verified repository SVGs and exact final values. |
| Supporting Information | ✓ Pass for analytical tables | Tables S1–S7 are committed, indexed and SHA-256 checked. |
| Data Accessibility | △ Needs revision | Repository and supporting-file index are identified; permanent archive and GBIF DOIs remain Not verified. |
| Separate title page | △ Needs revision | A structured template is added; author-controlled fields remain Not verified. |
| Acknowledgements | ✗ Missing | Author-controlled information. |
| Funding | ✗ Missing | Author-controlled information. |
| Author contributions | ✗ Missing | Final authors and CRediT roles require approval. |
| Conflict of interest | ✗ Missing | Requires confirmation from every author. |
| Biosketch | ✗ Missing | Requires verified final authorship. |
| Taxon image | ✗ Missing | Image, permission and caption are not identified. |
| Phylogenetic model | ✗ Missing | No matched dated phylogeny or phylogenetic comparative result is present. |

## Statistical interpretation

### Focal result

Preferred wording:

> Geographically structured variation was negatively associated with realised moisture niche breadth (odds ratio = 0.426, family-clustered 95% confidence interval = 0.184–0.985; clustered Wald p = 0.0460), whereas permutation support was borderline (two-sided p = 0.0556).

This wording reports the direction, magnitude, confidence interval and inferential disagreement without reducing the evidence to a binary significance label.

### Multiplicity context

Five metrics were evaluated at thresholds of 10, 20, 30 and 50 occupied climate cells. In the broader evidence set:

- all 20 models converged;
- all metric odds-ratio estimates were below one;
- odds ratios ranged from 0.532 to 0.836;
- all family-clustered 95% confidence intervals included one; and
- clustered Wald p-values ranged from 0.052 to 0.552.

Moisture breadth must therefore remain described as the focal reported association, not a preregistered primary endpoint.

### Evidence-set sensitivity

| Set | Metric | OR | 95% CI | Wald p | Permutation p |
|---|---|---:|---:|---:|---:|
| Baseline-unambiguous | Moisture breadth | 0.426 | 0.184–0.985 | 0.0460 | 0.0556 |
| All classified evidence | Moisture breadth | 0.563 | 0.292–1.085 | 0.0861 | 0.0944 |
| Baseline-unambiguous | PCA hull area | 0.632 | 0.332–1.202 | 0.1620 | 0.3184 |
| All classified evidence | PCA hull area | 0.684 | 0.432–1.082 | 0.1046 | 0.2997 |

The strict-versus-broad contrast supports a claim of classification sensitivity. It does not prove effect heterogeneity or that one evidence set is error-free.

### OR, confidence intervals, permutation and Wald reporting

- Report odds ratios together with their confidence intervals.
- Identify the interval and Wald p-value as family-clustered.
- Report the Wald and permutation results together.
- Do not call p = 0.0556 “marginally significant”.
- Use “borderline permutation support” or report the exact value without an adjective.
- State that family-deletion refits do not provide phylogenetic correction.
- Treat non-significant alternative models as failures to clearly account for the association, not proof of absence.

## Hedge calibration

| Evidence | Preferred wording |
|---|---|
| Direct descriptive output | “contained”, “included”, “ranged”, “converged” |
| Modelled relationship | “was associated with” |
| Stable direction with borderline uncertainty | “suggestive”, “evidence-sensitive”, “the estimated direction remained negative” |
| Mechanistic possibility | “may reflect”, “is consistent with”, “remains compatible with” |
| Unsupported mechanism | Do not state as an explanation; identify it as untested |
| Null or imprecise secondary model | “no clear evidence”, “did not clearly account for” |
| Causal inference | Prohibited for the present analysis |

Avoid `demonstrates`, `confirms`, `climate drives`, `local adaptation explains`, `robust evidence` and `no effect`.

## Round-two change history

| Before | After | Reason |
|---|---|---|
| Literature claims followed by `[CITATIONS REQUIRED]` or `Not verified` | Added verified citations to Rausher, Wessinger & Rausher, Trunschke et al., Narbona et al., Koski & Ashman and Dalrymple et al. | Established the biological basis, formal FCP definition and climate-comparison context. |
| “Literature-derived pipeline” without reproducible detail | Added all eight OpenAlex queries, pagination defaults, binomial matching, context windows, scoring and negative-title filters | Made initial automated discovery reproducible without claiming unverified manual procedures. |
| Language bias not discussed | Stated that no explicit language filter was coded but English queries and metadata availability may bias coverage | Converted an implicit limitation into an explicit sampling boundary. |
| Validation effort stated qualitatively | Added 664 candidates, 111 validated cases, ΔAIC = 10.46 and the quadratic effort-model estimates | Quantified ascertainment without interpreting it as prevalence. |
| GBIF records described as “cleaned” | Replaced with exact API parameters, 300-record first-page cap, coordinate rules, origin exclusion and 0.001° deduplication; stated filters that were not applied | Prevented overstatement of data cleaning and made limitations auditable. |
| No occurrence QC counts | Added 468 requested taxa, 356 with coordinates, 53,366 records, 264 with ≥20 records and zero failed requests | Connected Methods to the final workflow artifact. |
| “Occupied climate cells” ambiguous | Defined them as unique within-species nine-variable raster climate vectors | Aligned terminology with the implemented code rather than implying a regular spatial-grid count. |
| Climate processing summarized | Added 23,145 climate-linked rows, 354 taxa, 262 complete taxa and PCA variance explained | Reported final data flow and diagnostic context. |
| Software and permutation details Not verified | Added Python 3.12, statsmodels 0.14.6, seed 20260719, exact two-sided permutation formula, 9,999 valid permutations and missing-row rules | Removed a major reproducibility blocker. |
| Full multiplicity matrix absent | Added Table S1 with all 20 models and summarized the complete OR, CI and p-value ranges in Results | Prevented selective presentation of the focal metric. |
| Candidate-versus-control result only qualitative | Added exact same-genus 20-cell results and full Table S2 | Made the secondary comparison reviewable. |
| Negative controls only qualitative | Added fragmentation sample sizes, added-term OR range, turnover result and six-specification range; committed Table S5 | Avoided using non-significance as an unsupported general claim. |
| No main tables | Added classification, focal-model and matched-control tables | Made definitions and estimates accessible without searching prose. |
| Figure legends missing | Added legends for the two verified SVG figures | Aligned figure numbering with actual repository files and final values. |
| No Supporting Information numbering | Assigned Tables S1–S7 and created a source/hash index | Made every cross-reference verifiable. |
| References missing | Added a checked APA-style reference list for all current in-text citations | Removed citation placeholders while avoiding invented sources. |
| Declarations embedded in anonymized manuscript | Moved author-controlled material to a separate title-page template | Improved compliance with double-anonymous file separation. |
| Broad limitation “GBIF cleaning incomplete” | Specified exactly which checks were and were not performed | Gives reviewers a transparent basis for judging occurrence-data quality. |

## Editor Check

### Provisional editorial decision: **Major revision**

The manuscript now has a coherent Journal of Biogeography narrative, a reproducible analytical description, complete focal statistics, explicit multiplicity context, source-backed figures and tables, and calibrated interpretation. The core question—whether documented local coexistence and geographic differentiation differ in occupied climate—is within the journal’s scope and is potentially publishable.

Major revision remains the appropriate decision because the manuscript has not addressed phylogenetic non-independence, and several submission-critical items depend on author-controlled or externally archived information. A reviewer may also judge the deterministic first-page GBIF sampling and limited coordinate cleaning to be insufficient for a global comparative analysis. Those issues cannot be solved by prose alone.

### Reviewer concerns in priority order

1. **Phylogenetic non-independence.** Family clustering and deletion are sensitivity analyses, not phylogenetic comparative methods.
2. **Post-analysis focal selection.** The rationale for foregrounding moisture breadth must remain transparent beside the 20-model matrix.
3. **GBIF sampling design.** A first-page cap of 300 and limited cleaning could bias occupied-climate breadth.
4. **Classification error.** Among-population classification may reflect failure to document local coexistence.
5. **Scale mismatch.** Species-level realised niches cannot identify morph-specific climatic sorting.
6. **Non-random literature sample.** Research effort, English queries, taxonomic interest and metadata availability affect inclusion.
7. **Borderline method-dependent support.** Wald and permutation results differ at the conventional threshold.
8. **Small-sample clustered inference.** The focal model has 34 species and 25 family clusters.
9. **Absence of a prospective analysis record.** Moisture breadth is exploratory.
10. **Controls are not verified monomorphic.** The matched analysis addresses candidature, not absence of colour variation.
11. **Coarse alternatives are observational.** Point-cloud components cannot represent populations or morph distributions.
12. **Manual evidence review is incompletely documented.** Search dates, screeners and disagreement procedures remain missing.

## Remaining submission blockers

- Obtain or construct a species-level dated phylogeny and run a phylogenetically appropriate sensitivity model, or provide a strong editorial justification for its absence.
- Reassess whether the GBIF first-page sample and limited cleaning are adequate; rerun with a citable GBIF download and stronger quality control if feasible.
- Record the original search dates, manual screeners and adjudication procedure.
- Archive the exact submission release and add its persistent DOI.
- Create and cite the GBIF derived dataset DOI.
- Verify the final authors, affiliations, corresponding author, ORCIDs, funding and CRediT roles.
- Obtain conflict-of-interest confirmations.
- Complete Acknowledgements and the biosketch.
- Select and clear the required taxon image.
- Confirm that all final files and captions match the submission-system upload categories.
