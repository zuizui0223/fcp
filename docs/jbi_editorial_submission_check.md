# Journal of Biogeography editorial submission check

## Verified source state

- Source manuscript reviewed: `docs/jbi_manuscript_draft.md`, blob `91d061710613790091b0489e4035246bbdffd6ee`.
- Primary analysis: workflow run `29972327794`, artifact digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`.
- GBIF and Open Tree sensitivity: workflow run `30067762848`, artifact digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`.
- Fixed-seed dated-megaphylogeny sensitivity: workflow run `30076757379`, artifact `8590190840`, digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`.
- Automated literature chronology: repository commits document global discovery and recorded follow-up/enrichment activity from 16 to 19 July 2026; the reconstruction is preserved in `docs/jbi_literature_search_provenance.md`.
- Integrated submission and exact-GBIF validation: the persistent `PR JBI submission package validation` workflow checks journal structure, numerical guardrails, S1–S17, the exact occurrence bundle and the source-backed author prefill on every relevant change.
- Exact GBIF DOI-preparation bundle: 58,455 records, 34 species, 58,455 unique occurrence keys and 389 parent datasets; exact archive SHA-256 `f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094`.
- Deterministic release preview: the persistent `PR JBI release dry run` workflow records the PR head, package files, SHA-256 values, dirty-tree boundary and expected placeholder paths in its artifact manifest.
- Taxon-image candidate: the *Ipomoea purpurea* three-colour photograph is documented as CC0 on Wikimedia Commons; source, author, caption and interpretation boundary are recorded in `docs/jbi_taxon_image_candidate.md`.
- Author metadata prefill: current author-provided, university-linked and official institutional records support **ZHANG RUIQI**, **Division of Forest and Biomaterials Science**, **Graduate School of Agriculture, Kyoto University**, the general postal address **Kitashirakawa Oiwake-cho, Sakyo-ku, Kyoto 606-8502, Japan**, and institutional email **zhang.ruiqi.77h@st.kyoto-u.ac.jp**. The source record describes current status as research student.
- Final author order, any additional authors, laboratory or room, corresponding-author status and ORCID remain unverified. No ORCID linked to the institutional email was found in the searched mailbox records.
- No funding statement, archive DOI, GBIF DOI, human-screener identity, declaration, coauthor identity, laboratory assignment or molecular-phylogeny result has been invented.

## Journal-format and package check

| Item | Status | Verified finding |
|---|---|---|
| Title | ✓ Pass | 74 characters and non-causal. |
| Running title | ✓ Pass | 35 characters. |
| Structured abstract | ✓ Pass | Uses Aim, Location, Taxon, Methods, Results and Main conclusions; validated length is 245 words. |
| Keywords | ✓ Pass | Seven alphabetized keywords. |
| Main-text length | ✓ Pass | Introduction–Discussion contains 4,537 words in the validated branch state. |
| Introduction | ✓ Pass | Defines the spatial-organization question and separates polymorphism from geographic differentiation. |
| Automated literature chronology | ✓ Pass | Initial global output, deferred follow-up, evidence aggregation and ambiguous-case enrichment are dated from 16–19 July 2026 by repository commits. |
| Classification review state | △ Rule-derived; human review pending | Labels are source-traceable and rule-derived; the resolved queue remains `unreviewed`. S18 is a blank blinded 34-species review sheet, S19 is the separate key, and the evaluator computes agreement without changing analysis inputs. |
| Classification transparency | ✓ Pass | Four operational classes, deterministic rules, the frozen manifest, 34/34 rule reproduction, source matching, plural-expression audit and the human-review boundary are documented. |
| GBIF methods | ✓ Pass analytically | Primary and paginated sensitivity workflows are described exactly. |
| Exact GBIF citation bundle | ✓ Prepared and validated | The exact subset, parent-dataset counts, species counts, broad request, Derived Dataset metadata and manifest are fixed and cross-validated. |
| GBIF DOI issuance | △ Authenticated action required | A broad download DOI and exact-subset Derived Dataset DOI have not yet been issued. |
| Climate and metrics | ✓ Pass | WorldClim variables, resolution, PCA and five metrics are reproducible. |
| Statistical methods | ✓ Pass | Formula, covariance, confidence intervals, seeds, permutations, family deletion and diagnostics are reported. |
| Open Tree sensitivity | ✓ Pass with limitation | Two 100-replicate topology-based analyses completed; the model retained 30 species and used Grafen branch lengths. |
| Dated-megaphylogeny sensitivity | ✓ Pass with limitation | Six fixed-seed V.PhyloMaker2 models completed on time-scaled `GBOTB.extended.LCVP` trees; all retained 34 species, but six species were inserted. |
| Results | ✓ Pass | Primary, broader, matched-control, occurrence-sampling and both phylogenetic treatments are reported with uncertainty. |
| Discussion | ✓ Pass | Causal overstatement is removed; agreement in direction and uncertainty in inference are both explicit. |
| References | ✓ Pass for cited text | Fifteen references support the current text. Final GBIF data citations await issued identifiers. |
| Tables and figures | ✓ Pass | Four main tables, two figure legends and Tables S1–S19 are cross-referenced. |
| Supporting provenance | ✓ Pass | Analysis tables, literature chronology, S18/S19 review packet, classification protocol and rule audit, placement audit, phylogenies, exact GBIF bundle and manifests are indexed or directly referenced. |
| Automated package QA | ✓ Pass | File presence, figure links, required estimates, forbidden claims, S1–S19 row counts, all indexed SHA-256 values, exact GBIF counts and author-prefill boundaries are checked in CI. |
| Data Accessibility | △ Ready for DOI insertion | The exact citation boundary and archive files are documented; final repository, broad-download and Derived Dataset DOI citations remain. |
| Deterministic release package | ✓ Preview passes | The preview records its source commit and deterministic ZIP hash in the external workflow artifact, avoiding self-reference inside packaged documents. |
| Strict release package | △ Deliberately blocked | Four placeholder groups remain: cover letter, manuscript DOI/human-review wording, title page and Zenodo metadata. |
| Verified author identity | ✓ Partial | ZHANG RUIQI is prefilled from a current author-provided record. |
| Verified division and institution | ✓ Partial | Division of Forest and Biomaterials Science, Graduate School of Agriculture, Kyoto University is prefilled using official English names. |
| Verified postal address and email | ✓ Partial | The general Kyoto University Agriculture address and `zhang.ruiqi.77h@st.kyoto-u.ac.jp` are prefilled from official and university-linked records. |
| Final authorship and declarations | ✗ Missing author confirmation | Author order, coauthors, laboratory or room, corresponding-author designation, ORCIDs, funding, CRediT, conflicts, acknowledgements and biosketch remain unverified. |
| Author confirmation workflow | ✓ Partially prefilled | Verified identity, division, institution, postal address and email are entered; all remaining author-controlled fields and sign-offs are consolidated in `docs/jbi_author_confirmation_form.md`. |
| Cover letter | △ Partially prefilled | Scientific claims and the verified identity, division, institution, postal address and email are entered; corresponding-author designation, exclusivity, conflict and DOI declarations require approval. |
| Taxon image | △ Candidate verified | A high-resolution CC0 image and caption are identified; final author approval, download and submission upload remain. |
| Archive metadata | ✓ Template prepared | Release protocol and Zenodo metadata template are present; final creator list, licence, version and identifiers require approval. |

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

## Round-nine change history

| Before | After | Reason |
|---|---|---|
| Author division, address and email marked unresolved | Added the verified Division of Forest and Biomaterials Science, general institutional postal address and Kyoto University email | Reduced avoidable submission entry while keeping role and authorship decisions separate. |
| Cover-letter signature lacked verified contact details | Prefilled the division, address and institutional email behind a corresponding-author confirmation guard | Preserved the distinction between verified contact information and unconfirmed correspondence role. |
| ORCID field generic | Recorded that no ORCID linked to the institutional email was found in the searched mailbox records | Prevents Web of Science ResearcherID or another identifier from being misreported as ORCID. |
| Laboratory could be inferred from supervisor context | Kept laboratory unverified because available records list multiple fields and do not uniquely assign the author | Avoids a plausible but unsupported affiliation detail. |
| Author prefill protected only by prose | Extended `validate_jbi_author_prefill.py` to require the exact division, address and email and to reject any other email or ORCID-like value | Prevents metadata drift or silent invention. |
| GBIF DOI listed only as a missing item | Frozen and validated the exact 58,455-record subset, 389 parent-dataset contributions, broad request and Derived Dataset metadata | Reduced DOI work to authenticated submission and permanent hosting rather than data reconstruction. |
| No deterministic final package | Added a release builder with fixed ZIP metadata, file hashes and source-commit manifest | Makes the journal-facing package reproducible from a frozen commit. |

## Editor Check

### Provisional decision: **Major revision before submission**

The manuscript now addresses occurrence-sampling sensitivity and phylogenetic non-independence using two distinct phylogenetic constructions, documents the automated evidence-search chronology, freezes the exact GBIF occurrence subset and passes integrated package, citation, author-prefill and release-preview audits. Directional consistency is strong: the moisture-breadth estimate remained negative in every family-deletion refit, every Open Tree replicate and every dated-megaphylogeny scenario. Inferential certainty remains limited because all phylogenetic confidence intervals include one, moisture breadth was selected from a 20-specification matrix and the literature-derived response can contain classification error.

The remaining `Major revision` label no longer reflects missing search dates, a missing phylogenetic analysis, an unlocated image, an unreconstructed occurrence subset, unknown core author contact details or an internally inconsistent submission package. It reflects exploratory focal selection, observational scale mismatch, non-random evidence assembly, incomplete documentation of human review, unresolved final authorship and declarations, and DOI issuance that requires external accounts.

### Likely reviewer concerns, in priority order

1. **Post-analysis focal selection.** Moisture breadth was selected from five metrics and four thresholds.
2. **Classification error.** Among-population status may reflect under-documentation of local coexistence.
3. **Scale mismatch.** Species-level realised niches cannot identify morph-specific climatic sorting.
4. **Phylogenetic residual uncertainty.** Both phylogenetic treatments preserve the negative direction, but every interval includes one; six taxa were inserted into the dated backbone.
5. **Non-random literature sampling.** Research effort, English queries and metadata availability shape inclusion.
6. **Human review documentation.** The repository establishes that the current labels are deterministic screening labels and that the resolved queue remains unreviewed; independent human verification must still be completed or explicitly reported as not performed.
7. **Occurrence-data limits.** The paginated sample remains capped even though its exact contents and citation contributions are now frozen.
8. **Small comparative sample.** The strict comparison contains 34 species.
9. **Mechanistic non-identifiability.** Occurrences are not labelled by flower-colour morph.

## Remaining submission blockers

- Confirm final author order, any additional authors, laboratory or room when required, corresponding author and ORCID in `docs/jbi_author_confirmation_form.md`.
- Record the actual human screening arrangement and approve the resulting manuscript wording.
- Submit the prepared broad GBIF request through an authenticated account.
- Publish the exact gzip at a permanent URL and register the GBIF Derived Dataset.
- Freeze and archive the final repository release, then insert its DOI.
- Complete and approve CRediT, funding, conflicts, acknowledgements, biosketch and cover-letter declarations.
- Approve the CC0 *I. purpurea* image candidate, download the original-resolution file, retain the licence record and upload the final image.
- Run the release builder in strict mode and require `strict-pass` on the frozen submission commit.
- A bespoke species-level molecular tree remains optional editorial strengthening rather than an unaddressed analytical omission.
