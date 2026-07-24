# Journal of Biogeography submission completion checklist

This checklist contains only fields that still require verified author input or external archiving. Do not infer missing values.

## 1. Literature-review record

| Field | Status | Required action |
|---|---|---|
| Automated OpenAlex queries and filtering rules | ✓ Verified | Documented in the anonymized manuscript. |
| Automated search and follow-up dates | ✓ Verified | Repository commits record evidence-pipeline activity from 16 to 19 July 2026; see `docs/jbi_literature_search_provenance.md`. |
| Human screener names and number of screeners | ✗ Not verified | Record the people who actually screened or classified evidence. |
| Independent duplicate screening or agreement assessment | ✗ Not verified | State what was actually done; explicitly say when it was not performed. |
| Disagreement-resolution procedure | ✗ Not verified | Describe the actual procedure; do not create a retrospective protocol. |
| Language and database-coverage limitations | ✓ Verified | Already stated in Methods and Discussion. |

## 2. Authorship and declarations

| Field | Status | Verified value or action |
|---|---|---|
| Final author order | ✗ Not verified | Obtain approval from all authors. |
| Affiliations | ✗ Not verified | Enter official institution and department names. |
| Corresponding author | ✗ Not verified | Confirm name, postal address and email. |
| ORCID identifiers | ✗ Not verified | Enter only verified ORCIDs. |
| CRediT roles | ✗ Not verified | Agree roles for every author. |
| Funding and grant numbers | ✗ Not verified | Copy exactly from award records. |
| Conflict of interest | ✗ Not verified | Obtain a declaration from every author. |
| Acknowledgements | ✗ Not verified | Confirm contributors, permits, institutional help and material support. |
| Biosketch | ✗ Not verified | Prepare after author list is frozen. |

## 3. Data and code archiving

| Field | Status | Required action |
|---|---|---|
| Submission branch and manuscript files | ✓ Verified | Draft PR #5 contains the edited manuscript, Appendix S1 and Supporting Tables S1–S17. |
| Primary-analysis artifact | ✓ Verified | Workflow run `29972327794`; digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`. |
| GBIF/Open Tree sensitivity artifact | ✓ Verified | Workflow run `30067762848`; artifact `8586932030`; digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`. |
| Dated-megaphylogeny artifact | ✓ Verified | Fixed-seed workflow run `30076757379`; artifact `8590190840`; digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`. |
| Submission-package validation | ✓ Verified | CI run `30079894645` passed all structural, numerical, reference, file-presence and SHA checks. |
| Permanent code/data archive DOI | ✗ Not verified | Freeze the exact submission commit, archive it and insert the DOI. |
| GBIF derived-dataset DOI | ✗ Not verified | Create a citable GBIF download or derived-dataset record matching the submitted occurrence data. |
| Data Accessibility Statement | △ Drafted | Replace both `Not verified` DOI placeholders after archiving. |

## 4. Phylogenetic evidence

| Field | Status | Verified action or remaining limitation |
|---|---|---|
| Family-clustered and family-deletion analyses | ✓ Verified | Reported in the manuscript. |
| Open Tree topology sensitivity | ✓ Verified | 30 species; 100 completed fits per occurrence dataset. |
| Time-scaled megaphylogeny sensitivity | ✓ Verified with limitation | V.PhyloMaker2 S1–S3 models retained all 34 species; 28 were present in the backbone and six were inserted. |
| Species-specific molecular phylogeny | Optional strengthening | Not required to claim that phylogenetic sensitivity was assessed; it would refine placement and branch-length uncertainty. |
| Permitted claim | ✓ Fixed | “The association was stable to stronger occurrence sampling and remained negative under both topology-based and time-scaled phylogenetic treatments, but all phylogenetic confidence intervals included one.” |

## 5. Submission files

| File or field | Status |
|---|---|
| Anonymized manuscript | ✓ Prepared: `docs/jbi_manuscript_editorial_revision.md` |
| Separate title page | △ Template prepared; author-controlled fields remain blank |
| Literature-search provenance | ✓ Prepared: `docs/jbi_literature_search_provenance.md` |
| Supporting Information index | ✓ Prepared and SHA-256 checked |
| Tables S1–S17 | ✓ Prepared |
| Open Tree topology and dated S1–S3 trees | ✓ Prepared |
| Main tables and figure legends | ✓ Included in manuscript |
| Taxon image rights and caption candidate | ✓ Verified candidate | *Ipomoea purpurea* image is CC0; source, credit, caption and interpretation boundary are recorded in `docs/jbi_taxon_image_candidate.md`. |
| Final taxon image approval and upload | △ Author action required | Approve the candidate, download the original file, retain the licence page and upload the final image. |
| Cover letter | △ Template prepared | `docs/jbi_cover_letter_template.md` contains only verified novelty and result statements; author approval, exclusivity, conflict and DOI declarations remain bracketed. |
| Suggested or opposed reviewers | ✗ Not verified |

## 6. Final release gate

Do not mark the manuscript submission-ready until all of the following are true:

- final author list and declarations are approved;
- the actual human screening/classification arrangement is documented;
- the permanent archive DOI and GBIF DOI are inserted;
- the cover letter's author-controlled declarations are completed and approved;
- the taxon image candidate is approved and the final original-resolution file and licence record are retained;
- all manuscript, table, figure and Supporting Information cross-references pass the validator against the frozen submission commit.
