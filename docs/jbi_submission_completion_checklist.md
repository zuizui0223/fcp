# Journal of Biogeography submission completion checklist

This checklist contains only fields that still require verified author input, authenticated external actions or final archiving. Do not infer missing values.

## 1. Literature-review record

| Field | Status | Required action |
|---|---|---|
| Automated OpenAlex queries and filtering rules | ✓ Verified | Documented in the anonymized manuscript. |
| Automated search and follow-up dates | ✓ Verified | Repository commits record evidence-pipeline activity from 16 to 19 July 2026; see `docs/jbi_literature_search_provenance.md`. |
| Human screener names and number of screeners | ✗ Author confirmation required | Complete Section 1 of `docs/jbi_author_confirmation_form.md`. |
| Independent duplicate screening or agreement assessment | ✗ Author confirmation required | State what was actually done; explicitly say when it was not performed. |
| Disagreement-resolution procedure | ✗ Author confirmation required | Describe the actual procedure; do not create a retrospective protocol. |
| Language and database-coverage limitations | ✓ Verified | Already stated in Methods and Discussion. |

## 2. Authorship and declarations

| Field | Status | Verified value or action |
|---|---|---|
| Consolidated author-input form | ✓ Prepared | `docs/jbi_author_confirmation_form.md` collects all author-controlled decisions and sign-offs. |
| Final author order | ✗ Author confirmation required | Obtain approval from all authors. |
| Affiliations | ✗ Author confirmation required | Enter official institution and department names. |
| Corresponding author | ✗ Author confirmation required | Confirm name, postal address and email. |
| ORCID identifiers | ✗ Author confirmation required | Enter only verified ORCIDs. |
| CRediT roles | ✗ Author confirmation required | Agree roles for every author. |
| Funding and grant numbers | ✗ Author confirmation required | Copy exactly from award records. |
| Conflict of interest | ✗ Author confirmation required | Obtain a declaration from every author. |
| Acknowledgements | ✗ Author confirmation required | Confirm contributors, permits, institutional help and material support. |
| Biosketch | ✗ Author confirmation required | Prepare after author list is frozen. |

## 3. Data, DOI and release preparation

| Field | Status | Required action |
|---|---|---|
| Submission branch and manuscript files | ✓ Verified | Draft PR #5 contains the edited manuscript, Appendix S1 and Supporting Tables S1–S17. |
| Primary-analysis artifact | ✓ Verified | Workflow run `29972327794`; digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`. |
| GBIF/Open Tree sensitivity artifact | ✓ Verified | Workflow run `30067762848`; artifact `8586932030`; digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`. |
| Dated-megaphylogeny artifact | ✓ Verified | Fixed-seed workflow run `30076757379`; artifact `8590190840`; digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`. |
| Integrated submission and GBIF validation | ✓ Verified | `PR JBI submission package validation` checks journal structure, numerical guardrails, S1–S17 hashes and the exact GBIF bundle on every relevant change. |
| Exact GBIF occurrence subset | ✓ Frozen and validated | 58,455 rows, 34 species, 58,455 unique occurrence keys and 389 parent datasets; archive SHA-256 `f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094`. |
| Broad GBIF download request | ✓ Prepared | Replace the account email and submit with `submit_jbi_gbif_download.py` using authenticated GBIF credentials. |
| GBIF Derived Dataset registration bundle | ✓ Prepared | Parent-dataset counts, exact archive and metadata template are in `docs/supporting/jbi_gbif_doi_bundle/`. |
| Broad GBIF download DOI | △ Authenticated action required | Submit the prepared request and record the completed download DOI. |
| Exact-subset archive URL/DOI | △ External archive required | Publish the frozen gzip at a permanent public URL before Derived Dataset registration. |
| GBIF Derived Dataset DOI | △ Authenticated registration required | Register the exact subset using the prepared 389-dataset contribution table. |
| Permanent repository release DOI | △ External archive required | Freeze the final commit and archive the exact release through Zenodo or another approved service. |
| Archive protocol and Zenodo metadata | ✓ Prepared | See `docs/jbi_archive_release_protocol.md` and `docs/jbi_zenodo_metadata_template.json`. |
| Deterministic release preview | ✓ Verified | `PR JBI release dry run` validates both packages and builds a manifest with source commit, file sizes, SHA-256 values and the deterministic journal-facing ZIP. |
| Release-integrity boundary | ✓ Verified | Known generated reports are ignored; any substantive dirty working-tree change fails strict mode. |
| Data Accessibility Statement | △ Ready for identifier insertion | Exact subset, citation boundary and protocol are documented; insert final release, broad-download and Derived Dataset DOI citations. |

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
| Author confirmation form | ✓ Prepared: `docs/jbi_author_confirmation_form.md` |
| Literature-search provenance | ✓ Prepared: `docs/jbi_literature_search_provenance.md` |
| Supporting Information index | ✓ Prepared and SHA-256 checked |
| Tables S1–S17 | ✓ Prepared |
| Open Tree topology and dated S1–S3 trees | ✓ Prepared |
| Main tables and figure legends | ✓ Included in manuscript |
| Taxon image rights and caption candidate | ✓ Verified candidate | *Ipomoea purpurea* image is CC0; source, credit, caption and interpretation boundary are recorded in `docs/jbi_taxon_image_candidate.md`. |
| Final taxon image approval and upload | △ Author action required | Approve the candidate, download the original file, retain the licence page and upload the final image. |
| Cover letter | △ Template prepared | `docs/jbi_cover_letter_template.md` contains only verified novelty and result statements; author approval, exclusivity, conflict and DOI declarations remain bracketed. |
| Archive/release protocol | ✓ Prepared | Includes preview and strict release commands and Zenodo metadata workflow. |
| Suggested or opposed reviewers | ✗ Author confirmation required | Enter conflict-checked candidates in the author confirmation form. |

## 6. Strict release blockers detected by CI

The current preview is expected to report exactly four unresolved placeholder files:

1. `docs/jbi_cover_letter_template.md` — author declarations and contact details;
2. `docs/jbi_manuscript_editorial_revision.md` — human-review wording and final DOI citations;
3. `docs/jbi_title_page_template.md` — authorship and declarations;
4. `docs/jbi_zenodo_metadata_template.json` — creator, licence, version and identifiers.

Any additional placeholder file fails the release-preview workflow. Do not remove the four expected guards until verified values are available.

## 7. Final release gate

Do not mark the manuscript submission-ready until all of the following are true:

- final author list and declarations are approved;
- the actual human screening/classification arrangement is documented;
- the broad GBIF download, exact-subset archive, GBIF Derived Dataset and repository-release identifiers are inserted;
- the cover letter's author-controlled declarations are completed and approved;
- the taxon image candidate is approved and the final original-resolution file and licence record are retained;
- all manuscript, table, figure and Supporting Information cross-references pass the integrated validator;
- `python prepare_jbi_submission_release.py --strict` returns `strict-pass` on the frozen commit.
