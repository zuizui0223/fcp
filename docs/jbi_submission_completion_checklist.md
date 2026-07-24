# Journal of Biogeography submission completion checklist

This checklist contains only fields that still require verified author input, external archiving or an editorial decision. Do not infer missing values.

## 1. Literature-review record

| Field | Status | Required action |
|---|---|---|
| Automated OpenAlex queries and filtering rules | ✓ Verified | Already documented in the anonymized manuscript. |
| Original execution date for each discovery and follow-up search | ✗ Not verified | Recover from workflow metadata, logs or dated records and enter exact dates. |
| Manual screener names and number of screeners | ✗ Not verified | Record the people who actually screened or adjudicated evidence. |
| Disagreement-resolution procedure | ✗ Not verified | Describe the procedure actually used; do not create a retrospective protocol. |
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
| Submission branch and manuscript files | ✓ Verified | Draft PR #5 contains the edited manuscript and Supporting Tables S1–S15. |
| Primary-analysis artifact | ✓ Verified | Workflow run `29972327794`; digest `sha256:87d8c9ba89f27685e362abeffa0e077330adb1923652f7a7df73572c5e274ac8`. |
| GBIF/phylogenetic sensitivity artifact | ✓ Verified | Workflow run `30067762848`; artifact `8586932030`; digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`. |
| Permanent code/data archive DOI | ✗ Not verified | Freeze the exact submission commit, archive it and insert the DOI. |
| GBIF derived-dataset DOI | ✗ Not verified | Create a citable GBIF download or derived-dataset record matching the submitted occurrence data. |
| Data Accessibility Statement | △ Drafted | Replace both `Not verified` DOI placeholders after archiving. |

## 4. Phylogenetic decision

| Field | Status | Required action |
|---|---|---|
| Family-clustered and family-deletion analyses | ✓ Verified | Reported in the manuscript. |
| Open Tree topology sensitivity | ✓ Verified | 30 species; 100 completed fits per occurrence dataset. |
| Dated species phylogeny | ✗ Not present | Decide before submission whether the topology analysis is sufficient or a dated phylogeny is required. |
| Permitted claim | ✓ Fixed | “The negative association was stable to stronger occurrence sampling but remained statistically unresolved after topology-based phylogenetic correction.” |

## 5. Submission files

| File or field | Status |
|---|---|
| Anonymized manuscript | ✓ Prepared: `docs/jbi_manuscript_editorial_revision.md` |
| Separate title page | △ Template prepared; author-controlled fields remain blank |
| Supporting Information index | ✓ Prepared and SHA-256 checked |
| Tables S1–S15 | ✓ Prepared |
| Main tables and figure legends | ✓ Included in manuscript |
| Taxon image, permission and credit | ✗ Not verified |
| Cover letter | ✗ Not drafted from verified final author and archive details |
| Suggested or opposed reviewers | ✗ Not verified |

## 6. Final release gate

Do not mark the manuscript submission-ready until all of the following are true:

- final author list and declarations are approved;
- literature-search dates and manual review procedure are documented;
- the permanent archive DOI and GBIF DOI are inserted;
- the phylogenetic decision is recorded;
- the taxon image and rights are confirmed;
- all manuscript, table, figure and Supporting Information cross-references are checked against the frozen submission commit.
