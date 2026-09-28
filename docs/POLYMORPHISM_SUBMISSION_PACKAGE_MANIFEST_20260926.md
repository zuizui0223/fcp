# New Phytologist submission / Zenodo package manifest — 2026-09-26

## Role

This is the journal-facing packaging layer for the active FCP manuscript.

It does **not** alter or recompute any biological result. The authoritative self-contained numerical provenance remains the frozen release asset:

- release tag: `fcp-np-provenance-20260926`
- asset: `fcp-np-provenance-20260926.tar.gz`
- exact embedded SHA256: verified against the release-side `fcp-np-provenance-20260926.tar.gz.sha256` file at package-build time
- source commit: read from `PROVENANCE_SOURCE.txt` inside the verified tarball

## Generated release candidate

The packaging workflow produces three top-level assets:

1. `fcp-new-phytologist-submission-20260926.zip`
   - journal-facing manuscript, cover letter, five main figures plus Fig. S9, Supporting Information, current data-lineage map and journal-only checksums.
2. `fcp-zenodo-ready-20260926.tar.gz`
   - the journal-facing package plus the exact self-contained provenance tarball.
3. `fcp-release-sha256-20260926.txt`
   - SHA256 identities for both packages.

Candidate GitHub release tag:

`fcp-np-submission-candidate-20260926`

## Journal-facing contents

The journal ZIP is intentionally minimal:

- `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`
- `docs/POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md`
- `docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`
- `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`
- main Figure 1–5 PDF and PNG files, Fig. S9 PDF and PNG, plus the figure manifest
- journal-only SHA256 manifest

The following repository-management surfaces are **not** included in the journal ZIP: submission-readiness audit, claim ledger, secondary mechanism evidence ledger, provenance-release manifest, reproducibility audit, Zenodo metadata draft, package manifest, and the deep self-contained provenance tarball. They remain available in the Zenodo-ready archive layer.

## Publication gate

The package is **scientifically frozen but administratively incomplete**.

Before a public Zenodo deposition is published, author approval is required for:

- final author list and order;
- corresponding author;
- exact affiliations;
- ORCIDs;
- funding / grant numbers;
- acknowledgements;
- competing interests;
- CRediT roles;
- archive licence.

These fields must not be inferred from repository history.

## Zenodo rule

The candidate package may be uploaded to a draft Zenodo deposition for byte verification, but **do not publish/mint the final version DOI until the author metadata above are approved**.

The final Zenodo version should archive the exact `fcp-zenodo-ready-20260926.tar.gz` bytes or a later package built by the same workflow after author metadata are finalized. Pull-request dry runs verify the currently published provenance tarball against its own release-side SHA256. On main, the submission/Zenodo candidate is built only after the provenance workflow completes successfully, preventing the candidate from racing a provenance refresh.