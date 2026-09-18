# New Phytologist review-document packaging QA — 2026-09-18

Status: **PASS — review-ready DOCX built and independently rendered/inspected.**

This receipt records packaging only. It changes no frozen biological result, cohort, threshold, null, decision rule or manuscript claim boundary.

## Source and workflow identity

- source branch: `analysis/h2-third-cohort-preopening-20260916`
- polished manuscript source commit: `40aac30542313371225f6367e64d26752f77e0ee`
- document workflow: `.github/workflows/polymorphism-new-phytologist-document.yml`
- workflow run: `35319622806`
- workflow conclusion: **success**
- workflow artifact ID: `10537005224`
- artifact name: `fcp-new-phytologist-review-manuscript`
- artifact digest: `sha256:84c2b4f450fc0a9dbdd9e66e23b4c664f4410657df14906967353d90494e15d6`

The source commit made two editorial-only corrections before the final build:

1. the long machine verdict in Table 1 was replaced by the human-readable cell text `Prospective H2 confirmed`;
2. an internal repository-audit note was removed from the submitted References section.

No scientific value or interpretation changed.

## Automated guards

For commit `40aac30542313371225f6367e64d26752f77e0ee`:

- New Phytologist submission guard run `35319622809`: **success**;
- review-document build run `35319622806`: **success**.

## Independent rendering QA

The final DOCX artifact was extracted and rendered through the repository document QA path.

Rendered document:

- pages: **23**
- continuous line numbering: **present**
- page numbering: **present**
- 1.5-line review layout: **present**
- main Table 1: **readable**
- Figures 1–5: **all readable at review scale**
- figure legends: **present and paired with figures**
- References: **no clipping or internal audit note**
- Supporting Information legends S1–S9: **present**
- page clipping/overlap/material text loss: **none observed**

Every rendered page 1–23 was inspected visually after the final editorial corrections.

## PDF preflight

A PDF rendered from the final DOCX passed preflight:

- pages: **23**
- encrypted: **false**
- openable with PyMuPDF: **true**
- likely scanned: **false**
- XFA present: **false**

Export-byte checksums from this QA session:

- final DOCX SHA256: `017c055c45cf62e603db0ed53612bb3d1b70ce7c166c41ebf487d567279aa5e4`
- rendered PDF SHA256: `164bc6f2d995a2ee24458d32b4b5f0ea1c0c238534e7e8df108fcaf92a241a20`

These checksums describe the exported review files inspected in this QA session; the GitHub Actions artifact digest above remains the canonical workflow-artifact identity.

## Remaining submission blockers

The technical document-packaging gate is closed.

Formal submission still requires external/author input that must not be inferred from repository history:

- final author order;
- affiliations;
- corresponding-author name/email and ORCID;
- funding/acknowledgements;
- competing-interests declaration;
- author-contribution statement;
- permanent archive DOI/version.

After those values are inserted, the same document build and guard should be rerun once to produce the administrative-final submission file.
