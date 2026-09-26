# FCP — current flower-colour polymorphism paper

This repository is organized around the **active New Phytologist flower-colour polymorphism manuscript**. Historical JBI, 34-species comparative, RGFCA-development and `disttrait` work remain preserved for provenance, but they are not part of the active manuscript execution surface.

## Start here

1. **Submission manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md)
2. **Reproducibility contract:** [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md)
3. **Claim-to-input lineage:** [`docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`](docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md)
4. **Frozen claim ledger:** [`docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`](docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md)
5. **Supporting Information map:** [`docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`](docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md)
6. **Figure plan:** [`docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`](docs/POLYMORPHISM_FIGURE_PLAN_20260918.md)
7. **Canonical figures:** [`docs/figures/polymorphism_20260918/`](docs/figures/polymorphism_20260918/)
8. **Current provenance receipt:** [`archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`](archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md)

## Active mainline — global flower-colour polymorphism

The paper treats within-species flower-colour diversity as a species-level phenotype, validates its reproducibility, tests a prospectively frozen achromatic–chromatic colour-space alignment in a species-disjoint third cohort, and relates polymorphism amount to within-species geographic colour organization. The frozen prospective verdict is `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`. The third-cohort confirmation is species-disjoint within the same iNaturalist opportunity universe but is **not an independent-source replication**. The manuscript keeps that confirmation separate from later technical-validity and environmental follow-ups. Machine-readable frozen results and the lineage map control numerical interpretation when prose and artifacts differ.

## Reproducibility boundary

The current paper deliberately keeps its established paths stable. Manuscript-facing `docs/`, frozen `results/`, analysis `scripts/`, and `archive/fcp_submission_20260925/` are **not moved during repository cleanup** because those paths are referenced by protocols, tests, receipts and the provenance package.

Large measured tables that are intentionally absent from current `main` are recovered from immutable Git commits and verified by SHA256. Exact WorldClim bytes are mirrored in the repository release `fcp-worldclim-2.1-10m-20260925`. The self-contained manuscript snapshot is maintained under `fcp-np-provenance-20260926`; its Git-tracked receipt records the source commit, asset size, SHA256 and packaged-file count for the currently published asset.

## Repository layout

| Path | Role |
|---|---|
| `docs/POLYMORPHISM_*` | Active manuscript, protocols, claims and reader-facing audits |
| `results/polymorphism_*` | Frozen/current machine-readable paper results |
| `scripts/analysis/` | Study-specific analysis and figure code |
| `tests/test_polymorphism_*` | Claim and submission regression guards |
| `archive/fcp_submission_20260925/` | Permanent intermediate inputs, checksums and release receipts |
| `.github/workflows/` | **Active current-paper automation only** |
| `archive/workflows/` | Frozen historical workflow definitions; intentionally inactive |
| `packages/disttrait/` | Reusable methods package, separate from the frozen study-specific numerical pipeline |
| `docs/JBI_*`, `docs/jbi_*` | Retained historical/separate-paper documentation |

## Historical lanes

Two older inferential lanes remain recoverable but are not executed automatically from the active Actions directory:

- the Chapter-1 spatial-photograph/JBI lane;
- the frozen 34-species comparative climatic-niche lane.

Their former workflow definitions are retained byte-for-byte under `archive/workflows/jbi/`. Historical `disttrait` development workflows are under `archive/workflows/disttrait/`. The repository history and archived workflow files preserve provenance without presenting those pipelines as part of the current paper.
