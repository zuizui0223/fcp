# Current paper reproducibility contract

This file is the shortest supported route for reproducing or auditing the **active New Phytologist FCP manuscript** without navigating the repository's historical development surface.

## Authority order

When sources disagree, use this order:

1. frozen machine-readable result JSON/TSV/CSV;
2. exact input SHA256 values and release receipts;
3. frozen protocol/decision documents;
4. claim ledger and data-lineage map;
5. manuscript prose.

Repository cleanup must not alter a frozen numerical result, biological decision rule, immutable source commit, or checksum.

## Reader route

- manuscript: `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`
- canonical manuscript: `docs/POLYMORPHISM_MANUSCRIPT.md`
- claim ledger: `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`
- claim-to-input map: `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`
- Supporting Information map: `docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`
- submission-readiness audit: `docs/POLYMORPHISM_NEW_PHYTOLOGIST_SUBMISSION_READINESS_20260918.md`
- provenance package manifest: `docs/POLYMORPHISM_PROVENANCE_RELEASE_MANIFEST_20260926.md`
- current provenance receipt: `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`

## Self-contained snapshot

Current snapshot tag: `fcp-np-provenance-20260926`

Asset: `fcp-np-provenance-20260926.tar.gz`

The tag is the stable reader entry point. The **Git-tracked receipt** is the authority for the currently published asset's source commit, byte count, SHA256 and packaged-file count:

`archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`

The snapshot is rebuilt when the current manuscript/evidence surface changes. Refreshing the package does not recompute or upgrade biological results; it repackages the already frozen evidence at the new repository head. Previous receipt identities remain audit-visible in Git history.

The receipt is deliberately stored **outside** the tar.gz it identifies. It is written only after the final package byte count and SHA256 are known, so a refreshed package cannot contain a stale copy of its own receipt.

## Immutable large-input routes

| Input | Immutable source | SHA256 |
|---|---|---|
| Legacy discovery measured rows | commit `5142f7951af0dde5364bb047a566d67e8c479e51`, `data/derived/global_monte_carlo_measured_photos_v1.csv` | `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4` |
| Legacy reserve measured rows | commit `5142f7951af0dde5364bb047a566d67e8c479e51`, `data/derived/rgfca_reserve_replication_measured_photos_v1.csv` | `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6` |
| Prospective-H2 measured rows | commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`, `data/derived/polymorphism_h2_third_cohort_measured_photos_v1.csv` | `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186` |

Exact WorldClim 2.1 10-arc-minute inputs are mirrored under release tag `fcp-worldclim-2.1-10m-20260925`:

- BIO archive SHA256: `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`
- SRAD archive SHA256: `c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026`

Raster-level hashes are frozen in `archive/fcp_submission_20260925/worldclim_checksums.txt`.

## Active Actions surface

Only current-paper preservation, claim guards, figure regeneration and bounded follow-up workflows live in `.github/workflows/`. Historical workflows are stored under `archive/workflows/` so GitHub Actions does not execute them automatically.

The repository-layout guard fails if a legacy `jbi-*`, `disttrait-*`, P500 acquisition, or old manuscript-consistency workflow reappears in the active Actions directory without an explicit repository-layout change.

## Local verification

```bash
python -m pip install -e .
python -m pip install pytest
python -m pytest \
  tests/test_polymorphism_manuscript_claims.py \
  tests/test_polymorphism_new_phytologist_submission.py \
  tests/test_make_polymorphism_manuscript_figures.py \
  tests/test_repository_layout.py -q
```

The provenance package itself is produced by `.github/workflows/build-np-provenance-snapshot-20260926.yml`, which recovers immutable large inputs, verifies hashes, adds checksum-pinned WorldClim bytes, writes a deterministic archive, updates the release asset and records the resulting identity in the Git-tracked receipt.
