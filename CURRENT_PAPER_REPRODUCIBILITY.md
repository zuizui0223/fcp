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
- reproducibility audit: `docs/POLYMORPHISM_REPRODUCIBILITY_AUDIT_20260927.md`

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

## Reproducibility levels

The paper distinguishes four levels of reproducibility.

1. **Numerical inference replay from frozen measured tables — demonstrated.** Exact legacy discovery/reserve and prospective-H2 measured tables are byte-frozen. Workflow `.github/workflows/polymorphism-reproducibility-replay.yml` reruns the repeated observer-disjoint H1 test, deterministic H1 stress test, legacy white-axis target analysis and prospective third-cohort H2. Run `36324187508` reproduced all four frozen JSON results at absolute tolerance `1e-12`. Run `36373083779` additionally replayed the D–spatial Step 5/5b/6/8/9 chain from the original 999-permutation null arrays and passed the frozen-result comparison. The primary replay runtime is pinned by `requirements-np-replay-20260928.txt` under Python 3.12.14. Secondary analyses are independently replayed by `.github/workflows/polymorphism-secondary-reproducibility-replay.yml`: run `36379690691` reproduced the third-cohort BIO5/BIO14/SRAD result, observer sensitivities and legacy BIO5 transport (three JSON objects plus seven CSV tables), and replayed H3a/H3b under their recorded R/package versions. H3a/H3b JSONs matched at absolute `1e-9` / relative `1e-11`, while the archived H3 permutation and PGLS outputs reproduced their original SHA256 values exactly.
2. **Secondary ecological-analysis replay — demonstrated.** Workflow `.github/workflows/polymorphism-secondary-reproducibility-replay.yml` reruns the third-cohort BIO5/BIO14/SRAD analysis, its observer sensitivities, the legacy BIO5 transport test, and H3a/H3b. Run `36379690691` reproduced the BIO5 JSON/tables from the original successful artifacts and reproduced H3a/H3b under their recorded R/package versions; the H3 permutation/PGLS tables matched the archived SHA256 identities byte-for-byte. The BIO5 Python environment is frozen in `requirements-np-bio5-replay-20260928.txt` under Python 3.11.16.
3. **Outcome-blind sampling lineage and current evidence/package verification — supported.** The provenance tarball contains the exact 42,111-species opportunity frame, permanent P100/P500 selection bytes, deterministically reconstructed 3,230-species candidate frame, all four prior-ID exclusion sources, selected/candidate/authorized third-cohort metadata, the frozen measured tables, machine-readable results, historical execution-code snapshot, permanent H3/highlight intermediate archive, checksum-pinned WorldClim bytes, current figures and regression guards. Every file in the tarball is covered by `FILE_SHA256SUMS.txt`.
4. **Image-measurement implementation reconstruction — self-contained after byte completion.** The archive contains the exact ROI-v4 / fixed-palette measurement source, the trained ROI-v4 detector byte, the exact EfficientSAM encoder/decoder ONNX weights, and the frozen third-cohort source-photo metadata. Every model/input byte added at packaging is SHA256-verified.
5. **Raw source-image reconstruction — intentionally not guaranteed.** Candidate image pixels/masks were deleted after partition sealing. Frozen photo IDs/source URLs and per-row `image_sha256` values allow future reacquired bytes to be verified while the external provider still serves them, but the archive cannot reconstruct provider bytes that later disappear or change.

These levels should not be conflated: the paper's statistical results are directly replayable and the measurement implementation/model bytes are archived, but no claim is made that every original third-party image byte can be reconstructed indefinitely.

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

# GitHub Actions additionally performs historical-code numerical replays:
# .github/workflows/polymorphism-reproducibility-replay.yml
# .github/workflows/polymorphism-secondary-reproducibility-replay.yml
# .github/workflows/polymorphism-secondary-reproducibility-replay.yml
```

The provenance package itself is produced by `.github/workflows/build-np-provenance-snapshot-20260926.yml`, which recovers immutable large inputs and historical execution source, embeds and verifies the image-measurement detector/EfficientSAM bytes and frozen source-photo metadata, adds checksum-pinned WorldClim bytes, writes a deterministic archive, updates the release asset and records the resulting identity in the Git-tracked receipt.
