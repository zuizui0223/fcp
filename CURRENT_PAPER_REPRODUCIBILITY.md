# Current paper reproducibility contract

This file is the shortest supported route for reproducing or auditing the **active New Phytologist FCP manuscript**.

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
- D finite-sample sensitivity receipt: `results/polymorphism_D_finite_sample_sensitivity_20260928/result.json`
- mirrored frozen Step-8 source: `archive/fcp_submission_20260925/d_finite_sample/step8_frozen_result.json`
- distributed-polymorphism receipt: `results/polymorphism_distributed_polymorphism_posthoc_20261007/result.json`
- distributed-polymorphism robustness receipt: `results/polymorphism_distributed_polymorphism_robustness_20261007/result.json`
- phenotypic IBD/IBE-like receipt: `results/polymorphism_phenotypic_IBD_IBE_posthoc_20261007/result.json`
- BIO5 sorting robustness receipt: `results/polymorphism_BIO5_colour_sorting_robustness_20261007/result.json`

## Completed five-cohort original real-photo extension (separate from the frozen original H2 package)

The current five-cohort expanded study is **genuinely measured**, not a 447,500-photo source-metadata opportunity only. The original-method ROI-v4/EfficientSAM [run #38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400) completed all 256 partitions; the source-authenticated [real scientific join and H1/50-km/W replay #39049598596](https://github.com/zuizui0223/fcp/actions/runs/39049598596) completed all five actual jobs; the [permanent release #39050465077](https://github.com/zuizui0223/fcp/actions/runs/39050465077) independently authenticated actual pixel-result row identities and unchanged old source hashes/statistics.

**Canonical completed five-cohort real data:** [fcp-all4730-original-method-actual-colour-results-20261017](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017). Archive `fcp_all4730_original_method_actual_measured_data_and_inference.tar.gz`, SHA256 `99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a8de07edd38341b27889310bc25`, **108,115,413 bytes**. Full original source-joined measured photo CSV SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`; release receipt SHA256 `b0b152fcfaa1e0015a6e5524192265bb95c1663c10b5da3a4e885cb7be49747c`. The machine-readable final original 50-km/H1/W output and source-authenticated paper Results packets are included in this release, including JSON SHA256 `e4d5244fb7e9de0046d44cf01154573d12f2ea3418fa538022a8dd97b1d60dc6`.

**Reporting entry point:** `docs/FCP_FINAL_447500_REAL_FIVE_COHORT_RESULTS_20261018.md` and `results/fcp_original_447500_five_cohort_verified_reporting_20261018.json`. Exactly 4,475 fully measured independent 100-photo species, 447,500 terminals, 227,107 classifiable photos, 3,297 >=40-photo D-eligible species and **1,761 with >=30 local pairs at the post hoc 50-km radius**. These denominators are not interchangeable. Original discovery/validation/prospective source measurements and original unchanged spatial statistics were checked against immutable hashes/expected effects. Same iNaturalist provider is not independent-source replication. Five-cohort W is not a new untouched prospective experiment.

**Two distinct archive scopes:** the following longstanding `fcp-np-provenance-20260926` snapshot reconstructs the first **149,900-photo prospectively frozen study** and its original secondary evidence. It is **not** a standalone byte-complete archive of the new 447,500-photo extension. For the expanded paper, retain both the original prospective provenance snapshot **and the separate genuine completed five-cohort measured-results Release**. No underlying third-party JPEG source images are republished.

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
| Discovery-cohort measured rows | commit `5142f7951af0dde5364bb047a566d67e8c479e51`, `data/derived/global_monte_carlo_measured_photos_v1.csv` | `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4` |
| Validation-cohort measured rows | commit `5142f7951af0dde5364bb047a566d67e8c479e51`, `data/derived/rgfca_reserve_replication_measured_photos_v1.csv` | `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6` |
| Prospective-confirmation measured rows | commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`, `data/derived/polymorphism_h2_third_cohort_measured_photos_v1.csv` | `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186` |

Exact WorldClim 2.1 10-arc-minute inputs are mirrored under release tag `fcp-worldclim-2.1-10m-20260925`:

- BIO archive SHA256: `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`
- SRAD archive SHA256: `c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026`

Raster-level hashes are frozen in `archive/fcp_submission_20260925/worldclim_checksums.txt`.

## Reproducibility levels

The paper distinguishes five levels of reproducibility.

1. **Core numerical inference replay from frozen measured tables — demonstrated.** Exact discovery/validation and prospective-confirmation measured tables are byte-frozen. Workflow `.github/workflows/polymorphism-reproducibility-replay.yml` reruns repeated observer-disjoint H1, deterministic H1 stress, original-cohort white-axis localization, prospective confirmation H2, and the D–spatial Step 5/5b/6/8/9 chain. The current finite-sample reporting receipt aggregates the already frozen corrected-D sensitivities; raw-versus-corrected species ranks exceed 0.99997 in both high-depth cohorts, and corrected D–spatial p-values remain 0.007 in discovery and 0.025 in validation. The H1/H2 JSONs reproduce at absolute tolerance `1e-12`; the spatial chain passes absolute `1e-10` / relative `1e-12`. The runtime is pinned by `requirements-np-replay-20260928.txt` under Python 3.12.14.
2. **Secondary ecological-analysis replay — demonstrated.** Workflow `.github/workflows/polymorphism-secondary-reproducibility-replay.yml` reruns the prospective-cohort BIO5/BIO14/SRAD analysis, observer-paired/balanced sensitivities, BIO5 transport, and H3a/H3b. Stable cross-run verification is run `36380786248`. BIO5-family outputs pass absolute `1e-8` / relative `1e-10`; this tolerance is intentionally wider than machine epsilon because a separate runner exposed a conditional-logit CI difference of about `1.1e-10` despite identical software versions and decisions. H3a/H3b JSONs pass absolute `1e-9` / relative `1e-11`, and all six archived H3 permutation/summary/PGLS output files reproduce their original SHA256 values exactly. The BIO5 Python environment is pinned in `requirements-np-bio5-replay-20260928.txt`; H3 direct package versions are recorded in `docs/POLYMORPHISM_SECONDARY_REPLAY_ENVIRONMENT_20260928.md`.
3. **Post hoc spatial-evolution diagnostics — deterministic and source-frozen.** Distributed polymorphism, its observer/nonwhite/continuous/temporal falsification tests, the phenotypic IBD/IBE-like decomposition, and BIO5 sorting robustness all read the same immutable high-depth measured tables listed above. Their workflows preserve exact colour compositions or colour-vector multisets under matched vertex nulls and hard-code `confirmatory_decisions_changed=false`. These analyses are reproducible follow-up evidence but remain post-outcome and cannot be reclassified as prospective confirmation.

4. **Outcome-blind sampling lineage and evidence/package reconstruction — supported.** The provenance tarball contains the exact 42,111-species opportunity frame, permanent prior high-depth selection/exclusion bytes, deterministically reconstructed 3,230-species candidate frame, all four prior-ID exclusion sources, selected/candidate/authorized prospective-cohort metadata, the frozen measured tables, frozen execution source, permanent spatial null arrays, original BIO5/H3 artifacts, checksum-pinned WorldClim bytes, figures and regression guards. Every final payload file is covered by `FILE_SHA256SUMS.txt`.
5. **Image-measurement implementation reconstruction — self-contained after byte completion.** The archive contains the exact ROI-v4 / fixed-palette measurement source, trained ROI-v4 detector byte, exact EfficientSAM encoder/decoder ONNX weights, and frozen prospective-cohort source-photo metadata. Every model/input byte added at packaging is SHA256-verified.
6. **Raw source-image reconstruction — intentionally not guaranteed.** Candidate image pixels/masks were deleted after partition sealing. Frozen photo IDs/source URLs and per-row `image_sha256` values allow future reacquired bytes to be verified while the external provider still serves them, but the archive cannot reconstruct provider bytes that later disappear or change.

These levels should not be conflated. The named core and secondary statistical analyses above are directly replayed from archived evidence, while the measurement implementation/model bytes are archived but the original third-party image pixels are not. Supporting components outside the replay contract retain their own provenance status in the coverage matrix of `docs/POLYMORPHISM_REPRODUCIBILITY_AUDIT_20260927.md`.

## Active Actions surface

Only current-paper preservation, claim guards, figure regeneration and bounded follow-up workflows live in `.github/workflows/`. Provenance-only workflow definitions are retained under `archive/workflows/` and are not executed automatically.

The repository-layout guard fails if a provenance-only or superseded acquisition/method-development workflows reappears in the active Actions directory without an explicit repository-layout change.

## Local verification

```bash
python -m pip install -e .
python -m pip install pytest
python -m pytest \
  tests/test_polymorphism_manuscript_claims.py \
  tests/test_polymorphism_new_phytologist_submission.py \
  tests/test_make_polymorphism_manuscript_figures.py \
  tests/test_repository_layout.py -q

# GitHub Actions additionally performs frozen numerical replays:
# .github/workflows/polymorphism-reproducibility-replay.yml
# .github/workflows/polymorphism-secondary-reproducibility-replay.yml
```

The provenance package itself is produced by `.github/workflows/build-np-provenance-snapshot-20260926.yml`, which recovers immutable large inputs and frozen execution source, embeds and verifies the image-measurement detector/EfficientSAM bytes and frozen source-photo metadata, adds checksum-pinned WorldClim bytes, writes a deterministic archive, updates the release asset and records the resulting identity in the Git-tracked receipt.
