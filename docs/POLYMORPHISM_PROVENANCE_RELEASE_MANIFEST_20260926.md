# New Phytologist provenance snapshot manifest — 2026-09-26

## Role

This is a reproducibility/provenance snapshot for the active FCP flower-colour polymorphism paper.

It is **not** a declaration that the manuscript is finally submitted or author metadata is complete. Its purpose is to preserve the evidence chain behind the current numerical claims in one durable release asset.

Release tag:

`fcp-np-provenance-20260926`

## Snapshot refresh contract

The release tag is the stable reader entry point, while the Git-tracked receipt at `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md` defines the exact current asset identity. When the active manuscript, figures, current result surface, relevant analysis code/tests or reproducibility routing changes, the build workflow regenerates the deterministic package from that repository head, replaces the release asset and records its source commit, byte count, SHA256 and file count in the receipt.

Refreshing the package **does not recompute, alter or upgrade any frozen biological result**. The immutable raw-input commits and checksums below remain unchanged. Earlier receipt identities remain audit-visible in Git history.

The canonical receipt is intentionally **not embedded inside the tar.gz it identifies**. The workflow writes that receipt only after the package bytes and SHA256 are known, preventing a stale/self-referential receipt from being carried inside a newly refreshed package. The release tag is the current-snapshot asset container; the Git-tracked receipt defines the exact current bytes.

## Self-contained inputs added at package time

The release workflow recovers and checksum-verifies the large measured tables that are intentionally not carried on current main:

- legacy discovery measured rows
  - source commit: `5142f7951af0dde5364bb047a566d67e8c479e51`
  - source path: `data/derived/global_monte_carlo_measured_photos_v1.csv`
  - SHA256: `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`
- legacy reserve measured rows
  - source commit: `5142f7951af0dde5364bb047a566d67e8c479e51`
  - source path: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`
  - SHA256: `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`
- prospective H2 third-cohort measured rows
  - source commit: `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`
  - source path: `data/derived/polymorphism_h2_third_cohort_measured_photos_v1.csv`
  - SHA256: `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186`
- prospective H2 authorized source metadata
  - source commit: `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`
  - SHA256: `13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290`
  - contains the frozen photo IDs/source URLs used for the one-shot acquisition

The exact WorldClim inputs are copied from the already frozen release `fcp-worldclim-2.1-10m-20260925`:

- `wc2.1_10m_bio.zip`
  - bytes: 49,869,449
  - SHA256: `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`
- `wc2.1_10m_srad.zip`
  - bytes: 16,233,364
  - SHA256: `c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026`

## Repository material included

The package includes:

- active New Phytologist manuscript;
- canonical manuscript;
- Supporting Information evidence map;
- cover letter;
- claim ledger;
- paper architecture and figure plan;
- data-lineage map and data-lineage audit;
- RGFCA-to-polymorphism interpretation and 42,111-frame provenance;
- H1/H2/H3 protocols and frozen decision/claim notes;
- highlight-validity protocol/adjudication;
- BIO5/environment protocols;
- P500 terminal postmortem;
- all headline machine-readable result directories;
- publication figure source data and the five current manuscript figures;
- the exact highlight/H3a/H3b artifact freeze under `archive/fcp_submission_20260925`;
- relevant polymorphism/environment analysis scripts;
- relevant manuscript/provenance regression tests;
- relevant GitHub Actions workflow definitions;
- exact historical H1/H2 execution source;
- fuller frozen D–spatial Step 5/6/8/9 scripts and intermediate outputs;
- exact ROI-v4 / fixed-palette measurement implementation;
- trained ROI-v4 detector byte, SHA256 `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`;
- exact EfficientSAM encoder/decoder ONNX weights from revision `d525f622e6f640acf5a0fc37c7ca1f243da5bde0`, checksum-verified;
- frozen third-cohort source-photo metadata for byte-validated reacquisition.

## Packaging rules

The workflow:

1. checks out the exact source commit;
2. recovers the three measured tables from their immutable source commits;
3. verifies their frozen SHA256 values;
4. recovers exact historical H1/H2/spatial execution source and the frozen image-measurement code;
5. embeds and verifies the ROI-v4 detector byte, EfficientSAM ONNX weights and authorized source-photo metadata;
6. downloads the checksum-pinned WorldClim release assets and verifies them;
7. copies the repository evidence set into a staging tree;
8. writes a complete per-file SHA256 manifest;
9. creates a deterministic tar.gz with sorted paths, fixed mtime, numeric owner/group;
10. uploads the tar.gz and SHA256 manifest to the provenance release;
11. records a release receipt on main.

## Authority

Machine-readable frozen results and input SHA256 values remain authoritative over prose.

This snapshot does not upgrade any inferential claim. In particular:

- prospective H2 remains separate from post-H2 BIO5;
- BIO5 remains positive in the secondary third-cohort analysis but non-replicated under the fixed legacy transport rule;
- shared versus species-specific geographic mapping remains unresolved;
- exposure coupling of measured white remains a stated limitation.


## Raw-image reproducibility boundary

The archive does not contain the original iNaturalist image pixels because the prospective protocol explicitly deleted pixels and masks after partition sealing. It does preserve the frozen photo IDs/source URLs, per-row image SHA256 values, exact acquisition/measurement code, ROI-v4 detector byte and EfficientSAM weights. Reacquired bytes can therefore be validated while still available, but exact historical pixel reconstruction is not guaranteed if the external provider no longer serves those bytes.
