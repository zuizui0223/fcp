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

The release workflow also freezes the outcome-blind sampling lineage that precedes those measured tables:

- complete 42,111-species opportunity frame
  - source commit: `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`
  - SHA256: `5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`
- four prior-ID exclusion sources used before third-cohort selection
  - H9 exclusion ledger: `f9a6894740e9974399c055f92cba237be8ada41707d84e1807ba61b902c91b99`
  - H9 fresh metadata: `111d0f964618c0c3df749a6e4bd29f834214cda5d7a2d9bda7d43cdc9dbf4c6f`
  - 42,111 breadth-measurement panel: `38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605`
  - P500 candidate metadata: `a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595`
- third-cohort selected-species manifest
  - SHA256: `16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59`
- parent P100 / P500 outcome-blind selection artifact
  - original run: `34708044962`
  - artifact: `10302477571` (`polymorphism-h2-prospective-u100-selection-20260913`)
  - artifact digest: `sha256:7e43763e3f1d9fa3ff108154dc6ef9443430b70fdb62532e2f6845582d607abb`
  - P100 SHA256: `1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c`
  - P500 SHA256: `f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4`
- deterministically reconstructed third-cohort candidate frame
  - 3,230 species
  - builder: historical `build_polymorphism_h2_third_cohort_frame_20260916.py`
  - SHA256: `7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e`
- complete third-cohort candidate metadata
  - SHA256: `adec40e29e347b035872f2add95b67011906cb74e28510e6260ed1edbd075711`

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
- secondary mechanism evidence ledger separating study-derived BIO5/PAL-WAL analyses from source-derived biochemical/molecular context;
- PAL/WAL × FCP overlap-feasibility receipt showing that the direct bridge is not estimable (two PAL reserve overlaps, zero WAL and zero prospective overlap);
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
- publication figure source data, the five current main figures and Fig. S9 secondary-mechanism figure;
- the exact highlight/H3a/H3b artifact freeze under `archive/fcp_submission_20260925`;
- original complete BIO5 and legacy-BIO5 result artifacts under `external/secondary_original_artifacts`;
- relevant polymorphism/environment analysis scripts;
- relevant manuscript/provenance regression tests;
- relevant GitHub Actions workflow definitions;
- exact historical H1/H2 and third-cohort candidate-frame execution source;
- exact historical H3a and H3b execution source;
- pinned successful primary replay environment: Python 3.12.14 plus `requirements-np-replay-20260928.txt`;
- pinned successful BIO5 replay environment: Python 3.11.16 plus `requirements-np-bio5-replay-20260928.txt`;
- original complete BIO5/BIO14/SRAD and legacy-BIO5 result artifacts, copied before Actions expiry;
- fuller frozen D–spatial Step 5/6/8/9 scripts and intermediate outputs;
- permanent copies of the exact discovery/reserve 999-permutation spatial-null artifacts from Actions runs `34088925008` and `34178957447`;
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
6. permanently copies the original P100/P500 selection artifact and reconstructs/validates the exact 3,230-species candidate frame, with fallback to the previous permanent provenance release on later rebuilds;
7. copies the exact discovery/reserve spatial-null artifact families before Actions retention expiry, with fallback to the previous permanent provenance release on later rebuilds;
8. downloads the checksum-pinned WorldClim release assets and verifies them;
9. copies the repository evidence set into a staging tree;
10. writes a complete per-file SHA256 manifest;
11. creates a deterministic tar.gz with sorted paths, fixed mtime, numeric owner/group;
12. uploads the tar.gz and SHA256 manifest to the provenance release;
13. records a release receipt on main.

## Authority

Machine-readable frozen results and input SHA256 values remain authoritative over prose.

This snapshot does not upgrade any inferential claim. In particular:

- prospective H2 remains separate from post-H2 BIO5;
- BIO5 remains positive in the frozen primary secondary analysis but is weakened by post hoc observer sensitivities and non-replicated under the fixed legacy transport rule;
- shared versus species-specific geographic mapping remains unresolved;
- exposure coupling of measured white remains a stated limitation.


## Raw-image reproducibility boundary

The archive does not contain the original iNaturalist image pixels because the prospective protocol explicitly deleted pixels and masks after partition sealing. It does preserve the frozen photo IDs/source URLs, per-row image SHA256 values, exact acquisition/measurement code, ROI-v4 detector byte and EfficientSAM weights. Reacquired bytes can therefore be validated while still available, but exact historical pixel reconstruction is not guaranteed if the external provider no longer serves those bytes.
