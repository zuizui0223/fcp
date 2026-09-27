# FCP reproducibility audit — 2026-09-27

## Bottom line

The active New Phytologist paper has **strong result-level reproducibility**, but the original archive was not fully self-contained at the image-measurement layer.

This audit separates four reproducibility levels instead of treating them as one claim.

| Level | Question | Status after this revision |
|---|---|---|
| R1 | Can the exact frozen inputs and machine-readable outputs be recovered? | **PASS** |
| R2 | Can the headline H1/H2 results be recomputed from the frozen measured tables with the historical execution code? | **PASS, CI-verified** |
| R3 | Can the image-measurement machine itself be reconstructed without consulting historical Git branches? | **PASS for code/model/runtime bytes** |
| R4 | Can every original image pixel be reproduced bit-for-bit from the archive alone? | **NO — explicitly out of scope because source image pixels were not persisted** |

## R1 — frozen inputs and outputs

The provenance package contains or checksum-recovers:

- legacy discovery measured table — SHA256 `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`;
- legacy reserve measured table — SHA256 `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`;
- prospective third-cohort measured table — SHA256 `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186`;
- exact WorldClim BIO and SRAD archives and raster-level hashes;
- exact H3a/H3b/highlight intermediate files formerly stored only in expiring Actions artifacts;
- all headline result JSON/CSV/TSV files and the figure source freeze.

Machine-readable frozen results remain authoritative over manuscript prose.

## R2 — measured-table to headline-result recomputation

The exact historical execution tree at commit
`7e538e5c51c05a7cc47b2fcf53eea92634c8a863`
is now copied into the provenance package under
`historical/execution_source/third_cohort/`.

The workflow `.github/workflows/verify-np-core-reproduction-20260927.yml` independently creates a clean worktree at that historical commit, restores the frozen legacy measured tables, reruns:

1. the first-frozen 200-partition observer-disjoint H1 reliability analysis; and
2. the prospective third-cohort white-axis H2 analysis with its 999 structured-null worlds.

The rerun is compared with the frozen result objects and null/vector arrays. The verification is not allowed to change the manuscript or frozen result directories.

## R3 — image-measurement machine reconstruction

The provenance package now carries an exact measurement-source snapshot from commit
`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`, including:

- ROI-v4 source;
- ROI-v4 runtime;
- fixed palette measurement code;
- measurement firewall/reassembly code;
- acquisition / blind-measurement / sealing scripts;
- frozen contracts;
- `requirements-atlas-roi-v4.txt`;
- the ROI-v4 YOLO detector weight, checked against SHA256
  `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`.

The exact EfficientSAM revision and both ONNX files are also mirrored into the package:

- revision `d525f622e6f640acf5a0fc37c7ca1f243da5bde0`;
- encoder SHA256 `84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951`;
- decoder SHA256 `a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11`.

The exact prospective authorized metadata table is included with SHA256
`13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290`.
It contains the frozen photo IDs/source URLs used by the one-shot acquisition.

## R4 — raw-image byte limitation

The prospective protocol intentionally deleted image pixels and flower masks after each sealed partition. Therefore the archive cannot by itself regenerate the original source image bytes if iNaturalist later removes or changes those files.

This limitation is bounded rather than hidden:

- every measured row retains the iNaturalist photo ID;
- the frozen prospective source metadata retains the authorized source URL;
- every successfully measured row retains `image_sha256`;
- the exact acquisition and measurement code and all model weights are archived.

A future reacquisition can therefore be checked against the recorded image hash. If the provider no longer serves the original bytes, the exact pixel-to-measurement computation is **not recoverable** from this archive. The paper must not claim otherwise.

## Spatial-analysis status

The previous package carried reporting receipts for the D–spatial result. This revision additionally packages the exact historical spatial-analysis scripts and the frozen Step 5/6/8/9 result directories from source commit
`f14186590c11ac24c95e1985077908b732132e96`.

This improves auditability without upgrading the causal interpretation: the D–spatial association remains a structural correlate.

## Reproducibility ceiling

The correct claim is:

> The frozen measured tables, headline inferential analyses, null distributions, figures and environmental inputs are reproducible from the archived package. The image-measurement implementation and model weights are archived, and exact source-photo identities/hashes permit byte validation when the source files remain available. Because raw source image pixels were intentionally not persisted, the archive does not guarantee future bit-for-bit reconstruction of the image-acquisition stage.

That is the reproducibility statement to use in the manuscript / Data Availability / Zenodo record.
