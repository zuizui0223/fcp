# FCP reproducibility audit — 2026-09-27

## Conclusion

The active New Phytologist paper is reproducible at the level required to regenerate its frozen statistical claims from the archived measured tables.

A separate historical replay was executed in GitHub Actions (run `36324187508`) using the exact historical input tables and analysis code. It independently reproduced:

- H1 repeated observer-disjoint reliability;
- H1 deterministic stress test;
- legacy H2 white-axis target localization; and
- prospective third-cohort H2.

All four replayed JSON result objects matched the current frozen outputs with absolute floating-point tolerance `1e-12`.

This audit distinguishes that demonstrated numerical reproducibility from the stronger question of whether the original third-party image bytes can always be reconstructed.

## Four reproducibility levels

| Level | Question | Status |
|---|---|---|
| R1 | Can exact frozen inputs and outputs be recovered? | **PASS** |
| R2 | Can headline H1/H2 analyses be recomputed from frozen measured tables? | **PASS — independently replayed** |
| R3 | Can the image-measurement implementation be reconstructed from the archive? | **PASS after byte-completion revision** |
| R4 | Can every original source image pixel be regenerated from the archive alone? | **NO — intentionally not claimed** |

## R1 — exact evidence bytes

The provenance package contains or checksum-verifies:

- legacy discovery measured rows, SHA256 `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`;
- legacy reserve measured rows, SHA256 `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`;
- prospective third-cohort measured rows, SHA256 `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186`;
- exact highlight/H3 intermediate files formerly held only in expiring Actions artifacts;
- exact WorldClim BIO/SRAD archives and raster hashes;
- frozen result JSON/CSV/TSV files, figure sources, figures and claim guards.

Machine-readable frozen outputs remain authoritative over manuscript prose.

## R2 — demonstrated measured-table to result replay

The active workflow `.github/workflows/polymorphism-reproducibility-replay.yml` recovers the exact historical execution code and frozen measured tables and reruns the four headline analyses.

Successful replay run `36324187508` reported:

- `H1 repeated: EXACT_JSON_MATCH_WITH_1E-12_FLOAT_TOLERANCE`
- `H1 stress: EXACT_JSON_MATCH_WITH_1E-12_FLOAT_TOLERANCE`
- `legacy H2 target: EXACT_JSON_MATCH_WITH_1E-12_FLOAT_TOLERANCE`
- `prospective H2: EXACT_JSON_MATCH_WITH_1E-12_FLOAT_TOLERANCE`

This is direct evidence of computational reproducibility rather than only provenance bookkeeping.

## R3 — image-measurement implementation bytes

The provenance package preserves the historical measurement source from commit
`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`, including ROI-v4, fixed-palette measurement, acquisition, firewall, blind measurement and sealing code.

This byte-completion revision additionally embeds:

- ROI-v4 detector weight, SHA256
  `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`;
- EfficientSAM encoder, revision `d525f622e6f640acf5a0fc37c7ca1f243da5bde0`, SHA256
  `84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951`;
- EfficientSAM decoder, same revision, SHA256
  `a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11`;
- frozen third-cohort authorized metadata, SHA256
  `13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290`;
- frozen authorized-species table, SHA256
  `36a866b040b835e20539d318b0533bcb905cbe760d23df29e7da6587a054e593`.

The authorized metadata contain the exact photo IDs and source URLs used by the one-shot measurement run. The measured table preserves per-row `image_sha256` values.

The package also includes a fuller frozen D–spatial execution tree from commit
`f14186590c11ac24c95e1985077908b732132e96`, including Step 5/6/8/9 scripts and intermediate result directories.

## R4 — raw source image limitation

The prospective protocol deliberately deleted candidate image pixels and flower masks after each partition was sealed. Raw iNaturalist image bytes are therefore not stored in the archive.

Consequences:

- if iNaturalist still serves the original bytes, a future reacquisition can be checked against the recorded `image_sha256`;
- if the provider changes or removes those bytes, the exact historical pixel-to-measurement step cannot be rerun bit-for-bit from this archive alone;
- this does **not** affect replay of the paper's numerical analyses from the frozen measured tables.

The manuscript and archive must not describe the project as permanently self-contained at the raw-image-byte level.

## Reproducibility ceiling

The defensible statement is:

> The frozen measured tables, headline inferential analyses, null-model outputs, figures and environmental inputs are computationally reproducible from the archived evidence package. The exact image-measurement code and model weights are archived, and frozen source-photo identities plus per-row image hashes permit byte validation when original provider bytes remain available. Because raw source image pixels were intentionally not persisted, future bit-for-bit reconstruction of the original image-acquisition stage is not guaranteed.

No stronger raw-image claim is needed for the paper.
