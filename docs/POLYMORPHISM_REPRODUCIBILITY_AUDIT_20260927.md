# FCP reproducibility audit — 2026-09-27

## Conclusion

The active New Phytologist paper has demonstrated computational replay for its core H1/H2 and D–spatial analyses and for the named BIO5/H3 ecological filters. Supporting controls that are not part of those replay workflows are byte-frozen with explicit provenance boundaries. The outcome-blind sampling lineage is preserved from the 42,111-species opportunity frame through the prospective cohort.

Two historical replay checks were executed using exact historical input tables and analysis code.

Run `36324187508` independently reproduced:

- H1 repeated observer-disjoint reliability;
- H1 deterministic stress test;
- legacy H2 white-axis target localization; and
- prospective third-cohort H2.

All four JSON result objects matched the frozen outputs at absolute floating-point tolerance `1e-12`.

Run `36373083779` additionally replayed the ecological D–spatial chain from the original 999-permutation discovery/reserve null arrays:

- discovery Step 5;
- reserve Step 5b;
- sampled-span adjusted Step 6;
- technical-failure adjusted Step 8; and
- four-state ambiguity endpoint Step 9.

Step 5/5b reproduced exactly at `1e-12`; the complete chain passed a numerical identity gate of absolute `1e-10` / relative `1e-12`, chosen only to tolerate machine-level geodesic floating-point differences (the observed discrepancy that motivated it was ~3e-12 km).

Secondary replay was independently verified across runners. Run `36379690691` first reproduced the third-cohort BIO5/BIO14/SRAD models, observer-paired/balanced sensitivities, legacy BIO5 transport, and H3a/H3b. A later runner revealed only optimizer-level conditional-logit floating-point variation: one observed CI bound differed by about `1.14e-10`. The stable cross-run contract therefore uses absolute `1e-8` / relative `1e-10` for BIO5-family numerical fields. Run `36380786248` passes that contract. H3a/H3b rerun under the recorded R/package versions; their archived permutation-null, summary and PGLS tables reproduce byte-for-byte.

This audit distinguishes that demonstrated numerical reproducibility from the stronger question of whether the original third-party image bytes can always be reconstructed.

## Four reproducibility levels

| Level | Question | Status |
|---|---|---|
| R1 | Can exact frozen inputs and outputs be recovered? | **PASS** |
| R2 | Can the named core and secondary analysis chains be recomputed from frozen evidence? | **PASS — H1/H2, D–spatial, BIO5 sequence and H3 independently replayed** |
| R3 | Can the image-measurement implementation be reconstructed from the archive? | **PASS after byte-completion revision** |
| R4 | Can every original source image pixel be regenerated from the archive alone? | **NO — intentionally not claimed** |

## Replay coverage matrix

| Component | Permanent inputs/code | Independent replay in this audit | Current status |
|---|---|---|---|
| H1 repeated + deterministic stress | yes | yes | PASS |
| Legacy H2 white-axis localization | yes | yes | PASS |
| Prospective third-cohort H2 | yes | yes | PASS |
| D–spatial Steps 5/5b/6/8/9 | yes, including original 999-permutation null families | yes | PASS |
| Third-cohort BIO5/BIO14/SRAD + observer sensitivities | yes | yes | PASS at abs 1e-8 / rel 1e-10 |
| Legacy BIO5 transport | yes | yes | PASS at abs 1e-8 / rel 1e-10 |
| H3a phylogenetic signal | yes | yes | PASS; output SHA256 reproduced |
| H3b sampled-span replication | yes | yes | PASS; output SHA256 reproduced |
| Direct highlight-control reacquisition | technical table/join/high-clip outputs and analysis code are permanent | statistical artifact preserved; original pixel reacquisition is not archive-self-contained | BOUNDED |
| Fresh-image D transport imported into H1 | compact receipt, source result and transport artifact are permanent | not independently rerun inside this manuscript repository | FROZEN PROVENANCE |
| Discussion-only PAL/WAL maintenance reanalyses | scripts/results versioned | not part of the primary replay contract | SUPPORTING ONLY |

This matrix is the ceiling on the phrase "replayed": it applies only to rows marked **yes** in the independent-replay column.

## R1 — exact evidence bytes

The provenance package contains or checksum-verifies the full sampling-to-inference chain, including:

- complete 42,111-species outcome-blind opportunity frame, SHA256 `5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`;
- P100 3,730-species parent pool and frozen P500, permanently recovered from original artifact `10302477571`;
- reconstructed 3,230-species third-cohort candidate frame, SHA256 `7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e`;
- four exact prior-ID exclusion inputs, selected 500-species manifest, candidate 49,999-row metadata and authorized 49,900-row metadata;
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

This is direct evidence of computational reproducibility rather than only provenance bookkeeping. The replay workflow is now pinned to Python 3.12.14 with `requirements-np-replay-20260928.txt`, matching the package versions used in the successful replay.

## R2b — demonstrated secondary ecological replay

Run `36379690691` used the permanent provenance release as the evidence source.

### Environmental sequence

The replay used Python 3.11.16 and the exact successful pip environment frozen in `requirements-np-bio5-replay-20260928.txt`.

It reran:

- the third-cohort white/non-white BIO5, BIO14 and mean-SRAD analysis;
- observer-paired and observer-balanced sensitivity analyses;
- the legacy discovery/reserve BIO5 species-disjoint transport test.

The stable cross-run gate is absolute `1e-8` / relative `1e-10`, with exact non-numeric identity. A first runner passed a tighter `1e-10` gate, but a second runner exposed a harmless ~`1.14e-10` conditional-logit CI difference despite identical pinned Python/package versions; the wider gate prevents CPU/BLAS-level optimizer noise from being mistaken for scientific non-reproducibility.

### H3 alternative-explanation filters

The replay installed the package versions recorded inside the original H3 result JSONs:

- H3a: R 4.6.1, ape 5.8.1, phytools 2.5.2, jsonlite 2.0.0;
- H3b: R 4.6.1, ape 5.8.1, phylolm 2.6.5, jsonlite 2.0.0.

The H3 result JSONs matched at absolute `1e-9` / relative `1e-11`. More strongly, the H3a permutation-null/summary tables and H3b permutation-null/summary/PGLS tables reproduced their archived SHA256 values byte-for-byte.

Thus the H3 negative results and failed legacy BIO5 transport are computationally reproduced, not merely copied historical decisions.

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

> The named H1/H2, D–spatial, BIO5 and H3 analysis chains are computationally replayable from the archived evidence package under their recorded numerical tolerances, while supporting controls retain the provenance status shown in the coverage matrix. The exact image-measurement code and model weights are archived, and frozen source-photo identities plus per-row image hashes permit byte validation when original provider bytes remain available. Because raw source image pixels were intentionally not persisted, future bit-for-bit reconstruction of the original image-acquisition stage is not guaranteed.

No stronger raw-image claim is needed for the paper.
