# Polymorphism paper data-lineage map — 2026-09-25

## Purpose

This document gives one reader-facing route from every main inferential claim in the active New Phytologist manuscript to the protocol, machine-readable result, exact input source, immutable Git commit and/or GitHub Actions artifact needed to recover the analysis.

It does not recompute any biological result.

## How to read the paper in one pass

The paper uses four distinct inferential layers.

1. Legacy discovery/reserve resource — the same two high-depth species-disjoint measured-photo tables support H1, legacy H2 target localization, spatial organization, and H3.
2. Prospective H2 cohort — a new species-disjoint cohort selected after q_white was frozen provides the untouched H2 confirmation.
3. Secondary analysis of that prospective H2 cohort — only after H2 was terminalized, the same measured rows were reused for highlight validity and the pre-specified environmental/BIO5 follow-up. This is not part of the untouched H2 confirmation.
4. Historical/frozen diagnostics — P500, spatial Step-8/9 receipts, H3 trees/covariates and permutation artifacts constrain interpretation but have their own inferential roles.

The phrase third cohort therefore refers to one physical 499-species / 49,900-row measurement cohort, but two different chronological uses:
- first: untouched prospective H2 confirmation;
- later: post-H2 secondary validity/environmental analyses.

## Traceability grades

- A — protocol and result are on main; exact input is on main or recoverable from an immutable Git commit with a frozen hash.
- A- — all numerical results and exact source Git objects are frozen, but main carries a reporting receipt rather than the complete original analysis tree.
- B — exact result/protocol is traceable, but one or more essential intermediate inputs are stored only in a time-limited Actions artifact or an external versioned archive and should be copied to the permanent submission archive.

## Claim-to-data routing table

| Paper step | Physical cohort/data | Protocol / decision contract | Canonical result | Exact data route | Grade |
|---|---|---|---|---|---|
| Global opportunity frame | 42,111 metadata-discovered species | docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md plus upstream RGFCA protocols named there | reporting provenance note | frozen species CSV data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv, SHA256 5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc, retained in upstream Git history | A- |
| Legacy measured resource | discovery 500×100 + reserve 500×100 raw measured rows | acquisition lineage in manuscript/SI; H1 protocol verifies exact hashes | reused by H1/H2/spatial/H3 | immutable commit 5142f7951af0dde5364bb047a566d67e8c479e51: data/derived/global_monte_carlo_measured_photos_v1.csv SHA256 ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4; data/derived/rgfca_reserve_replication_measured_photos_v1.csv SHA256 0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6 | A |
| H1 observer-disjoint D reproducibility | legacy discovery/reserve | docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md | results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json | legacy measured resource above; result JSON repeats both SHA256 fingerprints | A |
| Legacy H2 white-axis localization | legacy discovery/reserve | docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md | results/polymorphism_white_axis_targeted_test_20260912/result.json | same frozen legacy measured tables; named axis is explicitly retrospective in these cohorts | A |
| Prospective H2 cohort selection | new species, no biological outcomes | docs/POLYMORPHISM_H2_THIRD_COHORT_SELECTION_PROTOCOL_20260916.md | results/polymorphism_h2_third_cohort_selection_20260916/result.json and selected_species_manifest.tsv | selected-manifest SHA256 16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59; selection-freeze commit 2768b2dd0f8baf4ed1185a1b64c1060768b36c00; candidate-universe CSV hash 7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e | A |
| Prospective H2 measurement | same new cohort | docs/POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md | results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json | immutable result commit 7e538e5c51c05a7cc47b2fcf53eea92634c8a863; measured table data/derived/polymorphism_h2_third_cohort_measured_photos_v1.csv, frozen SHA256 57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186; final artifact 10496492307, digest sha256:319c040aaffc30d3cd97dcbcc217409e75010ec0f67bcc86c6bb82d8275779c7 | A |
| Prospective H2 white-axis confirmation | same new cohort, before later mechanism work | same prospective protocol and result/claim freeze | results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json; docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md | reads the measured table above after support gate; result points back to the measurement receipt | A |
| Direct highlight validity | same prospective H2 cohort, reacquired post-H2 | docs/POLYMORPHISM_H2_THIRD_COHORT_HIGHLIGHT_VALIDITY_PROTOCOL_20260922.md | results/polymorphism_h2_third_cohort_highlight_validity_20260922/result.json | technical-seal workflow 35687421896; artifact 10709490106, digest sha256:adbca1112ba0635623bd056e21286a369d71f9a1cddbe38a920262d4630c64e3 | B |
| Secondary BIO5/BIO14/solar analysis | same prospective H2 cohort, after H2 terminalization | docs/POLYMORPHISM_WHITE_ENVIRONMENT_MECHANISM_PROTOCOL_20260925.md | results/polymorphism_white_environment_mechanism_20260925/result.json and result_receipt.json | biological artifact from run 35177668182; highlight seal from run 35687421896; WorldClim 2.1 10-arc-minute BIO and SRAD archives downloaded by .github/workflows/white-environment-mechanism-20260925.yml | B |
| BIO5 species-disjoint transport | legacy discovery/reserve | docs/POLYMORPHISM_LEGACY_WHITE_BIO5_REPLICATION_PROTOCOL_20260925.md | results/polymorphism_legacy_white_bio5_replication_20260925/result.json | exact legacy rows recovered with git show from source commit 5142f7951af0dde5364bb047a566d67e8c479e51 and hashes verified; WorldClim 2.1 BIO5 downloaded by frozen workflow | B |
| D–spatial organization | legacy discovery/reserve | frozen upstream analysis; current main carries reporting-only synthesis | results/polymorphism_spatial_organization_clue_20260918/result.json | source branch feat/polymorphism-paper-v0-1-post-step9, head f14186590c11ac24c95e1985077908b732132e96; source result paths and Git blob SHAs are enumerated inside the reporting receipt | A- |
| H3a phylogenetic signal | legacy D + frozen S1/S2/S3 trees | docs/POLYMORPHISM_H3A_PHYLOGENETIC_SIGNAL_PROTOCOL_20260912.md | results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json | source commit db2514f622468747a7d5896c2696fcbf13729ef6; run 34677042793 artifact 10292218669; tree preflight run 34676531258 artifact 10292662493; covariate preflight run 34676855989 artifact 10292767459 | B |
| H3b sampled-span replication | legacy D + frozen pre-outcome sampled-span panel | docs/POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md | results/polymorphism_h3b_reserve_span_20260912/result.json; docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md | execution head ff1eac2ec1378465d5027184d507cb597c67be27; run 34677468362; artifact 10292399238, digest sha256:34c5e646725f6865e313b17f1b70f2471db8169f443fb0649da83044d654386c; uses H3a pre-outcome covariate/tree artifacts and exact legacy row hashes | B |
| P500 measurement transport / failed H2 terminalization | separate prospective expansion | P500 frozen protocols and postmortem | docs/P500_PROSPECTIVE_TERMINAL_POSTMORTEM_20260916.md | complete measurement chain retained in postmortem; no durable biological H2 result exists and no replay may inherit prospective status | A- |

## What is and is not physically stored on main

The repository deliberately does not carry every large measured-photo CSV on current main.

That is acceptable only because the paper now records exact recovery routes:
- legacy discovery/reserve rows: immutable Git commit + file path + SHA256;
- prospective H2 rows: immutable biological result commit + file path + SHA256 + final workflow artifact;
- spatial results: immutable source commit + source result blob IDs;
- H3a/H3b preflight panels: workflow run/artifact IDs and hashes.

A filename in a protocol should therefore never be interpreted as meaning that the file is currently present on main.

## Remaining archival risks before submission

### 1. Time-limited GitHub Actions artifacts

Some exact intermediate objects are currently preserved only as Actions artifacts:
- third-cohort highlight technical seal;
- H3a S1/S2/S3 tree/preflight outputs;
- H3a pre-outcome sampling-opportunity panel;
- full H3b permutation-null and PGLS outputs.

These artifacts currently have retention expiries. Before final submission/Zenodo release, copy the exact artifact contents into the permanent archive or another non-expiring versioned store and record their SHA256 values.

### 2. WorldClim archive bit identity

The BIO5 analyses pin WorldClim 2.1, 10 arc-minute and the exact download URLs in the committed workflows. The original workflows did not record a SHA256 for the downloaded WorldClim zip archives. The scientific version is traceable, but a future bit-for-bit replay still depends on the provider serving the same v2.1 files.

Before archival release, retain the exact BIO and SRAD archives used for the frozen analyses, or record their archive/file checksums in the permanent release manifest.

### 3. H3a/H3b artifact-only intermediate panels

The compact decisions are on main, but the exact tree files and pre-outcome sampled-span panel are artifact-backed. These should be archived with the submission bundle.

## Reader shortcut

For a reader who wants to verify only the paper's headline logic:
1. Start at docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md.
2. Use Table 1 to identify the physical cohort.
3. Use this lineage map to identify the inferential execution.
4. Open the canonical result JSON.
5. Follow its protocol/source commit/artifact entry to the exact input.
6. If prose and a frozen machine-readable result disagree, the machine-readable result controls.

## Claim boundary

Traceability does not change inferential status.

In particular:
- the prospective H2 test remains untouched and precedes all later BIO5 work;
- BIO5 remains a pre-specified secondary post-H2 analysis, not part of the prospective H2 confirmation;
- legacy BIO5 transport remains failed under its frozen rule;
- H3a/H3b negative results do not establish equivalence or absence of all phylogenetic/geographic effects.
