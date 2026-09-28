# Polymorphism paper data-lineage map — 2026-09-25

## Purpose

This document gives one reader-facing route from every main inferential claim in the active New Phytologist manuscript to the protocol, machine-readable result, exact input source, immutable Git commit and/or GitHub Actions artifact needed to recover the analysis.

It does not recompute any biological result.

## Reader-facing data architecture

The paper is easiest to read as one inferential chain rather than as a history of individual executions:

1. **Global opportunity frame — 42,111 species.** Metadata-only iNaturalist discovery defines which species could have entered before flower colour was examined. This step is necessary to keep candidate-species assembly outcome-blind; it is not a flower-colour dataset or a prevalence denominator.
2. **Discovery + species-disjoint validation — 1,000 species / 100,000 photographs.** Fixed high-depth sampling estimates within-species distributions, while the species split separates pattern discovery from out-of-sample validation. This step is necessary because one or a few photographs cannot define a within-species distribution and because discovery on the same species is not confirmation.
3. **Fresh-image D transport — 136 overlapping evaluable species.** New photo IDs under the same measurement system test whether D transports beyond the original observer partitions. This step is necessary because observer splitting alone does not test a newly sampled image set.
4. **Prospective confirmation — 499 species / 49,900 new photographs.** After the white-versus-nonwhite axis is discovered and frozen, previously used high-depth species are excluded and a new outcome-blind species/photo cohort tests the fixed target. This step is necessary because a data-derived target cannot be prospectively confirmed in the data that generated it.
5. **Technical and explanatory annotations.** Highlight metrics, WorldClim variables, phylogenetic placements and sampled span are attached to existing image cohorts. These steps are necessary to test measurement coupling and alternative explanations; they are not additional independent flower-colour cohorts and cannot retroactively upgrade the prospective H2 test.

The prospective confirmation cohort is one physical 499-species / 49,900-row measurement dataset used in two chronological phases:
- first, untouched prospective H2 confirmation;
- only after H2 terminalization, post-confirmatory highlight-validity and environmental analyses.

## Traceability grades

- A — protocol and result are on main; exact input is on main or recoverable from an immutable Git commit with a frozen hash.
- B — one or more essential intermediate inputs remain available only through a time-limited artifact. No current headline analysis remains in this class after the 2026-09-25 artifact freeze.

## Claim-to-data routing table

| Paper step | Physical cohort/data | Protocol / decision contract | Canonical result | Exact data route | Grade |
|---|---|---|---|---|---|
| Global opportunity frame | 42,111 metadata-discovered species | docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md plus the upstream frozen protocols named there | reporting provenance note | frozen species CSV SHA256 5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc; exact CSV is copied into the permanent provenance package under raw_inputs/opportunity_frame | A |
| Discovery/validation high-depth resource | discovery 500×100 + validation 500×100 raw measured rows | acquisition lineage in manuscript/SI; H1 protocol verifies exact hashes | reused by H1/H2/spatial/H3 | immutable commit 5142f7951af0dde5364bb047a566d67e8c479e51; discovery SHA256 ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4; validation SHA256 0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6; exact repository paths are preserved in the provenance snapshot | A |
| H1 observer-disjoint D reproducibility | discovery/validation cohorts | docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md | results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json | discovery/validation resource above; result JSON repeats both SHA256 fingerprints | A |
| D finite-sample sensitivity | discovery/validation D plus previously frozen D–spatial/H3 sensitivities | raw Gini–Simpson D retained as primary; corrected D = D n/(n-1) | results/polymorphism_D_finite_sample_sensitivity_20260928/result.json | frozen Step-8 result mirrored under archive/fcp_submission_20260925/d_finite_sample; species-level D/corrected-D table and H3 summaries are checksum-stable archive inputs | A |
| Discovery/validation H2 white-axis localization | discovery/validation cohorts | docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md | results/polymorphism_white_axis_targeted_test_20260912/result.json | same frozen discovery/validation measured tables; named axis is explicitly retrospective in these cohorts | A |
| Prospective confirmation cohort selection | new species, no biological outcomes | frozen prospective-selection contract; exact repository identity in provenance snapshot | machine-readable selection receipt | prior high-depth exclusion artifact 10302477571 is copied permanently; component hashes 1473aad6… / f54d07fb…; 3,230-species candidate frame SHA256 7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e; selected manifest SHA256 16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59 | A |
| Prospective confirmation measurement | same prospective confirmation cohort | frozen prospective-measurement contract; exact repository identity in provenance snapshot | machine-readable measurement receipt | immutable result commit 7e538e5c51c05a7cc47b2fcf53eea92634c8a863; measured-table SHA256 57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186; final artifact 10496492307, digest sha256:319c040aaffc30d3cd97dcbcc217409e75010ec0f67bcc86c6bb82d8275779c7 | A |
| Prospective H2 white-axis confirmation | same prospective confirmation cohort, before later mechanism work | same frozen prospective contract and claim freeze | machine-readable H2 result | reads the measured table above after support gate; result points back to the measurement receipt | A |
| Direct highlight validity | same prospective confirmation cohort, reacquired post-H2 | frozen highlight-validity contract; exact repository identity in provenance snapshot | machine-readable validity result | run 35687421896 / artifact 10709490106; exact technical table, sealed join key, high-clip IDs and summary are permanently archived with per-file SHA256 | A |
| Secondary BIO5/BIO14/solar analysis | same prospective confirmation cohort, after H2 terminalization | docs/POLYMORPHISM_WHITE_ENVIRONMENT_MECHANISM_PROTOCOL_20260925.md | results/polymorphism_white_environment_mechanism_20260925/result.json and result_receipt.json; observer-sensitivity receipt at results/polymorphism_white_environment_observer_sensitivity_20260925/result.json | exact measured rows, highlight inputs and WorldClim bytes are in the permanent provenance package; original full result artifact is permanently copied; run 36379690691 replayed the models and observer sensitivities under the frozen Python environment and matched the original outputs | A |
| BIO5 species-disjoint transport | discovery/validation cohorts | frozen one-predictor transport contract; exact repository identity in provenance snapshot | machine-readable transport result | exact discovery/validation rows and WorldClim BIO bytes are permanent; run 36379690691 replayed the fixed test and reproduced the transport result/table | A |
| PAL/WAL persistence reanalysis | published Del Valle et al. (2019) Supplementary Tables S2 and S1 | descriptive table-reanalysis rules; no causal or meta-analytic upgrade | results/polymorphism_silene_decoupling_persistence_20260925/result.json; results/polymorphism_crossspecies_pal_wal_frequency_20260925/result.json | source tables and reporting scripts are preserved in the provenance snapshot; results are deterministic descriptive summaries | A |
| PAL/WAL × FCP bridge feasibility | 25 unique species from the Del Valle Table S1 registry × current discovery/reserve D-eligible species and 500 frozen prospective selections | exact-name overlap audit only; no inferential bridge test | results/polymorphism_pal_wal_h2_overlap_20260928/result.json | only two reserve overlaps, both PAL; discovery 0, prospective 0, WAL 0, so the direct bridge is explicitly not estimable | A |
| D–spatial organization | discovery/validation cohorts | frozen upstream analysis; current main carries reporting-only synthesis | results/polymorphism_spatial_organization_clue_20260918/result.json | source commit f14186590c11ac24c95e1985077908b732132e96; exact Step 5/6/8/9 scripts/intermediates and the original discovery/validation 999-permutation null artifact families are copied into the permanent provenance package | A |
| H3a phylogenetic signal | validation-cohort D + frozen S1/S2/S3 trees | frozen phylogenetic-signal protocol | machine-readable frozen result | signal/tree/covariate artifacts 10292218669 / 10292662493 / 10292767459 are permanently archived; exact source identity is in provenance; run 36379690691 replayed H3a under the recorded R/package versions and reproduced the permutation-null/summary SHA256 values byte-for-byte | A |
| H3b sampled-span replication | validation-cohort D + frozen pre-outcome sampled-span panel | frozen sampled-span replication protocol | machine-readable frozen result | artifact 10292399238 is permanently archived; exact source is frozen at commit ff1eac2ec1378465d5027184d507cb597c67be27; run 36379690691 replayed H3b and reproduced the 20,000-permutation null, summary and PGLS SHA256 values byte-for-byte | A |

## What is and is not physically stored on main

The repository deliberately does not carry every large measured-photo CSV on current main.

That is acceptable only because the paper now records exact recovery routes:
- discovery/validation cohort rows: immutable Git commit + SHA256 + bundled copy in the provenance snapshot;
- prospective H2 rows: immutable biological result commit + SHA256 + final workflow artifact + bundled copy;
- spatial results: immutable source commit + exact frozen scripts/intermediates + permanently copied discovery/validation 999-permutation null artifacts;
- H3a/H3b preflight panels: workflow run/artifact IDs and hashes.

A filename in a protocol should therefore never be interpreted as meaning that the file is currently present on main.

## External-byte archival status

### Resolved: expiring Actions artifacts

The exact highlight, H3a and H3b intermediate files were copied byte-for-byte from their original GitHub Actions artifacts into archive/fcp_submission_20260925. ARCHIVE_MANIFEST.md records the source run/artifact IDs and per-file SHA256 values. Their reproducibility no longer depends on Actions retention.

### Resolved: WorldClim provider dependence

The exact WorldClim 2.1 10-arc-minute input archives used by the environmental analyses are now mirrored as GitHub Release assets under tag:

fcp-worldclim-2.1-10m-20260925

Release assets:
- wc2.1_10m_bio.zip — 49,869,449 bytes — SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf;
- wc2.1_10m_srad.zip — 16,233,364 bytes — SHA256 c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026.

The release receipt is archive/fcp_submission_20260925/WORLDCLIM_RELEASE_RECEIPT.md. The committed checksum manifest also contains BIO5, BIO14 and all 12 monthly SRAD TIFF hashes, and the analysis workflows hard-fail on checksum mismatch.

The release itself is operationally mutable, but the canonical byte identity is immutable at the manuscript level because the required SHA256 values are committed in Git and checked before analysis. A replaced release asset would therefore be detected.

No headline analysis now depends solely on an expiring Actions artifact or on a remote provider serving unchanged bytes. The original discovery/validation spatial-null artifact families from runs 34088925008 and 34178957447 are also copied into the provenance snapshot before their December 2026 expiry, with later rebuilds falling back to the permanent release copy.

## Unified self-contained provenance snapshot

Readers who do not want to reconstruct the chain from multiple immutable Git objects can use the packaged provenance release:

- release tag: `fcp-np-provenance-20260926`;
- release URL: `https://github.com/zuizui0223/fcp/releases/tag/fcp-np-provenance-20260926`;
- asset: `fcp-np-provenance-20260926.tar.gz`;
- per-file checksum manifest inside the asset: `FILE_SHA256SUMS.txt`;
- canonical current-asset receipt: `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`.

The Git-tracked receipt is authoritative for the **currently published** package source commit, asset byte count, SHA256 and packaged-file count. The snapshot is regenerated when the active manuscript/evidence surface changes, so those package-instance values are intentionally not duplicated here. Earlier receipt revisions in Git history preserve prior package identities.

The package includes the exact 42,111-species opportunity frame, all four prior-ID exclusion sources used for the prospective cohort, the exact selected/candidate/authorized prospective-cohort metadata, discovery/validation cohorts measured tables, the prospective-H2 measured table, frozen spatial scripts/intermediates and original null artifacts, the fresh-D source receipt/artifact, permanent highlight/H3 intermediate archives, checksum-pinned WorldClim BIO/SRAD archives, image-measurement code/model bytes, and the active manuscript/SI/figures, protocols, machine-readable results and regression guards.

This package is a reproducibility snapshot of the current evidence chain, not a claim that author metadata or journal submission status is final.

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
- the BIO5 transport test remains unsupported under its frozen rule;
- H3a/H3b negative results do not establish equivalence or absence of all phylogenetic/geographic effects.
