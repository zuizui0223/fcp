# Fresh FCP scale-up — execution status (2026-10-15)

**Scientific endpoint:** cross-species commonality of *within-species* flower-colour variation, primarily the fixed white–chromatic axis `W_i=(u_i·q_white)^2`, in relation to environmental heterogeneity *within each species*. Mean flower colour among species is the separate `mois` question.

## Completely new biological sample, not an augmented old table

This branch is cut from `main` and contains only the fresh-data acquisition/re-measurement workflow, its provenance/biology guards, and tests. Historical measured photo CSVs, old 149,900 photo colour labels, old W species vectors and retrospective climate association results are **not** eligible inputs for this study. We use the historical source-code and pre-trained model *method* frozen at `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a` and its hashes. A previously photographed original photo ID may be acquired again; its old result is **never reused**, nor may it be claimed as an independent biological replication.

- Initial outcome-blind taxon opportunity frame: 42,111 taxa, fixed source table SHA256 `5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`.
- No `100 photographs per species` restriction and no `473,000-photo` global cap.
- All provider-frozen eligible, licensed photographs are traversed with an ID cursor. Every source photo ID receives one terminal status, including acquisition failures and ROI/segmentation uncertainty; no outcome-based replacement.
- Species with `n_classifiable < 40` remain explicit insufficient-support records. W eligibility follows the original paper's `0.10` coarse+continuous gates and a `0.20` sensitivity, without choosing a new colour direction.
- Source photos are remeasured **location blind** with original ROI-v4 + EfficientSAM + fixed palette. `scripts/acquisition/bootstrap_fcp_original_roi_v4.sh --full` restores and verifies the historical method source, detector and EfficientSAM model bytes before any fresh image is opened.
- Fresh metadata completion across every included taxon is mandatory before the full-sample pixel opening. The existing metadata pilot remains *incomplete*; do not turn it into a measured finding.

## Evidence completed before this branch

A bounded real **metadata-only** provider pilot was executed in [GitHub Actions run 38589322067](https://github.com/zuizui0223/fcp/actions/runs/38589322067), artifact ID `11906662449`, with frozen source date `2026-10-14`.

- One outcome-blind species: *Ribes californicum*, iNaturalist taxon ID 47129.
- Ten 200-observation pages attempted, ending at cursor `201872254`.
- **3,303 licensed source photo metadata rows**, in an **in-progress** source manifest, not 3,303 newly classified flower photos.
- Its source `photo_manifest.ndjson` checksum is `e3c5482c489081c4ea593f278f52166173049eec1e6c06d892a75983e0fb9c2b`.
- **0 fresh colour measurements** and **0 complete source species**. Other 42,110 taxa had not yet been accessed in that pilot.

The pilot demonstrated API paging and persistence, not completion of the global source universe. The artifact is checkpoint material, not scientific results. It may be advanced only with the same frame/date, exact checkpoint file bytes and cursor-monotonic API rules, respecting provider rate/usage policies.

## Technical execution command order

```bash
# On a trusted clone of this *fresh branch*; no biological outcomes are opened here.
bash scripts/acquisition/bootstrap_fcp_original_roi_v4.sh --code-only

# For actual future measurement, install historical source requirements and obtain
# authenticated original assets; a FULL run checks the original model hashes:
bash scripts/acquisition/bootstrap_fcp_original_roi_v4.sh --full

# Progress source metadata by deterministic shard; repeat until ALL species
# report complete. See full protocol for all partition and completion commands.
python scripts/acquisition/run_fcp_uncapped_shard.py \
  --frame /trusted/42111_outcome_blind_taxa.csv \
  --output-root /trusted/fcp_fresh_metadata \
  --shards 128 --shard 0 --created-before 2026-10-14 \
  --pages-per-species-per-run 10

# Freeze and remeasure ONLY once full metadata census has passed:
python -m fcp_pipeline.fcp_fresh_study freeze \
  --source-species /trusted/42111_outcome_blind_taxa.csv \
  --metadata-root /trusted/fcp_fresh_metadata \
  --shards 128 --worker-partitions 128 \
  --outdir /trusted/fcp_fresh_firewall
```

**Do not run this as a single unthrottled API crawl.** The full photograph count is unknown until the complete fixed-date source census is finished; provider quotas, suitable bulk extraction, compute/storage costs and provenance retention remain operational requirements. Do not silently treat a CI pilot or synthetic regression as the global dataset.

## Real 9-species source progress and switch to bulk

The follow-up metadata-only GitHub Actions [run 38639092021](https://github.com/zuizui0223/fcp/actions/runs/38639092021), artifact `11928883989`, continued the previously completed first taxon and **attempted another eight** under the same outcome-blind frozen date `2026-10-14`. The result is **4 complete taxa, 5 in-progress taxa, 42,102 not started**; **23,023 source photo-manifest rows** (metadata rows, not globally audited distinct flowers); and **zero fresh flower-colour measurements**. These numbers must not be reported as D, W or soil analysis outputs. All remaining partial taxa can resume from their provider ID checkpoints.

**Provider-rate corrective action:** The original iNaturalist observations API is intended for small/medium batches, about **one request per second and ~10,000 requests/day**; bulk access through the API is discouraged ([official recommended practices](https://www.inaturalist.org/pages/api+recommended+practices)). The recurring pull-request API collection job has therefore been **removed from CI**. This prevents repeated PR commits from automatically multiplying public provider requests. Historical API checkpoints stay as provenance/alternative limited route.

### Verified provider inventory and capacity (15 October 2026)

A **real byte-range/HEAD** GitHub Actions preflight ([run 38662030712](https://github.com/zuizui0223/fcp/actions/runs/38662030712), artifact `11939715787`) checked the complete original metadata headers at the provider's mutable root without downloading all rows. The observations (`13,420,521,411` bytes), photos (`20,668,099,997` bytes), and taxa (`39,999,332` bytes) required columns are present. All four data tables including photographers are now required by the new importer. The public bucket lists a **dated September 27, 2026 snapshot archive**, `metadata/inaturalist-open-data-20260927.tar.gz`, of **35,966,549,880 compressed bytes** (~36 GB). Dated snapshots are preferable for fixed-source analysis; mutable top-level CSV aliases are admissible for *schema/identity preflight only*.

**Do not download this giant snapshot onto a standard GitHub Actions runner.** Full production needs sufficiently large persistent disk for the compressed source, extracted metadata, SQL join and many fresh images, plus a separate documented data-ingest budget. Official provider media/API guidance also limits large repeated image downloads; coordinate an approved bulk/access plan before full-scale image acquisition. The current PR only proves format, source identity checks, bootstrap and complete-fail gates; it does not imply a 36-GB corpus or any photo pixels were imported.

### Real 42,111-taxon source identity: verified, not biologically sampled

In [GitHub Actions run 38662636329](https://github.com/zuizui0223/fcp/actions/runs/38662636329), the complete public `taxa.csv.gz` (39,999,332 bytes; SHA256 `a9440e4c1ad1e7b9a8e13d0e8c5b7d35e7a7c5a9f3ec53894f717656edf9348b`) was re-opened **from real provider bytes** against the frozen FCP 42,111-taxon frame. Among **1,668,467** actual taxon rows, **42,111/42,111** FCP taxon IDs were found, with **0 missing IDs, 0 name mismatches, 0 rank mismatches**. An additional 25,394 subordinate taxon nodes were mapped into the target species; the total source taxonomy nodes mapped to these taxa is **67,505**. This validates original taxon identity and ancestry mapping, NOT the availability, colour, genome, local co-occurrence or photo counts of those taxa. The complete 42,111-row source audit is stored in workflow artifact `11939771321`; the key numbers are also frozen in `results/fcp_fresh_source_audit_20261015/source_receipt.json`.

The original dated bundle `metadata/inaturalist-open-data-20260927.tar.gz` is **35,966,549,880 bytes**, confirmed by S3 inventory and bounded gzip/tar header inspection. Its first member is the directory `inaturalist-open-data-20260927`, not a root CSV file. The new tar-extraction helper allows that *verified* dated prefix and supports either plain .csv or pre-gzipped table members; it rejects mutable `latest` aliases, wrong SHA256 values, incomplete tables or reuse of an existing extraction folder. Downloading the full bundle and re-measuring its images have **not** occurred.

### How to freeze the complete source archive, before any new colour results

Run this only on a persistent storage/compute machine with adequate free disk, not the default small GitHub CI runner. The expected archive digest must be recorded *before* source admission; no value is invented here.

```bash
# The source object exists publicly but is a ~36-GB transfer.
aws s3 --no-sign-request --region us-east-1 cp \
  s3://inaturalist-open-data/metadata/inaturalist-open-data-20260927.tar.gz \
  /trusted/snapshots/inaturalist-open-data-20260927.tar.gz

sha256sum /trusted/snapshots/inaturalist-open-data-20260927.tar.gz
# Freeze the actual printed SHA256 in the source study ledger before proceeding.

python scripts/acquisition/pin_inat_open_data_snapshot.py \
  --archive /trusted/snapshots/inaturalist-open-data-20260927.tar.gz \
  --s3-key metadata/inaturalist-open-data-20260927.tar.gz \
  --archive-sha256 VERIFIED_FULL_DATED_ARCHIVE_SHA256 \
  --outdir /trusted/snapshots/20260927_extracted

# These extracted original four metadata tables and their verified SHA256 values
# are the ONLY permitted inputs to the full 42,111-taxon source photo census.
```

### Bulk snapshot importer: all eligible photos without a page cap

`scripts/acquisition/import_inat_open_data_snapshot.py` reads the official iNaturalist Open Dataset monthly metadata tables (`taxa.csv.gz`, `observations.csv.gz`, **`observers.csv.gz`**, `photos.csv.gz`). The fourth observer table is mandatory for each photo's copyright attribution (photographer may differ from observation owner). The [provider's metadata documentation](https://github.com/inaturalist/inaturalist-open-data/blob/main/Metadata/README.md) specifies the gzipped **tab-separated, unquoted** files and dated monthly snapshots.

This is the recommended scalable alternative. A user-selected **immutable dated snapshot**, its original source location, all three SHA256 values and the original 42,111-species frame identity must be recorded **before** import. A mutable `latest` alias is *not* a frozen scientific input. The importer maps included taxon ancestry without colour selection, retains all qualifying research-grade georeferenced licensed photos (including several photos per observation), and emits exactly the same complete per-species source manifests as the original FCP location-blind measurement firewall. Incomplete/zero-photo taxa remain in the census. Old measured flower-colour data are never read or used as a source.

Example using **already downloaded, checksum-frozen** snapshot files (the paths and SHA are placeholders for the actual pinned snapshot; no unverified hash is provided as genuine):

```bash
python scripts/acquisition/import_inat_open_data_snapshot.py \
  --species-frame /trusted/frozen_42111_species.csv \
  --taxa /trusted/snapshot/taxa.csv.gz \
  --observations /trusted/snapshot/observations.csv.gz \
  --observers /trusted/snapshot/observers.csv.gz \
  --photos /trusted/snapshot/photos.csv.gz \
  --snapshot-date YYYY-MM-DD \
  --taxa-sha256 VERIFIED_TAXA_SHA256 \
  --observations-sha256 VERIFIED_OBSERVATIONS_SHA256 \
  --observers-sha256 VERIFIED_OBSERVERS_SHA256 \
  --photos-sha256 VERIFIED_PHOTOS_SHA256 \
  --snapshot-pin /trusted/snapshots/20260927_extracted/fcp_source_archive_pin.json \
  --shards 128 --expected-species 42111 \
  --output-root /trusted/fcp_fresh_bulk_source
```

`fcp_pipeline.fcp_fresh_study.freeze` verifies that all 128 shards have a consistent provider, snapshot date, input SHA256 values and no modified photo ledger; mixing API and bulk-source manifests is prohibited. Its date field means the original **Open Dataset snapshot date**, not an API `created_d2` observation filter. An ingestion job may demand extensive disk and compute resources; the bulk metadata are **not yet downloaded or processed** and the full census remains incomplete.

## Model checksum and offline regression gate — PASS

The standalone fresh-code [GitHub Actions run 38663398717](https://github.com/zuizui0223/fcp/actions/runs/38663398717) passed **42 synthetic regression tests** plus the existing frozen New Phytologist manuscript claim checks. The original ROI-v4 detector and both EfficientSAM encoder/decoder ONNX weights were actually fetched and their SHA256 checks compared to the historical 1,499-species model identities; **all three matched**. This is recorded in artifact `11940305372` and the permanent `results/fcp_fresh_source_audit_20261015/source_receipt.json`. No source image or biological flower-colour value was opened. The full model-binary check is deliberately a one-shot provenance preflight, not a compulsory network download on every subsequent PR commit.

The dated source bundle importer now requires a preverified `fcp_source_archive_pin.json` for the actual 42,111-taxonomy study: (a) S3 key precisely corresponding to its declared snapshot date, (b) full original archive SHA256, (c) **all four** original metadata files' SHA256 values, and (d) the original 42,111-taxon frame's fixed checksum. The `fcp_pipeline.fcp_fresh_study.freeze` gate compares archive pin hashes on all 128 shards and will **not** allow mixing a mutable API sample with the dated bulk sample. Three fake tar fixtures (both raw and gzipped original table members, plus corrupted archive/alias cases) passed, but the 36-GB archive itself has NOT been fully downloaded yet.

## One-shot real ROI runtime and terminal-source validation — PASS

[GitHub Actions 38666471849](https://github.com/zuizui0223/fcp/actions/runs/38666471849), artifact `11941525917`, executed the *actual original frozen ROI-v4 YOLO + EfficientSAM estimator* on a single newly generated synthetic flower-like drawing, with the original measurement source and all 3 frozen model assets. The program produced exactly one source-blind terminal row with `morph=mixed_uncertain`, `measurement_status=not_evaluable_roi_or_flip_gate`; all 9 expected biological palette fraction columns were present and no identity, observer, source URLs, or locations were visible to the evaluator. The test successfully exercised original pixel decode, model load/inference and terminal output, **not evidence of a floral-colour morph**. It used **0 original or live provider photos**. The one-off heavyweight model-inference CI job was removed after the test passed.

**End-to-end source constraints now covered by 43 passing synthetic regression tests:** the official `observers.csv.gz` source table supplies the legal attribution for *photo creator* (not necessarily observation owner), and `attribution` passes through `photo_manifest.ndjson` → frozen `metadata_join_key.csv` → final `newly_remeasured_photos.csv`. The 3-column ROI worker packet remains location-, source- and attribution-blind. Terminal reassembly now verifies each photo belongs to its exact frozen worker partition, batch and within-partition row ordering, using `selected_source_worker_sha256`; a substituted or shifted terminal photo is rejected even when its CSV and self-declared terminal checksum agree. Incomplete batches still fail closed, not converted into a partial biological dataset.

## Real original iNaturalist image execution — 96/96, isolated technical pilot

**Actually executed, 15 October 2026:** [GitHub Actions run 38671341143](https://github.com/zuizui0223/fcp/actions/runs/38671341143) fetched original iNaturalist image bytes for 94 of the 96 taxon-equal, outcome-blind hash-selected photographs from the preexisting source metadata artifact `11928883989`. It ran the **unmodified frozen ROI-v4 + EfficientSAM source** and preserved per-image terminal palette/ROI output in artifact `11944359701`. The two apparent acquisition failures were NOT lost provider images: the custom technical-pilot validator mistakenly rejected existing original URLs ending with uppercase `large.JPG`.

**Corrective execution:** the original 96-row CSV was checked against its exact SHA256 `772dda9c51cf24fa4fdabf6408d16b2ba96fc965b532af82533c150cf0144c09`. A fresh, explicit two-ID repair on [run 38673995478](https://github.com/zuizui0223/fcp/actions/runs/38673995478) downloaded both `840402` and `3201197` with unchanged uppercase source URLs and passed them through the **same original** ROI-v4 model, preserving all other 94 measurement rows cell-for-cell. Independently opening both workflow ZIP artifacts confirmed exactly **94/94 prior rows unchanged**, exactly 96 unique photo IDs, both new original-image SHA256 values present, and the new full terminal CSV SHA256 **`911ec9edd10a7126bee28cfa67738cef0ff4d275b1f0ca1c06e3f26e7e290039`**. The final 96-row dataset is workflow artifact **`11945168406`**. The one-shot 2-photo download job was removed from recurring PR CI after success.

**Actual photo measurement outcomes:** **96** original images downloaded and decoded, **96** processed by frozen original ROI-v4 + SAM, **20** images classifiable by the original 4-state floral colour gate, **67** ROI/flip-gate failures, **9** accepted-ROI but ambiguous-palette failures, **0** image acquisition failures. The 66 original ROI/flip failures included flip colour instability in 62 and flip mask instability in 55 (reasons can co-occur); no quality threshold was relaxed. Strict QC success was **20/96 = 20.8%**: `Ribes californicum` 3/24, `Schoenoplectus californicus` 3/24, `Collinsia concolor` 9/24, `Allium falcifolium` 5/24. The 4-state outputs are *photographic machine classifications*, not experimentally verified flower morphs, and their species' classifiable support remains below the frozen **40-photo minimum**. The original photo metadata also included one *Allium* photo ID attached to three conflicting source observations: all its copies were quarantined by ID before blind selection; it was not counted as three independent plants.

**Quarantine and scope:** This small **24 photos/species** limit applies *only to the technical runtime pilot*, not the separate production protocol (which has **no per-species photo cap**). The pilot was selected by source photo-ID hash, not colour, temperature, spatial location or observer identity; the 4 source species are the 4 that happened to be complete in the earlier limited API metadata batch. It is **not** a representative census of the frozen 42,111 source taxa, not 96 observations of independently proven phenotypes, not the original 149,900-photo dataset, not a fitted W distribution or 999-world null, and not the basis for climate/soil effects. Its images were not persisted, and older image-derived flower-colour labels were **never opened**. End-to-end synthetic safeguards and frozen method verification passed **45 tests** in the repair run.

The current durable audit is `results/fcp_fresh_source_audit_20261015/source_receipt.json`. It distinguishes 96 real images measured in an *isolated method pilot* from **zero full-study images assembled** behind the mandatory all-42,111-species source census completion gate.

## What is still NOT accomplished

1. Source census finished across all 42,111 nominal taxa: **NO**.
2. Every frozen photo downloaded and newly classified with ROI-v4: **NO**.
3. New W cohort-wide 999-world structured-null result: **NO**.
4. New photo-coordinate-matched climate/SoilGrids environment completeness: **NO**.
5. Interspecific (species-equal) ecological inference about within-species W: **NOT YET TESTABLE**.

This project's existing 149,900-image W and prior environmental analyses remain published research-history/method validation, not new fresh-source evidence. No claim about pigment pathway, genotype-defined morphs, natural selection, or adaptation is implied by image-only environmental associations.
