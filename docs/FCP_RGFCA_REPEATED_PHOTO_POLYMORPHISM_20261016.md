# FCP × RGFCA: repeated random-photo robustness, not maximum photo volume (2026-10-16)

## Decision: change the primary sampling objective

**Research question**: Among plant species that have been selected **without looking at their flower colours**, which species repeatedly show at least two independently observable floral-display colour morphs when photographs are redrawn many times?

This adapts the **method** of the frozen Repeated Global Flower-Colour Atlas (RGFCA). The original global shared-hotspot conclusion is **not** carried forward: G1 concentration had `p=0.070` and the species-disjoint global commonness test had `p=0.856`. Those are different geographic questions, not evidence for or against the present photo-based polymorphism detection.

All stages use the dated **2026-09-27** provider source. Frozen count receipt: `ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9`, 42,111 taxon IDs and **93,209,351** qualifying source photo IDs. All those counts are PHOTO METADATA, not valid flower ROIs and not species known to be polymorphic.

## A. What has been executed and checked

The independently audited original counts have **24,740 eligible taxa with ≥128 source photographs**, split into exact opportunity tiers:

| Prior original source photo opportunity | Available taxa | Pre-colour sample quota |
|---|---:|---:|
| 128–255 | 4,990 | 334 |
| 256–1,023 | 9,187 | 333 |
| ≥1,024 | 10,563 | 333 |
| **Total** | **24,740** | **1,000** |

A SHA256-salted, taxon-name/ID-only source enrolment has been performed. Within each quota tier, priority comes from source identity and a fixed seed — **not** biological morphology, flower brightness, published polymorphism classifications, climate, soil or observer-reviewed colour.

Then RGFCA-balanced repeated species inclusion was computed on the **genuine 42,111-row source count CSV**:

- 200 RGFCA realizations, each selecting 250 distinct species;
- 4 realizations per cycle cover all 1,000 species exactly once;
- 50 randomized cycles → **every enrolled species appears exactly 50 times**;
- **50,000 total species × realization assignments**;
- convergence/coverage at R = 25, 50, 100, 150, 200 has exactly the original species-count balance property (at R=200: min=max=50);
- actual original provider source photo IDs across the 1,000 selected taxa: **3,430,502**, but these are only source opportunities;
- planned original ROI-v4/EfficientSAM unique image requests are **at most 1,000 × 128 = 128,000** for the first study pool, and can be fewer when observers, valid flowers or files are scarce;
- 200 × 250 × 120 = **6,000,000 repeated photo DRAW assignments** (if all selected pools support 80 each). These repeated assignments sample the **same fixed original image pool**: they do **not** require 4 million downloads and they do **not** create 4 million independent flowering plants.

Actual real-data CI: [RGFCA source-only execution #38828417100](https://github.com/zuizui0223/fcp/actions/runs/38828417100), **success**, artifact `12015683657`: `fcp-real-rgfca-photo-1000species-200x250-balanced-no-pixels`. Deterministic:
- 1,000-taxon frame SHA: `88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c`;
- exact 200×250 realization schedule SHA: `1f0d52f00d14d503c146c355934993eb7d4836fe4f5e523491ea18f08b49bb30`.

**These complete artifacts are selection schedules, not colour measurements.** No `P_poly`, H2, W, ecological or genetic effect has been computed from them.

## B. Photo-first outcome measurement, then repeated variation detection

Use the existing source-only [photo-first protocol](FCP_PHOTO_DEFINED_POLYMORPHISM_PROTOCOL_20261016.md), with the new optional `--rgfca-enrolment` and `--rgfca-receipt` to ensure the ORIGINAL SQLite index pulls only the selected **1,000** source taxa. That path still authenticates the full immutable source and checks each enrolled taxon's complete original metadata count, but it does not scan or select un-enrolled species' original image pixels. The full 42,111-frame count ledger remains available for source-opportunity/generalizability diagnostics.

Each selected taxon gets up to 64 original photos from each of two observer-disjoint folds, all selected **before opening pixels**; different localities/observers are preferentially represented under the fixed hashed plan. Only a real original ROI-v4 model + EfficientSAM evaluator may score those exact photo IDs. Failed or ambiguous images are retained as unclassifiable terminals.

The repeated robustness engine `fcp_pipeline/rgfca_photo_repeatability.py` then uses the *same* real sealed original photo results in each RGFCA inclusion to draw **60 photos per observer fold (120 per species/realization)**. Within each repeated draw:

- at least 40 photos must actually contain a classifiable original flower ROI;
- at least 20 of those must come from each observer-disjoint fold;
- >=2 coarse states must occur and the minority has >=10% of classifiable source photos;
- both morphs must have sufficient independent observations/observers, and the same pair must be replicated in each half (the existing original photo screen criteria);
- source ROI failures and colour-uncertain results are *not replaced with prettier or easier photographs*.

The species-level descriptive robustness estimand is

`P_poly(i) = n_positive_evaluable_realizations_i / n_evaluable_realizations_i`

with an eligibility ledger alongside it. A primary **robust photo-phenotype** flag requires at least **30 evaluable species realizations and P_poly >= 0.80**. Sensitivity thresholds 0.60 and 0.90 and photo-draw depths 40, 80 and 128 are separately recorded alongside the calibrated 120-primary draw to avoid a single arbitrary cutoff deciding conclusions. Non-detection is never monomorphism: `unevaluable`, `unstable`, and `single_photo_state_detected` remain separate. A 120-photo draw sampled from up to 128 frozen photos has substantial source-ID overlap across realizations and is repeated Monte Carlo *reweighting* of an existing photograph pool, not a biological replicate or an independent p-value.

For a true minority morph rate of 10% with independent perfectly scored photos, the probability of seeing at least two minority observations in **each of two 60-photo folds** is about 0.973. With failed original ROIs, dependent individuals or unequal observer contributions, power falls. A strict P_poly threshold must therefore be paired with sensitivity and missingness information.


### Explicit pre-image technical-power amendment (16 October 2026)

A separate **96-real-original-photo** SHA-frozen technical pilot (4 other taxa, not members of the RGFCA 1,000-species cohort) classified **39/96** original ROI-v4 images as flowers (40.625%). This is not an estimate of 1,000-species classifiability. If each photo's classification were independent with the same hypothetical probability, an **80-photo** draw would achieve 40 classifiable images only **5.64%** of the time, versus **95.85%** for a **120-photo** draw. The same 2×observer-fold minimum of 20 usable photos also matters. To avoid spending most of 200 realizations in the `unevaluable` state for a predictable sampling-depth reason, BEFORE opening any ROI photo of the newly sampled 1,000 taxa the primary was amended to **120 photos, 60 per observer fold**, with 40/80/128 sensitivities. This DOES NOT change the 1,000 selected species, 200×250 balanced inclusion list, or original 93,209,351 source opportunity denominator.

Using a deliberately simplistic shared-probability Beta(1,1)→Beta(40,58) update from that separate pilot yields a 95% interval of **31.3%–50.7%** for the pilot classifiability parameter; Beta-binomial predicted probability of reaching 40 usable photos is **12.45%** for an 80-photo draw and **88.25%** for 120. The two-fold 20-usable threshold pushes the true joint evaluability below these values. These are resource-planning sensitivity values only, *not* a validated general population photo-quality model or empirical polymorphism frequency.

Recalibrated genuine source-only sampling job: [#38832284491](https://github.com/zuizui0223/fcp/actions/runs/38832284491), **success**, artifact `12018230611`. Original species cohort and 200×250 species inclusion SHA digests are unchanged; the pre-image receipt now specifies 120/60 and comparison draw depths 40/80/128.

A separate actual-pixel calibration job [#38832829706](https://github.com/zuizui0223/fcp/actions/runs/38832829706) was launched for **six** new taxa drawn colour-blind from the frozen 1,000 species (two from each opportunity-depth tier). That job uses **live iNaturalist API observations created through 2026-09-27**, capped at two metadata pages per taxon; its photo IDs cannot be declared exact members of the *dated September 27 bulk archive*. Its source-only and ROI output (if successfully sealed) will be quarantined as technical QC. It does not produce scientific `P_poly`, even if real image masks are measured.

## C. Executable order once the true dated photo-ID stage is available

```bash
# Get the corrected pre-image 120-photo RGFCA source schedule from Actions #38832284491.
# The original 2026-09-27 all-photo metadata database must already be staged
# on persistent storage. The existing source-only 42,111 count CSV and
# provenance receipts must match their independently frozen SHA values.
python scripts/acquisition/plan_fcp_photo_first_screen.py \
  --staged-source /persistent/fcp/source_20260927/complete_42111_metadata_source \
  --original-frame /persistent/fcp/inputs/42111_species.csv \
  --source-counts /persistent/fcp/references/all_42111_species_photo_opportunity_counts.csv \
  --rgfca-enrolment /persistent/fcp/rgfca/rgfca_species_1000_original_source_only.csv \
  --rgfca-receipt /persistent/fcp/rgfca/rgfca_sample_receipt.json \
  --photos-per-fold 64 --max-per-observer 8 \
  --outdir /persistent/fcp/photo_first_rgfca_1000

# Freeze the resulting exact photo IDs into blinded original ROI-v4 workers,
# run each selected batch with the original source/model hashes and the
# authorized aggregate provider rate, then reassemble all selected ID
# results: see docs/FCP_PHOTO_DEFINED_POLYMORPHISM_PROTOCOL_20261016.md.

python -m fcp_pipeline.rgfca_photo_repeatability \
  --cohort /persistent/fcp/rgfca/rgfca_species_1000_original_source_only.csv \
  --balanced-schedule /persistent/fcp/rgfca/rgfca_200x250_balanced_species_schedule.csv \
  --schedule-receipt /persistent/fcp/rgfca/rgfca_sample_receipt.json \
  --new-original-photo-measurements /persistent/fcp/roi/newly_remeasured_photo_screen.csv \
  --measured-photo-screen-receipt /persistent/fcp/roi/sealed_photo_screen_terminal.json \
  --outdir /persistent/fcp/rgfca_photo_polymorphism_repeatability
```

**Operational dependency:** The 36GB dated provider archive and permanent original 93M photo-ID SQLite have not been retained by the source-count-only GitHub-hosted runner. The original source count was verified and the RGFCA 1,000-species selection was executed **without image pixels**. No all-source or 128,000-photo ROI measurement has been reported complete. A persistent source archive/SQLite snapshot plus appropriate provider-scale photo access remains the blocker for true biological colour results.

## D. Interpretation boundaries

The correct manuscript endpoint is *replicability of photographic within-species floral colour states under randomized, equal-species and independent-observer source selection*, followed by within-species geography/W/ecological covariance among confirmed photo-phenotype species. Selection based on available public photographs makes inference conditional on the sampled opportunity frame. There is **no** unbiased estimate of global polymorphism prevalence from the 1,000 species or a proof of hereditary polymorphism from pictures. The old RGFCA G1 global-hotspot null remains null: successful P_poly robustness does not turn it into a positive common-world-map result.
