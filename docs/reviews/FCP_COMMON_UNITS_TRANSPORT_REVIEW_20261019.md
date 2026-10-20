# Common physical units and transport: independent review

Review date: 2026-10-19. Scope: protocol sections A/B/C, native-unit implementation/tests and executed results, independent source and executed-output review of B/C, and the A/B/C numerical claims in `docs/FCP_COMMON_UNITS_TRANSPORT_MOLECULAR_20261019.md`. The reviewer implemented A and independently recomputed its key outputs; the B/C implementation is owned by a separate agent. This is an internal review, not an external validation. Molecular claims are outside this review's scope.

## A. Completed and verified

The common-unit analysis preserves all 1,574 species and 86,114 photos. It exactly transforms each joint slope and the full 3-by-3 covariance using the frozen within-species native SDs. No original result, trait join, admitted category, or photo fit was changed. Native results are in `results/fcp_environment_native_units_20261019/`.

| Increment | Equal-species mean | Genus-bootstrap 95% interval | REML mean | Modified-HKSJ 95% interval |
|---|---:|---:|---:|---:|
| BIO5, 1 degree C | 0.032329 | 0.017856 to 0.047277 | 0.017076 | 0.008974 to 0.025179 |
| BIO14, 10 mm | -0.077612 | -0.189997 to 0.031593 | 0.004859 | -0.005692 to 0.015411 |
| SRADmean, 1,000 kJ m^-2 day^-1 | -0.041635 | -0.115058 to 0.028889 | -0.003267 | -0.021067 to 0.014533 |

Source: `primary_summary.csv` in the native results directory. All nine GIFT nominal and approximate family9 intervals include zero, with exactly the prior coverage and categories (`slope_modifier_contrasts.csv`, `gift_scale_comparison.csv`). The positive average BIO5 association therefore persists under common physical increments; a GIFT modifier is still unresolved. These observations do not identify a causal mechanism or establish trait-group equivalence.

The most consequential scale diagnostic is BIO14. The 79 species in the lowest 5% of sampled BIO14 SDs contribute -0.086342 to the full equal-species mean of -0.077612 and 79.31% of the approximate sampling-variance sum. The 115 species whose sampled precipitation span is less than 10 mm remain included. These additive diagnostics explain why the physical-unit equal mean is imprecise; they are not response-driven trimming. Sources: `environment_range_diagnostics.csv`, `range_tail_contributions.csv`, `range_tail_species.csv`.

The means on the two scales average different quantities, and REML weights also change under species-specific unit conversion. Different numerical means, including different signs of imprecise rainfall point estimates, are not mutual refutations. Boundary tau estimates are not evidence that species are biologically identical.

### Verification evidence

- Six native-unit unittests pass, covering covariance conversion, invalid inputs, analytic REML/scale equivariance, zero-tau handling, unequal-genus-size arithmetic, and hash mismatch rejection.
- All 114 frozen input/code/cache hash references (112 distinct files) pass, including 72 GIFT source pages, the authentic source photos, both WorldClim archives, all 18 panel rasters, and the recorded source scripts/outputs.
- All 86,114 observed fitted probabilities are preserved; maximum absolute difference is 6.66e-16. The 1,574 joint covariance roundtrip, positive-definiteness, correlation, Wald, and observed slope-score variance checks pass. Covariance scope is slope-only because intercept-slope covariances were not saved by the original fit.
- All 270 real-data REML numerical-scale checks pass. The unchanged source REML function uses a common numerical normalization that is reversed in its dimensional outputs.
- Every original SD-scale GIFT result field is reproduced before changing units; all coverage, categories, original seeds and 5,000 shared draws per trait are preserved.
- After execution, all 17 receipt-listed output hashes and the script/test hashes were rechecked. Independent direct arithmetic reproduces the written native equal means and REML weighted means.

Detailed evidence is saved in `input_verification.csv`, `prediction_covariance_audit.csv`, `reml_numerical_scale_checks.csv`, and `execution_receipt.json`. The exact protocol snapshot is `protocol_at_execution.md`.

## B/C. Completed output review passed

The source-level review and all 16 final transport unittest checks found no consequential leakage or paired-scoring defect. B excludes the complete held-out genus before coefficient aggregation and evaluates all white/chromatic pairs. C uses train-only feature bases and longitude origin, great-circle buffers across the antimeridian, identical retained training rows for paired models, constant-response eligibility, and exact common successful held-out photo predictions. Its common smoothed intercept recalibration uses only training data and preserves slopes. Genus bootstrap arithmetic retains all sampled-genus member species and equal species weight.

All 16 receipt-listed transport artifacts passed hash checks, and every gzip/NPZ archive passed a complete integrity check. The original executed script snapshot matches the scientific execution hash; the current script and final test file match the separately recorded correction hashes. The source slopes and imported model/sampling code hashes match the frozen prior execution receipt, independently authenticated in A.

The saved map contains exactly 86,114 unique actual source photo IDs from 1,574 species. An independent join to the SHA-authenticated 447,500-photo source confirmed every species/taxon ID, observation ID, source cohort, latitude, longitude, classifiability and white/chromatic label. Independently recomputed SHA256 genus and cell assignments agree with every saved fold. The 764,081 within-species white/chromatic pairs are counted once per scoring model, not independent replicates.

For B, all 30 training-fold REML means agree with an independently normalized calculation to at most 4.32e-10. Every saved photo score and all 1,574 species' pair losses/AUCs for each of the two models reproduce exactly with round-trip CSV parsing. The small positive family2 loss gains are therefore supported by the saved predictions: BIO5 component +0.000699, approximately [0.000092, 0.001321], and joint three-component score +0.000634, approximately [0.000069, 0.001183]. Their AUCs are 0.516061 and 0.514402 (`genus_summary.csv`).

For C, the independent replay reconstructed all 15,740 species/buffer/fold splits and all 51,650 successful fits, including 42 successful constant-response fits. It checked retained training identities and counts, cell exclusion, great-circle buffers, train-only means/SDs/longitude origins, raw Firth coefficient score equations, the common training calibration target, raw and recalibrated test predictions, and corrected range flags. All 561,514 saved successful spatial prediction values reproduce exactly. Maximum achieved-training-mean error is 6.18e-14; maximum training Firth score is 4.40e-6; maximum independently calculated nearest-distance discrepancy is 3.64e-12 km. Every species loss/Brier gain, fold gain sum, reported coverage count and archived nominal/family2/family4 bootstrap endpoint also agrees.

| Buffer | Comparison | Species / photos | Mean loss gain | Approximate family4 interval | Corrected range-outside photo fraction |
|---|---|---:|---:|---:|---:|
| 0 km | Climate vs intercept | 1,570 / 83,890 | -0.059924 | -0.072569 to -0.049379 | 14.6549% |
| 0 km | Geography + climate vs geography | 1,237 / 61,394 | -0.149777 | -0.233317 to -0.081168 | 22.6358% |
| 50 km | Climate vs intercept | 1,566 / 78,866 | -0.073453 | -0.086930 to -0.060981 | 16.4799% |
| 50 km | Geography + climate vs geography | 1,234 / 56,519 | -0.237316 | -0.338889 to -0.146201 | 25.7931% |

Source: `results/fcp_environment_transport_20261019/spatial_summary.csv`, checked against `spatial_oof_predictions.npz`, `spatial_species_metrics.csv`, fold records and `spatial_bootstrap_draws.npz`. All four Brier gains are also negative with negative family4 endpoints.

### Output repairs and their scope

The original gzip photo map was truncated. Its owner reconstructed it from the authenticated source and frozen rasters, preserved the decoded prefix, and documented the original/new hashes and reconstruction in `execution_receipt.json` under `output_repairs`. The isolated truncation's root cause is not established. The complete replacement passed the independent source-photo join and every coefficient/prediction replay described above.

A separate regression test demonstrated that inverse standardization could mark an exact training-boundary value as outside the training range. The corrected code compares physical features directly with a documented 64-epsilon boundary margin. The independent replay exactly reproduces every corrected flag. Every nonrange spatial NPZ component has the same byte hash before and after the correction, including all raw/recalibrated scores, probabilities, train counts, distances, and standardized-distance diagnostics. The archived bootstrap arrays are unchanged. The receipt preserves the scientific execution snapshot, correction code/hash, before/after range counts and affected artifact hashes. Neither repair changes the scientific loss comparisons.

The previously raised covariance-interpretation audit item is resolved in the current profile and integrated report: `max_test_linear_predictor_variance` uses raw Firth covariance and must not be read as a complete uncertainty calculation for recalibrated predictions.

## Interpretation and manuscript-facing claims

The A/B/C numeric claims inspected in the integrated report agree with the corrected saved outputs. Section 3.3 correctly distinguishes B's common coefficients and within-species conditional pair loss from C's separately trained species models and binary photo loss. Their gains cannot be subtracted to quantify a single loss of generalization. C did not test transfer of B's common pooled coefficient into held-out locations; it tested the specified species-specific fitting and prediction procedure.

The reviewer recommended describing the 0 km comparisons as held-out cells rather than implying a guaranteed distance, and describing common units as aligning the physical increment rather than removing all effects of narrow sampled climate ranges. These are wording safeguards, not numerical defects.

The protocol correctly identifies the smoothed recalibration as a distinct fixed prediction procedure, not standard FLIC or guaranteed held-out calibration. The previously selected panel, residual cross-fold dependence, and bootstrap conditioning on fitted out-of-fold models remain inference limits. The four negative C results support lack of improvement by the tested procedure on its matched held-out predictions; they do not show absence of climate physiology, invalidate the positive pooled BIO5 association, or establish why a different model might succeed. No unresolved consequential implementation or output-integrity defect remains in the reviewed A/B/C analyses.

## Release packaging verification

The GitHub connector's 16 MiB request limit required recompressing `spatial_oof_predictions.npz` from ZIP DEFLATE to ZIP LZMA: 13,775,669 bytes became 10,388,989 bytes, with a 13,851,988-character base64 payload. All 10 ordered `.npy` members remain byte-identical; `np.load(allow_pickle=False)` confirms identical array dtypes, shapes, bytes and values. All 15 other previously recorded outputs, including every metric and bootstrap artifact, remain byte-identical. The generator, plotting script, test and protocol hashes are unchanged, and prior repair/correction history is preserved. The original archive is retained outside the repository under `fcp_inputs/transport_packaging_20261019/`; the receipt's `packaging_changes` records original/repacked hashes, member hashes and verification. `repack_spatial_predictions_20261019.py` reproduces the tested container exactly. A normal scientific rerun still produces DEFLATE and can be followed by this explicit packaging helper. This container-only change leaves every reviewed scientific result unchanged.
