# Fixed nonlinear environmental transfer test

## Question and provenance

Does a shared nonlinear response to all 36 measured environmental predictors improve photographic flower-colour prediction for entirely unseen genera in a different longitude sector, with a 500 km original-photo buffer? This tests a named unresolved model-form alternative in PR #169, not a causal mechanism or universal law.

Parent: `ae1136e7d43bc6c5711778a5a7c727b0b7ed8d7e` (PR #169). Newer PR #168 at `30a855c28ac060bed084337a85cd09ce1092ee56` already contains local training-count controls and analogue-support checks on a different photo-cell estimand. They are context and are not numerically pooled here. This design is retrospective and source-exposed; settings below are fixed before this new model's original-source outcomes are viewed. It is not preregistered or an independent biological replication.

## Frozen design

1. Reuse the exact original two ZIP SHA256 locks, 100,543 unique photo opportunities / 42,111 nominal taxa. Colour-unclassified photos stay unclassified and contribute to geographic exclusions. Main cohort: all36-complete classified 31,925 photos / 20,786 species. Existing strict-photo sensitivity: 6,457 species, using unchanged pixel/nuisance thresholds. No new taxa, raw photo relabelling, variable selection, or source-species resampling.
2. Reuse PR #169's `aggregate_species`, four 90-degree longitude sectors, deterministic five genus folds, and whole-species original-photo footprint exclusions. Evaluate the same target species-sector rows under genus-only, sector-purged, and sector-plus-500km training. Every species receives total evaluation weight one across its sectors, then genus-equal and species-equal summaries are both shown.
3. Fixed nonlinear learner: `ExtraTreesRegressor(n_estimators=128, max_depth=16, min_samples_leaf=20, max_features=1.0, bootstrap=False, criterion='squared_error', random_state=2026101161, n_jobs=2)`, with four continuous colour proportions jointly as outputs and equal training-species weights. All outputs are averages of observed training compositions and must stay on the simplex. No tuning, feature screening, target-label-dependent early stopping, or favourable model/fold selection. This permits nonlinearities and interactions without asserting that every such response is detectable.
4. Fit four predictors to each identical train/test split: original fixed-alpha40 ridge GEO; original ridge GEO+all36; nonlinear GEO; nonlinear GEO+all36. GEO is exactly `abs_lat, lon_sin, lon_cos`, inherited for direct comparison; this is a specified coordinate baseline, not comprehensive spatial adjustment. Tree training and preprocessing use training data only. Save all four-colour predictions as well as loss components and fold audits.

## Endpoints and decisions

Primary endpoint: genus-equal reduction in squared four-colour composition error from nonlinear GEO to nonlinear GEO+36 under the 500km regime, on the all36 main cohort. Positive values indicate predictive improvement. Use the inherited 1,999 paired genus bootstrap with fixed fits, and report 95% and 97.5% percentile intervals. No significance gate for source execution.

Necessary paired comparisons on the identical target rows: (a) nonlinear environmental gain minus ridge environmental gain; (b) ridge full-model loss minus nonlinear full-model loss; (c) ridge GEO loss minus nonlinear GEO loss. Report absolute baseline/full losses so a larger environmental increment cannot be mistaken for better prediction when its baseline worsens. A common nonlinear response is supported only at the predictive level represented by the primary contrast; model comparison and strict-quality sensitivity qualify robustness.

Descriptive checks: all three training regimes, strict-photo sensitivity, all four geographic sectors, white-axis and remaining-hue error components (total = white-axis + remaining-hue). Preserve negative and inconclusive cases. No new claims about pigment investment, selection, pollination, herbivory, fitness, individual genotypes, or universally absent/present environmental responses. Bootstrap conditions on existing fits and folds; four sectors are artificial domains, not four independent evolutionary replicates.

## Implementation and verification contract

Create `scripts/analysis/fcp_nonlinear_environment_transport_20261011.py`, focused tests `tests/test_nonlinear_environment_transport.py`, and a source-replay workflow. Keep parent source files unchanged. Reuse original helper functions instead of copying them. Tests must cover heldout-outcome invisibility, genus/species isolation, finite/simplex outputs, deterministic fits, paired-gain algebra, equal-species sector weighting, and rejection of duplicated or mismatched target rows. Original stable-fold and footprint tests also run. Save a machine-readable result, species-sector predictions/losses, all fold audits, fold-map hashes, source/script hashes, runtime versions, and exact model settings. Compare freshly recomputed ridge losses with the previously archived PR #169 outputs.

Execution sequence: freeze this protocol; tests fail for missing behaviour; implement; targeted tests pass; recover and hash original archives; execute the two fixed full36 cohorts; independently recalculate all point summaries and key paired intervals from saved predictions; review interpretation; save results and a draft PR. Final result receipt must state actual local/remote execution status. Main/H1/H2 and unused cohorts stay as they were.

## Pre-execution scientific review

The independent protocol review found no execution blocker. A positive nonlinear environmental increment alone must not be called proof of nonlinearity or restoration of the failed ridge transfer. That stronger predictive comparison additionally needs positive paired full-model improvement and nonlinear-minus-ridge environmental gain, with actual coordinate-baseline losses visible. The fixed trees average training compositions and do not generally extrapolate a smooth environmental trend; failure therefore cannot exclude a transferable nonlinear relation outside their learned support. These qualifications change no learner setting, source, fold, endpoint, or analysis scope.
