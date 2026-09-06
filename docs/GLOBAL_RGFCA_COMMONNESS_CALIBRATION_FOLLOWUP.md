# RGFCA commonness calibration follow-up

## Current decision

The upstream sampling target is defined, but this is not equivalent to full validation of a test specific to shared geographic boundaries. The 200 balanced outer maps reuse species and photographs; they are not 200 independent biological datasets. All completed observed support decisions remain unchanged.

This follow-up adds an executed synthetic proof-of-concept and a read-only Monte Carlo uncertainty audit. It does not fit new observed flower-colour outcomes.

## 1. Why concentration and commonness must be separated

For centered species fields X_s on a common uniformly weighted grid, with S fixed equally weighted species,

C_total = Var_cell(mean_s X_s)
        = sum_s Var_cell(X_s)/S^2
          + 2 sum_{s<t} Cov_cell(X_s,X_t)/S^2
        = C_self + C_cross.

Individual spatial structure contributes through C_self even without shared phases. Removing C_self is informative but is not sufficient by itself: the comparison distribution must also preserve the individual spatial structure. Cross-covariance is not causal evidence and can also reflect common observation structure.

For cell-dependent opportunity weights, the analogous algebra applies to the actual weighted species contributions on a common evaluation mask, not by substituting unweighted species fields into the formula above. That production transfer has not been executed here.

## 2. Executed structured-null specificity toy

Specification: `supporting/global_rgfca_structured_null_toy_contract_v1.json`, committed as `b86a5967fd53d998e0915070e107e5f870ce0ec9` before the new simulation outputs.

Runner: `scripts/analysis/run_global_rgfca_structured_null_toy.py`.
Tests: `tests/test_global_rgfca_structured_null_toy.py`.
All 64 rate rows, including both statistics and all calibration conditions, are retained in `supporting/global_rgfca_structured_null_toy_summary_v1.json`.

The toy has 40 synthetic species on 64 periodic cells, a fixed width-3 centered unit-SD Gaussian bump, independent Gaussian cell noise, amplitudes 0/0.5/1/2, and shared fractions 0/0.1/0.25/0.5/1. Calibration uses 2,000 replicates per amplitude; independent evaluation uses 1,000 per scenario. Seed: 2026090615. Zero-amplitude duplicates are collapsed.

Two null calibrations are compared: independent species phases retaining spatial shape; and within-species cell permutations preserving marginal values but destroying spatial shape. Both total and cross statistics use their own fixed empirical upper-95% thresholds. These are related methodological diagnostics, not independent biological discoveries.

Selected cross-statistic results (the full family is stored):

| Amplitude | Shared fraction | Shape-preserving calibration | Cell-scramble calibration |
|---|---:|---:|---:|
| 0 | 0 | 7.3% | 6.2% |
| 0.5 | 0 | 4.7% | 7.7% |
| 1 | 0 | 5.2% | 13.5% |
| 2 | 0 | 5.6% | 18.6% |
| 1 | 0.1 | 13.6% | 29.9% |
| 1 | 0.25 | 84.7% | 94.9% |
| 2 | 0.25 | 90.8% | 99.4% |

At amplitude 1 and zero sharing, the shape-preserving rate has a 95% Wilson Monte Carlo interval of 3.99-6.76%; the cell-scramble rate has an interval of 11.52-15.76%. At 25% sharing, shape-preserving power is 84.7% (82.34-86.80%). These intervals are conditional on the empirical thresholds and omit calibration-threshold uncertainty.

Important exception: at zero amplitude the shape-preserving cross rate is 7.3% (5.85-9.08%), not 5%. Therefore this pilot is not a certification of 5% size control across all conditions. No threshold, width, or seed was changed after seeing this result.

The exact toy expectation is E[C_cross] = q(q-1) A^2/S^2. Nine local unit tests passed. The maximum covariance-decomposition identity discrepancy in the executed calibration was 6.94e-17. The periodic ring is not the Earth, the observed clustered photo geometry, the tied categorical palette/JSD, or the balanced outer schedule. These rates do not estimate the actual fraction of species sharing a flower-colour boundary.

## 3. Read-only audit of the completed S1-S3 artifact

Runner: `scripts/analysis/audit_global_rgfca_design_simulation_mc.py`.
Source: design-recovery Actions run 34038925164. The downloaded result JSON was matched byte-for-byte by Git blob `87962d88f66f86bae534604d15e237a91b98cf6e` to `supporting/global_rgfca_design_power_identifiability_simulation_result_v1.json`. The 38 S1, 5 S2 and 105 S3 rows were checked, and S2 JSON/CSV values reconciled.

S1 independent null evaluation was 29/500 = 5.8% (Wilson 4.07-8.21%) for the classifiable frame and 18/500 = 3.6% (2.29-5.62%) for the complete synthetic frame. These are calibration diagnostics, not biological intervals. Frame contrasts also change species, locations and schedules; they do not isolate a causal missingness effect. The scalar absolute-difference proxy does not validate tied four-group JSD.

S2 at 200 draws: independently bootstrapped consensus pairs have median r=0.88360, with 632/2,000=31.6% reaching r>=0.9 (MC interval 29.60-33.67%). Correlation to the same frozen 200-map reference has median r=0.94505 and 92.75% reaching 0.9. The latter is not independent replication. Both bootstrap arms remain conditional on the same original maps.

S3 is the declared Gaussian shifted-null approximation. Across correlation scenarios 0/0.25/0.5, thermal target Holm power is 38.69-38.80% at an assumed effect of 0.0075, 65.04-65.37% at 0.01, and 95.85-96.00% at 0.015. Atmospheric energy/dryness is approximately 38%, 64-65%, and 95.6-95.7%, respectively. The first tested grid effect attaining 80% is 0.015 for those blocks; this is not an exact continuous MDE. Edaphic regime first reaches that level at 0.02. This does not establish that the observed effects are true or reopen their failed parent gates.

## 4. S4 interpretation retained

Sources: `supporting/global_rgfca_s4_low_information_weighting_execution_v1.json`, `supporting/global_rgfca_s4_low_information_weighting_result_v1.json`, and its runner. S4 independently samples each species' zero-centered empirical H6b permutation distribution; it is not a new H6 test.

Under this fixed centered-null simulation, equal-weight 188-species SD is 0.02739, pair-weighted SD 0.01442, inverse-null-variance SD 0.01336, and reliable-74 equal-weight SD 0.02294. Equal-188 variance is 4.2036 times inverse-variance-188 variance; means remain near zero. Heterogeneous information increases spread, not necessarily systematic mean bias. Weighted effective species counts are 40.12 and 29.54, respectively; under heterogeneous true effects those weights can target a different species average. These results do not prove that H6's observed sparse signal is entirely noise or justify changing the frozen primary weights.

## Reproduction and remaining task

Run the focused tests with `python -m unittest discover -s tests -p test_global_rgfca_structured_null_toy.py -v` and the toy with `OPENBLAS_NUM_THREADS=1 python scripts/analysis/run_global_rgfca_structured_null_toy.py --output-dir NEW_OUTPUT`. Outputs refuse overwrite. The executed NumPy version was 2.3.5.

The next methodological task is transfer to actual species contribution fields with fixed geographic opportunity, ranges, missingness and palette measurement. A ring shift must not be applied directly to the Earth or move species outside their available range. The already-authorized environmental heterogeneity run 34039210812 remained in progress at the last check; no outcome from it is asserted here. No additional observed-data commonness test or gate change is authorized by this toy.
