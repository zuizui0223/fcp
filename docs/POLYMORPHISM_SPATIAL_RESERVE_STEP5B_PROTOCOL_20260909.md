# Polymorphism spatial reserve Step 5b — frozen replication protocol

Date: 2026-09-09 JST

## Question

Does the Step-5 discovery result — stronger within-species geographic flower-colour structure among species with a non-trivial second morph — recur in the independent reserve photo tranche under the same fixed morph-frequency thresholds?

This is a replication of the **subset definition and aggregation rule**, not a new search for a threshold. The reserve whole-frame results are already known, but reserve polymorphism labels have not been joined to per-species spatial outcomes for this analysis.

Before this protocol was frozen, one reserve shard was opened only to verify artifact structure and array dimensions; its first displayed species rows were visible. No reserve polymorphism membership, subset aggregate, subset p-value, or threshold comparison was computed. The thresholds and decision rules below are inherited unchanged from discovery Step 5.

## Frozen upstream reserve evidence

Reuse the exact per-species outputs and 999 null arrays from GitHub Actions run `34178957447` (`rgfca-reserve-species-replication-v1`). Each of 363 eligible species has four fixed metrics in this order:

1. `primary`;
2. `observer_pair_exclusion`;
3. `calendar_quarter_stratification`;
4. `matched_background_differential`.

No reserve geographic distance, colour distance, observer filter, calendar filter, background contrast, seed, or permutation is regenerated.

The runner must first reconstruct the already known whole-reserve primary result (approximately mean rho `0.0254826`, upper-tail `p=0.001`) and verify all 363 species × 4 metrics × 999 null values are present before accepting subset results.

## Frozen reserve polymorphism frame

Use `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`.

For each species:

- retain rows with `global_classifiable == true` and `morph` in `white`, `yellow_orange`, `red_pink`, `blue_purple`;
- require at least 40 classifiable rows;
- calculate `D = 1 - sum_k p_k^2`;
- calculate `second_fraction` as the second-most-common morph frequency.

The frame must contain exactly 363 species and reproduce the already audited reserve fingerprint (D max approximately `0.680272`, second morph >=10% approximately `41.60%`, >=20% approximately `24.79%`).

## Subsets

Primary: `second_fraction >= 0.10`.

Sensitivity: `second_fraction >= 0.20`.

Complement for the dilution contrast: `second_fraction < 0.10`.

No alternate threshold is admissible.

## Primary reserve replication

For each metric, average the species-level observed rho across the primary subset. For each of the 999 frozen permutation indices, average the corresponding null rho over exactly the same species.

One-sided p-value:

`p = (1 + number(null_mean >= observed_mean)) / 1000`.

The primary reserve replication is the `primary` metric.

## Nuisance robustness

The observer-pair exclusion and calendar-quarter stratification metrics are prospectively retained robustness checks.

A strong reserve replication label requires:

- primary observed mean rho > 0 and p < 0.05;
- observer-pair exclusion observed mean rho > 0 and p < 0.05;
- calendar-quarter stratification observed mean rho > 0 and p < 0.05.

The matched-background differential is reported as a secondary flower-specific diagnostic. It is not required for Claim 3 because the earlier whole-reserve flower-specific conjunction did not pass. If it is positive with p < 0.05 in the fixed polymorphic subset, it may be described as new exploratory evidence that the subset signal exceeds matched background; otherwise the claim remains photo-derived flower-colour structure rather than a flower-specific biological spatial effect.

## Replication of dilution contrast

For the primary metric only, calculate the same >=10% subset minus <10% complement mean-rho contrast as discovery Step 5, and the corresponding contrast at each frozen permutation index.

A reserve dilution replication requires observed delta > 0 and one-sided p_delta < 0.05.

## Continuous diagnostic

For the primary metric only, correlate reserve D with the species-level observed primary rho across all 363 species. Compare with the 999 correlations obtained using the same fixed D and each permutation's species-level primary rho.

This threshold-free D gradient is secondary.

## Frozen interpretation

- `RESERVE_POLYMORPHIC_SUBSET_REPLICATED_ROBUSTLY` if primary, observer-pair, and quarter tests all pass in the >=10% subset.
- Add `WITH_DILUTION_REPLICATION` only if the reserve primary subset-minus-complement contrast also passes.
- Add `WITH_MATCHED_BACKGROUND_SUPPORT` only if the fixed matched-background differential is positive with p < 0.05.
- The >=20% sensitivity is always reported and is not a rescue threshold.

Failure of reserve Step 5b does not erase discovery Step 5; it limits Claim 3 to exploratory discovery evidence.

## Claim ceiling

Even a full replication supports only a repeated photo-derived within-species geographic colour association in species meeting the frozen morph-frequency threshold. It does not establish a common boundary, common climate threshold, adaptation, causal selection, or population-genetic differentiation.
