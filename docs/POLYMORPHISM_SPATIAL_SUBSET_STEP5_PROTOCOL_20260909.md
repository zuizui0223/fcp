# Polymorphism spatial subset Step 5 — frozen protocol

Date: 2026-09-09 JST

## Question

Do species that meet a prespecified discrete flower-colour polymorphism threshold retain the previously observed within-species geographic distance–colour dissimilarity association when analysed on their own?

This is the Step-5 test after the FCP mainline pivot to species-level polymorphism. The whole-frame 369-species spatial omnibus result is already known (`observed mean rho = 0.02702130374565584`, `p = 0.001`), so this subset analysis is explicitly **post-outcome exploratory**. It is not a new confirmatory test.

## Frozen upstream spatial evidence

Reuse, without recomputation or tuning, the exact per-species observed and 999 within-species colour-label permutation statistics from GitHub Actions run `34088925008`, artifact family `global-rgfca-within-species-omnibus-shard-*`.

The upstream statistic is the species-equal mean Spearman association between within-species geographic pair distance and flower-colour Jensen–Shannon dissimilarity. The upstream run contains one observed statistic (`permutation_index = -1`) and the same 999 permutation indices (`0..998`) for each of 369 species.

No geographic distances, colour distances, seeds, species-level rho estimates, or permutation labels may be regenerated for Step 5.

## Frozen polymorphism definition

Reconstruct the 369-species discovery frame from `data/derived/global_monte_carlo_measured_photos_v1.csv` exactly as in Steps 1–3:

- `global_classifiable == true`;
- `morph` in `white`, `yellow_orange`, `red_pink`, `blue_purple`;
- at least 40 classifiable rows per species;
- `D = 1 - sum_k p_k^2`;
- `second_fraction` = frequency of the second-most-common of the four morphs.

The reconstruction must reproduce 369 species, maximum D approximately `0.707645`, second-morph >=10% fraction approximately `0.4661`, and >=20% fraction approximately `0.2656` before spatial subset outcomes are accepted.

## Primary and sensitivity subsets

Primary polymorphic subset:

`second_fraction >= 0.10`.

Sensitivity subset:

`second_fraction >= 0.20`.

These thresholds were already part of the species-level polymorphism summary before Step 5. No alternative threshold may be selected after opening Step-5 spatial outcomes.

## Primary statistic

For the primary subset, calculate the equal-species mean observed upstream rho. For each of the 999 frozen permutation indices, calculate the equal-species mean rho across exactly the same subset.

One-sided exact-randomization p-value:

`p_primary = (1 + number(null_mean >= observed_mean)) / 1000`.

The same calculation is repeated for the >=20% sensitivity subset.

## Dilution contrast

To distinguish “the polymorphic subset itself is spatially structured” from the stronger claim “the whole-frame signal was diluted by low-polymorphism species”, partition the 369 species by the primary >=10% threshold.

Calculate:

`delta = mean_rho(primary_subset) - mean_rho(complement)`.

For each frozen permutation index, calculate the same difference from the corresponding per-species permutation rhos. The one-sided p-value is:

`p_delta = (1 + number(null_delta >= observed_delta)) / 1000`.

The dilution interpretation is allowed only if `observed_delta > 0` and `p_delta < 0.05`.

## Continuous secondary diagnostic

As a threshold-free diagnostic, calculate Spearman correlation across all 369 species between discrete Simpson D and the species-level observed spatial rho.

For each of the 999 frozen permutation indices, correlate the same fixed D values with that permutation's species-level rho values. The one-sided p-value is defined analogously.

This diagnostic is secondary and cannot replace failure of the primary >=10% subset test.

## Frozen interpretation rule

- `POLYMORPHIC_SUBSET_SPATIALLY_STRUCTURED` if the primary >=10% subset has observed mean rho > 0 and `p_primary < 0.05`.
- Add `WITH_DILUTION_CONTRAST` only if `observed_delta > 0` and `p_delta < 0.05`.
- The >=20% subset is a sensitivity result and must be reported regardless of direction.
- A positive continuous D–spatial-rho diagnostic may be described as a gradient, but only as exploratory support.

Failure of the primary subset test does not invalidate the already completed 369-species omnibus; it means Claim 3 cannot be sharpened specifically to discrete polymorphic species with this threshold.

## Claim ceiling

Even if Step 5 passes, the result establishes photo-derived within-species geographic colour structure among species meeting the frozen morph-frequency threshold. It does not establish adaptation, a common boundary, a common environmental threshold, population genetic differentiation, or biological prevalence of spatial structuring.

No new images, species filtering, threshold tuning, distance metric changes, or new permutations are allowed.
