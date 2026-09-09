# RGFCA South Florida M2 flower-side species-multipart stress test — result

Date: 2026-09-09 JST

Protocol frozen before outcome opening at `3222db6f96ac51a9581773abe890be8913ff0bbd`.

## Fixed target

Only the previously observer-consensus, flower-supported South Florida candidate was tested: reporting cell `(-16,6)`, mode `M2`, positive flower-only direction. Twenty fixed SHA256 species partitions (`salt=0..19`) were used; no salt, threshold, cell or mode was selected after outcome opening.

## Result

- salts tested: **20/20**
- salts where both species-disjoint halves passed: **11/20**
- pass fraction: **0.55**
- prespecified requirement: **>=16/20**
- `species_multipart_robust`: **false**
- weaker-half positive-fraction median: **0.75**
- weaker-half 2.5–97.5% range: **0.3725–1.00**
- weaker-half minimum: **0.35**

Thus the candidate is not stable to broad changes in which species constitute the local signal. Some partitions retain a strong positive direction in both halves, but others place the signal almost entirely in one half.

One salt had insufficient supported shifts in half A and remains in the denominator as a failed salt.

## Reproduction

A complete second execution produced byte-identical outputs for all four retained CSV files:

- `species_multipart_by_salt.csv`: `1b9a2b2a939a8dee0ab59d486bb3ff31cd3d91f55868de18c181d4147576198f`
- `species_multipart_global_summary.csv`: `bf76cad006b16d28ac24dbda9f47398e91381f35ee8e48430a6bc113757c9862`
- `species_multipart_realizations.csv`: `549f815fb24f41cbba368839f7e00926142e6f6bb7efa8e95b91222cd2cc833f`
- `species_multipart_species_scores.csv`: `8e804e1e45a3aa7e78a57b3c15bf8c02f8b2ed79bed7132e7e8faef2f0422922`

## Interpretation boundary

South Florida M2+ remains observer-consensus and flower-supported, but it is **not species-composition robust** under the prespecified multi-partition stress test. It should therefore not be promoted as a shared local biological rule from the current reserve. The stronger retained result is the observer-disjoint recurrence of the global M1–M3 signal axes themselves; local geographic placement remains more fragile.
