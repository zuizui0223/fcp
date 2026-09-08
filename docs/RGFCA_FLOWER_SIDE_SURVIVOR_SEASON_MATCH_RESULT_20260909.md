# RGFCA surviving flower-side candidate season/year matched nuisance audit — result

Date: 2026-09-09 JST

Protocol frozen before outcome opening at `0d448eba4676d01c30d5e2bbd6c006cb8542b77a`.

## Fixed target

Only the previously surviving flower-side candidate was audited: reporting cell `(-23,11)`, M2 positive (approximately 116.6°W, 51.7°N). No new cell, mode, date threshold, or ecological predictor was selected.

## Primary result

The prespecified season/year match required same species, comparison state outside the candidate cell, >=500 km separation, circular day-of-year distance <=30 days, and mean observation-year difference <=5 years.

- supported realizations: **200/200**
- matched species per realization: median **18**, range **12–27**
- species-equal candidate-minus-matched M2 contrast: median **+0.032355**
- 2.5–97.5% realization range: **-0.009769 to +0.073693**
- fraction positive: **0.93**
- prespecified season/year-match robustness gate: **PASS** (requires >=100 supported realizations and positive fraction >=0.70)

## Prespecified sensitivities

| scenario | supported reps | median matched species | median contrast | 2.5–97.5% | positive fraction |
|---|---:|---:|---:|---:|---:|
| <=30 d, <=5 y (primary) | 200 | 18 | +0.032355 | -0.009769 to +0.073693 | 0.93 |
| <=30 d, no year restriction | 200 | 18 | +0.034441 | -0.008414 to +0.074874 | 0.92 |
| <=30 d, <=2 y | 200 | 16 | +0.026756 | -0.024823 to +0.076158 | 0.85 |

The direction therefore persists under both prespecified date sensitivities, although the realization distribution still includes negative values and the stricter-year sensitivity is less stable.

## Reproduction and integrity

The recurrence engine reproduced the frozen upstream 200-mode realization ledger with maximum numerical discrepancy `4.440892e-16`.

A complete second execution produced byte-identical outputs for all four retained CSV files:

- `season_match_realizations.csv`: `c506d3e86cb10ef0a0cc015342adc769c6ed3ca4b0b89978da0aa803adfd31fa`
- `season_match_species_contrasts.csv`: `ad5c73370aabf9d3bd8bfb60f9ab2519995efa34b26518246da01db4cd636f93`
- `season_match_summary.csv`: `1b92f77d92e030e3b6c729d6a5b9c1cfe1f8c177a12e04d526423fc38a756ad2`
- `upstream_mode_reproduction_check.csv`: `541b45fd47e27ca3e22a91106f8c09f1e558bf3e29255d87ac2df6e82615fe1f`

## Interpretation boundary

Within the current reserve, the sole taxonomically robust M2+ local candidate is not readily removed by coarse matching on flowering season and observation year. This does **not** remove illumination, exposure, camera processing, local habitat/background, phenological stage within the 30-day window, or other photographic nuisance. It is not independent ecological replication and does not establish pollinator perception, adaptation, or a biological transition boundary.
