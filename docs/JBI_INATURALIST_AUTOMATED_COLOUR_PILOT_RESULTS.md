# Frozen automated-colour pilot result

## Decision

The location-blind measurement workflow completed, but the locked spatial confirmation did not support any of the three eligible species. The 50-species atlas therefore stops before opening its 20,200 candidate images until an independent flower-tissue localization benchmark and signal-recovery gate are frozen and passed.

## Frozen counts

- Development: 6 species, 480 encounters, 886 photographs; 3 species passed and 3 became `not_evaluable`.
- Locked: 3 species, 360 encounters, 717 photographs.
- Cache audit: 717 valid, 0 missing, 0 partial, 0 unexpected.
- Locked admission: 306 encounters and 508 photographs admitted.
- Coordinate firewall: opened only after all three species passed the locked completeness rule.

| Species | Admitted encounters | Primary rho | p | BH q | Flower-background rho | Contrast q | Decision |
|---|---:|---:|---:|---:|---:|---:|---|
| *Erythranthe lewisii* | 101/120 | 0.00774 | 0.3979 | 0.59685 | 0.04356 | 0.5106 | not detected |
| *Hesperis matronalis* | 97/120 | -0.05822 | 0.9637 | 0.9637 | -0.10022 | 0.9664 | not detected |
| *Orchis mascula* | 108/120 | 0.01996 | 0.3364 | 0.59685 | -0.08759 | 0.9664 | not detected |

All probabilities used 9,999 frozen within-species whole-vector permutations. Non-rejection is not proof of randomness.

## Immutable outputs

- `docs/supporting/jbi_inaturalist_automated_colour_locked_analysis_input_v2.csv`
- `docs/supporting/jbi_inaturalist_automated_colour_locked_analysis_input_manifest_v2.json`
- `docs/supporting/jbi_inaturalist_automated_colour_spatial_manifest_v2.json`
- `docs/supporting/jbi_inaturalist_automated_colour_spatial_v2/species_spatial_mark_results.csv`
- `docs/supporting/jbi_inaturalist_automated_colour_spatial_v2/descriptive_equal_pair_variogram.csv`

The six-species Chapter 1 Stage A/B results and the 34-species comparative analysis remain unchanged.
