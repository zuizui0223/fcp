# GloBI frozen-12 source-semantics audit — v0.1

## Decision

The schema-only census produced 12 fixed candidates. After reading only README, GloBI mapping configuration and citation/sampling descriptions, **two advance as fresh metadata candidates**:

1. **CarniDIET** — population-level carnivore diet records with explicit site/time, coordinates, protocol, quantification and sample-size metadata.
2. **CropPol / OBservData** — crop-pollinator data with field/study/site geometry plus explicit pollinator sampling effort.

No interaction row or partner identity was opened.

## Why the other ten do not advance as fresh evidence

- **BMPO** — integrated parasite-host association databases; spatial presence is not a standardized local-network rewiring denominator.
- **GATEWAy (Brose)** — valuable food-web compilation, but network architecture is already analyzed in the associated paper; calibration only.
- **Fricke 2020** — global plant-frugivore network homogenization is itself the published outcome; calibration only.
- **GlobalAMFungi / GlobalFungi / Limbu 2025** — occurrence/metabarcoding or dominant-plant context is not automatically a realized pairwise interaction network.
- **OSAL** — opportunistic museum specimen associations, not standardized local networks.
- **Chile pollination catalogue** — heterogeneous literature catalogue and deprecated GloBI configuration.
- **Russo 2022** — excellent effort-resolved plant-pollinator data, but network diversity is already analyzed; calibration only.
- **SCAR diet** — heterogeneous compiled diet studies without one common local-network sampling denominator.

The next gate is therefore bounded to CarniDIET and CropPol. No candidate can be swapped in because its biology looks appealing later.
