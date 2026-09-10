# RGFCA Step 7A — photo-derived local co-occurrence and spatial segregation

Date: 2026-09-10 JST
Status: prospectively frozen before opening C*/S* outcomes.

## Goal

Separate two independent properties of flower-colour polymorphism from the frozen RGFCA photos:

- **C\***: repeated local-scale co-occurrence of different admitted colour morphs;
- **S\***: whole-range spatial segregation of admitted colour morphs.

These are photo-derived organization descriptors. C* is **not** claimed to demonstrate within-population coexistence; S* is **not** claimed to demonstrate selection or genetic differentiation.

## Frozen inputs

Discovery: `data/derived/global_monte_carlo_measured_photos_v1.csv`.
Reserve: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`.
Literature validation only after photo-derived states are computed: `data/frozen/frozen_34species_coexistence_segregation_v22.csv`.

Only rows with `global_classifiable == True` and `morph` in `white`, `yellow_orange`, `red_pink`, `blue_purple` are admitted. Species require >=40 admitted rows. Organization analyses are restricted to species with second-morph fraction >=0.10; the >=0.20 subset is a sensitivity descriptor only.

## C* — repeated local co-occurrence

Primary local radius: **100 km**. Sensitivity radii: 250 and 500 km.

For all unordered photo pairs within a species, compute great-circle distance. C* support at 100 km requires all of:

1. at least 30 close pairs whose two photos come from different observers;
2. at least 5 close, observer-disjoint pairs with different morph labels;
3. those cross-morph pairs involve at least 4 unique photos and at least 3 unique observers.

Report the local cross-morph pair rate and the exact global cross-morph pair rate. No enrichment threshold is used to define C*: C* records repeated local co-occurrence, not over-mixing relative to a null.

## S* — spatial segregation

For each eligible species compute

`delta = median(distance | different morph) - median(distance | same morph)`.

Positive delta means same-morph observations are geographically closer than cross-morph observations.

Permute morph labels among the fixed coordinates **999 times**, preserving morph counts, and recompute delta. S* is supported when:

- observed delta > 0; and
- one-sided permutation p <= 0.05.

No distance threshold is tuned for S*.

## Four organization states

- `local_cooccurrence_only`: C*=1, S*=0
- `spatial_segregation_only`: C*=0, S*=1
- `cooccurrence_and_segregation`: C*=1, S*=1
- `unresolved_photo_state`: C*=0, S*=0

`unresolved_photo_state` is not biological absence.

## Literature validation firewall

The 34-species literature C/S table is not used to tune radii, pair-count gates, or the S* statistic. After all photo-derived states are frozen, intersect names and report agreement separately for literature C and S labels. Low overlap or poor agreement does not trigger threshold changes.

## Replication rule

Discovery and reserve are evaluated independently. Report state frequencies and continuous C*/S* metrics separately. If the same species occurs in both, state agreement is descriptive; no threshold is altered.

## Next-stage gate

Step 7B minimum-data rarefaction is allowed only after Step 7A is complete. Step 7B will ask how many species and photos/species are needed to recover full-data D, C*/S* composition, geographic organization, and the global atlas. Step 7A results may determine feasibility but not redefine Step 7A.