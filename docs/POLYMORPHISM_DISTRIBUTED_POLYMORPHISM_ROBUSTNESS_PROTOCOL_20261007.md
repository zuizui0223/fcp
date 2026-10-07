# FCP distributed-polymorphism robustness protocol — 2026-10-07

## Status

Post hoc falsification audit defined after the three-cohort local-depletion result was observed.

The established result to challenge is that, at 50 km, nearby conspecific observations contain less four-state flower-colour diversity than expected from each species' fixed overall colour composition.

## R1 — different-observer local pairs

Potential failure mode: observers are geographically clustered and observer/camera-specific colour bias could create local homogeneity.

For each species:

- retain all globally classifiable four-state rows;
- define local pairs as <=50 km;
- exclude every pair whose two photos have the same observer ID;
- require at least 30 retained local pairs.

Calculate the observed fraction of retained local pairs with different coarse morph labels.

Null: permute complete morph labels among fixed photo positions 199 times, leaving coordinates and observer identities fixed.

Species-wide pair diversity is unchanged.

Cohort support requires positive equal-species local depletion and matched-null p < 0.05.

## R2 — nonwhite-only three-state local depletion

Potential failure mode: the result is driven by the exposure-coupled white classifier or by a trivial achromatic/chromatic split.

For each species:

- remove all rows whose frozen coarse morph is white;
- retain yellow/orange, red/pink and blue/purple;
- require at least 30 retained nonwhite rows;
- define local pairs as <=50 km and require at least 30 such pairs.

Calculate finite-sample three-state species-wide pair diversity and local pair diversity.

Null: permute the three nonwhite morph labels among fixed positions 199 times.

Cohort support requires positive equal-species local depletion and matched-null p < 0.05.

## Interpretation

- R1 support in discovery, validation and third cohort reduces direct observer/camera-clustering explanations.
- R2 support in discovery, validation and third cohort shows that geographic partitioning is not restricted to white versus coloured states.

These are robustness diagnostics only. Neither establishes adaptation.

## R3 — continuous nine-colour local homogeneity

Potential failure mode: the distributed-polymorphism result is created by coarse morph classification rather than continuous flower-colour structure.

For each original D-evaluable species, use the nine biological flower-colour coordinates (legacy palette_count_*; third-cohort flower_fraction_*), normalized within rows.

Calculate mean pairwise Jensen–Shannon dissimilarity across all photo pairs and among pairs <=50 km. Require at least 30 local pairs.

Define continuous local depletion as species-wide mean JSD minus local mean JSD. Permute complete nine-colour rows among fixed coordinates 199 times.

Cohort support requires positive equal-species continuous depletion and matched-null p < 0.05.

R3 support in all three cohorts shows that geographic local homogeneity is not restricted to the coarse four-state representation.
