# FCP distributed-polymorphism post hoc protocol — 2026-10-07

## Question

Does species-wide flower-colour ITV reside disproportionately among geographic localities rather than as unrestricted local coexistence?

This directly tests the geographic pattern highlighted by the FCP literature: a polymorphic species can consist largely of locally homogeneous populations arranged in clines or mosaics.

## Data

Use the exact frozen measured tables for:

- discovery: 500 sampled species, 369 D-evaluable;
- validation: 500 species-disjoint sampled species, 363 D-evaluable;
- third cohort: 499 newly sampled species, 377 D-evaluable.

Eligibility is unchanged: globally classifiable rows in white, yellow/orange, red/pink or blue/purple, with at least 40 rows per species.

## Species-wide diversity

For each species, finite-sample pairwise diversity is

D_pair = proportion of all unordered observation pairs whose coarse morph labels differ.

This is exactly the finite-sample-corrected Gini–Simpson diversity for the observed sample.

## Local diversity

For radius r, local diversity is

D_local(r) = proportion of unordered observation pairs separated by <= r km whose coarse morph labels differ.

A species is evaluable at radius r only if it has at least 30 local pairs.

Radii are frozen before analysis:

- 25 km;
- 50 km — primary scale;
- 100 km;
- 250 km.

No radius may be selected after viewing results.

## Composition-preserving spatial null

Within each species, keep fixed:

- all coordinates;
- exact coarse-morph counts;
- sample size;
- D_pair.

Randomly permute complete morph labels among fixed coordinates 199 times.

For every radius and species calculate local diversity in each null world.

Define observed local depletion:

depletion(r) = D_pair - D_local(r).

Positive values mean nearby observations are more homogeneous than the species-wide composition predicts.

## Cohort-level tests

### Test 1 — local depletion

For each cohort/radius, use the equal-species mean depletion.

Matched null worlds are formed by averaging the same statistic across species for each permutation index.

Support requires:

- observed mean depletion > 0;
- upper-tail matched-null p < 0.05.

### Test 2 — does greater species-wide D correspond to stronger partitioning?

For each cohort/radius calculate Spearman(D_pair, depletion).

For every matched null world calculate the same cross-species correlation using the null depletion for each species.

Support requires:

- observed rho > 0;
- upper-tail matched-null p < 0.05.

This directly tests the opportunity explanation because the null preserves each species' D exactly.

## Replication criterion

The primary 50-km result is considered replicated only if discovery and validation both support the same test.

Third-cohort transport is reported separately and cannot retroactively change the original 500+500 evidence.

Multi-scale support is descriptive: consistency across 25, 100 and 250 km strengthens interpretation but cannot rescue a failed 50-km primary test.

## Interpretation

If Test 1 replicates, the allowed statement is:

> Species-wide flower-colour diversity is geographically partitioned: nearby conspecific observations contain less colour-state diversity than expected from each species' overall colour composition.

If Test 2 also replicates:

> Species with greater species-wide colour diversity allocate a larger share of that diversity among geographic localities rather than within local neighbourhoods.

This would provide a direct comparative analogue of the clinal/mosaic FCP structure described in population-level literature.

## Hard nonclaims

This analysis does not distinguish:

- divergent selection from drift;
- restricted gene flow from habitat filtering;
- genetic differentiation from plasticity;
- local adaptation from nonadaptive geographic history.

Those require environmental, fitness or genomic evidence.
