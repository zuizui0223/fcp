# RGFCA Step 7F — state-construction null diagnostic

Date: 2026-09-10 JST
Status: prospectively frozen after Step 7E and before opening Step 7F null-world outcomes.

## Aim

Test whether the Step-7A photo-state classifier itself can generate an apparent association between pure S* state and the Step-7D total climatic hypervolume, even when morph labels are spatially exchangeable within species.

This is a post-result diagnostic. It does not rescue or reopen Step-7D directional hypotheses.

## Fixed inputs

- discovery measured-photo table: `data/derived/global_monte_carlo_measured_photos_v1.csv`;
- reserve measured-photo table: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`;
- Step-7D species environmental metrics, including `hv3d_rarefied20`;
- Step-7A eligibility and exact C*/S* definitions.

## Null construction

Within each species independently, retain the observed coordinates, observer IDs, sample size, and observed morph counts, and shuffle morph labels across photos.

For each eligible species use two independent fixed-seed permutation banks, each of size 999:

1. **calibration bank**: used only to construct the null distribution of segregation delta;
2. **evaluation bank**: used to construct null worlds.

This avoids evaluating a permuted dataset against a null distribution that contains the same permutation.

Eligibility is unchanged by permutation because sample size and morph counts are preserved.

## Reproduce C* exactly

For every evaluation shuffle, use the Step-7A C* rule at 100 km:

- >=30 observer-disjoint close pairs;
- >=5 cross-morph close pairs;
- >=4 unique photos participating in those cross-morph close pairs;
- >=3 unique observers participating in those cross-morph close pairs.

## Reproduce S* exactly

For each species, calibration shuffles generate the null distribution of

`median(distance among different-morph pairs) - median(distance among same-morph pairs)`.

For each evaluation shuffle, compute the same delta and a one-sided empirical upper-tail p-value against the 999-member calibration bank:

`p = (1 + count(null_delta >= evaluation_delta)) / 1000`.

Set S*=TRUE only when evaluation delta > 0 and p<=0.05.

## Null worlds

Evaluation permutation index b=1..999 defines null world b by taking evaluation shuffle b for every species. Within each null world classify each species into:

- pure C*: C*=TRUE, S*=FALSE;
- pure S*: C*=FALSE, S*=TRUE;
- mixed;
- unresolved.

For every world and tranche, calculate the Step-7D H0 statistic on real frozen hypervolume values:

`median(hv3d_rarefied20 | pure S*) - median(hv3d_rarefied20 | pure C*)`.

A world is valid only when both pure C* and pure S* are represented.

## Observed comparison

Observed Step-7D H0 deltas are fixed at:

- discovery: 9.37182419530179;
- reserve: 9.661780091846277.

For each tranche report:

- valid null worlds;
- median, 95% interval, and 99% interval of null-world H0 deltas;
- one-sided empirical p for null H0 >= observed H0;
- null distributions of pure-C* and pure-S* counts.

## Diagnostic classification

- `construction_sufficient`: observed H0 is not above the 97.5th percentile of the classifier-null distribution in both tranches.
- `construction_insufficient`: observed H0 is above the 97.5th percentile and one-sided p<=0.05 in both tranches.
- otherwise `mixed_or_unresolved`.

## Boundary

If construction is insufficient, the residual association still does not imply selection or adaptation. Remaining candidates include real sampled climatic breadth, sampling geometry not captured by radius95, taxonomic/biogeographic composition, and other omitted structure.
