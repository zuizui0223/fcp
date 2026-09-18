# disttrait v0.6 package status — 2026-09-19

Status: validated candidate with external non-flower empirical transport.

Package:

`packages/disttrait/`

Version:

`0.6.0`

## 1. What v0.6 adds

v0.6 retains the v0.5 validation stack and adds two linked capabilities.

### Continuous scalar trait support

New public API:

- `absolute_pairwise`
- `continuous_spatial_rho`
- `continuous_spatial_permutation_null`

For a scalar individual-level trait `z`, the within-group spatial statistic is

`rho_i = Spearman(pairwise geographic distance, |z_j - z_k|)`.

The matched null fixes coordinates and the full multiset of observed trait values while permuting the assignment of complete trait values to positions.

This extends the package beyond categorical/compositional traits measured with Jensen-Shannon divergence.

### External empirical transport

The external transport uses a fixed historical snapshot of the San Francisco street-tree inventory distributed through TidyTuesday.

Frozen fixture:

- `packages/disttrait/fixtures/sf_street_trees_v1.csv`
- `packages/disttrait/fixtures/sf_street_trees_v1_metadata.json`

Source snapshot:

- repository: `rfordatascience/tidytuesday`;
- path: `data/2020/2020-01-28/sf_trees.csv`;
- Git blob: `bdc06c1297b7dd88bea0df77de4007eaab30198e`;
- upstream provider: San Francisco Public Works;
- upstream license: Open Data Commons Public Domain Dedication and License (PDDL) 1.0.

The fixture is a historical 2020 snapshot and is not used as a current census of San Francisco trees.

## 2. Outcome-blind transport design

Rows are admitted only when:

- DBH > 0;
- latitude and longitude are finite;
- tree ID is nonblank;
- the source taxon label contains at least two tokens;
- the generic `Tree(s)` label is excluded.

Taxon labels are inherited from the municipal inventory and are not taxonomically revalidated.

Among the 52 taxon labels with at least 500 usable rows, exactly 20 are selected by deterministic FNV-1a ordering of

`DISTTRAIT_SF_SPECIES_V1|taxon_label`.

Within each selected taxon, 80 observations are selected by deterministic FNV-1a ordering of

`DISTTRAIT_SF_TREE_V1|tree_id`.

The resulting fixed empirical frame contains:

- 20 taxon labels;
- 80 observations per taxon;
- 1,600 total observations.

Neither DBH values nor geographic coordinates enter the selection order after the usability gate.

## 3. Empirical estimand

Trait:

`natural log(DBH)`

Within-taxon dissimilarity:

`absolute difference in log(DBH)`

For every taxon:

1. calculate all pairwise great-circle distances;
2. calculate all pairwise absolute log-DBH differences;
3. compute Spearman rho between the two vectors;
4. generate 99 matched vertex-permutation nulls.

Across taxa, the primary empirical transport statistic is the equal-taxon mean spatial rho against the matched null worlds.

A secondary distribution-spatial diagnostic asks whether the standard deviation of log(DBH) across individuals predicts the taxon's spatial rho.

## 4. Frozen result

Canonical receipt:

`results/disttrait_sf_street_tree_empirical_v0_6_20260919/result.json`

Dedicated runner:

`packages/disttrait/examples/sf_street_tree_empirical.py`

### Species-conditioned result

- taxa = **20**;
- observations = **1,600**;
- mean within-taxon spatial rho = **0.0921946**;
- matched-null upper-tail p = **0.01**.

Thus the package detects positive geographic organization of individual DBH variation in this external continuous-trait dataset.

### Distribution amount versus spatial organization

- rho(SD log(DBH), spatial rho) = **-0.186466**;
- matched-null p = **0.83**.

There is no supported evidence here that taxa with greater DBH spread are more spatially organized.

### Naive pooled comparison

Ignoring taxon identity across all 1,600 selected trees gives:

- pooled rho = **0.0278272**;
- nominal p = **1.68e-217**.

This is not treated as evidence that the pooled method is biologically stronger. The pooled and species-conditioned analyses target different estimands, and the enormous pair count makes the ordinary pooled p-value unsuitable as an independent-replication count.

## 5. What v0.6 establishes

v0.6 closes the largest remaining transport gap identified after v0.5:

> the reusable spatial layer operates unchanged on an external, non-flower, continuous individual-level trait dataset with real geographic coordinates.

This external demonstration is independent of the FCP flower-colour measurement representation and does not use the RGFCA discovery/reserve data.

## 6. What v0.6 does not establish

The street-tree result does not identify why DBH is spatially organized.

Possible contributors include:

- planting history;
- neighborhood design;
- taxon/cultivar deployment;
- tree age;
- management;
- pruning/removal history;
- environmental conditions;
- municipal inventory processes.

Therefore v0.6 is an **empirical transport validation**, not an ecological causal analysis of urban tree size.

It also does not establish:

- universal external validity;
- taxonomic correctness of every municipal label;
- suitability of absolute trait difference for every continuous phenotype;
- robustness to arbitrary MNAR sampling;
- superiority over all model-based spatial methods.

## 7. Validation gates

Full package workflow:

`.github/workflows/disttrait-package.yml`

Dedicated empirical workflow:

`.github/workflows/disttrait-empirical-transport.yml`

Both workflows:

1. install `disttrait`;
2. run the relevant package/transport analysis;
3. regenerate the empirical result;
4. compare it against the frozen v0.6 receipt.

The first successful empirical transport run was:

- PR #46;
- workflow run `35381530119`;
- job `105718578832`;
- artifact `10563155571`;
- conclusion: success.

## 8. Methods-paper state after v0.6

The package now has:

1. mathematical unit tests;
2. synthetic integration tests;
3. exact compact FCP/RGFCA algorithm-equivalence fixtures;
4. optional raw-artifact replay hooks;
5. non-flower synthetic generalization;
6. confounding benchmark;
7. 24-cell performance surface;
8. alternative species-conditioned comparators;
9. a model-based fixed-effect logistic comparator;
10. **external non-flower empirical transport on a continuous trait**.

The next major methods-paper gap is no longer simple external transport. It is broader external replication and model-misspecification coverage, especially continuous/multivariate processes where the binary logistic comparator is intentionally wrong.
