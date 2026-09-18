# disttrait v0.9 package status — 2026-09-19

Status: validation candidate with a second external, non-flower empirical transport.

Package:

`packages/disttrait/`

Candidate version:

`0.9.0`

## 1. New in v0.9

v0.9 adds an external transport that is deliberately different from the
San Francisco street-tree example.

Source:

- ShareTraitDatabase;
- frozen source commit `9190652410f782de399ef2d9757a80beffd7a076`;
- ShareTrait release DOI `10.5281/zenodo.16537297`;
- source dataset `TRADAT039`;
- source dataset DOI `10.17605/OSF.IO/56TNH`;
- publication: Shokri et al. 2022, Journal of Experimental Biology,
  DOI `10.1242/jeb.244842`;
- ShareTrait database license: CC BY 4.0.

Fixture:

`packages/disttrait/fixtures/sharetrait_gammarus_metabolic_v1.csv`

Provenance:

`packages/disttrait/fixtures/sharetrait_gammarus_metabolic_v1_metadata.json`

## 2. Empirical frame

The fixture contains **Gammarus insensibilis** individuals collected from three
Adriatic populations:

- Quarantia: 123 individuals;
- Lesina: 130 individuals;
- Acquatina: 122 individuals.

Total:

- **375 individuals**;
- **3 populations / sites**;
- metabolic-rate unit = **Joule/day**;
- missing metabolic-rate values = **0**.

Individuals inherit the site coordinate of their sampled population. The data
therefore support between-population geographic organization, not fine-scale
within-site spatial inference.

## 3. Frozen analysis

Runner:

`packages/disttrait/examples/sharetrait_gammarus_empirical.py`

Primary trait:

`log(metabolic_rate)`

Primary dissimilarity:

`absolute pairwise difference in log(metabolic_rate)`

Primary statistic:

`Spearman(pairwise geographic distance, pairwise trait difference)`

Calibration:

- 999 within-species vertex permutations;
- site coordinates fixed;
- the complete observed metabolic-rate multiset fixed;
- metabolic-rate values reassigned among individuals.

Raw metabolic rate is retained as a sensitivity analysis using the same
permutation schedule.

## 4. Result

Frozen receipt:

`results/disttrait_gammarus_empirical_v0_9_20260919/result.json`

### Primary log-rate result

- observed rho = **-0.00682997**;
- matched-null p = **0.68**;
- null 95% interval = **[-0.02797, 0.03269]**.

### Raw-rate sensitivity

- observed rho = **0.0109306**;
- matched-null p = **0.198**;
- null 95% interval = **[-0.02427, 0.03016]**.

The frozen result is therefore **non-support** for geographic organization of
individual metabolic-rate differences across these three sampled populations
under this disttrait estimand.

## 5. Why a null transport is useful

The first external transport, San Francisco street-tree DBH, produced a positive
species-conditioned spatial signal.

The second external transport produces a null result.

This is methodologically preferable to selecting only external datasets that
produce significance. Together the examples show that:

- the same continuous-trait inference layer transports across unrelated
  biological systems;
- the package can return positive or non-supporting empirical outcomes;
- external transport is being used to demonstrate executable generality, not
  to manufacture a universal ecological pattern.

## 6. Claim boundary

The Gammarus result does not show that:

- latitude has no effect on metabolic physiology;
- the Shokri et al. experimental conclusions are unsupported;
- temperature, body mass or climate are unimportant;
- the three populations are biologically identical;
- metabolic rate lacks geographic structure under every model.

The source study analyzes metabolic physiology with predictors including body
mass, temperature and collection site. The disttrait transport asks a different,
deliberately generic question: whether **absolute individual metabolic-rate
differences increase with pairwise geographic separation** under the frozen
three-site observational geometry.

## 7. Validation

Dedicated workflow:

`.github/workflows/disttrait-gammarus-transport.yml`

The workflow:

1. installs disttrait;
2. re-runs the 375-individual empirical analysis;
3. compares all primary/null statistics with the frozen v0.9 receipt;
4. requires both the log-rate primary and raw-rate sensitivity to remain
   non-supporting under the frozen result contract.

Initial exploratory run:

- PR #49;
- run `35395871877`;
- job `105764444293`;
- artifact `10567594920`;
- conclusion: success.

## 8. Methods-paper state after v0.9

The external evidence now contains two qualitatively different transports:

1. **San Francisco street trees** — multi-taxon, individual continuous DBH,
   fine-grained coordinates, positive matched-null result;
2. **ShareTrait Gammarus** — one species, three populations, individual
   metabolic rate, population-level coordinates, null matched-null result.

The largest remaining strengthening gaps are now:

1. nonlinear within-species processes;
2. multivariate continuous traits;
3. MNAR observation processes;
4. standalone package repository/licence/release metadata.

A second external empirical dataset is no longer an unmet gate.
