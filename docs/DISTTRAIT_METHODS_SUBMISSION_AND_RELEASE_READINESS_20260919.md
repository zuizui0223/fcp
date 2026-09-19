# disttrait methods-paper and standalone-release readiness — 2026-09-19

Status: scientific benchmark stack complete for the tested families; manuscript and software distribution are technically reproducible; public release metadata remains intentionally unresolved.

## 1. Methods manuscript state

Canonical draft:

`docs/DISTTRAIT_METHODS_MANUSCRIPT_V0_1_20260919.md`

Current counts:

- Abstract: **225 words**;
- Introduction: **607 words**;
- Materials and Methods: **1,414 words**;
- Results: **799 words**;
- Discussion: **1,012 words**;
- main text, Introduction–Discussion: **3,832 words**;
- total working manuscript including references/legends: **5,447 words**.

Display package:

- **5 figures**;
- **1 method–estimand table**;
- reporting-only generator:
  `scripts/analysis/make_disttrait_methods_figures.py`;
- validated display artifact:
  workflow run `35417532675`, artifact `10576197828`;
- source manifest:
  `docs/figures/disttrait_methods_20260919/manifest.json`.

The manuscript is therefore no longer blocked by missing simulation classes, external transport or display architecture.

## 2. Scientific validation stack

Validated package version:

`disttrait 0.12.0`

The frozen evidence stack includes:

1. exact compact FCP/RGFCA implementation-equivalence fixtures;
2. observer-disjoint species-level phenotype reliability;
3. species-conditioned matched spatial nulls;
4. pooled-confounding demonstration;
5. effect × imbalance × MCAR power/calibration surface;
6. alternative species-conditioned aggregation comparators;
7. correctly specified binary logistic comparator;
8. continuous signed-direction heterogeneity;
9. permutation-calibrated species-specific slope/meta inference;
10. two external non-flower empirical transports;
11. explicit joint trait-by-location observation-process boundary;
12. nonlinear curvature benchmark;
13. multivariate continuous-trait orientation benchmark.

The paper claim is bounded to the tested estimands and observation-process assumptions. Universal superiority, universal type-I control and arbitrary MNAR robustness remain hard nonclaims.

## 3. Distribution validation

Existing distribution gate:

`.github/workflows/disttrait-release-check.yml`

Validated actions:

- source distribution build;
- wheel build;
- Twine metadata check;
- fresh-environment wheel installation;
- import and smoke test outside the monorepo working directory.

Recorded validated run:

- run `35412903843`;
- job `105815872877`;
- artifact `10573679286`;
- conclusion: success.

## 4. Standalone export dry-run

Exporter:

`scripts/release/export_disttrait_standalone.py`

Workflow:

`.github/workflows/disttrait-standalone-export.yml`

The exporter creates an internal standalone candidate containing:

- package metadata and README;
- runtime source;
- examples;
- benchmarks;
- scripts;
- fixtures;
- self-contained unit tests;
- source-commit provenance;
- a SHA256 manifest;
- an explicit record of public-release metadata that remains unresolved.

Monorepo regression tests that require frozen FCP/result receipts are deliberately not rewritten as standalone unit tests. Their relationship to the package is preserved in provenance rather than silently severed.

The standalone dry-run is not a public release. It deliberately fails the conceptual public-release gate until ownership metadata is resolved.

Validated dry-run record:

- PR #56;
- workflow run `35417843905`;
- job `105829812507`;
- artifact `10575954565` (`disttrait-standalone-0.12.0-candidate`);
- export: success;
- self-contained tests: success;
- sdist/wheel build: success;
- Twine validation: success;
- fresh-environment wheel install/import: success;
- candidate archive: success.

## 5. Remaining public-release blockers

Only decisions that should not be guessed remain:

- software licence;
- final package author/maintainer names and ordering;
- final `CITATION.cff`;
- standalone repository name and visibility;
- final repository URLs in `pyproject.toml`;
- first tagged release;
- archival DOI;
- optional TestPyPI/PyPI publication.

These are ownership/release decisions rather than missing analysis or software-build work.

## 6. Remaining manuscript-side tasks

The draft can now be revised as a paper rather than extended as a benchmark programme.

Remaining manuscript work is editorial:

- choose target journal and adapt article structure/length;
- freeze author list and affiliations;
- decide whether package and methods manuscript share identical authorship;
- add acknowledgements/funding/contributions/conflicts;
- replace monorepo code-availability language with the final standalone repository and release DOI once created;
- final reference and figure-caption style pass.

## 7. Current decision boundary

The next technically executable work after this audit is limited.

A public standalone release cannot be completed without choosing at least the software licence and ownership metadata. Automation should not infer those choices from the FCP repository owner or manuscript history.
