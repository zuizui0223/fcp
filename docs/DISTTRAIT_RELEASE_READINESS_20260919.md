# disttrait 0.12.0 release readiness — 2026-09-19

Status: distribution artifacts validated; ownership/release metadata remains intentionally unresolved.

Package:

`packages/disttrait/`

Validated version:

`0.12.0`

## 1. Distribution build validation

PR #54 validates the package as a non-editable distribution.

Release workflow:

`.github/workflows/disttrait-release-check.yml`

Run:

- workflow run: `35412903843`;
- job: `105815872877`;
- artifact: `10573679286`;
- artifact name: `disttrait-0.12.0-distributions`;
- conclusion: **success**.

Validated steps:

1. build source distribution;
2. build wheel;
3. `twine check` on both artifacts;
4. create a fresh Python 3.11 virtual environment;
5. install the built wheel, not the source tree;
6. change working directory outside the repository;
7. import the installed `disttrait` distribution;
8. verify installed package version = `0.12.0`;
9. run basic diversity, contrast and multivariate-distance smoke checks.

Thus the current package is technically buildable and installable independently
of editable monorepo installation.

## 2. Scientific validation state

The methods-paper architecture now includes:

- exact compact equivalence with frozen FCP/RGFCA algorithms;
- observer-disjoint species-phenotype reliability;
- structured trait-geometry nulls;
- species-conditioned spatial permutation inference;
- scalar and multivariate continuous traits;
- pooled-confounding benchmarks;
- observation imbalance and MCAR surfaces;
- model-based comparators;
- direction heterogeneity;
- permutation-calibrated species-slope meta-analysis;
- explicit joint-MNAR boundary testing;
- nonlinear curvature;
- multivariate orientation heterogeneity;
- two external non-flower empirical transports.

The validated scientific package version is therefore `0.12.0`.

## 3. Remaining blockers to an external public release

These are not technical failures.

They require explicit ownership/release decisions:

- software licence;
- final author/maintainer names and ordering;
- final `CITATION.cff`;
- standalone repository name and visibility;
- external tagged release location;
- archival DOI;
- optional TestPyPI/PyPI publication.

The current GitHub connector does not expose repository creation, so a
standalone repository cannot be created from this chat with the available
GitHub actions.

## 4. Safe release sequence

1. choose licence and author/maintainer metadata;
2. create the standalone repository;
3. copy the validated `packages/disttrait/` tree plus relevant CI and frozen
   equivalence provenance;
4. rerun the same distribution build/install gate;
5. tag `v0.12.0`;
6. archive the tag and record the DOI;
7. add package-index publishing only after the previous metadata is frozen.

No PyPI credentials or public-release action should precede those decisions.
