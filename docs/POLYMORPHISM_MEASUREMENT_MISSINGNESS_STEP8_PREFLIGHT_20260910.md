# Polymorphism measurement missingness — Step 8 preflight

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-measurement-missingness-step8`
Parent: `analysis/polymorphism-genus-adjustment-decomposition-step7b` at `b831bbf9b4770d926aeb5595b658a54867c669f9`

## Purpose

Step 7b showed that genus-level clustering of photo-derived polymorphism diversity survives sampled-span adjustment in the independent reserve but is strongly attenuated by adjustment for `n_classifiable`.

Because every discovery and reserve species entered the image-measurement stage with a fixed 100-photo denominator, `n_classifiable` is not raw sampling effort. It is a post-acquisition measurement outcome that can mix clearly technical failures with colour-composition ambiguity.

Before any further D/genus or D/spatial association is computed, Step 8 first audits the reserve measurement schema and status vocabulary without reading or calculating any D association.

## Frozen source

`data/derived/rgfca_reserve_replication_measured_photos_v1.csv`

The preflight must inspect the complete file header and row-level status/classification fields only.

## Outcome firewall

The preflight is forbidden to:

- read a species-level `D` table;
- compute D from morph frequencies;
- join genus, sampled span, spatial rho, or any Step-7 outcome;
- stratify status counts by species, genus, D, spatial result, or sampled span;
- change any classification rule.

It may report only whole-file schema and whole-file value counts needed to define a later measurement-process decomposition.

## Candidate fields

Record the header exactly. For every column whose name contains any of

`status, class, evalu, morph, palette, roi, flip, fail, error, result`

(case-insensitive), report:

- non-missing count;
- number of unique non-missing values;
- complete value counts if <=50 unique values;
- otherwise the 50 most common values only, explicitly flagged as truncated.

Also report exact row count, unique species count, and whether every species has exactly 100 rows. These are geometry/denominator checks, not biological outcomes.

## Next-stage admission rule

A later Step 8 association test may define a **technical-failure rate** only from status values whose semantics unambiguously indicate acquisition/model/ROI/flip/runtime failure. Ambiguous palette composition, no-biological-palette-mass, or mixed/intermediate colour states must not be silently pooled into the technical-failure control unless the frozen row-level semantics establish that they are purely operational failures.

The later decomposition protocol must be committed only after this preflight result is frozen.
