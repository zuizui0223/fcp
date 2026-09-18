# H2 third-cohort prospective transport — preopening freeze

Date: 2026-09-16 JST
Status: PREOPENING ONLY — biological outcomes must remain unopened
Parent evidence state: `analysis/p500-postmortem-no-rerun-20260916` at `ee59e48cbcb770a6a8d8df5e4f7e96003e26f8e4`

## Purpose

Test transport of the already-fixed white-versus-nonwhite H2 geometry in a genuinely new cohort. This study is not a replay, rescue, recovery, or reinterpretation of P500. P500 run `34953307374` remains terminal and closed.

## Independence firewall

The third cohort must be selected without using flower-colour measurements, morph labels, D, H2 vector geometry, W, structured-null p-values, or any recovered/recomputed P500 H2 outcome.

The cohort must exclude every species used in the original discovery cohort, reserve cohort, and P500 frozen cohort. P500 measured rows and its sealed measurement table are prohibited inputs.

Eligibility may use only outcome-blind opportunity information needed to establish that the frozen high-depth design is feasible. Species replacement after outcome opening is forbidden.

## Frozen biological estimand

No axis search or refit is permitted.

`q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`

Statistic:

`W = mean_i (u_i dot q_white)^2`

The palette, mode construction, orientation convention, structured-null construction, and construction-preserving constraints must be inherited unchanged from the frozen H2/P500 protocol. No non-white hue target may be introduced as a confirmatory endpoint.

## Frozen admission and analysis tiers

Primary tier: the existing 0.10 species-delta/admission rule.

Strict sensitivity tier: the existing 0.20 species-delta/admission rule.

The exact thresholds, seeds, structured-null replicate count, decision rule, and no-axis-refit rule must be copied mechanically from the P500 authorization/executor before cohort outcomes are opened. If any value cannot be recovered unambiguously from the frozen protocol, execution stops before biological opening; it must not be guessed or retuned.

## Technical qualification gate — mandatory before biological opening

The next prospective execution is not authorized until an outcome-free end-to-end qualification passes. The qualification must exercise the same production serialization and artifact path with synthetic/dummy data only and demonstrate all of the following:

1. H2 result object can be constructed through `H2_COMPLETE`.
2. Relative and absolute measurement-result paths serialize without exception.
3. Final JSON is written and can be read back with `allow_nan=False`-compatible finite numeric fields.
4. Terminal result validation passes on the written file, not only on an in-memory object.
5. Artifact packaging/upload path is exercised with a dummy artifact and its expected filename is verified.
6. The regression tests added after P500 remain green, including `tests/test_p500_h2_result_serialization.py`.
7. Qualification code has no access to third-cohort biological colour outcomes.

A unit test of `repo_display_path()` alone is necessary but not sufficient for this gate.

## Chronology requirement

The immutable order is:

1. freeze this protocol;
2. freeze third-cohort species/row identifiers and opportunity-only eligibility receipt;
3. pass and persist the synthetic end-to-end technical qualification;
4. freeze one bounded execution authorization with exact commit SHA, inputs, seeds, thresholds, artifact names, and no-rerun rule;
5. only then open/acquire biological colour measurements;
6. execute once;
7. persist result JSON, validation receipt, and artifact digest before interpreting H2.

If the one authorized biological run fails after outcome access, no rerun/replacement/recovery may be promoted to untouched prospective evidence. Technical repairs remain future-only.

## Decision classes

Only a durably serialized and validated result from the authorized new cohort can be assigned a prospective biological verdict. Technical/support failure is not a biological negative. Missing output is not a biological negative.

The manuscript claim ceiling changes only if the third cohort produces an admissible terminal H2 result under this contract.

## Immediate next gate

Build an outcome-blind candidate ledger from the 42,111-species frame, subtract discovery + reserve + P500 species, apply only frozen opportunity/support eligibility, and freeze the resulting candidate universe before any colour outcome is inspected.
