# FCP v2 post-Pass-B MV1–MV4 analysis freeze

Date frozen: 2026-09-24 JST  
Status: **FROZEN BEFORE ANY PASS-B BIOLOGICAL RESULT ARTIFACT EXISTS**

This document operationalizes the already frozen FCP v2 measurement-validity programme. It does not add a new ecological hypothesis, change any counterfactual, or authorize rescue analyses.

## Required handoff

The analysis may begin only after `fcp-v2-biological-seal-v1` exists and reports 40,000 terminal rows, 256/256 biological partitions, zero replacement rows, source-byte identity status for every row, and frozen technical table SHA256 `f757c90eddcb00f43d504a2ae6ccaeff4639d7ce468f231bc3b75675d66440e9`.

Rows with Pass-B acquisition failure or source-byte drift remain terminal rows and are never replaced.

## Frozen biological states and palette

Four coarse states: white, yellow_orange, red_pink, blue_purple. `mixed_uncertain` remains structural measurement missingness.

Nine-palette order: `[white, yellow, orange, red, pink, magenta, purple, blue, bronze]`.

The fixed H2 axis remains `q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`.

## MV1 — image-level invariance

All-row fixed-mask conditions compared with base are EV -1.0, -0.5, +0.5 and +1.0. EV 0 is an identity check.

For each condition report: base-classifiable rows; counterfactual-classifiable rows; paired-classifiable rows; classifiable→nonclassifiable rows; nonclassifiable→classifiable rows; coarse-state flip probability among paired-classifiable rows; the complete four-state transition matrix; Hellinger distance between normalized nine-palette vectors when both have positive palette mass; and the flower Lab delta-E already frozen in Pass T.

Every quantity is reported equal-image and equal-species. Equal-species means first summarizing within species, then taking an unweighted mean/median across species with at least one eligible row. No new effect-size threshold is introduced.

The same summaries are reported as explicitly secondary heavy-subset sensitivities for full-pipeline EV -1.0, -0.5, +0.5, +1.0, mid-grey background neutralization, and each of the seven deterministic prompt variants. Heavy-subset results never replace the all-row estimand.

## MV2 — species-level D invariance

For any condition, four-state diversity is `D = 1 - sum_k p_k^2` among classifiable rows.

For base versus each fixed-mask EV condition, a species enters the primary paired comparison only if it has at least 40 classifiable rows in both conditions. Report separately for Panel P, Panel N and pooled species: paired species count, Spearman rho, Lin-style CCC, mean absolute D difference, signed mean D difference and median absolute D difference. No pass/fail threshold is created for these effect-size summaries.

For heavy-subset counterfactuals, calculate D whenever a species has at least one classifiable row in the relevant 20-row subset and report the full classifiable-row denominator distribution. These are sensitivity summaries only.

Panel-P fresh-image transport uses the earlier frozen D for the same species without recomputation or refitting. The exact upstream D artifact and hash must be recorded before execution.

## MV3 — H2 geometry invariance

The H2 construction is unchanged: n_classifiable >= 40 for all-row conditions; coarse second-state fraction >= 0.10; normalized nine-palette rows; deterministic two-means in Hellinger coordinates; continuous minor-cluster fraction >= 0.10; nonzero Delta; `u_i = Delta_i / ||Delta_i||`; `W = mean_i (u_i^T q_white)^2`.

For base and each all-row fixed-mask EV condition report construction-gate counts, vector species count, W, W minus base W, vector-set retention relative to base, and species overlap with the base vector set. The axis is never refit and no new confirmatory structured-null p-value is required.

Also report base W after separately excluding each pre-biological technical stratum: high flower near-clip, high background near-clip, low flower-background Lab separation, and ROI unstable fixed. Strata are never recombined or retuned after biological opening.

Heavy-subset W-like summaries are secondary finite-sample sensitivities; their classifiable and vector denominators must be printed in full and they cannot replace the all-row MV3 estimand.

## MV4 — spatial measurement sensitivity

MV4 reuses the existing FCP species-specific spatial-organization statistic and geometry-preserving null without adding environmental, pollinator, phylogenetic, demographic or genetic predictors.

Primary comparisons are base biological measurement, ROI-stable-only measurement using the frozen Pass-T stratum, and the matched flower-minus-background response using already frozen technical flower/background measurements. Existing species-support rules and the existing spatial statistic implementation must be reused. Attrition is reported, not replaced.

Heavy 20-row spatial summaries are secondary and are not interpreted as the same finite-sample estimand as the 100-row base analysis.

## Ordering and no-rescue rules

Execution order after biological seal: `MV1 -> MV2 -> MV3 -> MV4`.

Forbidden after outcome opening: new EV values; new clipping thresholds; alternative neutral background; alternative ROI jitter; species/photo replacement; q_white refit; H2 null retuning; dropping an inconvenient technical channel; ecological-predictor rescue; or changing denominators after seeing results.
