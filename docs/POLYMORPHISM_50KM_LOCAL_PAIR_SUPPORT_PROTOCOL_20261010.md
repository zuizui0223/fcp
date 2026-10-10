# FCP 50-km local photo-pair support ladder: protocol (2026-10-10)

## Motivation

The current manuscript's original 50km visible flower-colour local-depletion result includes 166 discovery, 181 validation and 204 third-cohort species, conditional on >=40 classifiable photographs and >=30 local same-species photographic pairs. Previous composition-preserving adversarial synthetic relabelling shows 64 of 70 chosen photo-label exchanges were concentrated in species with fewer than 100 local pairs. This is a reason to audit photographic sampling leverage, NOT a finding of misclassification.

## Before real result: fixed source, thresholds, estimand

- Three exact historical high-depth measurement CSV sources, hardcoded original SHA256. No original photo labels, photos, geographic positions, or cohort memberships may be changed.
- Import the existing original distributed-polymorphism analysis and reuse the exact 50km site geodesic pair graph, 4-state mismatch statistic and 199 complete-vertex composition-preserving label shuffles.
- Original baseline must reproduce n=166/181/204 and species-equal absolute local depletion 0.020529254583812922/0.018672971642749295/0.01468491968437185, respectively.
- Six original local-pair minimum thresholds: **30, 50, 100, 200, 500, 1000**. No reporting of only a favourable threshold. Threshold 100 is the principal sparse-photo-pair sensitivity.
- For each threshold/cohort report species n, equal-species mean local depletion, relative depletion to observed species-wide photo pair discordance, observed positive-species fraction, 199-null unadjusted one-sided p and 1,999 species-bootstrap interval.
- A minimum of 30 species per cohort is required for a positive descriptive support decision. With fewer, retain the calculated sample diagnostic but mark **HOLD_INSUFFICIENT_SPECIES**.
- For each retained species, compute the maximum fraction of local colour-photo pair edges incident to any *single* photo. It is an edge-leverage diagnostic, not an independent population sample size.
- 4 synthetic numerical/falsification guards must pass before source files are recovered, and the workflow must check exact source immutable table SHA256 and original published baseline.

## Hard interpretation gates

Photo-pair edges share photographs, so 100 local pair edges are not 100 independent observations. Subsetting high pair-support species changes the eligible species pool, preventing causal comparisons between thresholds. Higher-level species-level geographic photographed colour clustering cannot distinguish adaptive evolution, genetically segregating morphs, individual-level plasticity, pollinator effects, or photographic morph-specific errors.

This is **post hoc**, not a new prospective experiment or confirmation of the mechanism. No H1/H2 or original source biological decisions are changed. The future selected 2,000+730 source taxa are untouched. Full 405-photo expert flower-organ/colour reannotation remains pending and cannot be replaced by sample-support sensitivity.
