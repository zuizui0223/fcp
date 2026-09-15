# FCP H1–H3 publication status

Status: IN DEVELOPMENT — not submission-ready.

## Current scope

The paper concerns reproducibility and geometry of photo-derived within-species
colour diversity. The 34-species literature comparison is retired from the main
paper and supplement. Six-species Stage A/B remains legacy. Neither is a source
of support for H1–H3, and historical frozen data/results must not be edited.

## Evidence and claim ceiling

- H1: direct observer-disjoint four-state D reproducibility supports a sampling
  reliability claim. The direct result from c34c5e8 is now retained unchanged in
  `results/polymorphism_h1_observer_disjoint_d_20260913/result.json`.
  The [active manuscript](FCP_H123_MANUSCRIPT.md) and
  [companion supplement](FCP_H123_SUPPLEMENT.md) record its source hierarchy.
  The older repeated-partition and stricter failed split remain visible; do not
  select a primary test according to its favorable result.
- H2: white-versus-nonwhite targeted geometry is supported in already opened
  cohorts, but the target was chosen after the broad geometry was inspected.
  Do not call it prospective confirmation or an identified biological mechanism.
- H3a/H3b: broad phylogenetic signal and sampled photographic span did not
  replicate. No post-hoc predictor search is authorized by this status document.
- Measurement: the ROI estimator pools flower regions without verified focal
  petal identity. Monarda localization qualification failed; repeatability does
  not cure this limitation. Highlight/exposure confounding remains unresolved.

## Required publication work

1. Reconcile the latest source-backed H1 hierarchy and make one active manuscript
   with a matching supplement and traceable numerical/figure inventory.
2. Audit the biological measurement claim, retaining failed external validation
   and distinguishing observed image states from verified biological morphs.
3. Qualify P500 technical-table recording, fixed high-clip membership, response
   join barrier and within-species coupling implementation on synthetic inputs.
   Freeze unresolved implementation choices before real responses.
4. Establish evidence-backed pre-opening chronology. The present acquisition
   JSON is not the strict execution record. Missing historical false flags
   cannot be invented. P500 opening remains unauthorized.
5. Only if all gates pass, execute the frozen prospective test with complete
   attrition accounting. Otherwise report the unresolved/negative outcome and
   assess the manuscript at its supported, narrower measurement scope.
6. Verify references, main/supplement consistency, figures and reproducibility
   before calling the package submission-ready. Journal acceptance is not
   guaranteed; submission itself requires separate user authorization.

## Development checkpoint

H3 Methods and Supplement Tables S3–S4 now preserve exact cohort/scenario
bindings for raw, corrected and opportunity-adjusted results. Tree coverage
matches the preflight manifest, including 22 unmatched reserve species that
are not coded as absent signal. Source inspection records permutation counts,
residualization order and secondary-statistic roles. The text explicitly
distinguishes non-support from equivalence or demonstrated power. Tree-model
reruns, full statistical references and submission-package closure remain open.

Main Figure 1 now presents frozen H1 reserve split comparisons, all four H2
targeted comparisons, all H3a placement scenarios and both H3b cohorts. It
retains the later strict H1 failure and H3 non-replication, distinguishes
bootstrap intervals from partition/null ranges, and invents no unavailable
H3 confidence intervals. PNG and the exported vector PDF were visually
inspected after correcting overlapping footer text. Tests check the committed
plot-data bindings and same-runtime duplicate rendering. Supplementary
diagnostic figures and a full manuscript rendering are still required; this
overview alone is not the complete submission package.

H2 methods now distinguish square-root coordinates for clustering from original
composition coordinates for displacement, fixed equal-species weighting, and
the 999-world plus-one Monte Carlo rule. Supplement Table S2 preserves the
structured-null medians and observed excesses, not merely p-values. Source
inspection also identified an important scope limit: null worlds retain the
observed selected species without reapplying the continuous minor-cluster gate,
and permutations do not preserve observer/season/geographic strata. The
manuscript reports this limitation without asserting a measured bias or changing
the historical null. Stronger biological calibration remains unresolved.

The active manuscript now reports the completed Monarda gate failure in its
abstract, Methods and Results, and Supplement Table S1 carries all three
operational criteria. Publication tests independently recompute its pooled and
median values from the unchanged count table and check exact source hashes,
complete image census and labeled table cells. This strengthens the measurement
limitations record; it neither repairs the failed gate nor establishes a
general atlas error rate. The reference-target mismatch and annotation
completeness remain unresolved rather than assigned a convenient explanation.

Technical snapshot implementation now exists in
`fcp_pipeline/p500_technical_snapshot.py`: it accepts only photo identity and
near-clip fraction, requires exact coverage of a supplied evaluable identity set,
applies the unchanged high-clip rule, and validates canonical bytes against a
separately supplied trusted digest. Synthetic tests cover missing/duplicate rows,
response-bearing columns, invalid values, row-order invariance and altered
membership. These are integrity tests, not historical non-access evidence.
An exclusive-create writer now refuses existing paths and validates the snapshot
before writing. A complete one-to-one response join verifies the snapshot before
iterating responses; missing, duplicate or unknown IDs are rejected rather than
silently reducing the cohort. Artificial-data tests only have exercised this
route. This is not tamper-proof storage or authenticated chronology; a failed
write may leave a partial file that must be retained as a failure, not repaired
silently. There is no acquisition or opening token.
The [coupling implementation specification](P500_COUPLING_IMPLEMENTATION_SPEC.md)
now records a pure species-conditioned fitter tested against analytic matched
pairs and an independently enumerated multi-photo likelihood. Separation,
numeric warnings and failed curvature checks produce no estimate. All fitting
so far is artificial-data only. This does not establish observer-robust coverage,
measurement validity, production-scale qualification or chronology.
A caller could still supply a wrong cohort or an untrusted digest; those
must be bound by the future execution recorder and chronology qualification.
Do not treat this partial implementation as a complete P500 execution gate.

Starting point: PR 33 head 9e87458185bba8b6df901280b7a7078ce3544de8.
Its dedicated 44 tests and manuscript checks passed; legacy reproduction also
passed but is not a publication advancement for this paper. The new work is on
a separate publication branch. Manual-only legacy reproduction is intentional
retirement, not a bypass for a failed check. Repository branch-protection rules
have not been changed and may require separate administrator reconciliation.
