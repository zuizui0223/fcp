# FCP H1–H3 companion supplement

Working evidence inventory — incomplete, not submission-ready.

## S1. H1 source hierarchy and chronology correction

The direct observer-disjoint protocol was committed as
743e13d3597c88ea7a725c596ea733b30d978024 at 2026-09-13 01:53:54 JST.
The repeated-partition protocol was committed as
92488ada5f4f535f6731bba6519101d8ed8a7e1a at 2026-09-13 02:03:49 JST.
Thus the older evidence-ledger description of the repeated-partition protocol
as the earliest protocol is not used in this manuscript. Commit chronology
alone is not authenticated evidence of every historical non-access event.

Direct result: `results/polymorphism_h1_observer_disjoint_d_20260913/result.json`,
copied without modification from c34c5e85504bee76aba7c5f780ef162863aef09e.
CI run 34706708406, artifact 10302770009. This is the direct current-D result,
not the historical approximate 0.971 summary.

Other H1 analyses remain retained:
`results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`
and `results/polymorphism_h1_observer_disjoint_D_reliability_20260913/result.json`.
Their different partition rules and the latter's failed threshold must remain
visible. No frozen result or threshold is changed by this correction.

### Table S5. Direct H1 value agreement and finite-sample sensitivity

| Cohort | Outcome | Species | Spearman rho | Lin CCC | Median absolute half difference | 90th percentile difference |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| discovery | D | 369 | 0.837015 | 0.866192 | 0.060547 | 0.176005 |
| discovery | D_unbiased | 369 | 0.836795 | 0.867324 | 0.063492 | 0.183023 |
| reserve | D | 363 | 0.810916 | 0.858294 | 0.062435 | 0.172113 |
| reserve | D_unbiased | 363 | 0.811614 | 0.859420 | 0.064615 | 0.178916 |

These are diagnostics from the direct split, not from the repeated or later
strict split. CCC is `2*cov(A,B)/(var(A)+var(B)+(mean(A)-mean(B))^2)`, using
population-divisor empirical variances/covariance in the retained implementation.
It penalizes location/scale disagreement as well as imperfect correlation;
Spearman rho alone does not do so. Absolute differences are in D units, not
error relative to verified biological ground truth. D_unbiased multiplies each
half's D by its own n/(n-1); it does not calibrate image classification.

The [direct protocol](POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260913.md),
[preflight implementation](../scripts/analysis/run_polymorphism_h1_observer_split_preflight_20260913.py)
and [direct analysis implementation](../scripts/analysis/run_polymorphism_h1_observer_disjoint_d_20260913.py)
were read, not executed, for this methods audit. The assignment and eligibility
panel hashes are verified by the runner before reading state labels. The
opportunity gate retained all 369/363 eligible species at 20 per half. This does
not mean every source photograph contributes: the source cohorts each contain
50,000 photos, but only 25,377 discovery and 24,885 reserve rows carry the
classifiability flag before the species/observer/half filters. Those are source
flag counts, not asserted final half-analysis totals.

The direct bootstrap samples paired species indices, recomputing ranks after
resampling. It accepts at least 99% finite draws; both raw-D cohorts report all
5,000 finite draws. The same random generator then permutes the fixed half-B
ranks 20,000 times, with a 1e-15 comparison tolerance. Seeds for raw D are
20260914 (discovery) and 20261014 (reserve), with one added for corrected D.
No bootstrap or permutation was rerun to create this table. The source hashes,
test counts and agreement statistics do not prove historical non-access or
independence across species. A globally observer-disjoint split or clustered
uncertainty would require separately specified designs and qualification, not
silent replacements of the frozen analysis.

## S2. H2 source and retrospective status

`results/polymorphism_white_axis_targeted_test_20260912/result.json` records the
targeted primary/strict results and its retrospective claim boundary.
CI run 34675697583, artifact 10292077657. The coarse-state-preserving null
must not be replaced with the easier isotropic comparison.

### Table S2. White-axis alignment against the structured construction null

| Threshold | Cohort | Species | Observed W | Null median | Excess over median | Upper-tail p |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Primary 0.10 | discovery | 152 | 0.514625 | 0.430808 | 0.083817 | 0.001 |
| Primary 0.10 | reserve | 129 | 0.514586 | 0.466546 | 0.048039 | 0.001 |
| Strict 0.20 | discovery | 75 | 0.542355 | 0.443626 | 0.098729 | 0.001 |
| Strict 0.20 | reserve | 65 | 0.510517 | 0.469943 | 0.040574 | 0.008 |

These medians describe 999 null worlds, not uncertainty intervals for observed
W. Excess is observed W minus null median. A p-value of 0.001 is the minimum
attainable with the plus-one rule, not zero probability. The selected species
are weighted equally regardless of photographic depth. The fixed contrast is
not a learned pigment axis, and W is not polymorphism prevalence.

### Operational source audit

The [target freeze](POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md),
[observed geometry code](../scripts/analysis/run_polymorphism_delta_geometry_validation_20260912.py),
[structured-null code](../scripts/analysis/run_polymorphism_delta_geometry_structured_null_20260912.py)
and [targeted runner](../scripts/analysis/run_polymorphism_white_axis_targeted_test_20260912.py)
were inspected, not rerun, for this methods description. The runner uses seed
20260912 for primary discovery, 20260913 for primary reserve, 20261012 for strict
discovery and 20261013 for strict reserve.

Two-means initializes one centre at the row farthest from the grand mean and
the other at the row farthest from that first centre, in square-root composition
space. Nearest-centre assignments and arithmetic centre updates repeat until
labels stabilize, with at most 200 iterations. Distance ties follow first-array
index order. If a cluster empties, a farthest row is reassigned; identical input
rows receive a deterministic half split, but zero displacement is excluded
from observed axes and errors in null axes. These implementation details are
not a biological assertion that each species has exactly two morphs. Exact
row order therefore belongs to reproducibility; it cannot be silently shuffled.

Null palettes are reassigned only within coarse-state pools within a cohort;
species-by-state counts and pooled palette rows remain fixed. The observed
selected species are retained in every null world. The continuous minor-cluster
fraction filter is **not reapplied** in null worlds, and no observer, season or
geographic strata constrain this permutation. Consequently this comparison
does not calibrate the entire data-dependent admission procedure or eliminate
photographic dependence. The audit identifies a scope limitation; it does not
establish the magnitude or direction of any resulting bias. No null was changed,
no additional result was generated and no historical p-value was replaced.
Residual-diagnostic tables and exact source/input bundle closure remain pending.

## S3. H3 explanatory boundaries

H3a source inventory:
`results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`.
H3b result note: `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`.
The manuscript figures must include reserve results and preserve all placement
scenarios. Missing trait coverage is not a biological negative.

The CI artifact summary tables were retrieved without rerunning either test:

- H3a: run 34677042793 / artifact 10292218669; retained as
  `docs/supporting/h123/h3a_signal_by_scenario.csv`. SHA256
  `8b37b0edcdf0823351aad82d443d2ff6686d3684294855f20900e1770614f986`
  matches the existing frozen manifest and retrieved artifact bytes.
- H3b: run 34677468362 / artifact 10292399238; retained as
  `docs/supporting/h123/h3b_span_summary.csv`. Retrieved SHA256
  `a8f58d2cf7282ec8aa55fd5741eb34e9440112fd319710654c4f03c12551460b`.

These copies use LF newlines; CI checks LF-canonicalized hashes to accommodate
Windows checkout conversion. This does not rewrite any original exact-artifact
contract. Manuscript H3 values are checked against these tables. This verifies
reporting, not independent reconstruction of tree models or permutation draws.

### Table S3. H3a: raw and opportunity-adjusted phylogenetic signal

| Cohort | Scenario | Tips | Raw K | Raw K p | Lambda | Lambda p | Residual K | Residual K p |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| discovery | S1 | 368 | 0.101549 | 0.0045 | 0.237590 | 0.00000003 | 0.054239 | 0.4832 |
| discovery | S2 | 368 | 0.012378 | 0.2209 | 0.235472 | 0.00000002 | 0.000843 | 0.9587 |
| discovery | S3 | 368 | 0.094884 | 0.0093 | 0.236230 | 0.00000003 | 0.053043 | 0.4912 |
| reserve | S1 | 341 | 0.071019 | 0.2716 | 0.048566 | 0.16569261 | 0.069291 | 0.3500 |
| reserve | S2 | 341 | 0.060148 | 0.4134 | 0.051129 | 0.15551420 | 0.065332 | 0.2464 |
| reserve | S3 | 341 | 0.070758 | 0.2674 | 0.047464 | 0.17023776 | 0.069156 | 0.3510 |

Raw K is primary; lambda likelihood-ratio p-values are secondary and are not
interchangeable with K randomization p-values. No reserve scenario passes raw
K p<0.05, and no opportunity-adjusted scenario passes. The all-scenario rule
does not choose the best placement after inspecting results. Tree coverage
is 368/369 in discovery and 341/363 in reserve; the missing 1 and 22 species
remain missing, not negative outcomes. The backbone/package revision and exact
tree identities are retained in the
[pre-outcome tree manifest](../results/polymorphism_h3a_phylogeny_preflight_20260912/frozen_tree_manifest.json).
Use of that phylogenetic infrastructure does not restore the retired 34-species
comparative analysis as evidence.

### Table S4. H3b: sampled-span association and opportunity sensitivities

| Cohort | Species | Analysis | Rank association | Two-sided p |
| --- | ---: | --- | ---: | ---: |
| discovery | 369 | Raw D | 0.179879 | 0.000900 |
| discovery | 369 | Corrected D | 0.179683 | 0.000800 |
| discovery | 369 | Partial ranks | 0.158508 | 0.002700 |
| reserve | 363 | Raw D | -0.002586 | 0.958602 |
| reserve | 363 | Corrected D | -0.002374 | 0.962902 |
| reserve | 363 | Partial ranks | 0.005519 | 0.916204 |

Raw/corrected rows are Spearman correlations. Partial ranks are Pearson
correlations of residualized centered ranks, not raw Spearman correlations.
All tests use 20,000 permutations, with a plus-one denominator of 20,001.
Their p-values do not test discovery–reserve differences directly. Rank-PGLS
on the 341 reserve tips remains a secondary diagnostic, reported in the
[frozen H3b result note](POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md);
it cannot override the failed primary test.

The rank-PGLS implementation first matches rows to the pruned tree, then
calculates centered average ranks within that matched set. Its formula is
`rank(D) ~ rank(log1p(span)) + rank(log1p(n_classifiable)) + rank(log1p(n_observers))`,
including an intercept, fitted with `phylolm(..., model="lambda")`. Thus these
rank values are not necessarily the ranks used in the full 363-species primary
test. The recorded span p-value is computed explicitly as
`2 * pnorm(-abs(beta_span / se_span))`: a normal-reference Wald calculation,
not the primary 20,000-permutation test or a package-summary t-test.
The coefficient is in centered-rank units, not change in D per kilometre.
The frozen note reports reserve lambda approximately 1e-7 in all scenarios;
this fitted value does not establish absence of phylogenetic dependence in
every model or sampling design. Runtime package versions and optimizer
diagnostics require the original artifact, not assumptions from current docs.

### Implementation and interpretation checks

The [H3a protocol](POLYMORPHISM_H3A_PHYLOGENETIC_SIGNAL_PROTOCOL_20260912.md),
[H3a runner](../scripts/analysis/run_polymorphism_h3a_phylogenetic_signal_20260912.R),
[H3b protocol](POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md) and
[H3b runner](../scripts/analysis/run_polymorphism_h3b_reserve_span_20260912.R)
were read to reconstruct Methods. In H3a, the implementation requests 10,000
K maps, checks that the first is the observed map, and uses the remaining 9,999
as randomizations. Both runners use a 1e-15 comparison tolerance in counting
upper-tail exceedances. Rank ties use average ranks. H3a fits opportunity
residuals before pruning to matched tree tips; H3b residualizes both D and span
against usable-image and observer counts. No alternative controls were fitted
for this manuscript. These code checks are not a rerun of the frozen analyses
or an independent validation of permutation exchangeability.

Non-support does not establish equivalence, absence of every phylogenetic effect,
or adequate power for every ecological effect size. No equivalence margin or
minimum detectable effect was established by these results. The source-derived
photographic span can reflect sampling and biology simultaneously; it cannot
identify environmental heterogeneity or pollinator mechanisms. Trait-coverage
failures described in the H3b protocol remain untested hypotheses, not negative
ecological results.

## S4. Measurement failures and unqualified gates

Retain `docs/RGFCA_MONARDA_REGION_AGREEMENT_RESULTS.md` and
`docs/RGFCA_ROI_QUALIFICATION_AUDIT.md` as measurement limitations.
P500 implementation qualification is documented in
`docs/P500_WHITE_CONTROL_IMPLEMENTATION_ERRATUM_20260915.md`.
The static gate does not authorize opening. No P500 biological results exist
in this package. CI success is not a measurement-accuracy result.

The [measurement-reference audit](FCP_H123_MEASUREMENT_REFERENCES.md) records the
primary sources supporting the main text's photography/calibration limitations.
Their validation does not transfer automatically to this estimator or cohort.

### Table S1. Completed Monarda operational measurement diagnostic

| Criterion | Observed | Fixed minimum | Decision |
| --- | ---: | ---: | --- |
| Pooled prediction precision | 0.56824250 | 0.70 | Failed |
| Pooled reference recall | 0.42608787 | 0.35 | Passed |
| Median annotated-image prediction precision | 0.51480059 | 0.70 | Failed |

All three conditions were required; the overall decision is **failed**.
The denominator is 109 annotated images, including 15 empty predictions, from
110 accounted-for images. One annotation-unknown image is not a verified
negative and does not enter overlap metrics. Pooled reference, prediction and
intersection counts are respectively 28,428,697, 21,316,820 and 12,113,123 pixels.
Pooled precision divides intersection by prediction; recall divides intersection
by reference. Per-image median precision retains empty predictions with zero
precision under the frozen convention. Larger regions contribute more weight
to pooled metrics than to the image median. These are descriptive operational
metrics, not independent pixel-level confidence estimates.

The original source files remain unchanged:

- [Saved result](../data/validation/monarda_region_agreement_v1/monarda_region_agreement_result_v1.json),
  exact SHA256 `49b2f017a6accb305867d62e7b4afe4e8ed96fcd03f8aafde86bda2da915fa4c`.
- [All-image count table](../data/validation/monarda_region_agreement_v1/monarda_region_agreement_rows_v1.csv),
  exact SHA256 `67cbeaa2cd0ea1f1d3f98516b4fe9f7bb1abc8648ac80268e6dc4ec3e754dc22`.
- [Original independent verification](../data/validation/monarda_region_agreement_v1/independent_verification.json).

Publication checks rederive Table S1 from retained integer pixel counts and
compare the saved result, exact byte identities, image census and table cells.
They do not verify annotation truth or rerun the estimator. Generic flower
polygons are not exhaustive focal-taxon petal labels. The observed disagreement
cannot separate localization error, incomplete annotation and target mismatch;
no such explanation is selected as a proven cause. Provider train/valid/test
labels are not FCP holdouts, and this opened export cannot serve as a new
untouched validation set. Neither these metrics nor the JRC box qualification
calibrates reflectance, ultraviolet signal or biological morph frequencies.

## S5. Package completion checklist

The main bibliography now includes the verified Lin, Phipson–Smyth, Blomberg
and V.PhyloMaker2 references. Their source access and bounded support are
recorded in the [statistical-reference audit](FCP_H123_STATISTICAL_REFERENCES.md).
This does not validate the particular FCP null or establish historical
non-access. Lambda, regression and remaining model citations are still pending.

Main Figure 1 is rendered by
`scripts/analysis/make_h123_evidence_figure.py` using only the six retained
summary files named in its plotted-data JSON. It does not open photographs or
fit models. PNG and vector PDF are provided, with LF-canonicalized source
hashes (not replacements for original exact-byte provenance). Tests check
source-to-plot data identity and repeated-render byte equality within one
runtime. Cross-platform byte identity is not asserted. H1's partition spread
and H2's null spread are explicitly distinguished from confidence intervals.
The figure is an overview, not the complete set of supplementary diagnostics.

- [ ] All source-to-number mappings checked automatically.
- [ ] Dedicated main/supplement figures generated and visually inspected.
- [ ] Methods and bibliographic references audited together.
- [ ] Measurement interpretation and P500 outcome/STOP finalized.
- [ ] Main text and supplement validated as one submission package.

The retired 34-species comparison and six-species spatial analysis are excluded
from this supplement, while their original files remain unchanged in history.
