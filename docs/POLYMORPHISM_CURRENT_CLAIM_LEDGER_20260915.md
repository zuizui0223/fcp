# Flower-colour polymorphism current claim ledger

Date: 2026-09-15 JST

This ledger consolidates the currently admissible claims on branch `analysis/polymorphism-42111-h1-h2-gates-20260912`. It does not run a new statistical test and does not replace the machine-readable receipts cited below.

## 1. Paper mainline

The current paper should be organized around two positive claims and two bounded negative results:

1. **H1 — measurement validity:** continuous four-state species-level flower-colour polymorphism `D = 1 - sum_k p_k^2` is reproducible across observer-disjoint photo sets in the high-depth validation design, with later split analyses showing modest design sensitivity rather than split-invariant reliability.
2. **H2 — geometry:** within-species colour polymorphism in the existing high-depth cohorts is disproportionately aligned with an achromatic–chromatic, white-versus-nonwhite axis. No recurrent non-white hue direction survives the construction-preserving audit.
3. **H3a — phylogeny:** broad tree-wide phylogenetic signal in D does not replicate in the species-disjoint reserve cohort.
4. **H3b — sampled geographic span:** the discovery association between sampled span and D does not replicate in reserve.

The old global shared-boundary/sharedness line is not part of this paper's positive biological evidence; its qualification failures remain an identifiability/detection-limit result.

## 2. H1 — positive measurement-validity claim

Canonical result freeze:

- `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RESULT_FREEZE_20260913.md`
- `results/polymorphism_h1_observer_disjoint_d_20260913/result.json`
- preflight: `docs/POLYMORPHISM_H1_OBSERVER_SPLIT_PREFLIGHT_FREEZE_20260913.md`
- chronology audit: `docs/POLYMORPHISM_H1_CHRONOLOGY_AUDIT_20260915.md`
- corrected reconciliation: `results/polymorphism_h1_evidence_reconciliation_20260914/result.json`

### Prospective chronology

The deterministic H1 test is the first frozen H1 outcome test, not a post-hoc rescue:

- outcome-blind preflight workflow `34706579235` completed successfully at 2026-09-12 16:53:41 UTC;
- protocol commit `743e13d3597c88ea7a725c596ea733b30d978024` froze the test at 16:53:54 UTC;
- preflight receipt commit `a77e23bd7a098307c4961d77271fb2e5a294ce88` followed at 16:54:57 UTC;
- primary outcome workflow `34706708406` ran successfully from 16:55:59 UTC and produced immutable artifact `10302770009`, digest `sha256:7a87432d91c257a1c010b9cf075eeb8a5f15e2ea016d55ed6d8cfaddc6f0b64e`.

The later 200-partition protocol was not frozen until 17:03:49 UTC. The previous 2026-09-14 reconciliation receipt incorrectly omitted this earlier deterministic run; that chronology has now been corrected without changing any statistic or threshold.

### Outcome firewall

- observer assignment and the per-half opportunity gate were fixed from `species`, `observer_id`, and `global_classifiable` only;
- `morph`, palette coordinates, D, and D_unbiased were not opened during preflight;
- all 369 discovery species and all 363 reserve species passed the frozen minimum of 20 classifiable photos per observer-disjoint half.

### Primary reserve result

- n = 363 species;
- Spearman `rho(D_A, D_B) = 0.8109164415`;
- bootstrap 95% percentile interval = `[0.7648383249, 0.8474976907]`;
- 20,000-permutation two-sided p = `4.99975e-05`;
- Lin CCC = `0.8582944314`;
- median `|D_A-D_B| = 0.0624349636`;
- finite-sample-corrected `D_unbiased` gives rho = `0.8116141836`.

Frozen decision gates were reserve rho >= 0.80, bootstrap lower bound > 0.70, and permutation p < 0.001. All pass.

Discovery is calibration/support only: n = 369, rho = 0.8370151412, CCC = 0.8661920927.

### Later robustness analyses

These analyses were frozen after the canonical primary outcome was already known and therefore cannot retroactively replace it, but their results remain part of the evidence record.

- repeated 200 observer-disjoint partitions: reserve median rho = **0.789103**, q05 = **0.765165**, q95 = **0.810940**, median CCC = **0.854831**;
- later strict deterministic stress test: reserve rho = **0.792728**, bootstrap 95% interval = **[0.741868, 0.832366]**, CCC = **0.847429**; its deliberately strict rho >= 0.80 criterion fails.

Thus H1 supports reproducible D measurement, but the exact correlation is somewhat split-design-sensitive around rho ~0.79–0.81. Near-perfect or split-invariant reliability is not supported.

**Allowed H1 claim:** under the frozen high-depth photo design, continuous four-state flower-colour polymorphism D is reproducible at the species level across disjoint observer sets in the species-disjoint reserve cohort (first-frozen reserve rho = 0.811, 95% bootstrap interval 0.765–0.847). Later alternative observer partitions yielded rho around 0.79, so near-perfect or split-invariant reliability is not claimed.

**Important boundary:** the historical project-summary observer reproducibility value near 0.971 is not current-D evidence and must not be used. H1 validates measurement/sampling stability in the high-depth cohorts; it does not estimate global prevalence, true population morph frequencies, genetic discreteness, or image-classification accuracy.

## 3. H2 — positive but targeted geometry claim

Canonical evidence ledger:

- `docs/POLYMORPHISM_H2_EVIDENCE_LEDGER_20260912.md`
- targeted result: `results/polymorphism_white_axis_targeted_test_20260912/result.json`
- target freeze: `docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`

H2 modes were reconstructed label-free from the normalized nine-colour palette by deterministic two-means after H1 admission. The broad leading-axis geometry first exceeded an isotropic null and then exceeded a coarse-state-preserving structured null.

The later audit showed that this geometry is dominated by a fixed white-versus-equal-nonwhite contrast

`q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`.

For `W = mean_i (u_i dot q_white)^2` under the same coarse-state-preserving structured null:

Primary 10% gate:

- discovery: N = 152, W = 0.514625, null median = 0.430808, p = 0.001;
- reserve: N = 129, W = 0.514586, null median = 0.466546, p = 0.001.

Strict 20% gate:

- discovery: N = 75, W = 0.542355, null median = 0.443626, p = 0.001;
- reserve: N = 65, W = 0.510517, null median = 0.469943, p = 0.008.

When the white contrast is projected out, residual directional concentration does not exceed the structured null. The non-white-only top-two subset also fails the discovery construction-controlled concentration test.

**Allowed H2 claim:** across the existing high-depth validation cohorts, within-species flower-colour polymorphism is disproportionately aligned with an achromatic–chromatic (white versus non-white) axis beyond expectation from the frozen coarse-state composition and the global state-to-palette mapping.

**Status boundary:** the specific white/non-white target was isolated after the initial H2 geometry had been opened. Existing-cohort evidence is therefore a targeted post-audit decomposition, not an untouched confirmatory test. `q_white`, W, the H1/H2 admission rules, and the structured null are now frozen as a prospective target for future high-depth expansion.

Do not claim a general recurrent non-white hue axis, a particular pigment pathway, direction of colour evolution, or an adaptive/pollinator/climate mechanism from H2.

## 4. H3a — closed negative result

Canonical frozen manifest:

- `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`

Reserve is the fresh primary cohort; each S1–S3 tree retains 341 tips.

- S1: K = 0.0710190, p = 0.2716; lambda = 0.04857, p(lambda=0) = 0.1657; opportunity-adjusted p = 0.3500.
- S2: K = 0.0601476, p = 0.4134; lambda = 0.05113, p(lambda=0) = 0.1555; opportunity-adjusted p = 0.2464.
- S3: K = 0.0707577, p = 0.2674; lambda = 0.04746, p(lambda=0) = 0.1702; opportunity-adjusted p = 0.3510.

Verdict: `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

This does not imply absence of all taxonomic structure or biological irrelevance of phylogeny. Discovery-only genus clustering and placement-sensitive discovery K cannot rescue the reserve failure.

## 5. H3b — closed negative result

Canonical result freeze:

- `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`

Discovery calibration:

- n = 369;
- rho(D, log1p sampled span) = 0.1798786;
- permutation p = 0.00089996.

Fresh reserve replication:

- n = 363;
- rho = -0.0025855;
- permutation p = 0.9586021;
- observer/classifiable-adjusted partial-rank rho = 0.0055187, p = 0.9162042;
- S1–S3 rank-PGLS all give beta about -0.00559, p about 0.9137.

Verdict: `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

Do not rescue this result with alternate span definitions, threshold changes, or new predictor hunting. The tested predictor is sampled span in the fixed photo design, not true biological range size.

## 6. Sampling frame and inference boundary

The scientific sampling frame is the 42,111-species global flower-colour universe, not the historical 369-species discovery cohort.

The high-depth cohorts are validation sets:

- discovery: 369 species with `n_classifiable >= 40`;
- reserve: 363 species, species-disjoint from discovery.

For the 42,111-species frame, observer-capped raw opportunity includes 4,730 species with at least 100 photos (`U100`). This is an opportunity ceiling, not the observed H1/H2 sample size and not a polymorphism prevalence estimate.

Failure to pass H1/H2 measurement gates in a future expansion must be coded as missing/underidentified unless the biological state is actually resolved; it must not be silently recoded as monomorphism.

## 7. Manuscript claim ceiling

The strongest defensible paper-level statement at the current evidence state is:

> Species-level flower-colour polymorphism can be measured reproducibly from high-depth citizen-science photographs, with observer-split correlations around 0.79–0.81 rather than near-perfect split invariance. In two species-disjoint high-depth cohorts, the dominant recurrent geometry of within-species colour variation lies along a white-versus-nonwhite axis rather than a general hue axis. Broad phylogenetic signal and a discovery association with sampled geographic span do not replicate.

This paper is therefore primarily about **measurement + geometry**, not about a discovered ecological or phylogenetic predictor of polymorphism.

## 8. Next prospective gate

The next result that can materially raise the claim ceiling is a response-blind high-depth expansion from the 42,111-species frame using the already frozen H1/H2 rules and fixed `q_white` target. The purpose is not to retune the axis; it is to ask whether the white/nonwhite geometry transports beyond the 369/363 validation cohorts.
