# H1 observer-disjoint D reliability — frozen result

Date frozen: 2026-09-13 JST

Protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md`

Canonical workflow run: `34707537360` (run #2; implementation-only preaggregation speedup, frozen estimand unchanged)

Artifact: `10302466831` (`polymorphism-h1-observer-disjoint-reliability-20260913`)

Artifact digest: `sha256:19185ad9160b24ba47ef8d99a33d952ab5f5b6b62243ff756d41917bfc587c53`

Canonical machine-readable result: `results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`

## Frozen verdict

**`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED+DISCOVERY_CONSISTENT`**

The reserve cohort is the prespecified primary decision cohort. All three frozen support criteria pass under the >=20-classifiable-per-half gate.

## Reserve primary result

Across 200 deterministic outcome-blind observer partitions:

- all 363 reserve species retained >=40 observer-known classifiable four-state observations before half-size filtering;
- median paired species per split = **329** (range 315–339), versus frozen minimum 100;
- median split-half Spearman rho = **0.789103**, versus frozen minimum 2/3;
- 5th percentile rho = **0.765165**, versus frozen minimum 0.50;
- 95th percentile rho = 0.810940;
- median Spearman-Brown projected full-estimate reliability = **0.882121**;
- 5th percentile projected reliability = 0.866961;
- median CCC = **0.854831**;
- median absolute split difference in D = **0.078981**;
- median signed bias A-B = **-0.000980**;
- 95th percentile absolute signed bias = 0.010581.

The result therefore clears the pre-outcome reliability rule by a substantial margin. The near-zero signed bias indicates no persistent A/B direction under the symmetric observer partition, while CCC and MAE are agreement diagnostics rather than decision criteria.

## Discovery transport/calibration diagnostic

The independent discovery species set is consistent with the reserve result:

- all 369 species retained observer-known eligibility;
- median paired species per split = **340**;
- median rho = **0.811545**;
- 5th percentile rho = **0.788380**;
- median Spearman-Brown projected reliability = **0.895970**;
- median CCC = **0.851426**.

Discovery consistency is diagnostic and was not allowed to rescue reserve failure.

## >=15-per-half sensitivity

The broader sensitivity also supports the same conclusion but is not used to rescue the primary gate:

- reserve median paired N = 363;
- reserve median rho = 0.806052;
- reserve rho 5th percentile = 0.780115;
- reserve projected full-estimate reliability = 0.892612.

## Interpretation allowed in the main paper

Within the existing high-depth design, the continuous four-state species-level polymorphism score D is reproducibly recoverable from observation sets contributed by completely disjoint observers. This supports using D as a measurable species attribute in the high-information cohorts rather than treating it as an unstable artifact of particular observers.

## Hard non-claims

This result does not establish:

- global prevalence of flower-colour polymorphism;
- representativeness of the 369/363 high-depth species;
- perfect correctness of image-level colour labels;
- absence of geographic morph-frequency structure;
- ecological or phylogenetic causes of D;
- the H2 white/non-white geometry result, which has a separate structured-null evidence chain.

## Manuscript consequence

H1 can now serve as the measurement-validity foundation of the polymorphism paper. H2 remains the positive biological/descriptive geometry result. H3a phylogenetic signal and H3b sampled-span prediction remain frozen negative results and must not be used to tune H1/H2.