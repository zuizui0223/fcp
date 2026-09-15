# H1 observer-disjoint four-state D — frozen result

Date frozen in repository: 2026-09-15 JST  
Analysis date: 2026-09-13 JST

Protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260913.md`  
Preflight freeze: `docs/POLYMORPHISM_H1_OBSERVER_SPLIT_PREFLIGHT_FREEZE_20260913.md`  
Script: `scripts/analysis/run_polymorphism_h1_observer_disjoint_d_20260913.py`

Workflow run: `34706708406`  
Successful job: `103587961180`  
Artifact: `10302770009` (`polymorphism-h1-observer-disjoint-d-20260913`)  
Artifact digest: `sha256:7a87432d91c257a1c010b9cf075eeb8a5f15e2ea016d55ed6d8cfaddc6f0b64e`  
Head SHA: `daabd3b63d540238a2cd27185330a59ca45e9d2b`

Exact machine-readable receipt: `results/polymorphism_h1_observer_disjoint_d_20260913/result.json`

## Frozen verdict

**`H1_OBSERVER_DISJOINT_D_REPRODUCIBLE`**

The current four-state polymorphism score

`D = 1 - sum_k p_k^2`

is reproducible across disjoint sets of observers in the species-disjoint reserve validation cohort under the frozen high-depth design.

## Outcome firewall

Observer assignment was frozen before morph-state opening. Preflight verification reports:

- `morph_opened_during_preflight_verification = false`
- `D_computed_during_preflight_verification = false`
- assignment SHA256: `f5e4333e596a69cdb032d840ec95eb02eadf1ec322658e2738e1b043fcde0cc1`
- panel SHA256: `bb30b6f73c6ecac06022c7ddd275cd01ea71d4bc40f54dd8889814819d1f92eb`
- selected minimum classifiable rows per observer-disjoint half: 20.

All 369 discovery species and all 363 reserve species passed the frozen split-depth gate. Thus the primary reserve result is not created by post-outcome species attrition within these high-depth validation cohorts.

## Primary reserve result

Reserve is the primary cohort; discovery cannot rescue it.

- n = **363 species**
- Spearman `rho(D_A, D_B)` = **0.8109164415**
- bootstrap 95% percentile interval for rho = **[0.7648383249, 0.8474976907]**
- two-sided 20,000-permutation p = **4.99975e-05**
- Lin CCC = **0.8582944314**
- mean D_A = 0.2297388470
- mean D_B = 0.2300604009
- median `|D_A-D_B|` = **0.0624349636**
- 90th percentile `|D_A-D_B|` = 0.1721134272

The frozen decision gates were:

- reserve rho >= 0.80;
- bootstrap lower bound > 0.70;
- permutation p < 0.001.

All three gates pass.

## Finite-sample sensitivity

For `D_unbiased = D*n/(n-1)` in reserve:

- n = 363
- rho = **0.8116141836**
- bootstrap 95% interval = **[0.7664741451, 0.8478134694]**
- permutation p = **4.99975e-05**
- Lin CCC = **0.8594198789**

The result is effectively unchanged by the finite-sample correction.

## Discovery calibration only

Discovery n = 369:

- rho = **0.8370151412**
- bootstrap 95% interval = **[0.7987951355, 0.8674412513]**
- permutation p = **4.99975e-05**
- Lin CCC = **0.8661920927**

These values are supportive calibration only; the claim is accepted because reserve independently passes the frozen rule.

## Historical 0.971 value

The historical project-summary observer reproducibility value near 0.971 is **not used as evidence** for H1. The canonical evidence is the direct current four-state-D analysis above.

## Allowed claim

> Under the frozen high-depth photo design, continuous four-state flower-colour polymorphism D is reproducible at the species level across disjoint observer sets, including an independent species-disjoint reserve cohort.

## Hard nonclaims

This result does not establish:

- unbiased global prevalence across all 42,111 species;
- biological independence of photographs;
- absence of image-classification error;
- genetic discreteness of the four coarse states;
- true population morph frequencies;
- ecological or evolutionary causes of polymorphism;
- direction of colour evolution;
- representativeness of the 369/363 high-depth validation cohorts for the 42,111-species universe.

## Manuscript consequence

H1 is closed as a positive **measurement-validity** result. It should precede H2 in the paper but should not be sold as a biological predictor result. H3a and H3b remain closed negative results and must not be used to rescue or reinterpret H1.
