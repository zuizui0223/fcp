# H1 observer-disjoint D reproducibility — frozen result

Date frozen: 2026-09-13 JST

Protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260913.md`
Preflight freeze: `docs/POLYMORPHISM_H1_OBSERVER_SPLIT_PREFLIGHT_FREEZE_20260913.md`

Workflow run: `34706708406`
Artifact: `10302770009` (`polymorphism-h1-observer-disjoint-d-20260913`)
Artifact digest: `sha256:7a87432d91c257a1c010b9cf075eeb8a5f15e2ea016d55ed6d8cfaddc6f0b64e`

## Status

**Frozen verdict: `H1_OBSERVER_DISJOINT_D_REPRODUCIBLE`.**

The historical reliability value near 0.971 was not used as evidence because its provenance could not be tied unambiguously to the current four-state D definition and split procedure.

## Outcome-blind gate

The observer split was fixed before morph labels were opened. The preregistered automatic opportunity rule selected a minimum of **20 classifiable photographs in each observer-disjoint half**.

- discovery: 369 / 369 eligible species retained
- reserve: 363 / 363 eligible species retained
- every observer is confined to one half within each species
- frozen preflight assignment SHA256: `f5e4333e596a69cdb032d840ec95eb02eadf1ec322658e2738e1b043fcde0cc1`
- frozen preflight species-panel SHA256: `bb30b6f73c6ecac06022c7ddd275cd01ea71d4bc40f54dd8889814819d1f92eb`

## Fresh reserve — primary decision

n = **363 species**.

Four-state D:

- Spearman rho(D_A, D_B) = **0.8109164**
- species-bootstrap 95% CI = **[0.7648383, 0.8474977]**
- two-sided 20,000-permutation p = **4.99975e-05**
- Lin CCC = **0.8582944**
- median |D_A - D_B| = **0.0624350**
- 90th percentile |D_A - D_B| = **0.1721134**
- mean D_A = 0.2297388
- mean D_B = 0.2300604

All frozen decision conditions passed:

1. rho >= 0.80: pass
2. bootstrap lower bound > 0.70: pass
3. permutation p < 0.001: pass

Finite-sample-corrected sensitivity:

- rho(D_unbiased,A, D_unbiased,B) = **0.8116142**
- bootstrap 95% CI = **[0.7664741, 0.8478135]**
- permutation p = **4.99975e-05**
- Lin CCC = **0.8594199**

## Discovery calibration

n = **369 species**.

- Spearman rho(D_A, D_B) = **0.8370151**
- bootstrap 95% CI = **[0.7987951, 0.8674413]**
- permutation p = **4.99975e-05**
- Lin CCC = **0.8661921**
- median |D_A - D_B| = **0.0605469**

Discovery agrees with the reserve result but is not needed to rescue or define the primary decision.

## Allowed interpretation

Under the fixed high-depth validation design, the continuous four-state flower-colour polymorphism score D is a reproducible **ranked species-level measurement** across completely observer-disjoint photo subsets.

The near-equality of reserve half means and CCC around 0.86 also argues against the result being only a rank-order artifact, but the frozen H1 decision is based on Spearman rho, its bootstrap lower bound, and the permutation test.

## Hard nonclaims

This result does not establish:

- unbiased global prevalence of flower-colour polymorphism;
- biological independence of photographs;
- absence of image-classification error;
- ecological or evolutionary causes of D;
- direction of colour evolution;
- representativeness of the 369/363 validation cohorts for all 42,111 species.
