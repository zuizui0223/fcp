# FCP: disjoint 50-km photographic pairing — prereplay audit (2026-10-10)

## Research objective

The original New Phytologist flower-colour ITV geographic allocation claim uses all photographed within-species pairs <=50 km. Multiple pairs share photographic endpoints; therefore increasing an arbitrary minimum number of pair *edges* does not by itself ensure that each observation contributes once. The prior source-frozen photo-pair opportunity stress test (PR #151) remains positive for 30, 50 and 100 local edges per species in discovery, validation and third source cohorts. This post-outcome falsification examines a different aspect of observational pseudoreplication: **does the local colour concentration remain when each original photograph appears in at most one local pair?**

## Source and frozen estimator

- Exactly original discovery and validation source commit `5142f7951af0dde5364bb047a566d67e8c479e51`, third source commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`. Byte-level SHA256 verification before source analysis; no new images/labels/coordinates or new source taxa.
- Species with >=40 classified four-state original flower-colour photographs; geographic within-species pair cutoff **50 km**.
- Reuse the original implementation's great-circle distance, four-state label coding, species-wide photographed pair discordance `D_pair`, and 199 complete-vertex colour-label permutations that conserve exactly each species' four-state frequency.
- **Photo-disjoint greedy maximal matching** of local graph edges, based ONLY on observed geographical distance and fixed photograph IDs. Primary matching considers nearest permissible source edges first (deterministic geodesic, row/photo-ID ties); a secondary outcome-blind fixed-seed random edge order asks how sensitive the result is to geographical matching choice.
- Primary policy permits all original photo pairs; secondary policy requires different nonmissing photographer identities across each selected pair, still without using the colour outcome during pair construction.
- Independent paired *source photograph identity* count is 2 times disjoint pairs. This is NOT independent genetic individuals or field localities.
- In each source cohort and matching scheme/policy, inspect minimum **5, 10, 15, 20 and 30 disjoint photographed pairs/species**. Main stress threshold is **10** disjoint pairs. Minimum **30 species per cohort** to state positive exploratory support; otherwise report explicit `HOLD_INSUFFICIENT_SPECIES`.
- Each species contributes equal weight to mean `D_pair - D_disjoint-local` rather than being weighted by pairs. One-sided 199 matched composition-preserving permutation p and 1,999 species bootstrap draws. For cautiously supportive exploratory wording require (i) n>=30, (ii) positive mean, (iii) p<0.05, (iv) positive bootstrap lower 95% endpoint.
- Program must recover the original 50km all-pair baseline exactly: n=166/181/204 and +0.020529254583812922/+0.018672971642749295/+0.01468491968437185, with high-depth source species 369/363/377. Fail closed otherwise.

## Inference firewall

- This is outcome-exposed post hoc sensitivity, not an untouched cohort or confirmation. Threshold and method effects are exploratory and p values are unadjusted across multiple chosen thresholds, policies and matching algorithms.
- Matching is **greedy maximal**, not guaranteed to be maximum-cardinality; different-order sensitivity explicitly tests some of that selection dependence.
- Disjoint photograph IDs remove within-matching photo endpoint reuse only. Observers may cover multiple observations, repeated photos may concern the same botanical plant, and local geographically proximate photos can still be clustered spatially.
- Real photographic white/pigmented morph error and true botanical organ classification are unmeasured without independently blinded human annotation. This diagnostic is not a replacement for 405 expert review, biological genotypes or fitness.
- Keep all original FCP H1/H2 frozen results and main manuscript unchanged; untouched prospective 2,000+730 future taxa stay unopened.
- All negative or underpowered cases are reported with equal visibility. Do not collapse 'not estimable' into biological absence.
