# FCP commonality among species — 4,475 genuinely measured taxa

**Analysis date:** 18 October 2026. **Research question:** Do photographed within-species flower-colour polymorphisms follow general patterns in *morph frequency skew* and *geographic allocation* that recur irrespective of the identity of the dominant hue? This is about FCP architecture, **not** whether white/yellow/red/blue flowers are more numerous globally.

**Machine provenance:** The exact 447,500-original-photo image measurements and original inference were successfully completed in [real ROI #38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400), [actual old-method 50-km/H1/W #39049598596](https://github.com/zuizui0223/fcp/actions/runs/39049598596), and permanently authenticated in [the original results release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017). The separate [genuine full-five-cohort commonality rerun #39091197295](https://github.com/zuizui0223/fcp/actions/runs/39091197295) passed unit tests and original numeric gates; its archived derived artifact contains both a source-verified JSON and species-level tables. Actual source-photo CSV SHA256 **8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6**, spatial species CSV SHA256 **960e8c4bae8a9ecb01fefbef06b6b88cd8fdc63d86af3e4c354cab4eb1aedd35**, original spatial JSON SHA256 **79aeab0555ec2efd865606a9bdc3e06574a5168c38029ee16e22b945e1a9eacf**. Derived JSON SHA256 **3a51d537fd99766e64846ad318c53d3e70f1d1fe0a1222ecfcae2c6854e551e2**. Main reproducer: \`scripts/analysis/fcp_five_cohort_polymorphism_commonality_20261018.py\`.

## Source, measurement and conditioning

The **4,730** species are an **outcome-blind photo-opportunity frame**, not 4,730 successfully measured taxa. Among them, **4,475 taxa × 100 genuine photos = 447,500 measured photo terminals**, 227,107 four-state classifiable photos, **3,297** taxa with >=40 classifiable colour states for species-wide D, and **1,761** of these with >=30 within-50-km local photo pairs. Of the other 255 taxa, 252 lacked a genuine full 100-photo source and three new taxa reused source photo/observation IDs and were conservatively excluded whole. Historical original three cohorts and their frozen prospective H2 test were not overwritten.

The four categories are *white, yellow/orange, red/pink, blue/purple*. They are photographed, algorithmically defined coarse phenotype states, **not a genotype-confirmed population-level count of inherited flower-colour morphs**. Classifiability (~50.8% of measurement terminals) is nonrandom and digital exposure/observer/location biases remain possible. Accordingly, all percentages below refer to the **retained photos**, never the fraction of natural plant populations maintaining an FCP.

## Pattern 1 — photographed FCP is frequently strongly skewed

Among **3,297** taxa with >=40 classifiable genuine original photos (1,229 inferred genus labels):

| Within-species photographed state descriptor | Observed |
|---|---:|
| Species Gini–Simpson D, median (IQR) | **0.17236** (0.06342–0.39294) |
| Most common state share, median | **0.90566** (90.6%) |
| Species with photographed majority state >=90% | **51.4%** |
| Species with 1 / 2 / 3 / 4 observed coarse states | **268 / 1,411 / 1,080 / 538** |

The species-equal tendency is toward a **dominant photographed type with relatively infrequent alternative states**. This does not demonstrate selection for one morph, a pigment-synthesis cost or stable wild morph frequency. The same descriptive quantity at >=60 classifiable photos is provided by the original code/derived source JSON for measurement-depth sensitivity.

## Pattern 2 — local colour sorting is not restricted to one named hue

For 1,761 genuinely eligible taxa the frozen species-global-minus-within-50-km mismatching pair fraction is **+0.0176149**, under the unchanged **199 complete-vertex permutations per species** (upper p=0.005). The new fifth cohort alone supplied 1,020 evaluable species (mean +0.0175987, p=0.005), closely matching the earlier four-cohort mean +0.0176372.

Instead of asking whether white or red flowers occupy particular climates, we **rename each species' most frequently photographed coarse morph the majority state** and merge the other 1–3 states as the alternative class. The **same original vertex permutations and species coordinates** then test local clustering of these identity-free majority/other states. The **species-equal 50-km depletion remains +0.0164851, p=0.005**. This is not an independent new confirmation; it is a *post hoc, hue-agnostic representation* of the same photographed geographic signal.

Original four-state depletion also remains positive against matched nulls within all **four dominant-state groups**:

| Species' photographically most common state | 50-km eligible species | Mean depletion | Original 199-null upper p |
|---|---:|---:|---:|
| White | 807 | +0.018552 | 0.005 |
| Yellow/orange | 521 | +0.010326 | 0.005 |
| Red/pink | 203 | +0.025927 | 0.005 |
| Blue/purple | 230 | +0.023502 | 0.005 |

Dominance class is itself inferred from the same photographic states and the subgroup tests are **post hoc**. These four upper p-values are not a comprehensive project-wide multiplicity-adjusted test. Similarly, across four sampled species-global D rank quartiles, the mean depletion ranges **+0.00290 (lowest D quartile)** to **+0.03675 (highest D quartile)**, with upper p=0.005 in each matched species-composition permutation comparison. Effect magnitude depends strongly on the amount of available colour diversity, so these strata do **not** imply that any one value of D causes local clustering.

Equal weighting across **832 sampled genera** (genus inferred from the first taxonomic name token) still yields **+0.0152512** in the original four-state analysis, and **+0.0146322** for identity-free majority/other (both upper p=0.005). This mitigates numerical domination by oversampled genera but does not correct complete phylogenetic nonindependence.

## Pattern 3 — photograph-level falsification and a critical exception

Reapplying the same composition-preserving null after discarding specific local photo pairs:

| Restriction at 50 km | Eligible species | Mean depletion | Original-style matched upper p |
|---|---:|---:|---:|
| Every genuine classifiable photo pair | 1,761 | +0.017615 | 0.005 |
| Local photos from different observers only | 1,745 | +0.017070 | 0.005 |
| Local photos observed in different calendar years only | 1,643 | +0.015845 | 0.005 |
| At least 60 four-state classifiable photos per species | 845 | +0.017066 | 0.005 |
| **Remove every white observation**, classify among three nonwhite states | **763** | **+0.006595** | **0.005** |

The nonwhite analysis is **not positive in all individual cohorts**: original prospective n=78 had upper p=0.14; fourth wave n=79 showed a negative mean −0.002277 and upper p=0.72. Consequently a pooled nonwhite effect cannot be represented as a universal, independently replicated chromatic-only sorting principle. The remaining white/background exposure sensitivity and geography-linked observer ascertainment require biological field checking.

## FCP commonality thesis and limits

The evidence favours a **two-layer, source-conditional FCP pattern**:

1. **Composition:** many species have uneven photographed morph-state distributions, with one state often dominant.
2. **Geographic allocation:** the observed alternatives are, on average, **more locally sorted than expected from each species' fixed overall state proportions**, including when state identity is reduced to *majority versus alternatives*, across genera and after several observer/time controls.

The historically prospective **W** test adds a separate, construction-controlled *direction-of-colour-space variation* result; its pooled fifth-cohort n=1,223 and W=0.51879 against null median 0.45537 (p=0.001) do **not** show that white increases evolutionarily or that pigment is selectively costly. Cross-species photomorph distribution is **not equivalent** to an evolved balancing-selection equilibrium. All analyses after the original freeze, including the dominance/stratum comparisons in this document, are exploratory. The 50-km scale was itself post hoc.

**Next biological discrimination, not demonstrated by the source:** systematically field-counted within-population morph frequencies, replicated local sites, genetic structure, reproductive success and pollinator data are needed to distinguish persistent balancing selection, spatially heterogeneous selection, history/drift, plasticity and camera/observer artifacts.

## Paper architecture

Keep the *New Phytologist* chronology intact: Stage 0 outcome-blind 4,730 photo-opportunity frame; Stage 1 observer-disjoint D validity; Stage 2 originally untouched species/photo-disjoint frozen W confirmation; Stage 3 bounded measurement, climatic and phylogenetic mechanism filters; Stage 4–5 genuine additional same-provider image cohorts tested with unchanged methods. **The new commonality layer is an exploratory *synthesis across the five measured cohorts*, not a retroactively frozen hypothesis.** It elevates the comparison of photographed *within-species FCP structure* above any simple global comparison of particular hues.
