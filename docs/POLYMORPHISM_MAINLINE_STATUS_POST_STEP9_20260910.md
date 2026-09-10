# FCP polymorphism mainline after Step 9

Date: 2026-09-10 JST
Branch: `analysis/polymorphism-ambiguity-bounds-step9`
Evidence endpoint: `f49d5488c94dc0a3f36e7c34b3f1bfce40cf286d`

## Mainline decision

Stop organizing this paper around a globally shared flower-colour boundary, a privileged recurrent colour direction, or a strong genus-level propensity claim.

The flagship result is now:

> **species with greater discrete flower-colour polymorphism show stronger within-species geographic organization of flower colour, even though the geography/direction of that organization is not shared among species.**

This relationship replicates in a species-disjoint reserve photo tranche and survives observer/time controls, matched-background subtraction, sampled-span adjustment, clear technical-failure adjustment, and both uniform ambiguity-completion endpoints.

## Claim 1 — flower-colour polymorphism is a measurable species-level property

### Retain

Discovery frame:

- 369 eligible species;
- four-state Gini-Simpson diversity `D = 1 - sum(p_k^2)`;
- D approximately 0 to 0.708;
- second morph >=10% in 46.61% of the admitted frame;
- second morph >=20% in 26.56%.

Reserve frame:

- 363 eligible species from a species-disjoint fixed-photo tranche;
- same four-state construction and eligibility logic.

### Boundary

These are properties of the admitted photo-derived frames, not prevalence estimates for flowering plants globally.

D is intentionally conservative: ambiguous-palette rows are excluded from the baseline four-state estimator and are not recoded as hidden morphs.

## Claim 2 — greater polymorphism is associated with stronger within-species geographic organization

### FLAGSHIP — retain as the main biological result

### Discovery

Existing 369-species within-species spatial exact randomization:

- whole-frame mean spatial rho > null, p=0.001;
- second morph >=10% subset mean rho = 0.034565, p=0.001;
- >=10% minus <10% contrast = 0.014130, p=0.014;
- continuous rho(D, species spatial rho) = 0.089213, p=0.034.

The subset/continuous re-aggregation in discovery was formulated after the whole-frame result and remains explicitly exploratory.

### Species-disjoint reserve replication

Fixed before reserve polymorphism membership was joined to reserve spatial outcomes:

- >=10% subset primary mean rho = 0.032299, p=0.001;
- observer-pair exclusion p=0.001;
- quarter-stratified p=0.001;
- matched-background differential p=0.005;
- >=10% minus <10% contrast p=0.041;
- continuous rho(D, primary spatial rho) = 0.101601, p=0.025.

### Sampled geographic-span robustness

Step 6/6b showed that sampled geographic opportunity is not a simple sufficient explanation.

- discovery partial rho(D, spatial | sampled span, n_classifiable) = 0.170198, geometry-preserving p=0.001;
- reserve primary partial rho = 0.078186, p=0.080;
- reserve matched flower-minus-background partial rho = 0.087180, p=0.037.

The reserve primary continuous gradient weakens after the broad `n_classifiable` control, but the flower-specific differential remains supported.

### Pure technical-failure robustness

Step 8 separated clear ROI/flip operational failures from ambiguous palette rows.

After controlling sampled span + pure technical-failure rate:

- discovery primary partial rho = 0.126637, geometry-preserving p=0.007;
- reserve primary partial rho = 0.099288, p=0.025;
- reserve matched flower-minus-background partial rho = 0.116241, p=0.010.

Thus the result is not explained by sampled span plus clear ROI/flip failures alone.

### Ambiguous-palette endpoint robustness

Step 9 considered the sensitivity model in which every ambiguous-palette row is an unobserved member of one of the same four morph classes.

For each species, `D_min4` is the exact minimum and `D_max4` the exact maximum possible four-state D under that assumption.

Reserve primary spatial relationship after sampled span + technical-failure control:

- D_min4 partial rho = 0.097078, p=0.029;
- D_max4 partial rho = 0.125286, p=0.008.

Reserve matched flower-minus-background:

- D_min4 partial rho = 0.116299, p=0.009;
- D_max4 partial rho = 0.132785, p=0.006.

Discovery also passes both endpoints.

### Allowed interpretation

Within these fixed photo-derived global frames, species with greater flower-colour polymorphism have more geographically organized flower-colour variation. The association recurs in an independent species-disjoint photo tranche and the flower-specific component survives the strongest implemented measurement/opportunity sensitivities.

Different species can be organized in different geographic locations and colour directions. No common global boundary is required.

### Not allowed

Do not translate this result into:

- adaptation;
- local selection;
- population-genetic differentiation;
- causal effects of range size;
- shared climatic thresholds;
- a universal flower-colour boundary.

## Claim 3 — genus-level taxonomic clustering

### Retain as secondary / qualified

Step 4 discovery:

- 61 repeated genera / 169 species;
- clustering gain = 0.202963;
- p=0.002300; Holm p=0.011499.

Species-disjoint reserve Step 7:

- 54 repeated genera / 146 species;
- raw gain = 0.153654, p=0.021699.

Cross-tranche comparison of the 23 genera represented by >=2 different species in both tranches:

- genus-mean D concordance rho = 0.530632, p=0.009900;
- D_unbiased concordance rho = 0.521739, p=0.010950.

### Opportunity / measurement decomposition

- sampled-span-only adjustment does not remove the reserve genus pattern: gain=0.141642, p=0.022999;
- finite-sample-corrected D_unbiased also replicates raw: gain=0.153907, p=0.020749;
- sampled span + pure technical ROI/flip failure still passes narrowly: gain=0.121211, p=0.044848;
- broad `n_classifiable` adjustment removes the reserve pattern, but n_classifiable is a post-acquisition classification-success variable under a fixed 100-photo denominator and is strongly coupled to ambiguous-palette missingness.

### Ambiguity endpoint qualification

Step 9:

- reserve D_min4: gain=0.110212, p=0.061647;
- reserve D_max4: gain=0.141546, p=0.023799.

Therefore genus clustering is **not uniform-endpoint robust** under the four-state ambiguity sensitivity model.

However, cross-tranche genus-mean concordance remains supported at both endpoints:

- D_min4 rho=0.559289, p=0.005900;
- D_max4 rho=0.571146, p=0.005200.

### Paper role

Use genus clustering as evidence that species-level polymorphism is not completely taxonomically exchangeable, but do not make it a headline claim of intrinsic genus propensity.

Call it `genus-level taxonomic clustering`, never formal phylogenetic signal.

## Species-attribute results — what remains and what does not

### Sampled geographic span

Discovery Step 4:

- rho(D, sampled span)=0.181031;
- raw p=0.000800;
- Holm p=0.004800.

But reserve:

- rho(D, sampled span)=-0.002567.

Therefore do **not** promote range/span as a general ecological predictor of polymorphism. Treat the discovery association primarily as an opportunity diagnostic and control variable.

### Latitude

Discovery absolute-latitude association unsupported:

- rho=-0.033906;
- Holm p=1.

Do not retain as a positive explanatory result.

### Family

Family-level clustering unsupported in discovery after the prospective family:

- raw p=0.300635;
- Holm p=1.

Do not imply broad family-level phylogenetic conservatism.

### Pollination mode and life form

Prospective coverage gates failed before D was opened. They remain unavailable rather than post-hoc backfilled.

This is a design strength, not a null biological result.

## Rejected former Claim — recurrent colour direction

Do not restore the old M1-M2 variance-direction claim.

After controlling total continuous colour spread, D does not preferentially expand along the frozen M1-M2 subspace, and M1-M2 does not beat random two-dimensional orientations.

The positive mean-M3 result remains auxiliary:

- it concerns average position in the measured colour representation;
- it is carried by the flower component rather than background;
- it does not mean polymorphism expands along M3;
- no flowering-stage interpretation is licensed.

## Global shared-boundary result

Keep as a conceptual negative contrast, not as a failed centerpiece.

The programme repeatedly found recoverable within-species spatial organization without convincing/reliably identifiable cross-species common boundary geography.

The central contrast is now:

> **repeatable species-level organization does not require repeatable boundary location or direction across species.**

This is more informative than simply saying that the shared-boundary test failed.

## Measurement limitations after Step 9

The remaining limitation is explicit and bounded rather than hidden.

- Every discovery/reserve species begins from 100 fixed acquired photos.
- Clear ROI/flip technical failures are separable and do not explain the main spatial result.
- Ambiguous-palette rows are strongly associated with baseline D and cannot be declared purely technical or purely biological.
- Step 9 four-state completion intervals have median width about 0.074 in both tranches and 95th-percentile width about 0.25-0.28.
- Despite that uncertainty, the flagship spatial relationship survives both uniform D_min4 and D_max4 endpoint rules in discovery and reserve, including matched-background differential.
- Genus clustering does not pass both reserve endpoints and stays secondary.

## Current paper spine

1. **Flower-colour polymorphism is a continuous species-level property in fixed photo-derived frames.**
2. **The amount of polymorphism predicts stronger within-species geographic organization.**
3. **That organization is species-specific rather than a universal shared boundary/direction.**
4. **The spatial result replicates and survives opportunity, background, clear technical-failure, and ambiguity-endpoint stress tests.**
5. **Genus-level taxonomic structure exists but is ambiguity-sensitive and remains secondary.**

## Connection back to Camellia / chun

Keep the two papers separate but conceptually adjacent.

- `chun`: repeated visible flower-colour transitions can reuse molecular modules without replaying one complete molecular package; generation through evolutionary time is hierarchical and transition-dependent.
- `fcp`: the amount of standing/observed flower-colour variation is geographically organized within species, but the location/direction of that organization is not universally shared across species.

Together they support a broader programme in which **generation, persistence, and spatial organization are distinct levels of flower-colour evolution**. Do not merge their datasets or claim that one empirically proves the mechanism of the other.

## Next active task

Move from hypothesis discovery to manuscript assembly.

Do not add another broad ecological predictor family unless it resolves a named reviewer-level vulnerability that is not already addressed above.

Immediate paper work:

1. freeze title / one-sentence claim;
2. build Results in the evidence order above;
3. make the discovery-reserve spatial replication + ambiguity robustness the central figure;
4. place genus clustering as a secondary result/figure or supplement depending journal length;
5. move shared-boundary and directionality failures into the estimand contrast that motivates the final interpretation.
