# FCP polymorphism mainline after Step 5b

Date: 2026-09-09 JST

## Mainline decision

The paper should no longer be organized around a globally shared flower-colour boundary or around a privileged recurrent M1–M2 variance direction. The strongest empirical spine is now a species-level polymorphism gradient linked to within-species geographic organization.

## Claim 1 — flower-colour polymorphism is a measurable species attribute

Retain.

Discovery contains 369 species with four-state Simpson diversity D spanning approximately 0 to 0.708. The second morph reaches >=10% in 46.61% and >=20% in 26.56% of the admitted 369-species frame. These are conditional properties of the sampled/classifiable frame, not prevalence estimates for all flowering plants.

The mixed-uncertain Step 2 does not license treating excluded rows as hidden morphs or numerically correcting D. Most `mixed_uncertain` rows are technical/structural missingness; the prospectively fixed independent-centroid bridge enrichment narrowly failed. D remains a conservative four-state summary of admitted rows.

## Former Claim 2 — polymorphism follows recurrent M1–M2 variance directions

Reject.

After total within-species continuous flower-colour variance is controlled, D is negatively associated with the fixed M1–M2 variance share and M1–M2 does not outperform random two-dimensional orientations. The earlier positive SD(M1)/SD(M2) associations track total continuous colour spread rather than privileged directionality.

Do not use the stable M1–M3 recurrence result to restore this claim.

## Auxiliary colour-space result — mean M3 position, not M3 variance

Retain as a species-attribute correlate, not as a major variance-direction claim.

In discovery, D is positively associated with mean frozen M3 position even after green/brown/black flower-palette nuisance coordinates are removed and after nuisance fraction is controlled. In the reserve component diagnostic, the D association is carried by flower M3 rather than background M3. SD(M3) is essentially unrelated to D.

This means polymorphic species occupy a somewhat different average region of the measured flower-colour representation, but polymorphism does not preferentially expand along M3. No flowering-stage interpretation is allowed because no direct phenological-stage label exists.

## Claim 2 (replacement) — greater polymorphism predicts stronger within-species geographic colour organization

This should replace the rejected recurrent-direction claim as the second major empirical claim.

### Discovery

The already completed 369-species exact-randomization spatial omnibus was re-aggregated without regenerating any spatial statistic or permutation.

- second morph >=10%: 172 species, mean spatial rho = 0.034565108, p = 0.001;
- second morph >=20%: 98 species, mean spatial rho = 0.027081099, p = 0.001;
- <10% complement: 197 species, mean rho = 0.020434835;
- >=10% minus <10% difference = 0.014130273, p = 0.014;
- continuous rho(D, species spatial rho) = 0.089213, randomization p = 0.034.

Because this subset idea was formulated after the whole-frame spatial result was known, discovery Step 5 is explicitly post-outcome exploratory.

### Independent reserve photo tranche

The same thresholds and aggregation rules were frozen before joining reserve polymorphism membership to reserve spatial outcomes. The 363-species reserve inference arrays were reused without new spatial permutations.

Primary >=10% subset (151 species):

- primary mean rho = 0.032299123, p = 0.001;
- observer-pair exclusion mean rho = 0.031774569, p = 0.001;
- quarter-stratified mean rho = 0.032299123, p = 0.001;
- matched-background differential mean rho = 0.013621726, p = 0.005.

Sensitivity >=20% subset (90 species):

- primary mean rho = 0.039077034, p = 0.001;
- observer-pair exclusion p = 0.001;
- quarter p = 0.001;
- matched-background differential p = 0.003.

Reserve >=10% minus <10% difference = 0.011671678, p = 0.041. The threshold-free reserve gradient also recurs: rho(D, species primary spatial rho) = 0.101601, p = 0.025.

Thus the direction, subset association, dilution contrast, observer robustness, calendar robustness, and continuous D gradient all recur in a disjoint reserve photo tranche. In the reserve, the fixed polymorphic subset additionally passes the matched-background differential that did not pass in the whole 363-species frame.

### Admissible interpretation

The evidence supports the statement that, within these photo-derived global frames, species with non-trivial discrete flower-colour polymorphism show stronger within-species geographic organization of flower-colour variation, and that this pattern recurs in an independent photo tranche.

The result does not imply a globally shared boundary. Different species may be geographically structured in different places and directions. That distinction is now central rather than a limitation: **species-level spatial organization can be repeatable even when cross-species boundary geography is not identifiable or shared.**

Do not translate this into adaptation, common climate thresholds, population-genetic differentiation, or causal selection.

## Claim 3 — what kinds of species carry polymorphism

Still open; this is now the main unfinished empirical layer.

The pre-outcome Step-4 covariate audit established:

- genus: derivable internally from binomial species names;
- range-size proxy: available from frozen colour-blind metadata;
- latitude centroid: derivable from frozen colour-blind sampled coordinates;
- family: not yet available for the current 369-species frame;
- life form: no tracked source yet;
- pollination mode: no tracked source yet.

No D association with these unopened Step-4 covariates has yet been used for selection.

## Current paper spine

1. **Polymorphism is a measurable species-level property** in a large photo-derived global frame and is reproducible across observer/time controls already completed upstream.
2. **More polymorphic species are more geographically organized internally.** This is the strongest new result: discovery + independent reserve, threshold and continuous versions, nuisance controls, and reserve matched-background support.
3. **Species attributes associated with polymorphism** are the next layer to estimate after the covariate panel/source rules are frozen.
4. **No shared global flower-colour boundary is required.** The failed sharedness qualifications and common-boundary analyses become an estimand/identifiability contrast: within-species organization is recoverable while a universal cross-species boundary is not supported/identifiable with the present geometry.

The stable M1–M3 modes remain measurement/representation diagnostics. They are not the paper's main biological claim.

## Next active task

Freeze Step 4 as two predeclared blocks before opening any new D associations:

- Block A: internally available taxonomy/geography (`genus`, range-size proxy, latitude centroid, plus the already established mean-M3 species attribute as a known positive descriptor rather than a newly tested predictor);
- Block B: externally sourced `family`, `life form`, and `pollination mode`, with source hierarchy, matching rules, missingness policy, and test family fixed before acquisition/association testing.

Block B must be acquired regardless of Block-A outcomes so that external traits are not selectively back-filled after seeing internal covariate results.
