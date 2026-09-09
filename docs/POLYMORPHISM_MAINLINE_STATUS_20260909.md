# FCP polymorphism mainline status after Steps 1–3

Date: 2026-09-09 JST

## Current paper spine

### Claim 1 — species-level polymorphism is measurable: RETAIN

The frozen discovery frame contains 369 species with continuous four-state Simpson diversity D. The previously reported fingerprint is reproduced exactly enough for lineage tracking: maximum D 0.707645, second morph >=10% in 46.61% of species, and >=20% in 26.56%.

The mixed-uncertain Step 2 does not strengthen this claim by treating missing rows as hidden morphs. Of 50,000 measured discovery rows, 24,623 are `mixed_uncertain`, but only 2,332 (9.47% of mixed rows) are chromatically evaluable ambiguous-palette cases; the remainder are technical/structural missingness. The frozen bridge-geometry enrichment misses its gate (p=0.051974), so no numerical D-underestimation claim is allowed.

### Former Claim 2 — polymorphism follows a recurrent small set of variance directions: REJECT

After controlling total within-species continuous colour variance, D is negatively associated with the fixed M1–M2 variance share: rho(D,R12)=-0.164727, label-permutation p=0.9991, bootstrap 95% CI [-0.273377,-0.054861]. M1–M2 also does not beat random 2D orientations (p_rotation=0.723928).

The earlier positive SD(M1)/SD(M2) associations are therefore scale-sensitive consequences of total colour spread: rho(D,total variance)=0.825757. They are not evidence that polymorphism preferentially follows M1–M2.

### M3 diagnosis — mean-position association survives as a separate species attribute

M3 cannot rescue former Claim 2 because the Step-1 variance-direction test is already frozen negative. However, the mean M3 position is independently associated with D in the discovery frame and is flower-biological rather than background-driven:

- discovery rho(D, mean M3 full-12)=0.255991, permutation p=4.99975e-05;
- after removing green/brown/black, rho=0.300195, p=4.99975e-05;
- partial rho controlling flower nuisance fraction=0.271723, p=4.99975e-05;
- reserve diagnostic rho(D, flower M3)=0.316098 versus rho(D, background M3)=0.007629;
- discovery rho(D, SD M3)=-0.015573.

Thus M3 belongs in the next species-attribute analysis as a mean-position correlate, not as a polymorphism direction or flowering-stage result.

### Claim 3 — polymorphic species show geographic structure, but not a shared global boundary: RETAIN, next test pending

The existing discovery and reserve spatial results remain valid within their prior scope. The next unrun test is the within-polymorphic-species subset analysis, after the species-attribute covariate set is frozen.

## Next ordered tasks

1. Freeze and audit the species-level covariate panel before testing D associations: taxonomy (family/genus), life form, pollination mode, range-size proxy, latitude centroid, and the now-separated mean M3 position.
2. Run the frozen species-attribute association model(s) without outcome-driven covariate selection.
3. Re-run within-species geographic structure on a prospectively defined polymorphic subset only.
4. Keep independent global-axis acquisition stopped; it is no longer part of the active mainline.
