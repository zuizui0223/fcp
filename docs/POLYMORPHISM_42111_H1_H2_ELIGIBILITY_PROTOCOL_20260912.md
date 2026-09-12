# 42,111-species polymorphism H1/H2 eligibility protocol

Date frozen: 2026-09-12 (JST)

## Purpose

This protocol separates the global **sampling frame** from the species sets that are actually eligible for flower-colour polymorphism and direction-of-variation inference.

The outer sampling universe is the frozen 42,111-species RGFCA discovery frame. The existing 369-species discovery cohort remains a high-depth validation cohort and is not redefined as the population denominator.

No species is called monomorphic merely because it fails a depth, image-measurement, or second-mode gate.

## Existing evidence that constrains the new design

1. The 42,111-species capacity census retains every discovered species, including structural zero-capacity cases. After the observer cap, 4,730 species have at least 100 eligible raw photos; 5,558 have at least 80; 6,770 at least 60; 7,632 at least 50; 8,753 at least 40; 10,330 at least 30; 12,985 at least 20; 18,301 at least 10.
2. The existing 369-species discovery polymorphism analysis required at least 40 classifiable flower observations per species. In that cohort, 46.612% (172/369) had a second admitted four-state colour at frequency >=0.10 and 26.558% (98/369) at frequency >=0.20.
3. The frozen mixed-image analysis did **not** qualify `mixed_uncertain` as an additional biological morph. Its global bridge-geometry test missed the frozen 0.05 threshold (`p = 0.051974`), so mixed/ambiguous rows remain structural measurement uncertainty and cannot be used to manufacture a second mode.
4. The previous recurrent-direction claim was rejected. The old statistic `rho(D, R12)` was negative and its label-permutation and random-orientation tests failed. H2 below therefore tests a different estimand: the geometry of species-specific primary-to-secondary colour displacement vectors.
5. One-photo breadth measurement classified 18,457 / 42,111 species anchors (43.829%). This one-anchor endpoint identifies only the species-equal distribution of observed coarse colour states; it does not identify species polymorphism, a modal colour, D, or C*/S*.

## Analysis populations

### U0 — global sampling universe

All 42,111 frozen discovered species.

U0 is the denominator for sampling opportunity and missingness. It is **not** the denominator for biological polymorphism prevalence unless a later measurement model explicitly identifies prevalence under missingness.

### U100 — high-depth opportunity frame

Species with `after_observer_cap >= 100` in the frozen capacity census.

Current size: **4,730 species**.

This is the prespecified expansion frame for H1/H2 because 100 raw photos is the closest large-scale design to the legacy 100-photo high-depth cohorts and gives a realistic opportunity to recover the legacy `n_classifiable >= 40` gate. The previous Step-8 maximum of 20 photos remains useful for coarse diversity/opportunity summaries but is not treated as sufficient by itself for stable two-mode geometry.

### V — legacy high-depth validation frame

The frozen discovery validation cohort with `n_classifiable >= 40`.

Current size: **369 species**.

V is used to validate measurement behaviour and to prototype H2 without changing the outer 42,111-species denominator.

## H1 — valid non-trivial secondary flower-colour mode

H1 is a measurement-validity gate, not a prevalence claim.

For species `i`, let admitted observations be only rows that pass `global_classifiable` and belong to one of the four frozen biological states:

- `white`
- `yellow_orange`
- `red_pink`
- `blue_purple`

A species is **H1-primary eligible** iff all of the following hold:

1. `n_classifiable >= 40`;
2. at least two admitted biological states are observed;
3. the second-most frequent admitted state has frequency `second_fraction >= 0.10`;
4. `mixed_uncertain`, ROI/flip failures, no-biological-palette rows, and unresolved rows are not promoted to biological states.

A prespecified strict sensitivity uses `second_fraction >= 0.20` with the same other gates.

The 0.10 and 0.20 thresholds are reused from the already frozen high-depth polymorphism analysis; they are not selected from the new 42,111-species outcomes.

Current validation-frame counts are therefore:

- H1 primary (>=10% second state): **172 / 369**;
- H1 strict (>=20% second state): **98 / 369**.

For U100, the biological H1 eligible N is **not known before high-depth image measurement**. The exact premeasurement ceiling is 4,730 species.

## H2 — geometry of the second-mode displacement

H2 is evaluated only in H1-eligible species.

For each H1-eligible species `i`:

1. call the most frequent admitted four-state morph the primary state and the second-most frequent admitted state the secondary state;
2. within each of those two states, normalize the nine biological flower-palette coordinates (`white`, `yellow`, `orange`, `red`, `pink`, `magenta`, `purple`, `blue`, `bronze`) to unit row sum;
3. let `mu_i,1` and `mu_i,2` be the mean normalized biological palette vectors for the primary and secondary states;
4. define the species-specific displacement

   `Delta_i = mu_i,2 - mu_i,1`.

Because both centroids sum to one, every `Delta_i` lies in the eight-dimensional zero-sum tangent subspace.

### Primary H2 estimand: sign-invariant directional concentration

The primary question is whether the set of species-specific displacement **axes** is more concentrated than expected under isotropic orientation in the zero-sum subspace.

For non-zero `Delta_i`, define `u_i = Delta_i / ||Delta_i||` and

`M = (1/N) sum_i u_i u_i^T`.

The primary statistic is the largest eigenvalue `lambda_1(M)`. This is sign-invariant: species varying along the same colour axis in opposite primary/secondary orientations still contribute to the same recurrent direction.

The frozen null independently applies Haar-random orthogonal rotations within the eight-dimensional zero-sum subspace to each observed `u_i`, preserving sample size and vector norms. Use 10,000 null realizations and the upper-tail Monte Carlo p-value `(1 + #null >= observed)/(10001)`.

### Secondary H2 diagnostics

Report, but do not substitute for the primary test:

- the signed mean resultant `||mean(u_i)||`;
- the leading H2 axis loadings in the nine biological palette coordinates;
- the fraction of squared displacement captured by the leading axis;
- the distribution of `||Delta_i||`;
- results under the strict H1 threshold (`second_fraction >= 0.20`).

### Required validation before 42,111-frame interpretation

H2 must first be calculated in the existing high-depth discovery frame and then repeated in the species-disjoint reserve using the same eligibility rule and no rotation of the learned coordinate system. A direction found only in discovery is not transportable evidence.

## Planning calculation for U100

The one-anchor 42,111-species measurement had a classifiable fraction of 0.438294. This is not itself an H1 estimate. It can be used only for resource planning.

Under an explicitly labelled independent-binomial planning approximation, 100 raw photos with classifiability probability 0.438294 gives an approximately 0.808 probability of reaching `n_classifiable >= 40`. The observed capacity-stratum classifiability range from the breadth missingness audit (about 0.399 to 0.478) gives a planning envelope of about 0.532 to 0.953.

Applied to the 4,730-species U100 capacity ceiling, this corresponds to roughly 2,515–4,508 species, central planning value about 3,823, reaching the **measurement-depth gate only**. These are not biological H1 or H2 eligible counts and must never be reported as polymorphism prevalence.

## Decision ledger

- `U0 = 42,111`: sampling denominator.
- `U100 = 4,730`: exact high-depth opportunity ceiling.
- `V = 369`: frozen high-depth validation cohort.
- `V_H1_10 = 172`: exact validation count at the frozen 10% secondary-state threshold.
- `V_H1_20 = 98`: exact strict validation count at the frozen 20% threshold.
- `U100_H1`: unknown until high-depth measurement; do not impute from V.
- `U100_H2`: unknown until H1 is evaluated and non-zero `Delta_i` vectors are constructed.

## Hard stops

1. Do not treat one-photo `breadth_classifiable` as evidence that a species is monomorphic.
2. Do not treat `mixed_uncertain` as a fifth morph or as a secondary mode.
3. Do not reuse the rejected M1/M2 `rho(D, R12)` directionality test as H2.
4. Do not report the binomial planning projection as an observed eligible N.
5. Do not open H3 ecological/phylogenetic predictors to select or tune the H1/H2 thresholds.
