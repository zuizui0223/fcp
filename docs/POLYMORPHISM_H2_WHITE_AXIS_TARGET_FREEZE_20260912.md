# H2 narrowed target freeze — achromatic/chromatic axis

Date: 2026-09-12 JST

## Status

This is a **targeted post-audit decomposition**, not a new preregistered discovery test.

The label-free H2 analysis found strong discovery/reserve directional concentration and survived a coarse-state-preserving structured null. A subsequent audit then showed that removing an a-priori white-versus-equal-nonwhite contrast eliminates the excess concentration under the same structured null, and the non-white top-two subset does not show excess discovery concentration under that null.

Therefore the broad claim “flower-colour polymorphism follows recurrent hue directions beyond white/non-white structure” is closed as unsupported. The only remaining H2 candidate is the fixed achromatic/chromatic contrast itself.

## Fixed contrast

In the normalized nine-colour flower-palette coordinates

`[white, yellow, orange, red, pink, magenta, purple, blue, bronze]`,

define

`q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`.

This is a zero-sum unit vector comparing the white coordinate with the equal mean of all eight non-white coordinates. Its sign is irrelevant.

No loading is fitted from the data for this targeted statistic.

## Analysis population

Use the already frozen label-free H2 species sets separately for:

- primary threshold: coarse H1 second-state fraction >=0.10 and continuous minor-cluster fraction >=0.10;
- strict threshold: corresponding >=0.20 gates.

Do not add species after viewing this targeted statistic.

## Statistic

For each unit species displacement vector `u_i = Delta_i / ||Delta_i||`, calculate

`W = mean_i (u_i^T q_white)^2`.

This is the mean squared alignment with the fixed achromatic/chromatic contrast.

Report discovery and species-disjoint reserve separately. The reserve does not refit or rotate `q_white`.

## Structured null

Use the same coarse-state-preserving construction null already frozen for H2:

1. within each cohort and each frozen four-state coarse morph, permute normalized nine-colour palette rows across the already selected H2 species;
2. preserve every species x coarse-morph row count;
3. refit the label-free Hellinger two-means within each species;
4. reconstruct `Delta_i` and unit axes;
5. recalculate `W` against the fixed `q_white`.

Use 999 null worlds and the upper-tail Monte Carlo p-value `(1 + #null >= observed)/1000`.

## Decision rule

Existing-cohort targeted support requires **both** discovery and reserve structured-null `p < 0.05` at the primary 0.10 threshold. The strict 0.20 result is a sensitivity requirement for a strong version of the claim.

Possible verdicts:

- `WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT`
- `WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_ONLY`
- `WHITE_AXIS_TARGETED_NOT_SUPPORTED`

Because the candidate axis was motivated by already opened H2 results, even a positive verdict is retrospective evidence. It may define a frozen prospective target for the 42,111-frame high-depth expansion, but must not be described as an untouched confirmatory test in the existing 369/363 cohorts.

## Allowed claim if positive

“Across the high-depth validation cohorts, within-species flower-colour polymorphism is disproportionately aligned with an achromatic–chromatic (white versus non-white) axis, beyond that expected from the frozen coarse-state composition and global state-to-palette mapping.”

Do not replace “achromatic–chromatic” with a specific pigment pathway, pigment-loss mechanism, pollinator mechanism, or adaptive interpretation without additional evidence.
