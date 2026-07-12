# Two-patch spatial extension

## Question

When can spatially heterogeneous selection maintain floral colour polymorphism despite gene flow and an overall directional bias?

## Model

Let `p1` and `p2` be the frequencies of morph A in two equal patches.

\[
\dot p_1=p_1(1-p_1)\{d+h+b(1-2p_1)\}+m(p_2-p_1),
\]

\[
\dot p_2=p_2(1-p_2)\{d-h+b(1-2p_2)\}+m(p_1-p_2).
\]

Parameters:

- `b`: within-patch balancing strength.
- `d`: net directional bias shared across patches.
- `h`: spatial contrast in local directional selection.
- `m`: symmetric gene flow.

The model deliberately separates three forces that are often conflated empirically: balancing selection, spatial heterogeneity, and homogenization by gene flow.

## Rare invasion of A into B

Near the B-monomorphic boundary, write the rare A frequencies as the vector `x=(p1,p2)^T`. The linearized system is

\[
\dot x = J_B x,
\]

with

\[
J_B=
\begin{pmatrix}
 d+h+b-m & m\\
 m & d-h+b-m
\end{pmatrix}.
\]

The dominant eigenvalue is

\[
\lambda_{A|B}=d+b-m+\sqrt{h^2+m^2}.
\]

A can invade B iff this quantity is positive.

## Rare invasion of B into A

Likewise, linearization around the A-monomorphic boundary gives

\[
\lambda_{B|A}=b-d-m+\sqrt{h^2+m^2}.
\]

B can invade A iff this quantity is positive.

## Exact protected-polymorphism condition

Mutual invasibility requires

\[
\lambda_{A|B}>0
\quad\text{and}\quad
\lambda_{B|A}>0.
\]

Equivalently,

\[
\boxed{
 b+\sqrt{h^2+m^2}-m>|d|
}
\]

Define the effective spatial contribution

\[
H_{\mathrm{eff}}(h,m)=\sqrt{h^2+m^2}-m.
\]

Then the phase boundary becomes

\[
\boxed{b+H_{\mathrm{eff}}>|d|}.
\]

This recovers the one-population result `b>|d|` when `h=0` or when migration becomes very large.

## Interpretation

The result gives a compact decomposition:

\[
\text{balancing force}
+
\text{migration-attenuated spatial heterogeneity}
>
\text{directional bias}.
\]

Important limits:

1. No gene flow (`m=0`):

\[
H_{\mathrm{eff}}=|h|.
\]

Spatial contrast contributes fully.

2. Strong gene flow (`m\to\infty`):

\[
H_{\mathrm{eff}}\sim \frac{h^2}{2m}\to0.
\]

Spatial contrast is homogenized away.

3. No spatial contrast (`h=0`):

\[
H_{\mathrm{eff}}=0.
\]

The model reduces exactly to the original one-population boundary.

## Empirical predictions

The theory predicts that the probability of species-level floral colour polymorphism should:

- increase with ecologically relevant environmental heterogeneity;
- increase with evidence for balancing or negative frequency-dependent selection;
- decrease with strong overall directional selection;
- show an interaction between heterogeneity and connectivity/gene flow, because heterogeneity matters less when homogenization is strong.

The interaction is not optional decoration. It is a direct consequence of the analytical boundary.

## What the meta-analysis should test

The first empirical target is not a literal estimate of `b`, `h`, `m`, and `d` for every species. Instead, map measurable proxies onto the theoretical decomposition:

- `b_proxy`: evidence for negative frequency dependence, pollinator switching, rare-morph advantage, or other balancing mechanisms;
- `h_proxy`: within-range climatic, elevational, habitat, or pollinator-community heterogeneity;
- `m_proxy`: dispersal/connectivity proxies, range continuity, population isolation, or direct gene-flow estimates where available;
- `d_proxy`: consistent morph-specific fitness advantage across environments.

The key comparative prediction is a positive heterogeneity effect that weakens as connectivity increases.

## Next theoretical extensions

Priority order:

1. Unequal patch sizes and asymmetric migration.
2. More than two patches; ask whether a spectral summary replaces `H_eff`.
3. Temporal environmental switching.
4. Mutation/origin and finite-population loss, separating origin, invasion, persistence, and observation.
5. Detection error in the empirical classification of polymorphic versus monomorphic species.
