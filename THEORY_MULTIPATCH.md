# Multi-patch spectral phase boundary

## 1. General landscape model

Consider `n` habitat patches. Let `p_i` be the frequency of floral colour morph A in patch `i`.

Local selection is decomposed into:

- `b`: within-patch balancing strength, shared across patches;
- `d`: global directional bias, shared across patches;
- `h_i`: patch-specific directional deviation;
- `m`: migration scale;
- `L`: weighted graph Laplacian describing the migration network.

The nonlinear frequency dynamics can be written schematically as

\[
\dot p_i = p_i(1-p_i)\{d+h_i+b(1-2p_i)\} - m(Lp)_i.
\]

The phase boundary is determined by linearization at the two monomorphic boundaries.

## 2. Rare A invading all-B

Near `p=0`, rare A obeys

\[
\dot x = \left[(b+d)I + \operatorname{diag}(h)-mL\right]x.
\]

The asymptotic invasion exponent is the dominant eigenvalue

\[
\lambda_{A\mid B}
= b+d+H_A,
\]

where

\[
H_A = \lambda_{\max}\left[\operatorname{diag}(h)-mL\right].
\]

## 3. Rare B invading all-A

Near `p=1`, rare B obeys

\[
\dot y = \left[(b-d)I + \operatorname{diag}(-h)-mL\right]y.
\]

Hence

\[
\lambda_{B\mid A}
= b-d+H_B,
\]

with

\[
H_B = \lambda_{\max}\left[\operatorname{diag}(-h)-mL\right].
\]

Protected polymorphism requires mutual invasibility:

\[
\lambda_{A\mid B}>0,
\qquad
\lambda_{B\mid A}>0.
\]

## 4. Spectral compression

Define

\[
H_{\mathrm{bal}}=\frac{H_A+H_B}{2},
\qquad
D_{\mathrm{asym}}=\frac{H_A-H_B}{2}.
\]

Then

\[
\lambda_{A\mid B}
= b+H_{\mathrm{bal}}+d+D_{\mathrm{asym}},
\]

\[
\lambda_{B\mid A}
= b+H_{\mathrm{bal}}-(d+D_{\mathrm{asym}}).
\]

Therefore the exact multi-patch phase boundary is

\[
\boxed{
 b+H_{\mathrm{bal}}
 >
 \left|d+D_{\mathrm{asym}}\right|
}
\]

for protected polymorphism.

This is the central compression result.

- `H_bal` is the landscape contribution that supports reciprocal invasion.
- `D_asym` is the landscape-induced directional shift caused by asymmetric arrangement of selective environments on the migration network.

Thus spatial heterogeneity does not enter only through its variance. Its effect depends on how selective environments are arranged relative to connectivity.

## 5. Recovery of the two-patch result

For two equally connected patches with

\[
h=(+q,-q),
\]

we obtain

\[
H_A=H_B=\sqrt{q^2+m^2}-m.
\]

Hence

\[
H_{\mathrm{bal}}=\sqrt{q^2+m^2}-m,
\qquad
D_{\mathrm{asym}}=0,
\]

and the general criterion reduces exactly to

\[
\boxed{
 b+\sqrt{q^2+m^2}-m>|d|.
}
\]

So the two-patch model is a special case of the spectral landscape model.

## 6. Biological interpretation

The same amount of environmental variation can produce different polymorphism outcomes depending on network placement.

For example, patches favoring one morph may be:

1. clustered together,
2. alternating across the landscape,
3. isolated in peripheral patches,
4. concentrated in highly connected hubs.

These landscapes can have identical variance in `h_i` but different `H_bal` and `D_asym`.

The theory therefore predicts that a meta-analysis based only on environmental variance is incomplete. The relevant quantity is selection heterogeneity filtered through connectivity.

## 7. Meta-analytic predictions

The theory generates a hierarchy of increasingly strong tests.

### Prediction 1: coarse heterogeneity effect

Species occupying more heterogeneous environments should be more likely to show floral colour polymorphism.

### Prediction 2: connectivity interaction

The positive effect of environmental heterogeneity should weaken as effective gene flow or landscape connectivity increases.

### Prediction 3: topology matters

Among species with similar environmental variance, polymorphism should be more common when contrasting selective environments are spatially arranged in ways that increase `H_bal`.

### Prediction 4: asymmetric landscapes bias fixation

If the spatial arrangement of environments produces nonzero `D_asym`, one morph should gain a net directional advantage even when the arithmetic mean of `h_i` is near zero.

## 8. What is new conceptually

The theory shifts the explanatory target from individual mechanisms to a phase boundary.

Different ecological mechanisms may alter:

- `b`: balancing selection,
- `d`: global directional selection,
- `h`: spatially varying selection,
- `L` and `m`: homogenization structure.

But the observed monomorphic-versus-polymorphic state depends only on whether the compressed inequality is crossed.

This creates a direct bridge from mathematical theory to comparative meta-analysis while preserving mechanism uncertainty below the compression layer.
