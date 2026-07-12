# Stage 4 — Temporal and spatiotemporal theory

## 1. Pure temporal fluctuation in a well-mixed population

Let

\[
\dot p=p(1-p)\{d(t)+b(1-2p)\}.
\]

Near the boundary \(p=0\),

\[
\frac{d}{dt}\log p=b+d(t),
\]

so the long-run invasion exponent of A into B is

\[
\gamma_{A|B}=b+\overline d.
\]

Near \(p=1\), the reciprocal invasion exponent is

\[
\gamma_{B|A}=b-\overline d.
\]

Hence protected polymorphism requires

\[
\boxed{b>|\overline d|}.
\]

The variance of additive temporal fluctuation does not appear. Therefore temporal fluctuation in directional selection is not, by itself, a universal balancing mechanism in this minimal continuous-time frequency model.

This negative result matters for comparative work: a species experiencing strongly variable environments should not automatically be predicted to be more polymorphic. Temporal variation must interact with another mechanism, such as generation overlap, storage, density dependence, spatial structure, or nonlinear genotype-by-environment responses.

## 2. General spatiotemporal rare-invasion problem

For a rare morph distributed over \(n\) patches, let \(x(t)\) be its patch-abundance vector and write

\[
\dot x=A(t)x.
\]

The matrix \(A(t)\) combines local rare-morph growth and movement among patches.

The invasion exponent is the top Lyapunov exponent

\[
\gamma=\lim_{t\to\infty}\frac{1}{t}\log \|\Phi(t)\|,
\]

where \(\Phi(t)\) is the fundamental solution matrix. In a periodic environment of period \(T\), this becomes the Floquet exponent

\[
\gamma=\frac{1}{T}\log \rho\{\Phi(T)\}.
\]

For reciprocal invasions, let \(\gamma_A\) be A into B and \(\gamma_B\) be B into A. Then

\[
\boxed{\gamma_A>0\quad\text{and}\quad\gamma_B>0}
\]

is the exact protected-polymorphism criterion.

Define

\[
B_{ST}=\frac{\gamma_A+\gamma_B}{2},\qquad
D_{ST}=\frac{\gamma_A-\gamma_B}{2}.
\]

Then the criterion compresses to

\[
\boxed{B_{ST}>|D_{ST}|}.
\]

This is the most general phase-boundary form currently recovered in FCP.

## 3. Why temporal variance is not enough

When the matrices \(A(t_1)\) and \(A(t_2)\) do not commute, the long-run growth rate depends on temporal order, not only on marginal means and variances. Therefore two species can experience the same mean and variance of environmental conditions but have different invasion exponents because of:

- autocorrelation,
- switching rate,
- covariance between local selection and dispersal,
- spatial synchrony,
- sequence/order of environmental states.

The comparative implication is that raw temporal climatic variance is only a proxy. A stronger empirical target is the interaction between temporal variability, spatial structure, and life-history buffering.

## 4. Nested recovery

The theory now forms a nested hierarchy:

1. One population, static:
   \[
   b>|d|.
   \]
2. Two patches, static:
   \[
   b+\sqrt{h^2+m^2}-m>|d|.
   \]
3. Many patches, static:
   principal-eigenvalue invasion criterion.
4. Many patches, time varying:
   top Lyapunov/Floquet invasion criterion.
5. All cases:
   \[
   \text{effective balancing}>|\text{effective directional asymmetry}|.
   \]

The last expression is a compression of exact reciprocal invasion rates, not a claim that all mechanisms are mechanistically identical.
