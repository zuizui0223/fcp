# Theory note: the minimal phase boundary

## Question

Why are some species polymorphic in floral colour whereas others are monomorphic?

The first model deliberately does **not** attempt to reproduce every known mechanism. Its role is to derive the smallest phase boundary separating monomorphic and protected-polymorphic regimes.

## Model

Let `p` be the frequency of colour morph A. The log-fitness difference is

\[
\Delta(p)=\log W_A-\log W_B=d+b(1-2p),
\]

where:

- `b` is the strength of frequency-dependent balancing selection. Positive values mean rare-morph advantage.
- `d` is net directional selection favouring A when positive and B when negative.

Frequency dynamics follow

\[
\dot p=p(1-p)\Delta(p).
\]

## Boundary invasion rates

A rare A morph invading a B-monomorphic population has invasion rate

\[
r_{A|B}=\Delta(0)=d+b.
\]

A rare B morph invading an A-monomorphic population has invasion rate

\[
r_{B|A}=-\Delta(1)=b-d.
\]

Protected polymorphism requires mutual invasibility:

\[
r_{A|B}>0,\qquad r_{B|A}>0,
\]

which is equivalent to the phase boundary

\[
\boxed{b>|d|}.
\]

The stable internal equilibrium is

\[
p^*=\frac{1}{2}\left(1+\frac{d}{b}\right).
\]

Thus the minimal prediction is not “frequency dependence causes polymorphism”. It is sharper: **balancing forces must exceed the magnitude of net directional bias**.

## Why this is useful

Later biological mechanisms can be mapped onto these two effective axes:

\[
b=b_{\mathrm{pollinator}}+b_{\mathrm{spatial}}+b_{\mathrm{temporal}}+\cdots
\]

and

\[
d=d_{\mathrm{abiotic}}+d_{\mathrm{herbivory}}+d_{\mathrm{pigment\ cost}}+\cdots.
\]

Different mechanisms can therefore occupy the same point in `(b, d)` space and produce the same observed phase. This immediately creates a mechanism-identification problem: observing polymorphism alone does not identify which mechanism generated the effective balancing strength.

## Next extensions

1. Split `b` into mechanistically distinct components and test when only their sum is identifiable.
2. Add spatial patches and migration to derive a gene-flow-dependent boundary.
3. Add mutation/origin and finite-population loss, separating origin, invasion, persistence, and observation.
4. Link the theoretical phase boundary to a phylogenetic meta-analysis of polymorphic versus monomorphic species.
