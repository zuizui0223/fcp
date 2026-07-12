# FCP — Floral Colour Polymorphism phase theory

A mathematical and comparative framework for the question:

> Why are some plant species polymorphic in floral colour while others are monomorphic?

The project starts from phase boundaries rather than from a single candidate mechanism, then uses the recovered theory to generate explicit comparative predictions for meta-analysis.

## Stage 1 — minimal one-population result

For two colour morphs A and B,

\[
\Delta(p)=d+b(1-2p), \qquad
\dot p=p(1-p)\Delta(p).
\]

The two boundary invasion rates are

\[
r_{A|B}=d+b,\qquad r_{B|A}=b-d.
\]

Protected polymorphism therefore exists exactly when

\[
\boxed{b>|d|}.
\]

Here `b` is effective balancing strength and `d` is net directional bias.

See [`THEORY.md`](THEORY.md).

## Stage 2 — spatial heterogeneity plus gene flow

For two equal patches with opposite local selection contrasts `+h` and `-h`, coupled by symmetric migration `m`, the dominant rare-invasion exponents yield the exact condition

\[
\boxed{
 b+\sqrt{h^2+m^2}-m>|d|
}
\]

or

\[
\boxed{b+H_{\mathrm{eff}}>|d|},
\qquad
H_{\mathrm{eff}}=\sqrt{h^2+m^2}-m.
\]

This gives a compact theoretical result:

> Spatial heterogeneity promotes polymorphism, but its contribution is attenuated by gene flow.

The immediate comparative prediction is therefore an interaction: environmental heterogeneity should be associated with more floral-colour polymorphism, but that association should weaken as connectivity or gene flow increases.

See [`THEORY_2PATCH.md`](THEORY_2PATCH.md) and [`META_ANALYSIS_PLAN.md`](META_ANALYSIS_PLAN.md).

## Run

No third-party dependencies are required.

```bash
python -m unittest discover -s tests -v
python scan_phase_grid.py --n 101 --out results/phase_grid.csv
```

## Repository structure

```text
model.py                    # one-population analytical classifier and simulation
two_patch.py                # spatial heterogeneity + gene-flow extension
scan_phase_grid.py          # one-population parameter-grid generator
tests/test_model.py         # stage-1 analytical/numerical tests
tests/test_two_patch.py     # stage-2 spatial phase-boundary tests
THEORY.md                   # minimal derivation
THEORY_2PATCH.md            # exact two-patch invasion analysis
META_ANALYSIS_PLAN.md       # theory-guided comparative design
```

## Research strategy

1. **Recover theory first.** Derive the smallest set of forces and interactions that separate monomorphic and polymorphic regimes.
2. **Translate theory into measurable predictions.** Avoid unconstrained correlate hunting.
3. **Show the broad comparative pattern by meta-analysis.** Use phylogenetically informed models and explicitly address ascertainment and detection bias.
4. **Return to mechanism only where the broad pattern fails or remains non-identifiable.**

## Next theoretical targets

- unequal patch sizes and asymmetric migration;
- many-patch / network formulation and a spectral phase boundary;
- temporal environmental switching;
- finite-population origin, invasion, persistence, and loss;
- imperfect observation of polymorphism.

The long-term goal is not merely to catalogue mechanisms that can maintain floral colour variation, but to explain the transition between evolutionary regimes and determine which parts of that transition are identifiable from comparative data.
