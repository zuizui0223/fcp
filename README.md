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

For two equal patches with opposite local selection contrasts `+h` and `-h`, coupled by symmetric migration `m`, the exact condition is

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

This yields the comparative prediction that environmental heterogeneity should be associated with more floral-colour polymorphism, but that association should weaken as effective gene flow or connectivity increases.

See [`THEORY_2PATCH.md`](THEORY_2PATCH.md).

## Stage 3 — arbitrary landscape networks

For `n` patches connected by an undirected weighted migration network with graph Laplacian `L`, let `h` be the vector of patch-specific directional selection deviations. Define

\[
H_A=\lambda_{\max}[\operatorname{diag}(h)-mL],
\]

\[
H_B=\lambda_{\max}[\operatorname{diag}(-h)-mL].
\]

Compress these into

\[
H_{\mathrm{bal}}=\frac{H_A+H_B}{2},
\qquad
D_{\mathrm{asym}}=\frac{H_A-H_B}{2}.
\]

Then protected polymorphism exists exactly when

\[
\boxed{
 b+H_{\mathrm{bal}}
 >
 \left|d+D_{\mathrm{asym}}\right|
}.
\]

This is the current general phase boundary.

The result shows that environmental variance alone is not sufficient: the same set of selective environments can have different evolutionary effects depending on how they are arranged on the connectivity network.

See [`THEORY_MULTIPATCH.md`](THEORY_MULTIPATCH.md) and [`META_ANALYSIS_PLAN.md`](META_ANALYSIS_PLAN.md).

## Run

No third-party dependencies are required.

```bash
python -m unittest discover -s tests -v
python scan_phase_grid.py --n 101 --out results/phase_grid.csv
```

## Repository structure

```text
model.py                    # one-population analytical classifier and simulation
two_patch.py                # exact two-patch spatial extension
multi_patch.py              # arbitrary network spectral phase boundary
scan_phase_grid.py          # one-population parameter-grid generator
tests/test_model.py         # stage-1 analytical/numerical tests
tests/test_two_patch.py     # stage-2 spatial phase-boundary tests
tests/test_multi_patch.py   # stage-3 spectral-reduction tests
THEORY.md                   # minimal derivation
THEORY_2PATCH.md            # exact two-patch invasion analysis
THEORY_MULTIPATCH.md        # arbitrary landscape network derivation
META_ANALYSIS_PLAN.md       # theory-guided comparative design
```

## Research strategy

1. **Recover theory first.** Derive the smallest set of forces and interactions that separate monomorphic and polymorphic regimes.
2. **Translate theory into measurable predictions.** Avoid unconstrained correlate hunting.
3. **Show the broad comparative pattern by meta-analysis.** Use phylogenetically informed models and explicitly address ascertainment and detection bias.
4. **Return to mechanism only where the broad pattern fails or remains non-identifiable.**

## Next theoretical targets

- temporal environmental switching and a space-time phase boundary;
- unequal patch sizes and asymmetric migration;
- finite-population origin, invasion, persistence, and loss;
- imperfect observation of polymorphism.

The long-term goal is not merely to catalogue mechanisms that can maintain floral colour variation, but to explain the transition between evolutionary regimes and determine which parts of that transition are identifiable from comparative data.
