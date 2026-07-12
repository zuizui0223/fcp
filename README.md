# FCP — Floral Colour Polymorphism phase theory

A minimal mathematical framework for the question:

> Why are some plant species polymorphic in floral colour while others are monomorphic?

The project starts from a phase-boundary approach rather than from a single candidate mechanism.

## Minimal result

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

Here `b` is effective balancing strength and `d` is net directional bias. This gives a compact monomorphism–polymorphism phase boundary.

See [`THEORY.md`](THEORY.md) for the derivation.

## Run

No third-party dependencies are required.

```bash
python -m unittest discover -s tests -v
python scan_phase_grid.py --n 101 --out results/phase_grid.csv
```

The CSV contains the full phase classification over a parameter grid:

- `protected_polymorphism`
- `A_fixation`
- `B_fixation`
- `bistability`

## Repository structure

```text
model.py                # analytical phase classifier and replicator simulation
scan_phase_grid.py      # reproducible parameter-grid generator
tests/test_model.py     # analytical/numerical consistency tests
THEORY.md               # derivation and research extensions
```

## Research direction

The minimal model is intentionally a compression layer. Subsequent models can decompose effective balancing and directional terms into pollinator-mediated frequency dependence, spatial heterogeneity, temporal variation, abiotic selection, herbivory, pigment costs, migration, mutation, drift, and detection.

The key next question is whether distinct ecological mechanisms remain identifiable after they collapse onto the same effective phase boundary.
