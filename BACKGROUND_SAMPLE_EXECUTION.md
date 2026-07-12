# Background sample execution

## Fixed target

The current 30-stratum design expands to **565 species-draw slots**.  These slots are defined before looking up flower-colour polymorphism status.

This is the denominator sample for estimating quantities such as

\[
P(\text{polymorphism}\mid\text{pollination mode, display opportunity, ecology}).
\]

It is intentionally separate from the literature-enriched positive-case set, which estimates where known colour-polymorphism systems occur but cannot by itself provide prevalence.

## Workflow

1. **Freeze taxonomy.** Build a species pool for each family from one chosen taxonomic backbone and record backbone/version/date.
2. **Random draw within stratum.** Sample without replacement using a stored random seed. Known colour-polymorphism status must not enter the draw.
3. **Code pollination from evidence.** Use direct pollinator or reproductive-biology evidence where available. Do not infer a species-level pollination mode from floral syndrome alone.
4. **Code display opportunity separately.** Record whether there is a human-visible floral display and which structure carries it (petal/tepal, bract, stamen, inflorescence, reduced/no conspicuous display).
5. **Search colour variation symmetrically.** Apply the same search protocol to every sampled species, including wind-pollinated taxa.
6. **Do not equate silence with monomorphism.** `monomorphic` requires adequate negative evidence; otherwise use `unknown`.
7. **Code spatial outcome.** When multiple colours occur, distinguish `geographic_mosaic`, `mixed`, and `local_coexistence`.
8. **Record observation effort.** Capture proxies such as number of floras/monographs checked, occurrence-record count, and number of colour-focused sources.

## Primary comparisons

### A. All sampled angiosperms

Tests total association of pollination mode with recorded flower-colour polymorphism, while explicitly retaining display opportunity as a mediator/confounder rather than hiding it.

### B. Visible-display subset

Restrict to species with a comparable human-visible floral display. This asks whether biotic versus abiotic pollination differs after the opportunity to express visible flower colour is held more constant.

### C. Petal/tepal-display subset

A stricter sensitivity analysis that avoids comparing petals with catkins, spathes, coloured bracts, or other display structures.

## Outcome hierarchy

- `monomorphic_confirmed`
- `geographic_mosaic`
- `mixed`
- `local_coexistence`
- `unknown`

The primary analysis must retain `unknown` rather than silently recoding it as monomorphic.

## Pollination hierarchy

- `abiotic_wind`
- `abiotic_water`
- `ambophilous_or_mixed`
- `biotic_insect`
- `biotic_bird`
- `biotic_bat`
- `biotic_other_vertebrate`
- `biotic_generalized`
- `unknown`

Secondary guild fields may split insect pollination into bee, fly, moth, butterfly, beetle, wasp, or generalized insect systems when supported by evidence.

## Immediate deliverable

Run:

```bash
python build_background_manifest.py
```

This materializes the 565 preregistered draw slots as `data/background_sample_manifest.csv`. The next empirical step is to attach a frozen taxonomic species pool to each stratum and execute the random draws with a stored seed.
