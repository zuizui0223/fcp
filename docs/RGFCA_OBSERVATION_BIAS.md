# RGFCA observation-bias control

Status: **frozen before any RGFCA global flower-colour field is opened**.

## Current conclusion

The repeated-atlas design can reduce several important forms of opportunistic-photo bias, but it cannot make public photographs equivalent to a probability sample of all flowering plants.

The inferential target is therefore conditional: flower-colour biogeography among species and photographs that pass the prospectively frozen public-photo sampling frame. Robustness analyses test whether the resulting field is dominated by known dimensions of observation availability.

## Bias already reduced by the primary design

1. **Prolific observers.** Candidate acquisition caps each observer at two retained photographs per species before measurement.
2. **Local spatial clustering.** Raw candidate photos are selected by deterministic geographic maximin rather than by simple abundance.
3. **Popular species.** Every included species has equal field contribution, and balanced repeated species schedules prevent photo-rich taxa from dominating by archive size.
4. **Unequal geographic edge opportunity.** The colour numerator is divided by the same species-normalized geographic edge opportunity. A region with many possible sampled edges is not automatically a colour hotspot.
5. **One lucky sample.** Two hundred balanced world-map realizations expose instability to species and photograph composition instead of treating one draw as the world.
6. **Broad taxonomic/geographic dependence.** Leave-one-realm-out and major-family deletion are first-class robustness outputs.
7. **Measurement-context leakage.** Flower-colour measurement remains location blind and cannot use species name, coordinates, climate or ecological context.

These controls address dominance and geometry. They do not reconstruct data that were never observed.

## Metadata-only sampling-availability surface

Before any RGFCA colour outcome is opened, the same 18 x 9 equal-area world grid receives three count-only iNaturalist queries per cell:

1. `all_research_photo_records`: all research-grade, georeferenced species-rank photo records in the focal plant clade under the positional-accuracy rule;
2. `license_eligible_photo_records`: the same records after the fixed reusable-photo licence filter;
3. `flowering_annotated_eligible_records`: the same records after both licence and flowering-annotation filters.

The census therefore uses 162 x 3 = 486 fixed requests. It reads metadata counts only and downloads no candidate image pixels.

The three surfaces separate, imperfectly, three stages of observation availability:

`field/platform observation -> reusable-photo availability -> flowering-annotation availability`.

They are called **sampling-availability proxies**, not direct sampling probabilities. Record density combines true biological availability with observer effort.

## Frozen postmeasurement diagnostics

### Observer-unique realizations

As a sensitivity analysis, each species contributes at most one classifiable photograph per observer within a realization. A species enters this sensitivity only if at least 20 distinct observers remain. The primary analysis is not replaced if this sensitivity has smaller coverage.

### High-effort deletion

Cells in the top decile of the pre-frozen `all_research_photo_records` surface are removed, and the complete observed/null RGFCA program is recomputed on the remaining support. A global band that exists only in the most heavily recorded cells is treated as observation-sensitive.

### Sampling-availability negative controls

For jointly evaluable cells, opportunity-weighted Spearman correlations are calculated between the RGFCA field and:

- `log1p(all_research_photo_records)`;
- the licence-eligible fraction;
- the flowering-annotation fraction conditional on licence eligibility.

The same correlations are calculated for the exact 999 species-conditioned null colour fields. No alternative favourable null is introduced.

### Observer concentration audit

The final report includes distinct observers per species, the largest observer share, repeated-observer contribution within realizations, and primary-versus-observer-unique field similarity.

## What remains fundamentally unresolved

The design cannot recover species or regions with nearly zero public-photo inclusion probability. It cannot directly remove preferential photographing of unusual colour morphs, colour-dependent identification/annotation errors, or all camera/illumination/post-processing effects. It also cannot decompose public-record density into pure biological abundance and pure human observation effort.

Therefore a robust positive result supports **recurring flower-colour geography that is not readily explained by the measured observation-availability structure of the public-photo frame**. It does not establish an unbiased census of global flower-colour frequencies.
