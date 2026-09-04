# RGFCA observation-bias control

Status: **frozen before any RGFCA global flower-colour field is opened; metadata-only availability census completed**.

## Current conclusion

The repeated-atlas design can substantially reduce several important forms of opportunistic-photo dominance, but it cannot make public photographs equivalent to a probability sample of all flowering plants.

The empirical availability audit shows why the correction architecture is necessary. Across the same 162 equal-area cells, all research-grade photo records have a Gini coefficient of **0.885**. The 17 cells in the pre-frozen top observation-effort decile contain **84.8%** of all target-group photo records. The corresponding shares are **84.2%** after the reusable-photo licence filter and **83.3%** after the flowering-annotation filter. Thus an all-at-once analysis of the public archive would be strongly exposed to geographic observation concentration.

The inferential target is therefore conditional: flower-colour biogeography among species and photographs that pass the prospectively frozen public-photo sampling frame. Robustness analyses test whether the resulting field is dominated by measured dimensions of observation availability.

## Bias already reduced by the primary design

1. **Prolific observers.** Candidate acquisition caps each observer at two retained photographs per species before measurement.
2. **Local spatial clustering.** Raw candidate photos are selected by deterministic geographic maximin rather than by simple abundance.
3. **Popular species.** Every included species has equal field contribution, and balanced repeated species schedules prevent photo-rich taxa from dominating by archive size.
4. **Unequal geographic edge opportunity.** The colour numerator is divided by the same species-normalized geographic edge opportunity. A region with many possible sampled edges is not automatically a colour hotspot.
5. **One lucky sample.** Two hundred balanced world-map realizations expose instability to species and photograph composition instead of treating one draw as the world.
6. **Broad taxonomic/geographic dependence.** Leave-one-realm-out and major-family deletion are first-class robustness outputs.
7. **Measurement-context leakage.** Flower-colour measurement remains location blind and cannot use species name, coordinates, climate or ecological context.

These controls address dominance and geometry. They do not reconstruct data that were never observed.

## Completed metadata-only sampling-availability surface

The same 18 x 9 equal-area world grid received three count-only iNaturalist queries per cell, for **486/486 successful requests and zero errors**:

1. `all_research_photo_records`: all research-grade, georeferenced species-rank photo records in the focal plant clade under the positional-accuracy rule;
2. `license_eligible_photo_records`: the same records after the fixed reusable-photo licence filter;
3. `flowering_annotated_eligible_records`: the same records after both licence and flowering-annotation filters.

No candidate image pixel or flower-colour value was opened.

Observed totals and concentration:

- all photo records: 53,270,601; Gini 0.885; top-effort 17 cells contain 84.8%;
- licence eligible: 40,554,393; Gini 0.882; top-effort cells contain 84.2%;
- flowering annotated and licence eligible: 4,049,721; Gini 0.885; top-effort cells contain 83.3%.

Across cells where the fractions are defined, median reusable-licence retention is **0.790**, whereas median flowering-annotation retention conditional on licence is only **0.095**. Spatial ranks remain highly concordant across filters (`rho=0.999` all vs licence; `rho=0.983` all vs flowering annotated), so the annotation filter strongly reduces volume without eliminating the underlying geographic concentration.

These are **sampling-availability proxies**, not direct sampling probabilities. Record density combines biological availability with observer/platform effort.

## Frozen postmeasurement robustness gate

A recurrent colour zone is not treated as robust merely because the primary species-conditioned permutation test is positive. The observation-bias robustness package requires reporting the following pre-frozen checks.

### Observer-unique realizations

Each species contributes at most one classifiable photograph per observer within the sensitivity realization. A species enters only if at least 20 distinct observers remain. The primary analysis is not replaced if this sensitivity has smaller coverage.

### High-effort deletion

The 17 cells in the pre-frozen top observation-effort decile are removed, and the complete observed/null RGFCA program is recomputed on the remaining support. A band that exists only in those highest-observation cells is classified as observation-sensitive.

### Sampling-availability negative controls

For jointly evaluable cells, opportunity-weighted Spearman correlations are calculated between the RGFCA field and:

- `log1p(all_research_photo_records)`;
- reusable-licence fraction;
- flowering-annotation fraction conditional on licence eligibility.

The same correlations are calculated for the exact 999 species-conditioned null colour fields. No alternative favourable null is introduced. Strong alignment of an observed colour field with the availability surfaces, without separation from the matched null, limits the biogeographic interpretation.

### Observer concentration audit

The final report includes distinct observers per species, the largest observer share, repeated-observer contribution within realizations, and primary-versus-observer-unique field similarity.

## What remains fundamentally unresolved

The design cannot recover species or regions with nearly zero public-photo inclusion probability. It cannot directly remove preferential photographing of unusual colour morphs, colour-dependent identification/annotation errors, or all camera/illumination/post-processing effects. It also cannot decompose public-record density into pure biological abundance and pure human observation effort.

Therefore a robust positive result supports **recurring flower-colour geography that is not readily explained by the measured observation-availability structure of the public-photo frame**. It does not establish an unbiased census of global flower-colour frequencies.
