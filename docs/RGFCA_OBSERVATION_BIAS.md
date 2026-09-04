# RGFCA observation-bias control

Status: **frozen before any RGFCA global flower-colour field is opened; target-group and platform-wide metadata-only availability censuses completed**.

## Current conclusion

The repeated-atlas design can substantially reduce several important forms of opportunistic-photo dominance, but it cannot make public photographs equivalent to a probability sample of all flowering plants.

The empirical availability audit shows why the correction architecture is necessary. Across the same 162 equal-area cells, target-group research-grade photo records have a Gini coefficient of **0.885**. The 17 cells in the pre-frozen top observation-effort decile contain **84.8%** of all target-group photo records. The corresponding shares are **84.2%** after the reusable-photo licence filter and **83.3%** after the flowering-annotation filter. Thus an all-at-once analysis of the public archive would be strongly exposed to geographic observation concentration.

A second, taxon-unrestricted negative-control census confirms that generic platform activity is a major component of that concentration. Across the same 162 cells there are **154,803,786** all-taxa research-grade photo records versus **53,270,601** target-group records. Target-group records therefore comprise **34.4%** of the all-taxa count in aggregate. Platform-wide counts are themselves highly concentrated (Gini **0.871**) and are strongly concordant with target-group counts (`Spearman rho=0.9844`, `Pearson r=0.9814` on log1p counts). A descriptive equal-area-cell regression of `log1p(target-group records)` on `log1p(all-taxa platform records)` has slope **1.0716** and `R^2=0.9632`. This is descriptive rather than causal or spatially corrected, but it shows that the geographic concentration of the plant-photo frame closely follows generic iNaturalist recording intensity.

The high-effort geography is also nearly the same: **15 of 17** pre-frozen target-group top-decile cells are among the 17 highest platform-wide cells (`Jaccard=0.789`). Platform-top-17 cells contain **83.8%** of target-group records, while the frozen target-group top-17 cells contain **81.0%** of all-taxa platform records. The remaining target-group residual geography cannot be interpreted as pure biology because it may contain target availability and unmeasured observation processes.

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

## Completed metadata-only sampling-availability surfaces

### Target-group three-stage availability census

The same 18 x 9 equal-area world grid received three count-only iNaturalist queries per cell, for **486/486 successful requests and zero errors**:

1. `all_research_photo_records`: all research-grade, georeferenced species-rank photo records in the focal plant clade under the positional-accuracy rule;
2. `license_eligible_photo_records`: the same records after the fixed reusable-photo licence filter;
3. `flowering_annotated_eligible_records`: the same records after both licence and flowering-annotation filters.

No candidate image pixel or flower-colour value was opened.

Observed totals and concentration:

- all target-group photo records: 53,270,601; Gini 0.885; top-effort 17 cells contain 84.8%;
- licence eligible: 40,554,393; Gini 0.882; top-effort cells contain 84.2%;
- flowering annotated and licence eligible: 4,049,721; Gini 0.885; top-effort cells contain 83.3%.

Across cells where the fractions are defined, median reusable-licence retention is **0.790**, whereas median flowering-annotation retention conditional on licence is only **0.095**. Spatial ranks remain highly concordant across filters (`rho=0.999` all vs licence; `rho=0.983` all vs flowering annotated), so the annotation filter strongly reduces volume without eliminating the underlying geographic concentration.

### Platform-wide negative-control census

A separate metadata-only census queried all research-grade georeferenced photo records without restricting the taxon to the focal plant clade. It completed **162/162 requests with zero errors**, opened no candidate image pixels, and used no flower colour.

- all-taxa platform records: 154,803,786;
- zero-count cells: 4/162;
- equal-area-cell Gini: 0.871;
- maximum single-cell share: 13.6%;
- log-count Spearman correlation with the target-group availability surface: 0.9844.

The platform-wide layer is a generic platform-activity proxy. It is not interpreted as a pure observer-effort probability because its counts also reflect the distribution and detectability of all photographed organisms.

These four pre-colour surfaces are therefore **negative-control / sampling-availability proxies**, not inverse-probability weights.

## Frozen postmeasurement robustness gate

A recurrent colour zone is not treated as robust merely because the primary species-conditioned permutation test is positive. The exact decision rules are frozen in `docs/supporting/rgfca_observation_bias_decision_contract_v1.json` before any RGFCA global flower-colour field is opened.

### Observer-unique realizations

Each species contributes at most one classifiable photograph per observer within the sensitivity realization. A species enters only if at least 20 distinct observers remain. If the already frozen global inferential gates cannot be met, the sensitivity is reported `not_evaluable`; the observer threshold is not lowered after seeing the colour result.

### High-effort deletion

The **17 cells** in the pre-frozen target-group top observation-effort decile are removed, and the complete observed plus 999-null RGFCA program is recomputed on the remaining support. The same primary recurrent-field support gate must pass. A field that exists only with the highest-observation cells is labelled `observation_sensitive_high_effort_dependence`. This sensitivity cannot rescue a primary null.

### Four sampling-availability negative controls

For jointly evaluable cells, the absolute opportunity-weighted Spearman correlation is calculated between the observed RGFCA field and each pre-frozen surface:

1. `log1p(target_group_all_research_photo_records)`;
2. target-group reusable-licence fraction;
3. target-group flowering-annotation fraction conditional on licence eligibility;
4. `log1p(all_taxa_platform_research_photo_records)`.

The identical absolute correlations are calculated for the exact **999 species-conditioned colour-null fields** with the sampling geometry and complete Monte Carlo schedule fixed. For each surface,

`p = (1 + number of null |rho| >= observed |rho|) / 1000`.

The four p-values form one pre-frozen **Holm-corrected family at family-wise alpha=0.05**. Direction is not used to select a favourable test. Any Holm-significant negative-control alignment raises an `observation_alignment_flag`. No negative-control result can create or rescue a primary colour zone.

### Overall labels

The pre-frozen decision contract permits only the following interpretations:

- `primary_not_supported`: the primary recurrent-field gate does not pass; bias sensitivities cannot rescue it;
- `primary_supported_observation_robust`: primary passes, high-effort deletion passes, observer-unique is evaluable and passes, and no Holm-significant negative-control alignment is present;
- `primary_supported_observation_sensitive`: primary passes but at least one evaluable observation-bias diagnostic flags sensitivity;
- `primary_supported_observation_robustness_incomplete`: primary passes but one or more required diagnostics are not evaluable under their frozen gates.

### Observer concentration audit

The final report includes distinct observers per species, the largest observer share, repeated-observer contribution within realizations, and primary-versus-observer-unique field similarity.

## What remains fundamentally unresolved

The design cannot recover species or regions with nearly zero public-photo inclusion probability. It cannot directly remove preferential photographing of unusual colour morphs, colour-dependent identification/annotation errors, or all camera/illumination/post-processing effects. It also cannot decompose public-record density into pure biological abundance and pure human observation effort.

Therefore even the strongest `primary_supported_observation_robust` result supports **recurring flower-colour geography that is not readily explained by the measured observation-availability structure of the public-photo frame**. It does not establish an unbiased census of global flower-colour frequencies.
