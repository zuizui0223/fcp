# FCP: heat–water interaction and molecular-mechanism decision protocol (2026-10-19)

## Scope and status
Continuation of PR #248. This document fixes the next inference and audit criteria; it **does not** claim new model runs, RNA-seq differential expression, or measured pigment costs. The 447,500 measured photographs are the source frame; previous climate analysis used 1,574 species and 86,114 eligible photographs. Previously reported weak BIO5-white association and failure of species-specific spatial transfer are retained.

## Primary ecological estimands
1. Compare within-species white vs chromatic photographs using temperature in the warmest month, **precipitation and radiation from that same warmest month** (not BIO14, which may refer to a different month). Center and scale using training-only statistics. Document climate normal years and units.
2. Fit temperature × warmest-month precipitation interaction, with corresponding main effects and radiation, species intercepts and partial pooling of slopes. Report uncertainty for interaction and simple slopes at prespecified precipitation quantiles. A nonsignificant interaction is not evidence of zero cost or zero biological effect.
3. Distinguish geographic and photographic association from fitness/selection. Add sensitivity to observer, observation date, cultivation status, regional environmental availability, geographic cell, and color classification confidence where available.
4. Use complete held-out 1-degree spatial cells and 50-km conspecific exclusion; compare intercept-only, shared slopes, partial pooling and separate species models on *identical* held-out rows. Compute log loss, Brier, calibration and uncertainty clustered by genus and spatial block. Tune regularization inside training folds only. Do not claim transfer from in-sample likelihood.
5. Report both equal-species and photograph-weighted summaries; stratify by availability of white/chromatic records. Do not call geographic color variation within-population genetic polymorphism without population-level evidence.

## Mechanism: three mutually nonexclusive explanations
- **Genetic morph differentiation**: segregating alleles or stable regulatory differences; test common-garden morph persistence and genotype-by-environment.
- **Plastic pigmentation**: the same genotype changes petal pigmentation under light, temperature or drought; test paired shade/light and recovery.
- **Spatial sampling and ancestry**: cultivation, observation bias, genetic drift, gene flow and habitat availability create photographed morph associations; compare matched cells and genomic structure.

## Molecular analysis decision gates
1. Freeze public accession/sample sheet: organism, genotype, morph before treatment, treatment, organ, time, biological replicate or pool, library and accession.
2. For the reported shade-induced white-flower Trifolium dataset, check raw count provenance, gene annotations, biological replicate independence and confounding of shade with morph before testing. If raw counts are not comparable or replicate identity is unclear, restrict to descriptive plots.
3. Fit a factorial expression model (morph + treatment + morph:treatment) **only if** independent biological replication and crossed design exist. For shade-induced color change within one genotype, analyze treatment, not a fabricated genetically independent morph contrast.
4. Distinguish early flavonoid genes (CHS/CHI/F3H), branch enzymes (DFR/ANS/UFGT), MYB–bHLH–WD40 regulators, transport, degradation and tissue specificity. Orthology and annotation confidence must be retained.
5. Cost-benefit requires actual resource or fitness outcomes: anthocyanin concentration/flux, growth, seed set, herbivory and pollinator visitation under controlled environments. Transcript fold changes are **not** metabolic cost measurements.

## Claim gates
- Recurrent white–chromatic spatial axis: supported descriptively in prior 1,761-species 50-km analysis; 83.81% is an effect partition, not fraction of species.
- Shared weak temperature association: supported as conditional photographic association, not adaptive mechanism.
- Heat × water interaction: preliminary continuation reports non-support; independently verify source numeric outputs before quoting effect sizes.
- Partial pooling transfer: preliminary continuation reports more stable predictions, but improvement over intercept baseline unresolved; do not claim validated predictive gain until paired loss intervals are checked.
- Molecular shade response: dataset candidate identified, not a validated differential-expression result.
- Evolutionary trade-off: remains an experimental hypothesis.

## Next executable checkpoint
Obtain the current heat–water and partial-pooling output tables and test logs; recompute paired held-out log-loss contrasts and intervals from stored predictions. Then verify the six-sample shade dataset sample identities and expression counts. Freeze numeric conclusions only after independent reconstruction.
