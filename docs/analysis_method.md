# Journal of Biogeography analysis protocol

## Central question

Do documented cases of intraspecific flower-colour variation differ in occupied climatic niche according to whether colour variants coexist within populations or are geographically structured among populations?

## Terminology

The umbrella term is **intraspecific flower-colour variation**.

- `within_population`: two or more discrete colour variants coexist within at least one natural population. These cases meet the strict operational definition of flower-colour polymorphism.
- `among_population`: colour variants are documented in different populations or regions, but local coexistence is not established. These cases are treated as geographically structured flower-colour variation, not automatically as within-population polymorphism.
- `mixed`: evidence supports both local coexistence and geographic structuring.
- `unclear`: available evidence is insufficient to assign spatial organization.

The within-versus-among classification describes the spatial organization documented in retained sources. It is not a direct measurement of gene flow, selection, physiological tolerance or morph-specific climatic niches.

## Study population and estimand

The repository contains a global literature-derived candidate sample, not a random sample of angiosperms. The primary comparative estimand is therefore:

> the association between species-level occupied climatic niche metrics and the documented spatial organization of intraspecific flower-colour variation among evidence-supported cases.

Results must not be described as worldwide prevalence estimates or universal laws for flowering plants.

## Evidence hierarchy

### Baseline-unambiguous analysis

The primary inferential set retains species whose spatial classification is directly supported by unambiguous baseline evidence. This minimizes classification error and supplies the main effect estimate.

### Evidence-enriched sensitivity

A broader classification incorporates additional retained literature. It is used to test whether the result survives expansion to less direct or less complete evidence. Attenuation after enrichment is reported as evidence-quality sensitivity rather than hidden through pooled analysis.

### Manual-review set

Species lacking a matching source-level justification remain explicitly flagged for manual review and are not silently treated as confirmed classifications.

## Climatic niche response

Climatic niche metrics are calculated from cleaned, coordinate-bearing GBIF occurrences and climate-cell occupancy. They represent species-level realised occupied climate, not fundamental tolerance and not climate assigned to individual colour morphs.

The primary result is the association between moisture niche breadth and the probability that documented variation is geographically structured among populations rather than maintained within populations.

## Model hierarchy

1. Primary family-adjusted logistic model on the baseline-unambiguous evidence set.
2. Permutation inference for the primary moisture-breadth coefficient.
3. Leave-one-family-out analysis to assess concentration in any represented family.
4. Evidence-enriched sensitivity analysis.
5. GBIF sampling-effort and observed-distribution checks.
6. Fragmentation and generic environmental-turnover analyses as secondary negative controls, not as direct tests of morph-environment sorting.

Models that do not meet the pre-specified minimum species count are reported as `not_estimable` with zero valid permutations.

## Taxonomic non-independence

Family-clustered uncertainty and leave-one-family-out analyses reduce dependence on individual clades but do not constitute a complete phylogenetic comparative analysis. Manuscript claims must acknowledge this limitation unless a species-level phylogeny is added.

## Interpretation rules

- Use **association**, not causation.
- Do not state that narrow moisture breadth causes geographic colour differentiation.
- Do not infer local adaptation without morph-labelled locality evidence.
- Do not call geographically separated variants within-population FCP.
- Do not equate GBIF point-cloud components with biological populations or barriers.
- Report full-classification attenuation alongside the strict-evidence estimate.
- Report effect sizes, confidence intervals, permutation support and family sensitivity together.

## Next decisive evidence

The next mechanistic step is not further clustering of unlabelled GBIF records. It is to connect documented colour states or named populations to localities and test whether morph-labelled sites differ environmentally. Until then, the paper remains a comparative biogeographic association study.
