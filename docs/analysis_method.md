# Journal of Biogeography analysis protocol

## Central question

Do documented cases of intraspecific flower-colour variation differ in occupied climatic niche according to whether colour variants coexist within populations or are geographically structured among populations?

## Terminology

The umbrella term is **intraspecific flower-colour variation**.

- `within_population`: two or more discrete colour variants coexist within at least one natural population. These cases meet the strict operational definition of flower-colour polymorphism.
- `among_population`: colour variants are documented in different populations or regions, but local coexistence is not established. These cases are treated as geographically structured flower-colour variation, not automatically as within-population polymorphism.
- `mixed`: evidence supports both local coexistence and geographic structuring.
- `unclear`: available evidence is insufficient to assign spatial organization.

The classification describes the spatial organization documented in retained sources. It is not a direct measurement of gene flow, selection, physiological tolerance or morph-specific climatic niches.

## Study population and estimand

The repository contains a global literature-derived candidate sample, not a random sample of angiosperms. The primary comparative estimand is:

> the association between species-level occupied climatic niche metrics and the documented spatial organization of intraspecific flower-colour variation among evidence-supported cases.

Results must not be described as worldwide prevalence estimates or universal laws for flowering plants.

## Analysis chronology and confirmatory status

The manuscript must distinguish between hypotheses specified before inspecting a given result and priorities adopted after exploratory analysis.

- The broad question comparing within-population coexistence with geographic differentiation was theory-led.
- Multiple occupied-climate metrics were evaluated in the comparative pipeline.
- Moisture niche breadth became the focal metric after the comparative signal was identified and should therefore be described as the **primary reported association**, not as a preregistered confirmatory endpoint unless a dated record demonstrates otherwise.
- The manuscript must not present the moisture direction as an a priori prediction without documentary evidence.
- All alternative niche metrics remain part of the reported multiplicity context and must be available in the supplement.

This wording prevents a post hoc focal result from being misrepresented as prospectively specified.

## Evidence hierarchy and classification freeze

### Baseline-unambiguous set

The primary inferential set retains species whose spatial classification is directly supported by unambiguous baseline evidence.

A species enters this set only when all of the following are true:

1. at least one retained source explicitly supports the assigned spatial category;
2. the supporting passage can be linked to the species and classification in the source-level audit;
3. the source does not simultaneously require a conflicting category unless the species is coded `mixed`;
4. no climatic result is used when assigning or retaining the classification; and
5. the inclusion list is frozen before the final model is rerun for manuscript reporting.

The repository must preserve the frozen species list, classification, supporting source and decision note. Any later change must be logged as a data correction rather than silently replacing the inferential set.

### Evidence-enriched set

A broader classification incorporates additional retained literature. It tests whether the result generalizes when less direct or less complete evidence is admitted. Attenuation is reported as evidence-quality sensitivity rather than hidden through pooled analysis.

### Manual-review set

Species lacking matching source-level justification remain explicitly flagged for manual review and are excluded from confirmed classifications.

## Climatic niche metrics

Climatic niche metrics are calculated from cleaned, coordinate-bearing GBIF occurrences and climate-cell occupancy. They represent species-level realised occupied climate, not fundamental tolerance and not climate assigned to individual colour morphs.

The focal reported result concerns moisture niche breadth. Multivariate climatic space, temperature breadth and other evaluated metrics must be reported as secondary or exploratory analyses rather than omitted.

## Model specification

The primary model uses a binary response:

- `1 = among_population`
- `0 = within_population`

The focal predictor is standardized moisture niche breadth. Any occurrence-effort or literature-effort covariates included in the fitted model must be listed explicitly in the manuscript and model table.

For every primary and sensitivity model, the repository and manuscript must report:

- exact model formula;
- estimation package and version;
- number of species;
- number of represented families;
- class counts;
- coefficient on the log-odds scale;
- odds ratio;
- 95% confidence interval and how it was computed;
- permutation procedure and number of valid permutations;
- convergence or separation diagnostics; and
- missing-data and exclusion rules.

The phrase `family-adjusted` must not be used by itself. The manuscript must distinguish clearly among family fixed effects, random effects and family-clustered standard errors. The current approach uses family-clustered uncertainty plus leave-one-family-out sensitivity unless the implementation is changed.

## Inference hierarchy

1. Baseline-unambiguous association estimate.
2. Confidence interval for the focal coefficient.
3. Permutation inference for the focal moisture-breadth coefficient.
4. Leave-one-family-out sensitivity to assess concentration in represented families.
5. Evidence-enriched generalization analysis.
6. GBIF sampling-effort and observed-distribution checks.
7. Fragmentation and generic environmental-turnover analyses as secondary negative controls.

Models that do not meet the pre-specified minimum species count are reported as `not_estimable` with zero valid permutations.

## Multiplicity and language

Because multiple niche metrics were evaluated, the manuscript must not define discovery solely by whether one p-value crosses 0.05. The focal result must be presented with its effect size, confidence interval, permutation support, family sensitivity and evidence-set attenuation.

Use these formulations:

- `we detected an association in the strict evidence set`;
- `the estimate attenuated in the broader evidence set`;
- `the examined coarse metrics did not account for the association`.

Avoid these formulations:

- `moisture breadth determines spatial organization`;
- `the null analysis proved no fragmentation effect`;
- `the primary hypothesis was confirmed` unless prospective specification is documented.

## Taxonomic non-independence

Family-clustered uncertainty and leave-one-family-out analyses reduce sensitivity to single-family dominance but do not constitute a complete phylogenetic comparative analysis. The current evidence supports a cross-species association robust to represented-family deletion, not a fully phylogenetically independent evolutionary relationship.

## Interpretation rules

- Use **association**, not causation.
- Do not state that narrow moisture breadth causes geographic colour differentiation.
- Do not infer local adaptation without morph-labelled locality evidence.
- Do not call geographically separated variants within-population FCP.
- Do not equate GBIF point-cloud components with biological populations or barriers.
- Do not treat non-significance as proof of absence.
- Report full-classification attenuation alongside the strict-evidence estimate.
- Report all evaluated niche metrics in the supplement.

## Next decisive evidence

The next mechanistic step is to connect documented colour states or named populations to localities and test whether morph-labelled sites differ environmentally. Further clustering of unlabelled GBIF occurrences is not an adequate substitute. Until then, the paper remains a comparative biogeographic association study.
