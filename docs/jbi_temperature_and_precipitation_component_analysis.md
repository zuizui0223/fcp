# Temperature decomposition and ecological interpretation of the BIO12 result

## Question

The previous temperature-breadth metric was the arithmetic mean of BIO1, BIO5, BIO6 and BIO7 occupied-climate breadths. The precipitation analysis showed that its analogous moisture composite was almost entirely driven by BIO12. We therefore analysed the four temperature components separately and compared two multiple-testing families:

1. four temperature components only;
2. all eight temperature and precipitation components.

In this context, a multiple-testing `family` means a set of related ecological hypotheses, not a plant taxonomic family.

## Methods

For BIO1, BIO5, BIO6 and BIO7, breadth was the difference between the 95th and 5th percentiles across occupied climate cells. Each variable was standardised separately before fitting:

`among ~ component_z + effort_z`

The models used the same 34 baseline-unambiguous species, occurrence-effort adjustment, family-clustered sandwich covariance, 9,999 common label permutations and leave-one-plant-family-out refits as the precipitation-component analysis. The primary and paginated occurrence samples were analysed separately. Open Tree and fixed dated-megaphylogeny sensitivity models were also fitted.

## Temperature results

### Primary occurrence sample

No temperature component showed statistical support:

- BIO1 annual mean temperature breadth: OR 1.248, 95% CI 0.509–3.058; Wald p 0.629; permutation p 0.594.
- BIO5 warmest-month maximum-temperature breadth: OR 1.008, 95% CI 0.467–2.177; Wald p 0.984; permutation p 0.983.
- BIO6 coldest-month minimum-temperature breadth: OR 0.782, 95% CI 0.404–1.511; Wald p 0.464; permutation p 0.535.
- BIO7 annual temperature-range breadth: OR 0.732, 95% CI 0.333–1.606; Wald p 0.436; permutation p 0.424.

All four Holm-adjusted Wald and permutation p-values were 1.0. The legacy arithmetic-mean temperature breadth was also unsupported (OR 0.862, 95% CI 0.410–1.809).

### Paginated quality-filtered occurrence sample

BIO6 produced the strongest negative estimate, but remained unresolved:

- BIO1: OR 0.809, 95% CI 0.435–1.507; Wald p 0.505; permutation p 0.596.
- BIO5: OR 0.896, 95% CI 0.380–2.113; Wald p 0.801; permutation p 0.780.
- BIO6: OR 0.551, 95% CI 0.289–1.050; Wald p 0.0698; permutation p 0.187.
- BIO7: OR 0.773, 95% CI 0.358–1.671; Wald p 0.513; permutation p 0.558.

For BIO6, four-temperature Holm p-values were 0.279 (Wald) and 0.749 (permutation). The legacy temperature composite was unsupported (OR 0.671, 95% CI 0.322–1.398).

### Phylogenetic sensitivity

All temperature-component confidence intervals included one. The paginated BIO6 estimate remained negative under both phylogenetic treatments:

- Open Tree: OR 0.580, 95% CI 0.227–1.484, p 0.256;
- dated megaphylogeny: OR 0.605–0.619, CI envelope 0.242–1.522, p 0.283–0.296.

The temperature decomposition therefore provides no evidence that occupied thermal breadth generally distinguishes within- from among-population flower-colour organization.

## Multiple-testing family interpretation

### Four-component domain families

A four-temperature correction controls the probability of at least one false-positive claim among BIO1, BIO5, BIO6 and BIO7. Ecologically, it corresponds to the question:

> Is any major dimension of occupied temperature breadth associated with flower-colour spatial organization?

A separate four-precipitation correction corresponds to:

> Is any major dimension of occupied precipitation breadth associated with flower-colour spatial organization?

This separation is defensible when temperature and precipitation were specified as distinct ecological hypothesis domains before inspecting the results.

### Eight-component climate family

An eight-component correction treats BIO1, BIO5, BIO6, BIO7, BIO12, BIO14, BIO15 and BIO17 as one exploratory climate screen. Ecologically, it corresponds to the broader question:

> Is any tested occupied-climate breadth component associated with flower-colour spatial organization?

Under this stricter family, paginated BIO12 remained supported:

- Holm-adjusted Wald p = 0.0231;
- Holm-adjusted permutation p = 0.0440.

Primary BIO12 did not remain supported under the eight-component family:

- Holm-adjusted Wald p = 0.231;
- Holm-adjusted permutation p = 0.160.

Because the component tests are correlated, Holm correction is conservative but valid. The correction does not measure ecological independence; it protects against selecting the smallest p-value from a declared set of related tests.

## Ecological interpretation of broader BIO12 breadth and lower among-population odds

BIO12 breadth is the 95th–5th percentile range of annual precipitation across the species' sampled occupied climate cells. A negative coefficient means that species occupying a wider range of annual-precipitation regimes were more often classified as having local coexistence of colour variants, whereas species occupying a narrower annual-precipitation range were more often classified as geographically structured.

The most defensible interpretation is a spatial-organization hypothesis, not a direct moisture-selection mechanism:

1. **Broad precipitation occupancy may be associated with persistence of multiple colour variants across shared local contexts.** Species spanning many annual-precipitation regimes may possess demographic connectivity, ecological generality or repeated local conditions that allow colour variants to coexist within populations.
2. **Narrow precipitation occupancy may make regional separation more visible or persistent.** In climatically restricted species, colour variants may be confined to different regions because of dispersal limitation, demographic history, range fragmentation or correlated environmental differences.
3. **BIO12 may be a proxy for range structure rather than the causal agent.** Annual-precipitation breadth can covary with geographic extent, biome transitions, elevation, soil, pollinator communities, population connectivity and research coverage.
4. **The response is literature-derived.** Local coexistence may be under-recorded in studies focused on regional differentiation, so the association can partly reflect how species are studied and described.

The result does not show that precipitation breadth causes within-population polymorphism, that narrow-ranged species are locally adapted by flower colour, or that colour morphs themselves occupy different precipitation niches. GBIF records are not labelled by morph, and the predictor is species-level realised occupied climate.

## Recommended manuscript framing

Preferred wording:

> Species with broader occupied annual-precipitation ranges were less likely to be classified as showing geographically structured rather than within-population flower-colour variation. This association was concentrated in BIO12 and was not mirrored by temperature-breadth components. It remained supported in the paginated non-phylogenetic analysis after correction across eight climate components, but phylogenetic confidence intervals included one.

Avoid:

- broad precipitation tolerance maintains flower-colour polymorphism;
- annual precipitation causes the transition between within- and among-population variation;
- within-population species are precipitation generalists;
- among-population species are locally adapted to narrow precipitation niches.

## Reproducibility

- temperature workflow run: `30668225933`;
- artifact: `8807825673`;
- artifact digest: `sha256:a8d564229d06f38048237aadedc386c88dd2494e93bc9b67c3c89d0117bdbe76`;
- temperature non-phylogenetic table: `docs/supporting/jbi_table_s22_temperature_component_models.csv`;
- temperature phylogenetic table: `docs/supporting/jbi_table_s23_temperature_component_phylogenetic_summary.csv`;
- eight-component adjustment table: `docs/supporting/jbi_table_s24_climate_component_familywise_adjustment.csv`.
