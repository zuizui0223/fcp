# Flower-colour polymorphism and climatic niche breadth

## Literature-based novelty assessment

### What is already known

1. Flower-colour polymorphism (FCP) is conventionally defined as two or more discrete flower-colour variants co-occurring within a population. Reviews emphasize pollinator selection, abiotic selection, gene flow, drift, mating system and clonal reproduction, but do not provide a global comparative test of climatic niche breadth across flowering plants (Gómez et al. 2017, doi:10.1111/plb.12575).
2. Macroecological flower-colour studies have related species-level colour traits to biotic and abiotic environments, but they analyze mean or dominant colour traits rather than verified within-population discrete polymorphism (Dalrymple et al. 2020, doi:10.1111/nph.16737).
3. Within Protea, polymorphic species span broader elevational gradients than monomorphic species, suggesting environmental heterogeneity can predict polymorphism at a genus-wide scale (Carlson et al. 2015, doi:10.1098/rspb.2015.0583).
4. Individual FCP systems show morph-environment associations and distinct climatic niches, especially in anthocyanin polymorphisms (Arista et al. 2013, doi:10.1111/1365-2745.12151; Warren & Mackenzie 2001, doi:10.1046/j.1469-8137.2001.00159.x; Jiménez-López et al. 2024, doi:10.1186/s12870-024-05481-y).
5. Across animal clades, colour-polymorphic species have sometimes shown broader climatic niches and lower extinction risk (Takahashi & Noriyuki 2019, doi:10.1098/rsbl.2019.0228), but a phylogenetically controlled review of web-building spiders found no association with niche breadth or range size (Oxford et al. 2024, Journal of Evolutionary Biology 37:1345-1358). This establishes a comparative hypothesis, not a universal law.
6. In plants more generally, climatic niche breadth is associated with geographic range size and can differ between sister species with different reproductive systems, demonstrating that occurrence-derived niche breadth is a meaningful comparative response (Sheth et al. 2014, doi:10.1111/evo.12494; Grant & Kalisz 2020, doi:10.1111/evo.13870).

### Novelty claim supported by the search

We found no global, taxonomically matched comparative analysis testing whether verified natural within-population FCP in flowering plants is associated with broader or more heterogeneous climatic niches. Existing studies are predominantly single species, single genera, regional reviews, or analyses of average interspecific flower colour.

The defensible novelty is therefore:

> The first global screening analysis of the association between verified natural flower-colour polymorphism and occupied climatic niche breadth and heterogeneity in flowering plants, with taxonomically matched candidate controls and explicit correction for occurrence and research effort.

This is a pattern-and-association claim. It is not a demonstration that FCP causes niche expansion.

## Ecological question

**Are species with verified within-population flower-colour polymorphism associated with broader and more environmentally heterogeneous occupied climatic niches than taxonomically matched species without documented FCP?**

Because absence of evidence is not evidence of monomorphism, controls remain `species outside the verified candidate set`, not asserted monomorphic species.

## Competing hypotheses

### H1. Polymorphism–niche breadth association

Verified FCP species occupy broader climatic niches than matched controls.

Predictions:

- greater robust multivariate climatic dispersion;
- larger occupied PCA-space area;
- larger temperature and precipitation breadth;
- association persists after adjusting for GBIF record count and taxonomic matching.

### H2. Environmental heterogeneity association

Verified FCP species occur across more heterogeneous climates than matched controls.

Predictions:

- larger within-species climatic variance and interquantile ranges;
- strongest association for variables representing moisture, temperature extremes and seasonality;
- potentially stronger for anthocyanin-based FCP, when pigment information becomes available.

### H0. Detection and geography explanation

The apparent association is explained by observation effort, geographic range, taxonomy, or research culture.

Predictions:

- effects attenuate after matching on GBIF effort;
- results are sensitive to minimum occurrence thresholds;
- effects disappear in same-genus comparisons or leave-one-family-out analyses.

## Interpretation boundary

A positive association is compatible with two causal directions:

1. **ecological expansion:** polymorphism increases population-level response diversity and permits occupation of more climates;
2. **heterogeneity-generated polymorphism:** exposure to diverse climates produces spatially varying selection and maintains or generates colour variants.

Cross-sectional occurrence data cannot distinguish these directions. The paper should report the macroecological association and present these as competing explanations.

## Primary analysis design

### Units

- focal: species with high-confidence natural FCP;
- controls: taxonomically matched species outside the candidate set;
- strict sensitivity: high-confidence `within_population` FCP only;
- secondary comparison: `among_population` colour differentiation analyzed separately.

### Climate source

Use WorldClim 2.1 bioclimatic variables at 10 arc-minute resolution for the first reproducible CI analysis. It provides 19 standard temperature and precipitation variables for 1970-2000. A CHELSA 2.1 sensitivity analysis can follow because CHELSA provides terrain-informed approximately 1-km climate surfaces for 1981-2010.

### Occurrence cleaning

- accepted species name where possible;
- coordinate-bearing present records;
- remove invalid coordinates and obvious zeros;
- deduplicate rounded coordinates;
- cap records per species using reproducible spatial thinning or deterministic sampling;
- require at least 20 unique occupied climate cells for the primary analysis;
- sensitivity thresholds at 10, 30 and 50 cells.

### Niche metrics

Primary:

- robust multivariate dispersion in globally standardized BIO space;
- mean Euclidean distance to species centroid in global climate PCA (PC1-PC3);
- 95% range on PC1 and PC2;
- convex-hull area in PC1-PC2 where sample size permits.

Secondary:

- temperature breadth: robust ranges for BIO1, BIO5, BIO6 and BIO7;
- moisture breadth: robust ranges for BIO12, BIO14, BIO15 and BIO17;
- climatic heterogeneity: mean standardized within-species SD across selected variables.

### Models

1. Conditional logistic regression for verified FCP versus taxonomically matched controls.
2. Predictors entered separately first to avoid circularity and collinearity:
   - multivariate niche breadth;
   - climatic heterogeneity;
   - geographic range;
   - log GBIF occurrence count.
3. Combined model only after variance-inflation and correlation checks.
4. Same-genus-only, effort-balanced, minimum-record, leave-one-family-out and within-FCP-only sensitivities.
5. Report effect sizes and uncertainty; do not convert unknown controls into confirmed monomorphic species.

## Decision rule

The niche hypothesis becomes the paper's main macroecological result only if the direction is stable across taxonomic and effort-balanced sensitivities. If the signal is weak or absent, the null result remains informative because analogous colour-polymorphism studies disagree across animal clades, and no global plant test currently exists.
