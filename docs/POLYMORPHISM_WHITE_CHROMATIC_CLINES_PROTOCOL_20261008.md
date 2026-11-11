# FCP: directional white/chromatic clines and local coexistence — exploratory protocol (2026-10-08)

## Biological questions

The current New Phytologist main paper has strong evidence for geographically *partitioned* species-wide visible flower-colour variation. That does not establish a signed elevation/latitude cline or adaptive maintenance. This separate post-outcome analysis asks:

1. Within the same species, do photographed non-white flowers occur at higher absolute latitude and/or higher elevation than photographed white flowers?
2. Are white + coloured visible states common in the measured high-depth species, as opposed to two different *non-white* hues, under identical minimum-photo rules?
3. When a species has both white and coloured flowers, do the two states actually appear together within 50 km, or do they segregate among localities?

The outcome categories are deliberately **visible photographic states**; they are not anthocyanin measures. The same pipeline classifies yellow/orange as coloured, but that state may represent carotenoid rather than anthocyanin pigmentation. A secondary red/pink/blue/purple vs white comparison avoids this obvious mechanistic conflation without claiming biochemical specificity.

## Frozen inputs and eligibility

Reopen **only** the immutable 100-photo-per-species discovery, validation and prospectively sampled third-cohort measurement tables at their published SHAs. Data are used for a new exploratory outcome; nothing alters previously frozen H1/H2 or ecological claims. Altitude is sampled from WorldClim v2.1 10-arc-minute elevation, source ZIP SHA256:

6a9f4e9a37f289d594ed1567b63be344a4dacf247f0b69756e81df031087e3bd

The altitude hash was recovered independently from the successful download in the 2026-10-07 full-gradient Actions logs; it is not estimated from GPS photo metadata. Elevation is too coarse to represent individual plant microhabitats.

Require 40 classifiable photos per species for each cohort, then at least 5 photos per side of the white/nonwhite contrast. Require at least one degree of sampled *absolute* latitude span for latitude tests, or at least 100 metres of sampled WorldClim elevation span for altitude tests. No species, observation or colour-state substitutions are allowed after seeing outcomes. Each eligible species has equal weight.

## Estimands, nulls and uncertainty

For predictor x (|latitude| in degrees, elevation in metres), define for species i:

delta_i(x) = [mean(x | nonwhite) - mean(x | white)] / sd_i(x).

Primary contrast is mean(delta_i) across species in each of three disjoint cohorts, with a 1999-species bootstrap 95% interval. Fixed-count *within-species vertex permutations* (499) preserve the sampled coordinates, elevation, classifiable rows and number of each binary state; plus-one P values address the one-sided direction "coloured at higher |latitude|/altitude". Report two-sided P values and Holm adjustment across the two directional primary axes per cohort, whether the signs are favourable or not. The genus-balanced mean checks domination by speciose genera, without claiming full phylogenetic independence.

Secondary analyses: (i) repeat visible contrasts using red/pink/blue/purple only versus white, excluding yellow/orange; (ii) regress elevation on absolute latitude *within species* and describe the residual elevation contrast (not a separately powered predeclared test). Do not select only species with positive contrasts.

The primary high-latitude prediction is **not presumed true**: high UV often occurs at low latitudes, whereas low temperature, pollinator turnover, developmental plasticity and historical distribution can move in other directions. Any one-cohort sign is insufficient for a generality claim.

## Coexistence versus spatial differentiation

A species is described as **image-sampled white + colour** when >=5 classifiable photos show white and >=5 show nonwhite. It is **image-sampled nonwhite hue-variable** when >=5 photos occur in each of at least two nonwhite coarse hue classes. Report overlap rather than declaring the states mutually exclusive. These counts are from selected photo species and are NOT global incidence or actual frequencies in natural populations.

Using within-species photographs <=50 km apart, report observed local white-colour discordant pair counts divided by the expected number under exact species-wide fixed-morph counts: p(W,C)=2 n_W n_C/[n(n-1)]. Compute a parallel nonwhite-hue ratio and record missing/low-opportunity species. This is local *photo neighbourhood* co-occurrence, not proof of polymorphism in one mating population. Existing distributed polymorphism is the context, not a new adaptive selection test.

## Mechanisms and discriminating biological predictions

- **Abiotic fitness benefit of pigmentation** predicts repeatable, within-species visible colour enrichment in genuinely stressful UV/cold/dry microhabitats *after photo and spatial confounding control*, and eventually a fitness-by-morph interaction with the same sign.
- **Pollinator-mediated choice** predicts morph-specific visitation or pollen transfer conditional on pollinator guild and local frequencies; latitude alone is a weak surrogate.
- **Spatially varying selection + gene flow** can produce predominantly locally monomorphic patches with a few white+colour transition-zone populations. Mere clustering is also possible under dispersal and drift.
- **Balancing/frequency-dependent selection** requires persistent local coexistence and direct morph-specific fitness varying with frequency or temporal state. Photographs cannot establish this.
- **Developmental/genetic accessibility and pleiotropy** may explain white+colour variants being observed repeatedly without implying that white morphs are fitter. Published petal-limited anthocyanin-loss (PAL) polymorphisms and whole-plant-loss (WAL) variants must be distinguished; the current repository's literature-derived PAL/WAL ledger is not a population-unbiased cross-species frequency estimate.

A biologically strong final argument would jointly explain the *origin of white states*, the *spatial placement of morph frequencies* and their *persistence within mating populations*. None of those three steps follows solely from a macroecological latitudinal effect.

## Matched-composition robustness of the white/pigmented incidence comparison

The initial contrast between white+any-nonwhite and two different nonwhite hues is asymmetric: the former can pool three nonwhite colour categories to exceed the minimum. We additionally require the same >=5 *per category*: white >=5 and at least one individual nonwhite hue >=5 versus two distinct nonwhite hues each >=5. These counts may overlap. A further descriptive category asks whether white is among the two most abundant colour labels, with the second label >=5. These are photo-based and have no wild-population prevalence interpretation.

Within the **same species with estimable opportunities for both discordance types**, compare the observed/expected white-vs-nonwhite local pair ratio against the observed/expected different-nonwhite-hue local pair ratio. Paired mean difference, bootstrap uncertainty and a retrospective paired signflip test are reported, without promoting the test as independent or prospective. Genus-blocked signed nulls are similarly added as secondary sensitivities for the signed latitude/elevation contrasts. None of these overcome exposure or phylogenetic confounding.

## Hard inferential limits

The image-based coarse white classifier is **known to be exposure-coupled** in the third cohort (source: current H2 highlight-audit result), and this workflow does not have a complete direct highlight/reflectance control in all three cohorts. The sampled species are not globally representative; neither the WorldClim 10-minute elevation nor recorded mean climate equals the plant's immediate environment. Colour/photo associations cannot determine genotype, UV absorbance, pigment chemistry, selection coefficient, fitness or fitness-mediated local adaptation.

The existing FCP manuscript is not rewritten; output is a separate exploratory result and per-species receipt. A favourable association can motivate an **independent species/population level** test, not validate adaptation retroactively.
