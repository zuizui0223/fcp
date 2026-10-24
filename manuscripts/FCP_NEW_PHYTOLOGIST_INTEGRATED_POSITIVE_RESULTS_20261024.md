# Recurrent floral colour directions, geographic assembly and climatic associations across 4,475 angiosperm species

**New Phytologist — Full Paper | Integrated positive-results manuscript | 24 October 2026**

**Authors:** [Author list and order to be confirmed]

**Affiliations:** [Institutional affiliations to be confirmed]

**Corresponding author:** [Name and contact to be confirmed]

**Running title:** Geography and climate of intraspecific floral colour

## Summary

- Intraspecific floral colour variation may be expressed as alternative colours, gradual shifts in chromatic composition, or geographic replacement of visually contrasting states. We asked whether the direction of colour variation, its geographic assembly and its climatic associations show recurring patterns across flowering plants.
- We analysed 447,500 original photographs representing 4,475 angiosperm species. A fixed, location-blind measurement pipeline yielded 227,107 classifiable photographs with nine-component colour palettes and four photographic colour categories. Species-conditional null models, observer-disjoint analyses and whole-genus resampling were used to examine colour geometry, geographic structuring and climate.
- Colour variation was disproportionately oriented along a white–chromatic axis (mean alignment 0.4789 versus matched-null median 0.4320; P = 0.005). Nearby conspecific photographs were more similar than expected from species-wide colour composition (50-km mismatch depletion +0.0161 across 3,600 species; P = 0.005); the result was reproduced using continuous nine-colour distances.
- A shared thermal-climate axis was positively associated with the photographed white fraction (mean conditional coefficient +0.00821 across 2,837 eligible species), and 85.4% of the squared cross-species mean climatic colour vector aligned with the white–chromatic direction. The mean association persisted after observer thinning, removal of any one of five species cohorts, and independent-observer photo cross-fitting.

**Synthesis:** The direction of colour variation, its spatial arrangement and its climatic orientation are distinguishable, recurrent features of intraspecific photographic floral colour diversity. Together they provide a comparative framework for studying how floral phenotypic variation is assembled across geographic and environmental gradients.

**Key words:** angiosperms; climate; colour polymorphism; flower colour; geographic differentiation; intraspecific variation; macrophenotypic ecology; white–chromatic contrast.

---

## Introduction

Flower colour is among the most conspicuous traits through which flowering plants interact with their environments. Colour may contribute to pollinator detection and preference, but the expression and persistence of pigmentation are also associated with developmental regulation, abiotic stress, biogeographic history and reproductive ecology (Wessinger & Rausher, 2012; Narbona et al., 2018; Dalrymple et al., 2020; Sapir et al., 2021). One striking feature of flower colour is that alternative appearances often occur within a single named species. Yet the ecological meaning of such variation depends not only on which colours exist, but also on **how the colours differ and where they occur**.

Most broad comparative studies quantify differences among species or among regional floras, whereas the mechanisms of flower-colour polymorphism have often been examined in selected populations of individual species. Both approaches have yielded valuable insights, but neither by itself reveals whether the same directions of within-species phenotypic variation recur across thousands of taxa, or whether those directions are organized geographically in comparable ways. Distinguishing within-species trait geometry from geographic allocation is therefore essential for connecting studies of floral variation with macroecology (Albert et al., 2010; Carmona et al., 2019; Palacio et al., 2025).

Three aspects of intraspecific colour variation can be separated conceptually. The first is **direction**: does colour vary predominantly between white and chromatic appearances, or among distinct nonwhite colours? A recurring white–chromatic direction could emerge from broadly available routes of changes in floral pigmentation or development, although visually similar outcomes need not share particular alleles (Wessinger & Rausher, 2012; Larter et al., 2018). The second aspect is **geographic assembly**: do similar colour appearances cluster within particular parts of a species' photographed range, or are alternative appearances distributed without detectable structure? This is a spatial-statistical property, rather than a direct estimate of local coexistence or selection. The third aspect is **climatic orientation**: are the recorded colours associated with thermal, precipitation or radiation gradients, and does the climatic response share the same colour-space direction that recurs across species?

These dimensions need not be interchangeable. A species can exhibit a recurrent visual contrast but only weak spatial structuring; another can show strong geographic sorting of colours without a dominant white–chromatic contrast. The same geographic arrangements could arise from historically restricted dispersal, environmentally structured development, local selection, range expansion or sampling processes. Accordingly, geographic colour structure should be measured separately from any explanation for its origin or maintenance (Schemske & Bierzychudek, 2007; Wang & Bradburd, 2014).

Digital observations make a geographically extensive comparison possible. Citizen-science photographs now sample flowers from many habitats, climates and regions, including species rarely examined for colour polymorphism. Automated colour extraction can convert these images into reproducible multicomponent traits, provided that the photographic source, classifier, observation effort and the distinction between images and biological morphs remain explicit (Di Cecco et al., 2021; Laitly et al., 2021; Luong et al., 2023). This approach is particularly suited to detecting comparative regularities that can subsequently be tested with botanical phenotyping and genetics.

Here, we assembled a source-verified atlas comprising **447,500 original photographs of 4,475 angiosperm species**, sampled in five species-disjoint acquisition cohorts. We analysed four connected propositions. First, the geometry of within-species photographic colour variation should show a repeated white–chromatic direction. Second, the geographic arrangement of conspecific colours should be nonrandom beyond each species' overall colour composition. Third, the colour vector associated with a common thermal gradient should be preferentially oriented toward the white–chromatic axis, after accounting for other broad climate dimensions. Fourth, if the processes defining a colour difference and those arranging it geographically are not equivalent, the shared geometric direction should not necessarily predict the strength of geographic or climatic sorting. We integrated earlier categorical analyses with continuous nine-colour spatial metrics, fixed shared thermal and precipitation axes, genus-level inference and observer-disjoint source-photograph checks to obtain one species-level ecological synthesis.

## Materials and Methods

### Source photographs and the common taxonomic register

The study used a single frozen observational source comprising **4,475 named angiosperm species, each with 100 originally sampled photographs** (447,500 unique image records). Five species-disjoint acquisition cohorts contributed 500 discovery, 500 validation, 499 prospective, 498 fourth-wave and 2,478 additional species. All cohorts used the same iNaturalist-based sampling source and a fixed location-blind colour measurement procedure. The complete species register, rather than only taxa meeting a later statistical threshold, was retained throughout the synthesis.

A flower-region image pipeline combined detection and segmentation with a normalized nine-colour photographic palette. The components were white, yellow, orange, red, pink, magenta, purple, blue and bronze. For categorical analyses, an originally fixed classifier assigned images to white, yellow/orange, red/pink or blue/purple. Exactly **227,107 original photos (50.75%)** yielded classifiable colour measurements. The remaining 220,393 records were retained as unclassified rather than assigned a colour. The original image and climate archives, frozen source hashes, taxonomic ledgers and source-workflow records are available in the project's reproducibility materials.

We use the term **photographic floral colour variation** for all source-derived traits. A photographed colour class is an observed appearance, not by itself a proven hereditary morph or a single biological floral organ. Species with insufficient classified photographs, geographic pairs or sampled climate contrasts are represented by missing estimates in the relevant analysis rather than inferred zero effects.

### Quantifying species-level colour diversity and recurring colour directions

For the four frozen photographic colour categories, species-level sampled colour diversity was summarized by the Gini–Simpson index, D = 1 − Σ p_k², where p_k is the observed species-wide proportion of classifiable source photos in category k. The finite-sample-corrected form, 1 − Σ n_k(n_k − 1)/[n(n − 1)], was retained as a sampling sensitivity. This index represents diversity among the photographed records, including differences between locations, and is not equivalent to the frequency of multiple inherited morphs within a breeding population.

To characterize the **direction of continuous colour variation**, we normalized each nine-palette vector, applied a Hellinger transformation and identified a deterministic two-mode displacement axis within each species. The sign of the axis is arbitrary, so agreement with the fixed white–chromatic contrast was measured by the squared projection W_i = (u_i · q)², where u_i is the fitted species-level unit displacement vector and q is a fixed unit vector contrasting white with the eight nonwhite palette components. The species-equal mean W estimates the tendency for photographic colour displacement to align with the white–chromatic direction.

A matched source null preserved species composition, photo-colour category counts and acquisition cohorts while shuffling complete nine-colour photo vectors within the frozen source constraints. The complete 4,475-taxon white-axis analysis was refitted in each of 199 null worlds. Monte Carlo probabilities used the plus-one correction, so P = 0.005 is the minimum reportable value from 199 permutations (Phipson & Smyth, 2010). Previous smaller, species-disjoint confirmation cohorts were preserved as supporting historical analyses rather than silently pooled with the complete-source test.

### Species-conditional geographic structure

We estimated geographic colour organization using the same-species difference between **range-wide colour mismatch** and mismatch among conspecific photographic records within a chosen radius. For the four categorical appearances, this quantity is positive when local conspecific records are more similar than the species-wide photographed colour composition predicts. The principal comparison used a 50-km radius; additional thresholds and distance scales were treated as supporting analyses.

For each species, photo-colour labels were permuted among fixed sampled geographic positions, preserving the exact category counts, photo locations and source sampling opportunity. A species contributed one estimate to the across-species mean. We assessed 50-km photo-pair opportunity using the complete 4,475-taxon roster; species without eligible nearby classified photo pairs were not assigned a spatial effect. Observer thinning retained at most one photo per original observer within a named species.

To test whether spatial organization depended on the coarse four-colour discretization, we repeated spatial comparisons using **base-2 Jensen–Shannon divergence** (JSD) between the full original normalized nine-palette vectors. We compared 5-, 25-, 50-, 100- and 250-km distances, with month-matched colour-vector permutations preserving the photographic composition within each named species and calendar month. Genus and species sampling opportunities were retained at their observed values. The 50-km continuous analysis used 1,761 species under the strict photo-pair design; its observer-thinned counterpart contained 1,655 species.

In a separate earlier analysis, we compared species-wide sampled colour diversity D with geographic organization in the discovery and validation cohorts using observer and photographic effort covariates and species-conditional matched nulls. These are distinct estimands from the average range-wide minus local mismatch.

### Climate sampling and species-level photographic colour associations

Original source locations were linked to the archived WorldClim 2.1 10-arcminute climate normals for BIO1 (annual mean temperature), BIO4 (temperature seasonality), BIO5 (maximum temperature of the warmest month), BIO12 (annual precipitation), BIO14 (precipitation of the driest month), BIO15 (precipitation seasonality) and mean annual solar radiation calculated from 12 monthly fields. These climate values describe the photographed locations, not experimental exposure or necessarily weather at flowering.

We first quantified the signed within-species association between white-category photographic observations and BIO5, including a source-wide univariate Firth logistic estimator expressed in log odds per 1°C. The full source yielded 4,112 species with white and nonwhite classified photos suitable for this climatic comparison. As a separate continuous-colour analysis, we related the nine photo-palette fractions to climate variables, retaining both signed white-fraction covariance and nine-dimensional colour-vector direction.

Because BIO1 and BIO5 represent correlated aspects of a common thermal gradient, our integrated climatic synthesis used **five fixed, colour-outcome-blind environmental domains**. Within each species, the seven original climate variables were standardized by their photographed standard deviations. The shared thermal level was (zBIO1 + zBIO5)/2; the shared precipitation amount was (zBIO12 + zBIO14)/2. BIO4 temperature seasonality, BIO15 precipitation seasonality and mean solar radiation were retained separately.

For each named species, the five domain variables were first adjusted for photographic calendar month and a broad geographic basis using Earth-centred latitude–longitude coordinates. The residual of the target climate domain was then projected against the other four domains. A domain was estimable only when at least 5% of its original sampled variance remained after adjustment. The source analysis additionally required at least 40 classified, climate-complete original photos, 20 distinct observers and five occupied WorldClim raster grid cells.

The climatic response was the signed covariance between the uniquely residualized domain and each photographic nine-palette fraction, with the exact expectation of within-species/calendar-month permutation removed. Dividing that covariance by the domain's residual variance yielded a **conditional, slope-like photographic association** on each species' standardized climate scale. This is distinct from the per-°C Firth log-odds coefficient and from a causal thermal response. We calculated species-equal mean vectors, their squared alignment with q, and whole-genus bootstrap intervals using 999 resamples. White-direction intervals were evaluated simultaneously across the five fixed environmental domains. All five original taxonomic cohorts were also compared descriptively, and each cohort was omitted in turn in a whole-genus sensitivity analysis.

### Decomposing colour-category composition and colour shifts within categories

The original four photo-colour labels make it possible to express the climatic covariance of the full nine-palette vector as the exact sum of covariance **between photographic categories** and covariance **within photographic categories**. We estimated these two components on the same taxa and the same model. For the between-category component, a species required at least five classified photos in each of at least two original photographic colour categories. Because covariance vectors can cancel, the ratio of a signed category-between mean to its signed total is described strictly as an algebraic contribution rather than as an explained biological variance share.

The temperature and precipitation directions were compared again on **exactly the same species** with identifiable conditional photo-colour relationships to both climate domains. For each climatic domain, we computed the squared white-axis projection of the cross-species mean nine-colour coefficient vector; whole genera were resampled to quantify uncertainty in the *difference in direction* between domains. This is a geometric comparison of recorded photographic appearances, not a quantification of pigment pathways.

### Independent-observer photographic replication

Two complementary observer analyses were performed. In the full-source climate model, we deterministically retained one original image per observer within each species and repeated the shared-thermal estimation.

We also tested whether the positive mean could be reproduced without sharing observers across separate photograph sets. Within each species, original observers were assigned to two colour-outcome-blind groups before photo-colour eligibility or climate filtering. The thermal–white coefficient was fitted using only one observer group and then refitted after exchanging the groups. Held-out models preserved the same five environmental domains, geographic/month adjustments and climatic residual requirements, with eligibility recalculated from the held-out photographs. A separate version retained one image per held-out observer. All source taxa remained in the four panel registers. Whole-genus bootstrap intervals used 999 resamples.

### Contrasting colour-space geometry with ecological allocation

Finally, we compared the species-level geometric score W with geographic mismatch depletion and the original species BIO5 association in a common source-linked registry. Rank associations were adjusted for source photographic classifiability, photographed colour diversity, baseline white-label abundance, observer/sample opportunities and acquisition cohort, as appropriate. Genus bootstrap uncertainty was calculated without treating individual photograph pairs as independent species replicates. The analysis asks whether the **direction** of colour variation and its **geographic or climatic arrangement** can be summarized by one species-level metric, rather than inferring independence of evolutionary mechanisms.

## Results

### A recurring white–chromatic direction is detectable across the full source

The mean four-category, finite-sample-corrected photographed colour diversity across the 4,475 species was **D = 0.29187**. Of these species, 4,189 had at least two observed photographic colour categories, whereas **2,400 had at least two categories each supported by five or more photos**. These are observations of species-wide photo-colour heterogeneity.

The complete-source mean geometric alignment with the white–chromatic axis was **W = 0.47887**, exceeding the matched-null median **0.43204** (199-world upper P = 0.005; Fig. 1). The genus-bootstrap 95% interval for the observed mean was **0.46874–0.48978**. The direction exceeded the matched null in all five species-disjoint acquisition cohorts under the original protocol. Thus a particular orientation of colour variation recurred well beyond one or a few intensively studied taxa.

The geometric result concerns continuous nine-colour configurations: white and the eight nonwhite colours were preserved as distinct fractions for every classified photo. It therefore reflects a repeated colour-space contrast rather than simply counting the number of species whose photographs were labelled white.

### Conspecific photographs are geographically assembled rather than randomly intermixed

Among **3,600 species** with at least one classified conspecific photo pair within 50 km, the species-equal range-wide minus local mismatch was **+0.01612**, exceeding the within-species category-composition null (P = 0.005). When at least 30 nearby photo pairs and 40 classified source photographs were required, **1,761 species** remained and the difference was **+0.01761** (P = 0.005; Fig. 2). Deterministic one-photo-per-observer thinning retained **3,592 species** with a mean difference **+0.01590** (P = 0.005).

The spatial result was not limited to four-category labels. For the strict 1,761-species source panel, month-matched continuous nine-palette JSD showed **50-km geographic similarity excess +0.01146** relative to month-preserving permutations (P = 0.005). In the one-photo-per-observer continuous panel, **1,655 species** yielded **+0.01123** (P = 0.005). Additional geographic distances of 5, 25, 100 and 250 km revealed a multi-scale arrangement of conspecific colour composition.

An earlier independent-cohort comparison likewise linked species-wide photographic colour diversity to geographic organization, with adjusted rank associations **+0.12664** in the discovery cohort (369 species, matched P = 0.007) and **+0.09929** in validation (363 species, P = 0.025). Thus the abundance and geographic arrangement of photographed colour states are related in sampled species, while spatial organization remains a separate measurable property from overall colour richness.

### Thermal habitat gradients repeatedly point toward white-rich photographic colours

In the complete photo-category climate comparison, **4,112 species** had source locations with both white and nonwhite photographed colour classes. The species-equal univariate BIO5 association averaged **+0.02189 log odds per 1°C**, with whole-genus bootstrap 95% interval **+0.00846 to +0.03489**. This result establishes a small, average direction in the source photos, rather than a common per-species response magnitude.

The shared thermal-domain analysis provided a direct continuous nine-colour result. After conditioning the temperature domain on shared precipitation amount, temperature and precipitation seasonality, solar radiation, broad geographic coordinates and observation month, **2,837 species** had identifiable climatic colour associations. Their mean conditional coefficient for the photographed white-palette fraction was **+0.00821**. The squared alignment of the cross-species mean nine-colour response vector with the fixed white–chromatic axis was **85.4%**. The five-domain simultaneous whole-genus 95% interval for the scaled white-axis projection was **+0.00424 to +0.01318** (Fig. 3).

Using one original photograph per observer preserved this pattern (**2,771 species; mean white-fraction coefficient +0.00800; white-axis alignment 87.1%; scaled white-axis simultaneous interval +0.00408 to +0.01289**). Although only **54.2%** of individual eligible species had a positive photographed white/thermal coefficient, the species-equal mean was positive in each of the five original species-disjoint acquisition cohorts.

The result was also robust to the contribution of any single acquisition cohort. Removing each of the five cohorts in turn yielded a positive mean thermal–white coefficient with a whole-genus 95% interval above zero. Even after removing the largest additional cohort, 1,273 eligible species remained and the mean was **+0.00937** (95% interval **+0.00423 to +0.01410**). These comparisons show that the recurring mean photo-colour direction was not generated solely by one collection wave.

### Warmth and precipitation are associated with different mean colour-space directions

For the **same 2,689 species** with identifiable shared thermal and shared precipitation gradients, **83.7%** of the squared cross-species mean thermal colour vector followed the white–chromatic axis. The corresponding alignment for the shared precipitation domain was **24.3%**. The observed difference was **59.4 percentage points**, with whole-genus bootstrap 95% interval **6.0–86.8 percentage points**. In the observer-thinned same-taxon set (**2,625 species**), the respective alignments were **84.8%** and **12.7%**; the paired difference was **72.0 points** (95% interval **10.0–88.3**; Fig. 3).

These comparisons identify differences in the **directions** of the mean photographic colour response to two distinct climatic domains. They do not require that the dominant sampled colour be the same in all species or that the climatic predictor cause a particular developmental change. In this grouped analysis, the precipitation-domain white projection and the individual nonwhite hue components were less consistently resolved than the shared thermal–white direction.

### Climatic colour associations contain between-category and within-category components

For **1,376 identical species** containing at least two adequately supported photographic colour categories, the shared thermal-axis mean white-fraction coefficient was **+0.01243**. The exact photographic covariance partition attributed **+0.00843** to differences between the original categorical photo groups and **+0.00400** to continuous variation within the categories (Fig. 4). The between-category component represented **67.8%** of the signed species-equal total (whole-genus 95% interval **44.2–85.4%**).

The observer-thinned common-species comparison retained **1,299 species**, with an algebraic between-category share of **66.1%** (95% interval **39.2–84.3%**). Accordingly, the recorded climatic colour relationship was not restricted to small within-bin colour differences, nor could it be represented solely by changes among the four coarse photo labels. Both aspects contributed to the observed thermal geometry.

### Independent observers reproduce the across-species mean

To test the source photo-temperature association separately from photographed baseline colour-state assignment, the observer-disjoint experiment refitted the shared thermal/white slope from one set of photographers at a time. With one fold held out, **2,403 species** had an average white/thermal coefficient of **+0.01288** (genus-bootstrap 95% **+0.00733 to +0.01874**). Reversing the observer folds yielded **2,430 species**, mean **+0.00575** (95% **+0.00012 to +0.01128**). Restricting each held-out observer to one source photograph preserved positive means in both directions (**+0.01101** across 2,324 species and **+0.00616** across 2,364 species; Fig. 4).

This supports the repeated cross-species mean photo-colour association beyond reuse of the identical photographs or observers to compute baseline state and response. The crossfold estimates for *individual species* were more variable: among 1,919 species with both estimates, Pearson correlation of the two per-species thermal coefficients was **+0.018**. Thus the comparative mean is a more stable result than the ranking of individual species by response strength using only half their source photographs.

### Colour direction and ecological assembly are related but distinguishable

The recurring white–chromatic geometric score was not a strong predictor of how the same species' photographed colours were spatially allocated. The adjusted species-level association between white-axis alignment W and 50-km local mismatch depletion was **r = +0.03375** across 3,600 species. The corresponding association between W and signed BIO5 photo-climate slope was **r = −0.00927** across 4,112 species. These small associations indicate that repeatedly reaching a similar visual contrast does not, on its own, determine where its appearances are photographed.

Together, the results identify three complementary dimensions of intraspecific floral photographic variation: **a recurrent phenotype-space direction, a measurable geography of colour allocation, and a thermal orientation of the species-mean colour vector**.

## Discussion

### A repeated phenotypic direction across flowering plants

A consistent result of this study is that the species-level directions of photographic colour variation are not isotropic. Across the entire 4,475-species register, the white–chromatic direction occurs more prominently than expected when the original species-level photographic colour categories, opportunities and cohort structure are preserved. This extends the logic of comparative intraspecific trait analysis from trait magnitude or categorical richness to the **orientation of variation in trait space**.

The recurrence of white–chromatic differences invites developmental hypotheses. Anthocyanin regulation, alternative pigment production, tissue-specific gene expression and differing floral organ structure can all generate pale versus chromatic appearances (Warren & Mackenzie, 2001; Wessinger & Rausher, 2012; Casimiro-Soriguer et al., 2016; Larter et al., 2018). Yet repeated appearances need not represent repeated genetic substitutions. The contribution of our global comparison is to identify a reproducible *visual direction* for mechanistic studies rather than to assign a single molecular route to all 4,475 species.

### Geographic colour assembly is a separate ecological property

The stronger resemblance of conspecific colours at geographic distances within 50 km establishes an additional regularity: species-wide photographed colour compositions are spatially arranged rather than completely intermixed. This structure remains observable in continuous colour metrics and after excluding repeated observations by the same photographer, making it a general feature of the sampled photographic atlas rather than a consequence of choosing exactly four categorical labels.

The magnitude of spatial mismatch depletion is modest and highly heterogeneous among taxa. Its ecological interpretation is therefore comparative rather than categorical: **the atlas detects repeated excess geographic organization**, not a universal division into monomorphic versus polymorphic local populations. Spatially structured climate, dispersal limitation, demography, genetic history and developmental plasticity are distinct plausible contributors (Schemske & Bierzychudek, 2007; Wang & Bradburd, 2014). A shared phenotypic displacement direction and spatial assembly were only weakly associated, emphasizing that traits can converge in apparent colour without converging in geographic organization.

### Shared thermal environments orient the mean photographic colour response

The climatic analyses provide a third dimension. Temperature-linked colour associations have been reported at the level of dominant flower colours among species and at the level of different colours within selected populations (Dalrymple et al., 2020; Narbona et al., 2018). Our comparison evaluates how that association is *oriented* in a continuous colour space across thousands of taxa.

The positive thermal–white relationship survived grouping BIO1 and BIO5 into a fixed shared thermal coordinate and conditioning on precipitation amount, climatic seasonality and solar radiation. This formulation estimates a thermal-domain pattern, rather than forcing highly correlated temperature normals to yield artificially separable individual effects. The thermally associated mean colour vector was mainly aligned with the white–chromatic axis, and that alignment differed from the precipitation-associated vector on the same species.

The result was repeated across the five sampled taxon cohorts, under cohort-deletion resampling and using observers disjoint from those supplying baseline colour classifications. These converging analyses support **a recurrent species-mean photographic association between thermal habitat context and white-rich appearances**. The fact that the group mean is more consistent than individual species' slopes suggests that the comparative signal is expressed against substantial species-level heterogeneity, potentially including differences in colour composition, sampled climate range and photographic opportunity.

### Colour-state composition and continuous expression both matter

The exact between/within photograph-category partition shows why the ecological response should not be reduced to either categorical morph frequencies or continuous hue shifts alone. About two thirds of the signed thermal white coefficient within the eligible multicolour photo panel appeared between the original broad photographic colour categories, while the remaining third arose within them. The two contributions describe different statistical parts of the same nine-palette relationship.

This organization provides a useful language for subsequent flower-colour field studies. In a genotyped natural population, a between-morph contribution might reflect differences in genetic morph frequencies, whereas a within-morph contribution might reflect variation in pigment quantity, environmental response or flower development. In our photographic analysis the categories have not been assigned those biological meanings. Nevertheless, separating category composition from within-category hue change clarifies the observations that future field experiments need to explain.

### Towards an evolutionary framework of direction, assembly and maintenance

These results suggest organizing the evolution of within-species floral colour into three linked, but nonidentical, questions. **Which colour differences can be generated?** Repeated white–chromatic displacement identifies a common visual contrast. **Where are those differences assembled?** Nonrandom conspecific colour geography describes one outcome of history and environmental context. **Under which ecological gradients do frequencies or intensities change?** The shared thermal association connects photographic colour composition with an environmental direction.

A further question concerns which resulting variants persist as inherited morphs within natural breeding populations. Pollinator preferences, opposing physiological benefits, developmental pathways, dispersal and demographic processes may jointly determine that persistence (Delph & Kelly, 2014; Narbona et al., 2018; Sapir et al., 2021). Pigment maintenance costs and benefits are therefore plausible mechanisms motivated by this framework, not empirical parameters derived from the photographic atlas. The most decisive next mechanistic studies would pair validated fresh-petal reflectance and pigment assays with independent morphological or genetic morph identification, repeated local morph frequencies and component-specific fitness across thermal contexts.

Our photographic method also has a defined scope. About half the original source photographs did not receive a valid colour classification; the source cohorts share one platform and image-analysis pipeline. Photographs may contain flower-age, floral-organ and background differences, and historical site-climate normals need not equal the conditions during flower development. A separate background-based climatic photo-colour analysis was not confirmatory. Consequently, the patterns should be interpreted as **reproducible comparative photographic phenotypes**, without treating them as measurements of genetic colour polymorphism or direct climatic selection.

### Conclusions

Across 4,475 angiosperm species, intraspecific photographic floral colour variation exhibited three recurring features: **a strong white–chromatic component in the direction of variation, nonrandom geographic assembly of conspecific colours, and a shared thermal orientation of mean colour composition**. These features are separable rather than alternative names for one ecological process. Their joint analysis shifts comparative flower-colour research from asking merely how much colour diversity a species possesses to asking **which differences recur, how they are distributed and how they relate to environmental gradients**. This organization offers a general basis for connecting broad floral phenotype surveys with targeted tests of pigment function, population history and the maintenance of inherited colour morphs.

## Acknowledgements

[Funding, data providers, field collaborators and computational contributors to be added after author verification.] We acknowledge the contributors to the photographic biodiversity observations on which the source atlas is based.

## Author contributions

[To be completed and approved by all authors.]

## Competing interests

[To be completed by the authors.]

## Data availability

The source photographs, fixed computational provenance, derived species-level measurements, analysis scripts and input checksums are documented in the FCP repository: https://github.com/zuizui0223/fcp. The immutable original 447,500-photograph source is identified by release tag **fcp-all4730-original-method-actual-colour-results-20261017** and original compressed photo CSV SHA256 **8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6**. The original WorldClim BIO and solar archives, frozen observer controls and species-level output receipts are listed in the Supporting Information.

The complete-source white-axis, geographic and original BIO5 estimates are maintained in the existing New Phytologist hypothesis manuscript and its source-locked results. Later source analyses are traceable through FCP pull requests **#260–#271** and the respective authenticated GitHub Actions runs, including [nine-colour geographic structure](https://github.com/zuizui0223/fcp/actions/runs/39918785955), [shared thermal/precipitation photo-climate geometry](https://github.com/zuizui0223/fcp/actions/runs/40066224301), [five-cohort thermal holdout](https://github.com/zuizui0223/fcp/actions/runs/40066704796), [observer-disjoint source-photo climate test](https://github.com/zuizui0223/fcp/actions/runs/40080001765), and [independent crossfold source audit](https://github.com/zuizui0223/fcp/actions/runs/40080546906).

## References

Albert CH, Thuiller W, Yoccoz NG, Soudant A, Boucher F, Saccone P, Lavorel S. 2010. Intraspecific functional variability: extent, structure and sources of variation. *Journal of Ecology* **98**: 604–613. https://doi.org/10.1111/j.1365-2745.2010.01651.x

Carmona CP, de Bello F, Mason NWH, Lepš J. 2019. Trait probability density (TPD): measuring functional diversity across scales based on TPD with R. *Ecology* **100**: e02876. https://doi.org/10.1002/ecy.2876

Casimiro-Soriguer I, Narbona E, Buide ML, del Valle JC, Whittall JB. 2016. Transcriptome and biochemical analysis of a flower color polymorphism in *Silene littorea* (Caryophyllaceae). *Frontiers in Plant Science* **7**: 204. https://doi.org/10.3389/fpls.2016.00204

Dalrymple RL, Kemp DJ, Flores-Moreno H, Laffan SW, White TE, Hemmings FA, Moles AT. 2020. Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist* **228**: 1972–1985. https://doi.org/10.1111/nph.16737

Delph LF, Kelly JK. 2014. On the importance of balancing selection in plants. *New Phytologist* **201**: 45–56. https://doi.org/10.1111/nph.12441

Di Cecco GJ, Barve V, Belitz MW, Stucky BJ, Guralnick RP, Hurlbert AH. 2021. Observing the observers: how participants contribute data to iNaturalist and implications for biodiversity science. *BioScience* **71**: 1179–1188. https://doi.org/10.1093/biosci/biab093

Laitly A, Callaghan CT, Delhey K, Cornwell WK. 2021. Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution* **11**: 4071–4083. https://doi.org/10.1002/ece3.7307

Larter M, Dunbar-Wallis A, Berardi AE, Smith SD. 2018. Convergent evolution at the pathway level: predictable regulatory changes during flower color transitions. *Molecular Biology and Evolution* **35**: 2159–2169. https://doi.org/10.1093/molbev/msy117

Luong Y, Gasca-Herrera A, Misiewicz TM, Carter BE. 2023. A pipeline for the rapid collection of color data from photographs. *Applications in Plant Sciences* **11**: e11546. https://doi.org/10.1002/aps3.11546

Narbona E, Wang H, Ortiz PL, Arista M, Imbert E. 2018. Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology* **20** (Suppl. 1): 8–20. https://doi.org/10.1111/plb.12575

Palacio FX, Graco-Roza C, de Bello F, Carmona CP. 2025. Integrating intraspecific trait variability in functional diversity: an overview of methods and a guide for ecologists. *Ecological Monographs* **95**: e70024. https://doi.org/10.1002/ecm.70024

Phipson B, Smyth GK. 2010. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology* **9**: Article 39. https://doi.org/10.2202/1544-6115.1585

Sapir Y, Gallagher MK, Senden E. 2021. What maintains flower colour variation within populations? *Trends in Ecology & Evolution* **36**: 507–519. https://doi.org/10.1016/j.tree.2021.01.011

Schemske DW, Bierzychudek P. 2007. Spatial differentiation for flower color in the desert annual *Linanthus parryae*: was Wright right? *Evolution* **61**: 2528–2541. https://doi.org/10.1111/j.1558-5646.2007.00219.x

Wang IJ, Bradburd GS. 2014. Isolation by environment. *Molecular Ecology* **23**: 5649–5662. https://doi.org/10.1111/mec.12938

Warren J, Mackenzie S. 2001. Why are all colour combinations not equally represented as flower-colour polymorphisms? *New Phytologist* **151**: 237–241. https://doi.org/10.1046/j.1469-8137.2001.00159.x

Wessinger CA, Rausher MD. 2012. Lessons from flower colour evolution on targets of selection. *Journal of Experimental Botany* **63**: 5741–5749. https://doi.org/10.1093/jxb/ers267

## Tables

### Table 1. Frozen photographic resource and analytic opportunities

| Source or analysis | Number of species | Source photographs or eligibility |
|---|---:|---|
| Original source sampling universe | 4,475 | 447,500 photographs, 100 per species |
| Nine-palette colour-classifiable records | 4,475 | 227,107 source photographs; species enter separate analyses according to observation opportunity |
| White–chromatic geometric direction | 4,475 | Hellinger two-mode colour displacement, composition-preserving null |
| 50-km categorical geographic assembly | 3,600 | At least one valid conspecific nearby photo pair |
| Strict 50-km, categorical and continuous spatial comparisons | 1,761 | Species meeting strict original photo/pair opportunities |
| Shared thermal conditional nine-colour analysis | 2,837 | Original photos, five fixed climate domains |
| Observer-thinned thermal conditional analysis | 2,771 | At most one original photo per observer |
| Matched temperature versus precipitation alignment | 2,689 | Exactly the same source taxa with both climatic axes |
| Independent-observer held-out temperature comparison | 2,403 / 2,430 | Forward / reversed source-observer folds |

### Table 2. Principal positive species-equal results

| Hypothesis or estimand | Supported observation | Result and uncertainty | Denominator |
|---|---|---|---:|
| Recurrent visual direction | White–chromatic displacement W exceeds matched source null | W = 0.47887 vs null 0.43204; P = 0.005 | 4,475 species |
| Geographic colour assembly | Conspecific records within 50 km are more similar than the species-wide photographic composition | Mismatch depletion +0.01612; P = 0.005 | 3,600 species |
| Continuous geographic assembly | Nine-palette 50-km JSD is lower than its month-preserving source null | Similarity excess +0.01146; P = 0.005 | 1,761 species |
| Signed photo white versus BIO5 | Average white photo-category occurrence is associated with warm-season temperature | Mean log-odds coefficient +0.02189 °C⁻¹; genus 95% +0.00846 to +0.03489 | 4,112 species |
| Shared thermal colour orientation | Temperature-related mean nine-colour vector aligns with white | 85.4% of squared mean vector; white coefficient +0.00821; q-projection simultaneous 95% +0.00424 to +0.01318 | 2,837 species |
| Different climatic directions | Shared thermal and precipitation mean colour vectors have different white-axis alignments | 83.7% vs 24.3%; paired difference 59.4 points, genus 95% 6.0–86.8 | 2,689 identical species |
| Photo-category composition | Both original photo category differences and variation within categories contribute to thermal association | Total +0.01243 = between +0.00843 + within +0.00400 | 1,376 identical species |
| Observer-disjoint thermal reproduction | Mean white/thermal association remains positive when observations from baseline observers are excluded | Forward +0.01288; reverse +0.00575; both genus intervals above zero | 2,403 / 2,430 species |

## Figure legends and assembly plan

**Fig. 1 | The recurring white–chromatic direction of within-species floral photographic variation.** Full 4,475-species W geometry against source cohort- and photo-colour-category-preserving null expectations. Original verified source figure: [white-axis geometry](https://github.com/zuizui0223/fcp/blob/main/figures/fcp_np_all4475_white_axis_five_cohorts_20261021.svg).

**Fig. 2 | Geographic assembly of intraspecific photographic colours.** Species-equal range-wide minus local (≤50-km) photo-category mismatch, with original categorical null comparisons and separate continuous nine-palette spatial similarity. Source-verified original figure for the principal categorical result: [50-km geographic structure](https://github.com/zuizui0223/fcp/blob/main/figures/fcp_np_all4475_figure2_geo_50km_20261021.svg). The continuous nine-colour result uses the source-verified analysis associated with Actions run 39918785955.

**Fig. 3 | Shared climatic directions of continuous colour variation.** Main synthesis panels: (a) thermal versus precipitation colour-axis alignment on 2,689 identical taxa (83.7% versus 24.3%); (b) positive shared thermal white-fraction conditional coefficients in original and observer-thinned photos (2,837 and 2,771 taxa); (c) five acquisition-cohort thermal means, each positive; and (d) cohort-leave-one-out means. All quantities are frozen in the PR #269 source machine result and the five-cohort audit. [Source result](https://github.com/zuizui0223/fcp/blob/analysis/fcp4475-shared-thermal-moisture-axes-20261024/results/fcp4475_shared_thermal_moisture_photo_hues_verified_20261024.json). A new consolidated multipanel figure should be assembled from those existing numbers; no additional biological hypothesis test is required.

**Fig. 4 | Composition and observer replication of the thermal photographic colour association.** (a) Exact between-category and within-category partition for the same 1,376 source species; (b) positive species-equal means for two independent-observer folds, with whole-genus 95% intervals; and (c) a conceptual diagram separating colour direction, geographic assembly and ecological context. Source results: [thermal partition](https://github.com/zuizui0223/fcp/blob/analysis/fcp4475-shared-thermal-moisture-axes-20261024/results/fcp4475_shared_thermal_moisture_photo_hues_verified_20261024.json) and [observer folds](https://github.com/zuizui0223/fcp/blob/analysis/fcp-crossfitted-photo-colour-thermal-20261024/results/fcp4475_observer_disjoint_thermal_white_verified_20261024.json).

## Supporting Information

**Methods S1.** Source sampling, taxonomic inclusion and photograph-level classifiability.

**Methods S2.** Image detection, flower-ROI segmentation, four-colour labels, nine-palette measurement and archive checksums.

**Methods S3.** Species-conditional spatial permutations across 5, 25, 50, 100 and 250 km, calendar-month sensitivity and observer thinning.

**Methods S4.** Fixed WorldClim domain construction, unique environmental variance, matched-month covariance expectations and whole-genus bootstrap.

**Methods S5.** Photographic colour-category covariance decomposition and common-support climatic direction comparison.

**Methods S6.** Independent-observer split, held-out thermal modelling and reproducibility receipts.

**Table S1.** Complete 4,475-species photographic measurement and eligibility register.

**Table S2.** Five-cohort geometric and 50-km geographic comparison.

**Table S3.** Nine-hue, seven-climate and fixed macroclimate source summaries with observer thinning.

**Table S4.** Source-cohort omission and independent-observer thermal-climate replication.

**Figure S1.** Photographic geographic and colour-classification opportunities.

**Figure S2.** Multiscale continuous nine-hue geographical similarity.

**Figure S3.** Source photo missingness, camera and floral-organ interpretive sensitivity.

**Figure S4.** Evolutionary model distinguishing phenotypic generation, geographic assembly and morph maintenance.
