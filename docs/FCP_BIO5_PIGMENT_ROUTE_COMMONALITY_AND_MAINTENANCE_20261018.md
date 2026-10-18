# FCP common principle: environment × pigment-route modularity × population-genetic maintenance filter

**18 October 2026 — actual five-cohort BIO5 reanalysis and independent primary molecular evidence.** This document tests whether flower-colour polymorphism is constrained by a common evolutionary process, not whether a particular colour is globally common. **It is not evidence of a universal BIO5-mediated selective mechanism**, nor a gene-level analysis of the photographed iNaturalist individuals.

## Evidence provenance and nonnegotiable scales

This analysis replayed the immutable measured **4,475 species ×100 = 447,500 original photographs** (227,107 four-colour-classifiable) from [the authentic five-cohort original-pixel release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017). Frozen merged photo CSV SHA256: \`8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6\`. The predictor was original **WorldClim 2.1 BIO5** (maximum temperature of warmest month, **1970–2000 climatological mean**, 10 arc-minute resolution, not 2023–2025 photo-date temperatures). The source TIFF SHA256 was \`32e02d85868734c32547da70cf9c620f4739a3cf5dbc4f2b751e16ebb2338ae1\` from archive \`wc2.1_10m_bio.zip\` SHA256 \`00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf\`.

[Source-authenticated full replay #39100344149](https://github.com/zuizui0223/fcp/actions/runs/39100344149): 11 regression checks passed. Full output \`fcp_real_447500_BIO5_unfrozen_transport.json\` SHA256 **\`93267d3c492e8bf0392f2faf89d3338f54210d80ebb856d4c7ab710e562e0ae6\`**, independent molecular/photo resolution JSON **\`e7d747079645610d3002366ec66430621bd1c8962bb0fe77547479a1c7b9d64b\`**, phase-model JSON **\`be6146c317c8a80fabfa1fdbbef828ab2173df5f3186ef09e9445a6bd2519267\`**. The full species-level BIO5 effect table contains **taxon names and aggregate effects only**, no individual photo IDs/observer IDs/exact locality coordinates. Actual execution artifact is \`fcp-real-447500-bio5-transport-and-molecular-trait-resolution\` from the same run.

**Important chronology:** The earlier [*New Phytologist* original prospective environmental-filter result](https://github.com/zuizui0223/fcp/blob/main/results/polymorphism_white_environment_mechanism_20260925/result.json) used a distinct, response-blind **high-clip image exposure exclusion**, present only in the original third-cohort validity study. Its 281-species white-vs-colour BIO5 Wilcoxon-Holm p=0.0354487 and conditional near-clip-adjusted OR=1.07347 (95% CI 1.02940–1.11944, p=0.0009188) **remain unchanged**. Earlier observer-paired p=0.2023 and frozen discovery/validation joint transport **NOT SUPPORTED** remain unchanged. The five-cohort reanalysis below is **new, post hoc**, using the same original photo classification but without claiming identical high-clip corrections for later new-wave photos. Do not call it a re-confirmation of the original experimental pipeline.

## Actual test 1 — is a signed white–BIO5 contrast shared across species?

Each classifiable true original photo was assigned its SHA-verified fixed WorldClim cell BIO5. Within each species, climate was z-standardized (species mean and SD). A species was included only when >=5 white and >=5 nonwhite photographs occurred in finite BIO5 cells, with within-species BIO5 variance >0. **Species, not photo records, were equal independent summary units.** The response is the observed *photograph-location contrast* \`mean z(BIO5 | white) − mean z(BIO5 | nonwhite)\`, not an environmental effect on germination, plant fitness or inherited morph frequencies.

| Authenticated photo cohort | Eligible species | White-minus-nonwhite BIO5 median (species SD) | Species Wilcoxon two-sided p |
|---|---:|---:|---:|
| Discovery | 271 | **−0.006777** | 0.7433 |
| Validation | 260 | **+0.066165** | 0.0541 |
| Prospective photos (new no equivalent high-clip screen) | 286 | **+0.053601** | 0.0247 |
| Fourth wave | 297 | **+0.051032** | 0.0718 |
| **Fifth wave new U100** | **1,295** | **+0.080270** | **3.06 × 10⁻¹⁰** |
| **Full five-cohort pooled species-equal** | **2,409** | **+0.066010** | **1.20 × 10⁻¹¹** |

Among all 2,409 species, 55.83% have positive signed contrasts; pooled genus-equal summary across **1,075 genera** is median **+0.05436 SD**, Wilcoxon p=**2.13 × 10⁻⁶**. This means a weak, repeated **PHOTO sampling association** exists in the original measured photographic source. It is not monotonic across all separate species cohorts, and significance from 2,409 species does not eliminate sampling or observer selection bias.

### Same observer and same month/year are DIFFERENT tests

A more direct challenge is to compare *white versus nonwhite photos from the SAME identified photographer inside the SAME species*, first computing within-species × observer contrast then equal-weighting one summary per species. This restricts 2,409 species to **731** eligible species and 1,067 observer strata.

| Exactly the same 731 species | Median signed BIO5 (species SD) | Wilcoxon two-sided p |
|---|---:|---:|
| Ignore observer identity but restrict to exactly these species | **+0.054593** | **0.0004266** |
| Contrast within each species × SAME observer first | **0.000000** | **0.93348** |

The loss of support is therefore **not explained merely by losing species**: the unpaired contrast is supported within the identical 731-species subset. All individual observer-paired cohort tests are unsupported: discovery 71 species p=0.744, validation 67 p=0.140, prospective 89 p=0.667, fourth wave 84 p=0.512 and **new fifth wave 420 species p=0.930**. However, a same-observer restriction reduces spatial climate coverage, so the outcome also cannot strictly *prove* observer bias causes the whole original association.

Matching **species × the same calendar YEAR-MONTH** instead retains 2,360 species and produces median +0.07355 SD (Wilcoxon p=3.15 × 10⁻⁸). That year-month comparison does not resolve photographer/site confounding, and climatological BIO5 does not change with the year of observation. In contrast, a same-species positive BIO5 contrast among photo sites could be entirely generated by different observers taking differently coloured flower photographs in different parts of a species range.

**BIO5 conclusion:** an association in broad photos persists, but a common, observer-controlled thermal sorting/selection mechanism is **not supported**. Do not infer species-wide direction of anthocyanin gain/loss from this signed white-label effect. The underlying high-clip/exposure coupling of the white classifier and the earlier failed cross-cohort transport remain relevant.

## Actual test 2 — does one molecular route explain recurrent FCP?

A separate [five-cohort photo × independent primary molecular source registry](https://github.com/zuizui0223/fcp/actions/runs/39100344149) explicitly distinguishes pathway **association, causal genotype, morph fitness, reproductive isolation and taxon-name overlap**. These are independent literary sources, **not allele-phenotype matches to the photo individuals**.

| Source biological system | Independent molecular/mechanistic evidence | Relationship to photographed FCP and its maintenance |
|---|---|---|
| *Silene littorea* (Casimiro-Soriguer et al. 2016; doi:10.3389/fpls.2016.00204) | White-bud *F3h1* expression down >42-fold; *Myb1a* regulatory candidate; flavonoid/HPLC differences. No consistently colour-diagnostic SNP identified. | Flower-restricted molecular-expression **candidate** for inherited local FCP. **Not present in the 4,475 species photo table**. The regulatory causation remains untested. |
| *Ipomoea purpurea* (Zufall & Rausher 2003; doi:10.1093/jhered/esg098; Coberly & Rausher 2008; doi:10.1111/j.1558-5646.2008.00355.x) | A documented *F3′H* **coding** loss-of-function creates pink vs purple; a separate white *A* allele increases selfing transmission but reduces survival through pleiotropy. | **47** classified photos: blue/purple 28, red/pink 11, white 7, yellow/orange 1. **No genotyping of these photographed plants**. Demonstrates that not all colour polymorphisms are MYB regulatory lesions or pure pollinator selection. |
| *Mimulus guttatus* / *Erythranthe guttata* (anthocyanin polymorphism; doi:10.1086/694853 and 2012 PLA1 study) | Natural unpigmented allele frequency >20% in one population; lower MYB5 expression and PLA1/MYB cluster association; no fully established pleiotropic/fitness cause. | **60 classified photos, all yellow/orange.** Its actual anthocyanin variation occurs in floral/vegetative **spots/organs** invisible to the coarse whole-flower colour states. A colour-class "0% polymorphism" assertion would be biologically false. |
| *Mimulus aurantiacus* / *Diplacus aurantiacus* (Streisfeld et al. 2013; doi:10.1371/journal.pgen.1003385; premating isolation 2014) | Causal *MaMyb2* **cis-regulatory** change in red–yellow ecotype divergence; divergent pollinator preferences produce premating isolation while neutral loci can exchange across a hybrid zone. | **70** classified photo rows: yellow/orange 66, white 3, blue/purple 1. The photograph data do **not** resolve source ecotypes, *MaMyb2* genotype or genetic barriers; this is a **speciation-side example**, not protected within-locality FCP. |
| PAL versus WAL (*Silene* and literature, Del Valle et al. 2019; doi:10.1186/s12870-019-2082-6) | Flower-only anthocyanin loss (PAL) preserves photosynthetic-tissue pigment; whole-plant loss (WAL) eliminates more functions. In the compiled 13 PAL/13 WAL system literature and *S. littorea*, PAL morphs attain higher reported frequencies. | Suggests a **pleiotropic maintenance filter**, NOT a causal unbiased cross-species estimate. The published source set has ascertainment, varying sampling efforts and no direct genetic genotype×full-photo morph linkage. |

**Do not pool gene p-values across these studies.** Some provide a molecular *candidate*, others experimentally demonstrated alleles or reproductive fitness outcomes. Their biological response variables differ (petal hue, petal/leaf anthocyanin, nectar-guide spots, cross-ecotype isolation). The evidence instead supports a general *architecture distinction*: a conserved flavonoid/anthocyanin biosynthetic pathway permits different **structural, regulatory and tissue-specific** changes, with very different organism-wide costs and mating ecological effects.

The 2019 source-described *Silene littorea* PAL positive frequencies ranged **8–21%**, vs WAL **0.05–0.86%**; the literature table is 13 PAL and 13 WAL nonrepresentative reported systems. This is consistent with the idea that avoiding whole-plant pigment-function losses can improve persistence, but frequency ascertainment, age/mutation rate and field selection are unresolved. *M. guttatus* is an especially instructive counterexample: systemic anthocyanin-loss alleles **can** reach appreciable local frequency, so petal-only loss is not a necessary condition for polymorphism.

## An explicit genotype-level "maintenance corridor" instead of a new climate colour law

To separate phenotype formation from maintenance, let \`p_j\` be the **true inherited frequency of an alternative morph among interbreeding individuals**, NOT photo detections. A toy population-genetic approximation is

\[
\dot p_j = p_j(1-p_j)\{s_j(T,\mathrm{pollinators},\mathrm{habitat})-c_{\mathrm{route}}+k(1-2p_j)\}+
m(\bar p-p_j).
\]

- \`s_j\`: a **locally measured relative fitness advantage**, potentially determined by ecological conditions including actual thermal environment or pollinators.
- \`c_route\`: differential viability/pleiotropic cost of **the true inherited genetic alteration** (different for petal-restricted, vegetative+petal, coding vs regulatory variants). Tissue restriction is not synonymous with no cost.
- \`k\`: negative frequency-dependent protection of the rare morph, to be estimated by **frequency-manipulated field reproduction**, not the number of differently coloured photos.
- \`m\`: gene flow, to be measured from genetic markers, pollination and parentage.

With **no migration and k>0**, invasion of BOTH alternative genetic morphs when rare requires exactly

\[
s_j-c_{\mathrm{route}}+k>0,\qquad k-(s_j-c_{\mathrm{route}})>0
\]
\[
\Longleftrightarrow\quad |s_j-c_{\mathrm{route}}|<k.
\]

Thus a locally protected inherited FCP is possible only inside a **maintenance corridor** around the equilibrium fitness difference. The same pigment phenotype could be locally maintained in one population but fixed/lost in another if environment alters \`s_j\`; a globally observed geographic morph mosaic therefore does NOT by itself prove local frequency dependence. A tissue-specific regulatory alteration changes the *cost side of this inequality* without specifying whether warmer BIO5 favours white or pigmented flowers.

This corridor is a **testable mathematical prediction under the stated toy assumptions, not a fit to FCP photographs**. In \`scripts/analysis/fcp_bio5_pigment_route_maintenance_phase_20261018.py\`, explicitly arbitrary numerical costs \`0.02\` versus \`0.24\`, negative frequency-dependence \`k=0.20\` and a hypothetical signed climate-selection gradient are used to demonstrate different invasion/fixation regimes. They do not estimate actual plant costs or BIO5-specific fitness slopes.

**Non-identifiability:** even if a field population's true genetic morph frequency \`0<p_j<1\` were perfectly measured, for arbitrary \`k\` and \`m\` an appropriate \`s_j-c_route\` can reproduce the *same equilibrium*. The number of distinct photographed colours has less information still. Therefore static photo D, overall hue prevalence and 50-km spatial autocorrelation cannot identify selection, gene flow or speciation.

The degree of pre-mating/post-mating reproductive isolation is a **separate biological quantity absent from this model**. It must be measured before evolving phenotype differences are reclassified as ecological speciation.

## Specific falsifiable cross-species predictions and minimum data needs

1. **Thermal sorting beyond photographer allocation:** independently sampled marked individuals of the same inherited colour genotype should associate with local measured temperature after site, observer, date, background, habitat and genetic ancestry adjustment. A common sign across genetically validated species would add environmental generality; the current broad photo-BIO5 trend with observer-paired p=.933 does NOT meet this gate.
2. **Pathway modularity versus systemic pleiotropy:** genotype/sequence + petal/leaves HPLC/anthocyanin expression, common garden and segregating crosses should distinguish structural loss, petal-specific regulatory loss and whole-plant loss **within independent species**; test whether the same relative-cost difference predicts survival/morph persistence after controlling mutation origin and demography.
3. **Local true polymorphism versus geographic replacement:** marked, genetically confirmed colour-morph plants in the same breeding population followed for multiple flowering seasons, with manipulations of rare-versus-common morph frequency and measurement of morph-specific seed, pollen and germination fitness. The [recent strict 5-km/photo-date null](https://github.com/zuizui0223/fcp/blob/main/docs/FCP_LOCAL_POLYMORPHISM_MAINTENANCE_VS_SPECIATION_20261018.md) was non-support (194 photo-eligible species, p=.300), and the [interannual candidate audit](https://github.com/zuizui0223/fcp/blob/main/docs/FCP_RECURRENT_LOCAL_PHOTO_MORPH_QA_AND_MAINTENANCE_IDENTIFIABILITY_20261018.md) retained only eight *unconfirmed* strict 1-km species; the FCP photo atlas alone cannot supply this field gate.
4. **Reproductive-isolation/speciation boundary:** pollinator assortativity, compatibility, cross fitness, and neutral versus colour-associated genetic differentiation must be established; \`MaMyb2\` red/yellow *Diplacus* illustrates why pollinator-mediated divergence is not the same as two morphs coexisting in one interbreeding population.

## Paper scope and decision

**Recommended research-level central hypothesis (not yet empirically confirmed):** a conserved pigment network creates multiple phenotypic routes, but **fitness effects of tissue-specific regulatory/coding changes, local selection and gene flow together determine which variants persist as within-population FCP, become local geographic mosaics, or contribute to reproductive isolation**.

**Supported now:** (1) authentic photographed FCP spatial sorting, (2) positive broad photo-BIO5 association with strong observer-conditional failure, (3) molecularly distinct flower-colour mechanisms in independent primary systems, (4) strong limitations to trait matching and local genetic morph identification.

**Not supported now:** universal BIO5 selection, uniform gene-level genetic convergence, a pooled mutation-to-FCP causal relation across 4,475 taxa, a universally adaptive white phenotype, or a general reproductive-isolation mechanism.

Maintain the [still-open manuscript PR #243](https://github.com/zuizui0223/fcp/pull/243) on *biological-scope hold* until its central claims distinguish photographic variation from genetically maintained local FCP. This document and the new BIO5/molecular bridge do not reinterpret the original species-disjoint untouched H2 confirmation, and they keep the separate MOIS interspecific seasonal climate programme out of the FCP paper.
