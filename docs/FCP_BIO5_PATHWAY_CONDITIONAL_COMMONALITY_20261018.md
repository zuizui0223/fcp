# FCP: BIO5, pathway modularity and conditional commonality of within-species differentiation

Date: 18 October 2026. Status: **new source-overlap audit plus a falsifiable mechanistic research design**, not an estimated universal mechanism, a new molecular experiment, or a new climate-effect fit. Preserve the frozen original W/50-km/BIO5 results and the later negative local-month results. This note does not release the biological-scope hold on paper PR #243.

## 1. Main question and hypothesis

**When can a flower-colour change evolve without disrupting the physiological functions of the rest of the plant, and when does that accessible variation persist as polymorphism rather than disappear or fix?**

A testable candidate for cross-species commonality is not a universal warm-to-white slope or one universally mutated gene. It is **environment-dependent selection acting on alternative molecular routes with different tissue specificity and pathway side effects**. Altering a floral regulatory module may change a reproductive signal while leaving vegetative pigment functions relatively intact. Blocking an early shared biochemical step may change both visible pigmentation and other flavonoids. Identical-looking white flowers can therefore have different physiological costs, thermal responses and evolutionary fates.

This proposed principle primarily concerns biochemically verified anthocyanin-related morph contrasts. White/coloured classifications in the four-state photograph atlas do not diagnose anthocyanins, carotenoids, betalains or a causal locus. Keep pigment chemistry unknown where unmeasured.

## 2. Restore BIO5 as bounded evidence, not discard it or call it universal

Repository source: `docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md` and its frozen result JSONs.

- The original within-cohort BIO5 analysis included 281 species: median white-minus-nonwhite contrast +0.0690113 within-species SD, species-level Wilcoxon p=0.0118162, three-predictor Holm p=0.0354487. The conditional estimate was OR=1.073475 (95% CI 1.029395–1.119441), p=0.0009188. Retain this positive finding.
- The observer-paired analysis had OR=0.786636 (95% CI 0.543940–1.137618), p=0.202310; discovery/validation species-level transport p=0.743252/0.054067. These designs change both conditioning and eligible observations, so neither automatically overwrites the other. They do prevent a universal signed BIO5 claim.
- The unsigned continuous-colour BIO5-distance analysis is a different estimand from warm-to-white replacement. Its original positive components do not overcome the failed flower-minus-background and strict same-observer controls.

WorldClim BIO5 is the maximum temperature of the warmest month, derived from monthly climate summaries; WorldClim 2.1 represents 1970–2000. It is not the temperature of the photographed flower on the observation day, an experimental heat treatment, or the magnitude of cold/UV/water stress. Source: https://www.worldclim.org/data/bioclim.html and https://worldclim.org/data/worldclim21.html . Never substitute BIO5 for growing-season temperature in an external study without documenting the distinction.

The proposed analysis asks whether a **molecular-route by environment interaction** explains heterogeneity in within-species morph-environment associations beyond geographic, observer and developmental alternatives. An unmeasured route cannot be imputed from the sign of the BIO5 slope. Failed transport itself is not positive evidence for this interaction.

## 3. Independent primary evidence: multiple routes, different endpoints

| System | Molecular/biochemical evidence | Environmental or maintenance evidence | Ceiling |
|---|---|---|---|
| Parrya nudicaulis | White petals show roughly 24-fold lower CHS expression, with no corresponding expression change in leaves/sepals; substantial petal flavonoid loss | White-morph frequency rises with mean July temperature across 14 populations (rho=0.659, p=0.016) | Geographic association, not a temperature-manipulated fitness test; causal variant unidentified |
| Silene littorea | F3h1 downregulation and candidate Myb1a regulation in one molecular study; causal mutation not established | Independent tissue profiling distinguishes flower-restricted vs whole-plant anthocyanin loss; floral-restricted whites reach 8–21% in two repeatedly surveyed populations, whole-plant whites remain below 1% where observed | Frequency contrast is not a randomized causal effect of tissue restriction; two reports of the same system are not independent species |
| Ipomoea purpurea | Different A/CHS and W/regulatory routes can yield white flowers within the same species | In an experimental CHS-null comparison, mutant male and female fertilization success is reduced at high, but not low, temperature; a field study finds reduced survival opposing a transmission advantage | Fertility, survival and natural morph frequency are different endpoints; do not pool their signs as one BIO5 effect |
| Gymnadenia rhellicani | GrMYB1 protein-truncating alleles, downstream ANS expression, pigments and RNAi evidence are linked | Opposing bee/fly preferences and higher seed set of the heterozygous red morph support a particular overdominance mechanism | One population/system does not establish universal maintenance, no universal BIO5 response is tested, and transcription-factor coding mutations are not synonymous with cis-regulatory mutations |

Primary sources:
- Dick et al. 2011, DOI 10.1371/journal.pone.0018230.
- Casimiro-Soriguer et al. 2016, DOI 10.3389/fpls.2016.00204.
- Del Valle et al. 2019, DOI 10.1186/s12870-019-2082-6.
- Tiffin et al. 1998, DOI 10.1266/ggs.73.105.
- Coberly & Rausher 2003, DOI 10.1046/j.1365-294x.2003.01786.x; 2008, DOI 10.1111/j.1558-5646.2008.00355.x.
- Kellenberger et al. 2019, DOI 10.1038/s41467-018-07936-x.

The positive warm-region white association in Parrya and the high-temperature fertility disadvantage of an Ipomoea mutant are **not directly comparable effect sizes**, but demonstrate why visible colour alone cannot specify the relevant physiological mechanism. SlMyb1a and GrMYB1 are not automatically orthologs merely because both have MYB names. Compare sequence-supported orthogroups and regulatory targets, not spelling.

## 4. New executed bridge audit: 4 original overlapping species become 14, but no molecular label transfer

The original 26-row/25-source-species PAL/WAL registry was recovered from GitHub artifact 10855596712. Its CSV SHA256 is `3a19d4717cefbc69b13022982dff9e354b155a898316d84c11016129ed6dbb22`. It is the published Table S1-derived registry, not a newly systematic survey of all molecular cases. Preserve its intervals, qualitative rarity, unpublished-source flags and original taxonomy strings.

The actual measured photo table from artifact 12106604882 has SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`: 4,475 species, 447,500 terminals, 227,107 classifiable records, 3,297 species with at least 40 classifiable photos. The audit recomputed species support and D from those actual rows. It reused only four historically audited aliases: Delphinium nelsonii to D. nuttallianum; Mimulus guttatus to Erythranthe guttata; Mimulis lewisii to Erythranthe lewisii; Silene dioca to S. dioica. No fuzzy matching or new taxonomic acceptance claim was made.

| Natural-context registry records | Present among measured taxa | At least 40 classifiable photos |
|---|---:|---:|
| Flower-restricted anthocyanin-loss (PAL) | 9 | 7 |
| Whole-plant anthocyanin-loss (WAL) | 10 | 8 |
| Unique taxa after cross-class deduplication | 18 | **14** |

The older natural-context bridge comprised 4 taxa. The new 14 include **Ipomoea purpurea in BOTH classes**: the 7+8 rows are not 15 independent species. The greenhouse-only Erythranthe guttata record remains excluded. Parrya nudicaulis and Silene littorea are absent from this 4,475-species table under the documented matching policy.

Most importantly, a published PAL/WAL allele or population in a species does NOT identify the molecular route of arbitrary iNaturalist records of that species. Therefore this audit performs **no PAL-versus-WAL D test and no BIO5-by-route regression**. The correct future unit is a study/population/allele contrast with measured tissues, not a species assigned one permanent pathway class. Use explicit `floral_restricted_loss` and `systemic_loss` labels: the phenotype abbreviation PAL must not be confused with the enzyme phenylalanine ammonia-lyase.

## 5. A minimal mathematical distinction: accessibility is not maintenance

As a proposed two-heritable-type approximation, not a fitted Mendelian/diploid model, write

`dp_j/dt = p_j(1-p_j)[a_j(E,P) - lambda*c_v(E) - delta*c_f(E) + k(1-2p_j)] + m*(p_bar-p_j)`.

Here p is the true frequency of the depigmented type, a is its net floral/reproductive advantage (not assumed positive), lambda is coupling to vegetative costs, delta is loss of additional floral pathway functions, c_v/c_f are environmental costs, k>0 is a phenomenological stabilizing frequency effect, and m is migration. None of these parameters is estimated from photographed colour proportions. Heterozygote advantage requires an explicit genotype model before being equated to k.

At m=0 and k>0, let r=a-lambda*c_v-delta*c_f. Protected interior coexistence requires **-k < r < k**, with p*=(1+r/k)/2. At that root the derivative is -2*k*p*(1-p)<0. With k=0 and constant r there is no restoring force maintaining a protected local polymorphism. Spatially heterogeneous selection with migration can generate range-wide variation without local protection.

Reducing pleiotropic costs does **not** monotonically increase polymorphism: with arbitrary illustrative a=0.5, c_v=1, delta=0, k=0.2, m=0, lambda=0.8 loses white, lambda=0.5 gives p*=0.5, and lambda=0 fixes white. Thus functional decoupling can allow a variant to survive selection while a separate process determines coexistence, fixation or loss. This identity and counterexample were numerically checked, not biologically fitted.

A temperature slope can differ across routes because `dr/dE = da/dE - lambda*dc_v/dE - delta*dc_f/dE`. BIO5 need not have the same sign as a stress cost in all climatic regimes. This gives a falsifiable conditional prediction rather than an unfalsifiable claim that any colour pattern is adaptive.

## 6. Analysis sequence and stopping rules

1. **Retain measurement and prospective evidence.** Keep H1, the frozen W result and original spatial result as distinct phenotype-level observations. They are not genetic pathway tests. Continue floral-organ/developmental validation in parallel; do not make the eight local photo candidates the universe for molecular inference.
2. **Build the mechanism registry independent of effect significance.** Unit: species, population, allele or segregating contrast, tissue, developmental stage, gene/orthogroup, pathway position, coding/cis/trans evidence, pigment chemistry, vegetative retention, fertility/survival outcome and study design. Include null gene tests, unknown tissues and non-persistent mutants. Avoid selecting only published p<0.05 loci. Keep cultivars and interspecific contrasts separate from natural intraspecific contrasts.
3. **Estimate environment-by-route effects only after record-level linkage.** Match actual field populations/genotypes to climatic context. Model species/study-specific slopes and route interactions; do not convert a common mean into universal signed responses. Include predetermined moisture/radiation, geographic and observer/developmental controls. For RNA data use within-study, same-stage contrasts and petal-minus-leaf expression contrasts where paired tissues exist; never pool raw FPKM across species or treat missing leaf data as zero change.
4. **Test maintenance separately.** Use repeated field frequencies and genotype-specific reproductive/survival performance, frequency manipulations where available, and neutral population structure. A route must predict an independently observed outcome beyond neutral geography and sampling, not just reproduce the known colour contrast. Retain unsupported results. If route interactions fail to improve held-out study/species prediction, the proposed generalization is not supported.

Public-data starting points verified in primary-paper data statements: Silene RNA-seq SRP033277 and ABP sequences KT954895–KT954923; Gymnadenia BioProject PRJNA504609, GrMYB1 MK163677–MK163681, GrANS1 MK163682–MK163685, figshare 7314731 and 7321880; Parrya GenBank HQ215216–HQ215513. PRJNA504609 includes more than the focal species, so species/morph/tissue metadata filtering is mandatory. The Silene study has nine morph-by-stage libraries sequenced twice; repeat sequencing is not biological replication. In this session direct ENA/figshare metadata downloads failed due network access; these accessions are **literature-verified entry points**, not downloaded/reanalysed RNA-seq.

### Current end state

Completed: original full-photo/registry hash validation, 26-row crosswalk, 14-species deduplicated support audit, primary-source mechanism comparison, analytic conditions and falsifiable design. Not completed or claimed: new climate regression, raw-read differential expression, causal pathway assignment to atlas photos, universal FCP mechanism, or final revision/merge of the held New Phytologist manuscript.
