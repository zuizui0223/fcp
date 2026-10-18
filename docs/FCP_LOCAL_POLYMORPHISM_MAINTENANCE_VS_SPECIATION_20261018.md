# FCP research pivot: distinguish local polymorphism maintenance from spatial photo-colour structure

**Evidence date:** 18 October 2026. **Status:** source-verifiable new **post hoc falsification**, not a new confirmatory adaptation claim. **No original 199-null/50-km estimate or prospective H2 axis is overwritten.**

## Core biological distinction

"Same species" does not imply identical ecological niches across individuals/populations, and a flower's visible colour is not proof of any ecological niche divergence. Distinct flower-colour labels across a species' photographed range are not sufficient to establish genetically maintained *within-population* flower-colour polymorphism (FCP). A geographic colour cline can emerge through drift, limited migration, variable environments, plasticity, cultivation/observer concentration, or image exposure, with or without selective maintenance. Speciation specifically requires a barrier to reproduction/gene exchange, **not simply different phenotype frequencies**. A photographed morphological difference cannot identify reproductive isolation.

**The prior source-verified FCP result is real and bounded:** In 1,761 50-km spatially eligible species from 447,500 genuine original photo measurements, mean global-minus-local four-state colour dissimilarity was **0.0176149022** (original 199 species-vertex permutations, upper p=0.005). The unchanging colour composition and fixed sampling coordinates make that a test of photographic spatial sorting under global species-wise label exchangeability. They do not distinguish drift from locally divergent selection, show that two states coexist in a breeding population, or estimate morph fitness.

**Stronger question for future FCP:** *When do genetically inheritable alternative flower morphs coexist within an interbreeding population for multiple generations, and when do they instead sort geographically or become reproductively isolated?* This is a causal eco-evolutionary question distinct from the descriptive 50-km atlas.

## Completed within-site/season feasibility audit on the same REAL data

Input: SHA256-verified 447,500-photo measured CSV from the [immutable final original-method five-cohort release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017); original measured CSV SHA256 **8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6**. Only **227,107** classifiable four-state photo labels enter the within-site opportunity audit.

For every species separately, calculate actual geodesic pair distances. A local pair must have **distinct original iNaturalist observation IDs, different identified observers, same calendar month and year, <=5 km recorded separation and positional-accuracy fields <=1 km for each observation**. This is a conservative *photograph-pair opportunity* for local morph-coexistence research, **not proof that two different plants share one breeding population**. Same-day results are stricter in time but not individually audited morph-genotype photographs.

| Eligibility / photographed colour-state evidence | Species |
|---|---:|
| >=1 eligible same-month, different-observer photo pair | **801** |
| >=1 eligible same-month, different-observer different-colour pair | **300** |
| >=5 eligible local same-month pairs | **194** |
| >=30 eligible local same-month pairs | **14** |
| >=1 eligible same-day local pair | **258** |
| >=1 eligible same-day different-colour pair | **57** |

Observed same-month eligible pairs total **3,485**, including **726** pairs labelled as different coarse colours. These pair counts are dependent photographic comparisons, **not independent individual plants**, true frequencies or an experimental measure of stable morph coexistence. Simple morph-pair presence can result from photographic misclassification; photo revisits of one plant or nearby cultivars must be excluded during field validation.

## New falsification: does extra **fine-local** photographic sorting survive conditioning on season/year?

Condition each species' photo labels on the **observed full year-month-specific colour-state multiset** (rather than a single global species-wide colour composition). Shuffle only the photo labels **within each species and calendar year-month** 199 times while leaving all photo positions, observer IDs, timestamps and the month-constrained local-pair graph unchanged. The data-driven support gate is based solely on local pair geometry/date/observer identifiers, independent of colour label.

Statistic: equally weighted species mean of \`E_month-null(local dissimilarity) - observed local dissimilarity\`; one-sided upper p uses the 199 matching monthly null worlds. This is a NEW, explicitly post hoc stricter null, not the originally frozen 50-km test. It addresses whether 5-km local colour similarity remains after controlling species-wide month-by-year composition, but **not the spatially autocorrelated neutral/genetic null**.

| Strict support | Eligible species | Observed local pair dissimilarity | Matched month-null mean | Mean additional local similarity | p (199 month-conditioned shuffles) |
|---|---:|---:|---:|---:|---:|
| >=5 strict pairs | **194** | **0.18927** | **0.19314** | **+0.003870** | **0.300** |
| >=30 strict pairs | **14** | **0.27342** | **0.28287** | **+0.009450** | **0.275** |

**Result: the new small-scale month-controlled colour-sorting contrast was not supported at either support gate.** This is not evidence of no within-population polymorphism or proof that season explains the entire original 50-km geographic effect. Only 194 species have >=5 independent nearby photo pairs; the stricter test lacks broad geographic/replicate representation. Calendar month does not remove within-month flowering, camera exposure, repeat-plant or observer-specific detection biases; 5 km is not a scientifically delimited population.

### Status of the original New Phytologist paper

The existing 447,500 real-photo original null/W/H1 results remain valid statements about **photographed species-level composition**. Neither the strongly skewed most common photographed state nor the original 50-km local homogeneity test alone establishes biologically maintained FCP. A synthesis that calls such effects "maintenance of balanced polymorphism" or "ecological speciation" is overstated. The original prospective W result identifies an achromatic–chromatic *geometry of photographed variation* conditional on frozen coarse states, not a direction of genetic loss/reacquisition or reproductive isolation.

## Mainline theory: three processes and an actual branching point

Use a minimal two-morph, multi-locality evolutionary model as a falsifiable mechanism **in a separate simulation**, not to fit unidentifiable selection/genetic parameters directly from photos:

\[
\frac{dp_j}{dt}
 = p_j(1-p_j)[s_j+k(1-2p_j)] + m(\bar{p}-p_j)+\eta_j(t).
\]

- \`p_j\`: the true **genetic** morph frequency in population j (not the photo-colour proportion).
- \`s_j\`: the signed local fitness advantage; may vary with water, soil, pollinators or microhabitat.
- \`k>0\`: a simple negative frequency-dependence parameter favouring rarity, which can maintain **local sympatric polymorphism** without a differentiated niche (this is a toy model, not a claim about FCP).
- \`m\`: gene flow among populations. Spatially opposite selection plus gene flow can maintain **species-wide mosaic variation** even where individual populations are nearly fixed.
- \`eta_j(t)\`: drift/environmental stochasticity. Strong spatial differences can appear even without adaptive selection.

An additional reproductive-isolation process must be measured before using a **speciation** label. Flower-morph divergence alone, even between different ecological habitats, is insufficient.

Predictive observation/experiment distinctions:

| Hypothesis | Genetic local coexistence | Morph-specific habitat/fitness | Gene flow/reproductive barriers |
|---|---|---|---|
| Photo artefact / phenotypic plasticity | Not established; same genotype may change colour | May correlate with camera/condition or environment | No required barrier |
| Neutral spatial drift/colonisation | Randomly persistent or replaced; local frequencies drift | Colour clines no more distinctive than neutral loci | General spatial IBD, no morph-specific barrier |
| Divergent ecological selection with gene flow | Often local dominance, contact-zone mixtures | Replicated reciprocal morph advantage by habitat, beyond neutral spatial structure | Morphs still interbreed |
| Protected local FCP (frequency dependence / heterozygote advantage) | Multiple inherited morphs recur across local generations | Rare morph advantage/heterozygote fitness in replicated sites | Interbreeding maintained, possibly disassortative crossing |
| Ecological speciation / incipient species | Morphs may coexist or segregate | Divergent ecological response plus reproductive linkage | Direct assortative pollination, reduced cross-fertilization/hybrid fitness and genomic isolation beyond neutral geography |

**Decisive tests beyond the iNaturalist atlas:** replicated individual-marked flowering censuses in at least two years; same-population inherited morph/genotype concordance in a common garden; morph-specific seed/pollen fitness as a function of morph frequency and environment; neutral MIG-seq/SNPs versus colour-associated genetic clines; pollen flow and hand-cross compatibility between colour states. The specific point at which colour differentiation becomes *speciation* is where sustained reproductive isolation of morph-associated lineages is demonstrated, not where the photograph colours cross a spatial-distance threshold.

## Sources that motivate, but do not substitute for, these tests

- Narbona et al. (2018), *Plant Biology*, DOI 10.1111/plb.12575: distinguish within-population flower-colour polymorphism from spatially varying, locally monomorphic morph mosaics.
- Schemske & Bierzychudek (2007), *Evolution*, DOI 10.1111/j.1558-5646.2007.00219.x: *Linanthus parryae* colour clines supported local selection **only** after neutral markers and reciprocal-transplant contrasts.
- Bolnick et al. (2003), *American Naturalist*, DOI 10.1086/343878: individual niche variation can occur within species.
- Rundle & Nosil (2005), *Ecology Letters*, DOI 10.1111/j.1461-0248.2004.00715.x: ecological speciation requires ecological divergent selection **linked to reproductive isolation**.
- Gigord et al. (2001), *PNAS*: experimental rare-morph advantage in *Dactylorhiza sambucina* exemplifies a direct frequency-dependent maintenance test.

**Research directive:** preserve the negative fine-scale/month-conditioned null and the limited candidate sampling coverage. Do not optimize thresholds until it becomes significant. The FCP follow-up should predefine local sampling populations, observer/botanical review, power and gene-flow/fitness evidence **before** any upgrade from photographed spatial structure to genetic polymorphism maintenance.
