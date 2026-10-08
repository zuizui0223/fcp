# Global FCP mechanism maintenance: primary-source to photo-cohort feasibility
## 2026-10-08; retrospective, purposive source panel

### The biological inference problem

FCP's 149,900 photos of 1,499 selected species and 369/363/377 species clearing the primary colour measurement gate demonstrate a recurring geographic allocation of intraspecific photographed colours. That result is NOT yet the same as a particular evolutionary mechanism maintaining genetically stable within-population morphs.

A primary-study source audit distinguishes five experimentally established or well-supported routes:

1. Correlated-trait selection plus reproductive assurance — Silene littorea: Rodríguez-Castañeda et al. 2020 (DOI 10.3389/fpls.2020.588383) found an indirect flower-number-mediated corolla-colour effect on seed fitness +0.418 (p=0.022) in Barra and -0.257 (p=0.064) in Melide. The second estimate is borderline, not a significant population contrast. Buide et al. 2021 (DOI 10.1111/plb.13209) found pollinators preferred pink flowers while white flowers compensated via autonomous selfing (white vs pink autonomous fruit percentages 27 vs 20; seed percentages 20 vs 12). Del Valle et al. 2019 (DOI 10.1186/s12870-019-2082-6) independently classified tissue-restricted versus whole-plant anthocyanin loss.
2. Pollinator-mediated negative frequency dependence — Dactylorhiza sambucina: Gigord et al. 2001 (DOI 10.1073/pnas.111162598) used artificial morph-frequency arrays and found rare-colour fitness advantage in male/female reproductive components. This is not a white-flower metabolic production cost.
3. Neutral history and migration — Iris lutescens: Wang et al. 2016 (DOI 10.1093/aob/mcw036) recorded 1,120 plants from 41 populations with eight microsatellite loci, implicating drift/limited gene flow in geographically monomorphic Spanish populations, and gene flow among mixed French populations. Selection was not excluded.
4. Stress-mediated correlated traits — Boechera stricta: Vaidya et al. 2018 (DOI 10.1111/nph.14998) measured reduced herbivore damage in purple forms, environmentally inducible colour and a white reproductive advantage under well-watered glasshouse conditions. This is not a demonstration of universal anthocyanin synthesis resource costs.
5. Context-dependent habitat/fitness and plasticity — Castilleja coccinea (DOI 10.1002/ajb2.70094), Abronia fragrans (DOI 10.1002/ajb2.70142), Moricandia arvensis (DOI 10.1093/evlett/qrae017), and Leptosiphon parviflorus (DOI 10.1002/ajb2.70018) collectively distinguish habitat-dependent red/yellow bract variation, white/pink betalain-related display and antagonism, nonadaptive seasonal colour plasticity tied to flowering time, and soil/climate phenotypic correlation without direct fitness. These effects cannot be combined as independent estimates of the same selection coefficient.

The panel also includes Raphanus sativus and Ipomoea purpurea as historically known experimental floral-colour systems. It is PURPOSIVE, not a random or complete systematic sample; numerical class counts are not the worldwide prevalence of maintaining mechanisms.

### The empirical crosswalk question

Which of the source-verified organism systems are actually represented by the original FCP photo cohorts at >=40 classifiable photos? Among those, how many have at least five photographically white records and five nonwhite, and which have independent evidence mapping the **same photos** to an authenticated pigment genotype, a natural population, or an experimentally measured fitness outcome?

Source: data/curated/fcp_morph_maintenance_primary_studies_20261008.json.
Three measured photo SHA256 checks:
- Discovery ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4.
- Validation 0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6.
- Third 57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186.

Species names are compared only by their EXACT genus+species binomial. Aliases, near matches and spelling repairs require a separately audited taxonomic crosswalk, not the current task. Species appearing in more than one nominally disjoint photo cohort cause the test to fail closed. Every retained match reports photographic counts, photo-coarse white/chromatic counts, >=40 gate, >=5/arm gate, local <=50km photo pair information and literature evidence type.

### Inference ceiling and decision rule

Even a species present in two evidence sources cannot link the photo-colour state to the paper's genetic pigment phenotype and fitness simply from a species name. Published examples include distinct petals versus bracts, white/pink betalain versus anthocyanin, and seasonal colour from a single genotype. Individual photos are not associated with a field-verified mating population or survival/seed fitness measurement.

A general causal crosswalk requires at least six independent literature systems eligible for photographic white/chromatic variation AND at least three independent true-image-to-chemical/genotype/fitness links. None of the original image tables contains a validated linkage; therefore the second requirement is currently not satisfied. The audit reports direct biological bridge NOT ESTIMABLE without converting absence of sampling into absence of tradeoffs.

Do not sum evidence across noncommensurate fitness variables, infer the rate of balancing selection, relabel photographic white as anthocyanin loss, or rewrite the paper's frozen 50-km, H2, latitude/elevation or photographic-source results.

### Delivered analysis

The audit script scripts/analysis/run_polymorphism_mechanism_photo_overlap_20261008.py checks primary source identity and source-photo SHA256, summarizes selected systems in the three original photo cohorts, and writes result.json and study_system_photo_overlap.csv. Synthetic tests check taxon-name misidentification, photo vs pigment/genotype separation and cohort disjointness. This is a causal feasibility diagnosis, not another post hoc ecological fitness test.

### Source-morph semantic negative control (added before final execution)

Even exact species-name matches are not comparable if the original study concerns purple-vs-yellow forms, red-vs-yellow *bracts*, or within-genotype summer white plasticity. The registry therefore carries a separate, source-derived boolean for whether the focal published morph contrast includes white, whether genetically based white variation was demonstrated in that study, and whether white is expressed as a seasonal single-individual state. Every matched source photo row publishes these indicators alongside unreviewed photo-visible white counts. **A published white form plus a white-photo label STILL does not establish the photographed individual's biochemical or genetic identity.** Fail closed on absent source-phenotype declarations. Neither missing white in the focal study contrast nor a photo-white label can prove a source population lacks a rare white morph; direct reinspection would be needed.
