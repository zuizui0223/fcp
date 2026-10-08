# FCP cross-species geographic organization after flowering-calendar adjustment
## 2026-10-08 post-outcome analysis contract — main ecological question

### Why this route, rather than more one-species anecdotes

The FCP paper concerns repeated, geographically allocated flower-colour ITV in 149,900 photographs / 1,499 sampled species across three species-disjoint high-depth resources (500+500+499), not merely the identity of some photo-multimorph plants. Its strongest robust ecological quantity is fixed-species-composition local colour depletion at 50 km (discovery/validation/third +0.02053/+0.01867/+0.01468, matched 199-permutation P=0.005 each). This geographic pattern already survives excluding white, excluding same-observer pairs, using continuous nine-colour distance, quarter-conditioned calendar nulls, and cross-year photo pairs.

The primary remaining *comparative ecology* question is:

**Does geographic colour sorting remain detectable across many different plant lineages when calendar-season photo turnover is preserved, or can the world-scale photographic pattern be explained merely by photographing different visible morphs in different seasons?**

Prior art: Narbona et al. (2018, Plant Biology, DOI 10.1111/plb.12575) described geographically mosaic morphs in Mediterranean flower-colour-polymorphic systems, so we do not claim that discovering local monomorphic vs polymorphic sites is new. Sapir, Gallagher & Senden (2021, TREE, DOI 10.1016/j.tree.2021.01.011) review different balancing/nonselective processes. Our distinct contribution is a standardized, repeated, species-comparative, measurement-conditional and opportunity-controlled assessment of spatial partition, not proof of local adaptation.

### Cohorts and biological estimand

- Original 100,000 photos across 1,000 sampled species, split into discovery and validation groups, plus a fresh 49,900-photo / 499-species prospective image resource.
- Reopen only source-frozen measured tables at exact SHA256; never replace missing photos or change the established New Phytologist confirmatory H1/H2 verdicts.
- Eligible species: >=40 classifiable photographs with parseable dates (1990–2026), >=30 conspecific photo pairs within 50 km. Stratum-conditional modes require at least 10 image positions belonging to within-stratum photo sets that have >1 photo and >1 colour state (ensures a nondegenerate restricted null).
- Main *species-level* estimand for species i: Delta_i = D_all,i - D_local,i, where D is pairwise coarse visible-colour discordance. Compare Delta_i to the conditional expectation from within-species, within-time-stratum exchangeability. Compute equal-species mean excess over this restricted null for each of the three cohorts separately. A prospective probability sample of all flowering plants is NOT claimed.
- Every photo colour category is held fixed in total within each species; no white/pink/blue observations may be introduced, substituted or inferred.

### Null tiers, stronger phenology tests, sampling rules

1. Unconditional species-fixed morphology label shuffle (historical main comparator), included as a numerical calibration.
2. Calendar quarter (January–March, April–June, etc.) label shuffle within species, *pooling observation years*. This is an older, previously observed positive-robustness route and therefore not a new confirmatory test.
3. Calendar month within species across years (**primary added retrospective environmental calendar test**): retains each species' exact observed monthly colour-state counts and which photographs fall in each month, while randomizing only geography/observation-specific assignment *within* those month groups. If species-wide spatial organization were purely due to phenological month structure, no extra geographic colour depletion should appear beyond this null when it is adequately exchangeable.
4. Year × month within species: a deliberately stricter check. Much less within-block label exchangeability is expected with only 100 images/species. Report coverage-limited or NOT ESTIMABLE when there are insufficient exchangeable strata; do not convert lack of power to evidence for any null.

For *every* null tier, also analyze (A) all local conspecific photo pairs and (B) local pairs where both observer IDs are nonempty and distinct. Keep the same 50 km maximum local pair distance and >=30 local-pair minimum. Date filtering must happen before calculating both the species overall composition and its local diversity.

Use 199 fixed-seed label permutations per stratum/cohort/species, preserving complete species colour-class frequencies and exact month or month×year colour-class composition, local geographic positions and observer IDs. Aggregate null statistics across species with equal weights. 1,999 species bootstraps yield uncertainty for the mean excess. Genus-balanced mean and leave-one-genus-out range are secondary broad lineage sensitivity, not a phylogenetic correction.

### Explicit generality gates

No cross-species rule can be claimed if any cohort has <80 exchangeable and geographically evaluable species for the primary month-conditioned mode. An observational monthly-partition generality candidate requires, in **all three species-disjoint cohorts**, both local-pair policies to show:
- >=80 estimable species,
- positive mean local depletion beyond the exact month-preserving null,
- permutation upper-tail P<0.05.

Failure in one cohort or one policy must be reported as nontransport/limited generality. Passing does not establish a universal evolutionary driver or percentage of *all* angiosperm species. The three cohorts share the same provider and photographic measurement algorithm; this is within-system replication only, not independent external-source validation.

Report the conditional null expectation, mean additional observed spatial depletion, species-bootstrap CI, genus-balanced effect, leave-one-genus-out range, positive species fraction, number of species/genera, date coverage, effective within-month exchangeability, and mode status for *every cohort and mode*. No favourable subgroup selection after opening results, no family-based post hoc mechanistic interpretation without a distinct test.

### What the result can and cannot decide

If supported: repeated photographs show a general **spatial organization of visible flower-colour ITV beyond observed calendar-month composition** across many plant species within the measured opportunity frame. This would distinguish geographic allocation from one simple phenomenological confounder, but does not establish local adaptation, balancing selection, morph genotype, pollinator choice, pigment protection, flower-sex/organ identity, or morph emergence rates.

If unsupported: the data cannot separate photo-level spatial from seasonal organization under this permutation design. A zero effect is not proven and does not erase the pre-existing 50 km unconditional/quarter-conditioned result. A year-month constraint that leaves no exchangeable colour/locations is a coverage failure, not a biological refutation.

*Moricandia arvensis* remains a biologically verified prior-art example of nonadaptive, summer-induced within-individual flower-colour plasticity (Gómez et al. 2024, Evolution Letters, DOI 10.1093/evlett/qrae017). It warns that temporal vs spatial co-observations do not identify selection. *Silene littorea* tissue-limited vs whole-plant anthocyanin loss remains a source-derived biochemical benchmark, not a directly measured molecular covariate for all FCP species.

The mechanism question for a subsequent independent-source comparative study is not "which single latitude direction causes pigment?" but "how are colour morphs retained locally despite costs, versus repeatedly segregated among sites, and when is colour itself the target of selection rather than correlated plasticity?" Test this with verified wild populations, direct chemistry/genetic segregation, and reproductive fitness, not just colour-associated climate.
