# Flower-colour polymorphism is geographically organized within species without a shared global boundary

**Running title:** Geography of flower-colour polymorphism  
**Article type:** Research Article  
**Provisional target:** *Journal of Biogeography*  
**Version:** v0.1-post-Step-9+acquisition-provenance, 2026-09-11
**Numerical authority:** `paper_numbers.json` / `PAPER_NUMBERS.md`

## Abstract

### Aim

Flower-colour research commonly asks which colour state evolves, which environments favour particular colours, or whether similar transitions recur among lineages. A different property has received less comparative attention: how much flower-colour variation a species carries, and whether that variation is geographically organized within the species. We tested whether species-level flower-colour polymorphism predicts within-species geographic colour organization without requiring species to share the same colour-transition direction or global boundary.

### Location

Global community-photograph sample.

### Taxon

Angiosperm species represented in fixed-photo discovery and species-disjoint reserve frames.

### Methods

Each candidate species entered image measurement with a fixed 100-photo raw denominator. We represented admitted flower colours as four broad states—white, yellow-orange, red-pink and blue-purple—and quantified species-level diversity with the Gini–Simpson index, `D = 1 - sum(p_k^2)`. The discovery frame contained 369 species and the species-disjoint reserve contained 363 species meeting the frozen classifiable-observation threshold. Within each species, geographic colour organization was quantified from the association between geographic separation and flower-colour dissimilarity. Spatial nulls permuted colour structure while retaining each species' observed sampling geometry. We then tested replication, observer and calendar sensitivities, matched flower-minus-background structure, sampled geographic span, clearly technical ROI/flip failures, and exact minimum- and maximum-diversity completions of ambiguous-palette observations under a four-state sensitivity model.

### Results

Flower-colour diversity varied continuously in the discovery frame (`D = 0–0.7076`); a second colour state contributed at least 10% of admitted observations in 46.6% of species. Greater polymorphism was associated with stronger within-species geographic colour organization in discovery (`rho = 0.0892`, `P = 0.034`) and independently in the reserve (`rho = 0.1016`, `P = 0.025`). In the reserve, the relationship remained after adjustment for sampled geographic span and clear technical-failure rate (`partial rho = 0.0993`, `P = 0.025`) and was stronger for the matched flower-minus-background differential (`partial rho = 0.1162`, `P = 0.010`). The spatial result also survived the exact diversity-minimizing and diversity-maximizing uniform completions of ambiguous-palette observations: reserve primary `P = 0.029` and `0.008`, respectively, and reserve flower-minus-background `P = 0.009` and `0.006`. Genus-level taxonomic clustering recurred across tranches but was less robust to ambiguity completion.

### Main conclusions

Species carrying more flower-colour polymorphism are, on average, more geographically organized internally, even though species need not place their colour transitions in the same geographic regions or vary along one privileged colour direction. The repeatable feature is therefore the **degree of species-specific geographic organization**, not a universal global flower-colour boundary. Explicit separation of technical failure, background spatial structure and unresolved colour ambiguity shows that this result is not an artefact of any one of these measured processes.

**Keywords:** citizen science; colour polymorphism; flower colour; geographic variation; intraspecific variation; random labelling; spatial organization

---

## Introduction

Intraspecific phenotypic variation is often treated as residual variation around a species mean, but it can itself be a structured biological property. Alternative phenotypes may coexist locally, replace one another among populations, form clines, or occur in combinations of these spatial states. These configurations matter because the processes capable of maintaining variation within populations are not identical to those that generate differentiation among populations. Balancing selection, negative frequency dependence and opposing selective agents can maintain local diversity, whereas spatially varying selection, restricted migration, demographic history and drift can contribute to geographic differentiation (Delph & Kelly, 2014; Narbona et al., 2018). A comparative analysis of intraspecific variation therefore needs to distinguish not only *which* phenotype occurs, but also *how much* variation exists and *how that variation is arranged in space*.

Flower colour is particularly useful for making this distinction. Colour variants are conspicuous, are recorded extensively in community-photograph platforms, and can reflect multiple interacting ecological and evolutionary processes. Pollinator behaviour can impose colour-dependent selection, but evidence for pollinator-mediated selection is heterogeneous and cannot explain all flower-colour variation (Narbona et al., 2018; Trunschke et al., 2021). Abiotic conditions, herbivory, correlated traits, gene flow, drift and demographic history can also affect pigmentation or morph frequencies. Classic and recent studies show that colour morphs can differ geographically and ecologically within individual species, including plant systems such as *Linanthus parryae* and *Boechera stricta* and, more broadly, colour-polymorphic animals analysed from large photographic datasets (Schemske & Bierzychudek, 2001, 2007; Vaidya et al., 2018; Farquhar et al., 2023; Jansen et al., 2025). At macroecological scales, flower-colour distributions are also associated with both biotic and abiotic conditions (Dalrymple et al., 2020). These literatures establish that colour can be spatially structured, but they do not imply that independent species should share one spatial boundary or one colour-transition direction.

Three levels of variation are easily conflated. The first is the **amount** of within-species colour diversity: a species can be nearly invariant or can contain multiple common colour states. The second is the **direction** of variation in colour space: different species can reach similar levels of diversity through different combinations of colour states. The third is the **geographic organization** of those states: colour differences can be spatially random, locally mixed, or increasingly different with geographic separation. Repeatability at one level does not require repeatability at the others. In particular, two species may both show strong geographic organization while placing transitions in different regions and varying along different colour axes. A universal boundary is therefore a much stronger hypothesis than repeated species-specific organization.

Large community-photograph collections make it possible to compare these levels across hundreds of species, but they create their own inferential problems. Observation geometry differs among species; photographs contain spatially structured backgrounds; automated flower localization and colour measurement can fail; and some flower regions contain colour compositions that cannot be safely assigned to a single discrete morph. Treating all unclassified images as generic missingness can be especially misleading because classification success occurs *after* image acquisition and may itself depend on the visual complexity of the flower. A credible comparative test therefore requires fixed photo budgets, species-conditioned spatial nulls, an independent reserve, explicit background controls and a measurement model that keeps clearly technical failures separate from biologically ambiguous colour observations.

Here we use two fixed-photo, species-disjoint global frames to ask three questions. First, does four-state flower-colour polymorphism behave as a continuous species-level property rather than a rare binary exception? Second, are more polymorphic species more strongly geographically organized internally? Third, is species-level polymorphism taxonomically structured, and is that pattern as robust to colour ambiguity as the spatial result? We deliberately treat common boundary geography, recurrent colour direction, latitude, pollination and life form as distinct or secondary questions rather than folding them into one omnibus explanation. This design allows a sharper test of a simple proposition: **species may differ consistently in how much colour variation they carry and how geographically organized that variation is, even when the location and direction of differentiation are not shared among species.**

---

## Materials and Methods

### Fixed-photo discovery and reserve frames

The analysis was built from georeferenced community photographs assembled under a fixed-photo design. Both analytical tranches trace to the same frozen iNaturalist v1 observation candidate pool. The pre-pixel acquisition contract required `quality_grade=research`, photographs and georeferences, species-rank taxa, a frozen flowering-annotation gate (term 12, value 13), positional accuracy no worse than 5 km, unobscured coordinates and one of five allowed photo licences (`cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc` or `cc-by-nc-sa`). Candidate-page and candidate-species selection did not use flower colour, and candidate image pixels remained unopened during acquisition. Within each species, an observer contributed at most two photographs and the final 100 photographs were chosen by deterministic geographic maximin selection.

The frozen candidate pool contained exactly 1,000 species with 100 raw photographs per species. The original fixed hash-ranked 500-taxon measurement budget defined the discovery tranche; the reserve was the entire species-disjoint complement of the other 500 taxa from the same frozen candidate pool. The reserve constructor replays the original discovery budget against the historical Git object and verifies the candidate-pool SHA-256 before taking that complement. Thus discovery and reserve differ in species membership but inherit the same frozen acquisition frame. Each candidate species therefore entered image measurement with exactly 100 raw photographs. The fixed denominator was chosen to prevent common, heavily photographed species from dominating species-level comparisons simply because more images were available. Species inclusion in the final analytical frames was conditional on the upstream source, acquisition and measurement gates; the resulting species sets are not random samples of all angiosperms and are not used to estimate angiosperm-wide prevalence.

The frozen acquisition query did **not** apply a native-range restriction or explicit `captive=false`/`wild=true` parameter, and the local metadata parser did not add such a filter. We therefore use Research Grade only as the iNaturalist quality-grade criterion and do not interpret it as proof of native-range-only or exclusively wild-population sampling. Geographic estimands in this paper describe organization among the observed community-photograph records.

The discovery polymorphism frame contained 369 species with at least 40 observations admitted to the four-state colour representation. A separate reserve campaign supplied 363 eligible species and was species-disjoint from discovery. The reserve was used to test transfer of the already defined polymorphism and spatial patterns rather than to search for a new threshold, colour axis or predictor set.

### Flower-region measurement and four-state colour representation

Image measurement followed the frozen flower-region pipeline used in the global photo campaign. The operational representation reduced admitted flower observations to four broad visible colour states: white, yellow-orange, red-pink and blue-purple. These classes are image-derived visible-colour categories rather than spectrophotometric pigment measurements, biochemical states or claims about pollinator visual systems.

Every raw photograph was assigned one measurement outcome. In the 50,000-photo reserve candidate frame, 24,885 observations (49.77%) were admitted to a four-state morph, 22,681 (45.36%) failed the automated ROI/flip stability gate, 2,415 (4.83%) were retained as ambiguous palette compositions and 19 (0.038%) contained no admissible biological palette mass. The same conceptual distinction was preserved in later robustness analyses: ROI/flip-gate failure was treated as clearly technical; ambiguous-palette observations were not silently recoded as technical failures or hidden biological morphs.

### Species-level flower-colour polymorphism

For each species, let `p_k` be the frequency of admitted observations in colour state `k`, for `k = 1,...,4`. Species-level flower-colour diversity was quantified with the Gini–Simpson index

`D = 1 - sum(p_k^2)`.

`D = 0` indicates that all admitted observations occupy one colour state, whereas larger values indicate more even representation of multiple states. We also retained the conventional finite-sample correction

`D_unbiased = n/(n-1) * D`,

where `n` is the number of classifiable observations. The continuous `D` value is the main species-level response. For descriptive interpretation and pre-frozen subset comparisons, we additionally recorded the frequency of the second-most common state and used thresholds of 10% and 20%. These thresholds summarize the presence of a non-trivial secondary state but do not replace the continuous estimand.

### Test of a recurrent colour direction

Before emphasizing geography, we tested whether increasing polymorphism was preferentially expressed along a common learned colour direction. A frozen reference ordination defined two leading axes, M1 and M2, and the primary statistic asked whether species with larger `D` concentrated more variance in the corresponding recurrent two-dimensional subspace. Direction specificity was tested against random two-dimensional orientations. The purpose was to distinguish a general increase in colour spread from expansion along one privileged cross-species colour axis. Because the preregistered directionality claim was not supported, these analyses are used as a negative contrast rather than as a main positive result.

### Within-species geographic colour organization

Each species received a spatial statistic describing whether flower-colour dissimilarity increased with geographic separation among its observations. The analysis retained each species' own observed locations and sampling geometry. For the frozen spatial null, colour structure was randomized within species while coordinates and the species-specific geometry remained fixed. Thus the null asks whether the observed arrangement of colour carries more geographic structure than expected after conditioning on where that species was sampled; it does not compare a species with a spatially uniform world or with other species' ranges.

The global spatial statistic aggregates species equally so that a species with more spatial pairs does not automatically dominate the result. A positive species-level spatial `rho` means that observations farther apart tend to be more dissimilar in flower colour within that species. This statistic measures the *degree* of within-species geographic organization. It does not require species to share a boundary location.

The discovery spatial campaign preceded the present polymorphism-focused reanalysis. We therefore distinguish the original whole-frame spatial test from the post-outcome exploratory comparison of spatial strength against `D`. The independent reserve provides the stronger replication test of the polymorphism–spatial relationship.

### Species-disjoint reserve replication and nuisance controls

In the reserve, we reused the frozen 10% second-morph threshold and continuous `D` gradient. Four spatial responses were retained: the primary statistic, an observer-pair exclusion sensitivity, a calendar-quarter stratification sensitivity and a matched flower-minus-background differential. The background differential asks whether spatial structure in the flower exceeds spatial structure recoverable from the paired image background.

The reserve spatial nulls were not regenerated for the later robustness analyses. Instead, the same 999 within-species spatial randomizations were propagated through the subsequent partial-rank statistics. This preserves each species' observed coordinate geometry while testing whether the relationship between polymorphism and spatial organization exceeds the geometry-preserving null.

### Sampled geographic opportunity and measurement-process controls

We used sampled geographic span as an opportunity variable, not as an estimate of biological range size. In discovery, the frozen primary proxy was `log1p(maximum_span_km_after_observer_cap)`. In the reserve, sampled span was the log-transformed maximum pairwise great-circle distance among all 100 fixed photo coordinates for each species, calculated before filtering by morph or classifiability.

An initial robustness analysis also conditioned on `n_classifiable`. Subsequent measurement auditing showed why this quantity is not a clean sampling-effort covariate: every species started from the same 100-photo denominator, whereas `n_classifiable` is generated after ROI and colour-composition gates. We therefore separated the clearly operational ROI/flip failure rate from ambiguous-palette observations and repeated the partial-rank analyses controlling sampled span and **technical-failure rate** rather than broad classification yield. This is the primary measurement-process adjustment reported here.

### Ambiguous-palette endpoint sensitivity

Ambiguous-palette rows were not assigned biological morph labels in the primary analysis. To assess how strongly these unresolved rows could alter species rankings, we constructed exact species-level diversity bounds under a deliberately limited sensitivity model: every ambiguous row was assumed, counterfactually, to belong to one of the same four existing colour states.

For a species with observed class counts `c_k` and `A` ambiguous rows, the diversity-minimizing endpoint `D_min4` assigns all `A` rows to the currently most abundant state. The diversity-maximizing endpoint `D_max4` distributes the `A` rows among the four states so as to equalize their final counts as far as integer constraints allow. These endpoints are exact for species-level `D` under the four-state completion assumption.

We then applied each endpoint rule uniformly to all species and repeated the span-plus-technical-failure spatial tests using the original geometry-preserving spatial nulls. These are **uniform endpoint stress tests**. They do not exhaust arbitrary species-specific adversarial allocations and do not establish that ambiguous observations are true four-state morphs.

### Species attributes and genus-level taxonomic clustering

A prospectively frozen six-family species-attribute analysis tested sampled geographic span, absolute latitude, family-level clustering, genus-level clustering, pollination guild and life form. Pollination and life form did not pass their pre-outcome coverage gates and were assigned `P = 1` in the six-slot Holm family rather than being replaced after outcomes were opened. Of the active tests, sampled span and genus-level clustering were Holm-supported in discovery; latitude and family-level clustering were not.

For genus-level clustering, only genera represented by at least two species contributed. For each repeated genus, all within-genus species pairs were enumerated, each genus received equal total weight, and the statistic `W` was the weighted mean absolute difference in `D` between congeneric species. Smaller `W` means greater within-genus similarity. Species labels were permuted to obtain a lower-tail null. We report clustering gain as `1 - W_observed/W_null_mean`. This is a taxonomic clustering analysis, not a formal phylogenetic comparative test.

We repeated the same statistic in the species-disjoint reserve, compared mean `D` among the 23 genera represented by at least two different species in both tranches, and assessed sensitivity to sampled span, technical failure and ambiguity endpoints.

### Inferential hierarchy

The analytical sequence contains tests with different inferential roles. The original within-species spatial test and the Step-4 species-attribute family had prospective contracts. The discovery polymorphism–spatial gradient was a post-outcome exploratory reaggregation. The reserve is species-disjoint replication. Later sampled-opportunity, technical-failure and ambiguity analyses are explicitly labelled robustness or sensitivity analyses rather than retroactively being called preregistered primary tests. Throughout, Monte Carlo probabilities use the frozen randomization schemes of their corresponding analysis layer.

---

## Results

### Flower-colour polymorphism varied continuously among species

The discovery frame contained 369 species spanning the full observed range from `D = 0` to `D = 0.707645` (Figure 1). The distribution was continuous rather than separating naturally into invariant and polymorphic categories. A second colour state accounted for at least 10% of admitted observations in 46.61% of species and at least 20% in 26.56%. These percentages describe the admitted fixed-photo frame and should not be interpreted as estimates of the prevalence of flower-colour polymorphism across angiosperms.

Greater `D` did not correspond to expansion along one privileged recurrent colour direction. The primary relationship between `D` and variance in the frozen M1–M2 recurrent subspace was negative (`rho = -0.1647`) rather than positive, and the predeclared directionality claim failed. Direction specificity was also unsupported relative to random two-dimensional orientations (`P = 0.724`). In contrast, `D` was strongly correlated with total colour-space variance (`rho = 0.826`). Thus higher `D` primarily represented greater colour diversity, not preferential expansion along one common global axis.

### More polymorphic species were more geographically organized in discovery

The original 369-species whole-frame spatial analysis showed positive within-species geographic colour organization: equal-species mean spatial `rho = 0.027021`, with geometry-preserving `P = 0.001` (Figure 2). The subsequent polymorphism-focused reaggregation showed that this spatial signal was stronger among more polymorphic species. Species with a second colour state contributing at least 10% of observations (`n = 172`) had mean spatial `rho = 0.034565`, compared with `0.020435` in the 197-species complement. The difference was `0.014130` (`P = 0.014`). Across all discovery species, continuous `D` was positively associated with species spatial `rho` (`rho = 0.089213`, `P = 0.034`). Because this comparison was defined after the original spatial outcome was available, we treat it as exploratory discovery rather than independent confirmation.

### The polymorphism–spatial relationship replicated in species-disjoint reserve data

The reserve reproduced the relationship using different species. Among 363 eligible reserve species, continuous `D` correlated positively with primary species spatial `rho` (`rho = 0.101601`, `P = 0.025`; Figure 2). The 151 species meeting the frozen 10% second-morph threshold had primary mean `rho = 0.032299` (`P = 0.001`), compared with `0.020627` in the 212-species complement. Their difference was `0.011672` (`P = 0.041`).

The threshold-defined reserve subset remained spatially structured after observer-pair exclusion (`mean rho = 0.031775`, `P = 0.001`) and calendar-quarter stratification (`mean rho = 0.032299`, `P = 0.001`). Importantly, the matched flower-minus-background differential was also positive (`mean rho = 0.013622`, `P = 0.005`). The reserve therefore reproduced not only the direction of the continuous polymorphism–spatial gradient but also the stronger spatial organization of species with a non-trivial secondary colour state.

### Sampled span and clear technical failure did not explain the spatial relationship

Sampled geographic opportunity alone did not provide a general explanation for the pattern. In discovery, `D` was positively associated with sampled span (`rho = 0.181031`, Holm-adjusted `P = 0.00480`), but this relationship did not recur in the reserve (`rho = -0.00257`). More importantly, the polymorphism–spatial relationship remained positive after explicit opportunity and measurement-process controls.

The measurement audit confirmed that the raw photo denominator was exactly 100 per candidate species. Unclassified observations therefore could not be treated simply as lower acquisition effort. In the 50,000-photo reserve candidate frame, 45.36% of photographs failed the ROI/flip stability gate, whereas 4.83% were admitted by the ROI pipeline but left unresolved because of ambiguous palette composition (Figure 3). Species-level `D` was moderately associated with clear technical-failure rate in both discovery (`rho = 0.233`) and reserve (`rho = 0.217`), but much more strongly associated with ambiguous-palette rate (`rho = 0.728` and `0.704`, respectively). This difference motivated retaining ambiguity as its own measurement state rather than absorbing it into a generic classifiability control.

After controlling sampled geographic span and the **pure technical ROI/flip failure rate**, the discovery relationship remained positive (`partial rho = 0.126637`, geometry-preserving `P = 0.007`). The reserve primary relationship also remained (`partial rho = 0.099288`, `P = 0.025`), as did the reserve matched flower-minus-background relationship (`partial rho = 0.116241`, `P = 0.010`). Clear technical image-processing failure plus sampled geographic opportunity was therefore insufficient to reproduce the observed association.

### The spatial result survived exact ambiguity endpoints under the four-state completion model

Ambiguous-palette observations were common enough to matter but produced bounded uncertainty in species-level `D` under the four-state completion model. The median `D_max4 - D_min4` interval width was 0.0742 in discovery and 0.0732 in reserve; the 95th percentiles were 0.2533 and 0.2790, respectively. Observed `D` remained highly rank-correlated with both completion endpoints: in discovery, `rho(D,D_min4) = 0.9963` and `rho(D,D_max4) = 0.9730`; in reserve, the corresponding correlations were 0.9969 and 0.9610 (Figure 3).

The spatial conclusion was stable at both endpoints. After controlling sampled span and technical-failure rate, discovery primary partial `rho` was 0.133508 at `D_min4` (`P = 0.006`) and 0.117438 at `D_max4` (`P = 0.010`). In the reserve, primary partial `rho` was 0.097078 at `D_min4` (`P = 0.029`) and 0.125286 at `D_max4` (`P = 0.008`). The matched flower-minus-background relationship was likewise supported at both reserve endpoints: `rho = 0.116299`, `P = 0.009` at `D_min4`, and `rho = 0.132785`, `P = 0.006` at `D_max4`.

Thus, within the stated four-state sensitivity model, the spatial result did not depend on completing ambiguous observations in either the diversity-minimizing or diversity-maximizing direction. This does not identify the true latent state of ambiguous observations; it shows that two mathematically extreme uniform completions lead to the same qualitative spatial inference.

### Genus-level taxonomic structure replicated but was more ambiguity-sensitive

The prospectively defined discovery species-attribute analysis identified genus-level taxonomic clustering after multiplicity correction. Sixty-one repeated genera containing 169 species produced a clustering gain of 0.202963 (raw `P = 0.00230`, Holm-adjusted `P = 0.01150`; Figure 4). Family-level clustering was not supported (Holm `P = 1`), nor was absolute latitude. Pollination and life form were not tested because their prospective coverage gates failed.

Raw genus clustering recurred in the species-disjoint reserve: 54 repeated genera containing 146 species yielded a gain of 0.153654 (`P = 0.02170`). The 23 genera represented by at least two different species in both tranches also showed concordant genus means (`Spearman rho = 0.530632`, `P = 0.00990`). Sampled span alone did not remove reserve clustering (gain `= 0.141642`, `P = 0.02300`), and adjustment for sampled span plus pure technical-failure rate retained weak support (gain `= 0.121211`, `P = 0.04485`).

Unlike the spatial result, however, genus clustering was not uniformly supported across ambiguity endpoints. Under `D_min4`, reserve clustering gain declined to 0.110212 (`P = 0.06165`), whereas under `D_max4` it was 0.141546 (`P = 0.02380`). Cross-tranche genus-mean concordance remained positive at both endpoints, but the clustering statistic itself was sensitive to how ambiguous observations were completed. We therefore treat genus-level structure as a qualified secondary result rather than as evidence for formal phylogenetic signal.

### Repeatable organization did not require a common global boundary

The positive result above concerns the *strength of within-species organization*. Earlier global-boundary analyses asked a stronger question: whether independent species concentrate their strongest colour transitions in the same geographic regions. Those analyses did not provide confirmatory evidence for one shared global boundary, and the recurrent colour-direction test was also unsupported. The empirical picture is therefore asymmetric: within-species geographic organization is repeatedly recoverable, and its strength increases with polymorphism, but the geography and colour direction of differentiation need not be shared among species.

---

## Discussion

### The repeatable feature is species-specific organization, not a universal boundary

The main result is a replicated association between the amount of flower-colour polymorphism a species carries and the strength of its internal geographic colour organization. This relationship was first recovered in a 369-species discovery frame and then repeated in a species-disjoint 363-species reserve. In the reserve it survived observer and calendar sensitivities, subtraction of matched background structure, explicit sampled-span and technical-failure controls, and both diversity endpoints of a four-state ambiguity-completion model. The convergence of these tests is more informative than any single nominal probability: the same directional relationship persists when several distinct observational explanations are removed or stressed.

This result does **not** mean that highly polymorphic species share one geographic transition map. The distinction between organization and coincidence is central. A species can show strong geographic segregation of colour states along an east–west gradient, another along elevation, a third across historical barriers, and a fourth through several regional replacements. All can yield positive within-species geographic organization without producing a common global boundary. The failure of a universal-boundary claim is therefore not a contradiction of the positive spatial result; it identifies the level at which repeatability occurs.

### Amount, direction and geography are separable properties of colour variation

The failed recurrent-direction analysis sharpens this interpretation. Species with greater `D` occupied more colour-space variance overall, but that extra variation was not concentrated along the same learned M1–M2 direction. In other words, the scalar *amount* of polymorphism can be comparable across species even when its phenotypic direction differs. The same logic applies to geography. A repeatable relationship between `D` and spatial organization requires neither identical colour transitions nor identical locations.

Separating these levels may help reconcile apparently heterogeneous flower-colour systems. Reviews of flower-colour polymorphism emphasize multiple mechanisms—pollinator selection, non-pollinator selection, gene flow, drift and genetic architecture—that can act in different combinations among species (Narbona et al., 2018; Trunschke et al., 2021). Such heterogeneity makes a universal mechanistic boundary unlikely, but it is compatible with a higher-level regularity: when multiple colour states are maintained at appreciable frequencies, their distribution is more likely to become geographically structured than in species dominated by one state.

### The present data identify organization, not its ecological cause

Several biological processes could generate the observed association. Spatially varying abiotic conditions can favour different pigmentation states; pollinator assemblages can differ regionally; dispersal limitation and demographic history can preserve geographic differentiation; and drift can structure morph frequencies among partially isolated populations. Flower-colour systems provide examples consistent with several of these routes. In *Linanthus parryae*, colour morph fitness and spatial differentiation have been linked to precipitation and local environmental variation (Schemske & Bierzychudek, 2001, 2007). In *Boechera stricta*, floral pigmentation covaries with drought, elevation and herbivory (Vaidya et al., 2018). More broadly, macroecological analyses show associations between flower colour and biotic or abiotic conditions (Dalrymple et al., 2020).

None of those mechanisms is identified by our spatial statistic. Positive geographic organization is observational evidence that colour differences increase with spatial separation relative to a species-conditioned null; it is not a test of adaptation or local selection. The reserve result that `D` is essentially uncorrelated with sampled span also argues against turning the discovery span association into a general biological range-size claim. Future work can use the present result to prioritize mechanistic comparisons, but environmental and pollinator drivers should be tested directly rather than inferred from the existence of spatial structure.

### Measurement uncertainty is part of the inference, not a nuisance footnote

Community photographs offer scale but force explicit decisions about what a failed or ambiguous image means. A useful feature of the present design is that raw acquisition effort was fixed at 100 photographs per species. This made it possible to see that `n_classifiable` was not merely effort: it combined clearly technical failures with observations whose flower regions were measurable but whose palette composition resisted safe assignment to one discrete state.

Separating these processes materially changed interpretation. Controlling the clearly technical ROI/flip failure rate did not remove either the reserve spatial gradient or the flower-minus-background signal. Ambiguous-palette rate, by contrast, was strongly associated with `D`, so conditioning on total classifiability would risk adjusting away part of the visual complexity that defines the measurement problem itself. Rather than choosing whether ambiguity was “error” or “biology,” we bounded its possible effect under an explicit four-state model. The main spatial conclusion survived both species-level diversity endpoints. This does not solve ambiguity in a latent-state sense, but it shows that the flagship result is not balanced on one arbitrary treatment of the ambiguous rows.

The matched-background result supplies a complementary safeguard. Spatial autocorrelation can enter photographs through landscape, lighting, seasonality and photographic practice. By subtracting paired background spatial structure and retaining a positive polymorphism association, the reserve analysis shows that the signal is not simply an image-wide geographic pattern. Together, fixed denominators, technical-failure decomposition, background matching and ambiguity bounds provide a more informative robustness argument than a single generic “image quality” covariate.

### Genus structure is suggestive, but not yet a phylogenetic result

Species-level `D` was more similar within genera than expected under taxonomic exchangeability, and this pattern replicated in a species-disjoint reserve. Cross-tranche genus means were also positively concordant even though the contributing species differed. These findings suggest that the propensity captured by the photo-derived polymorphism metric is not distributed randomly across taxonomy.

However, the genus result is clearly weaker than the spatial result. The reserve clustering statistic crossed the conventional 0.05 threshold at the `D_min4` ambiguity endpoint, whereas the spatial gradient remained supported at both endpoints. We therefore avoid language such as “phylogenetic signal” or “phylogenetic conservatism.” Genus membership is only a coarse taxonomic grouping, and ambiguous-colour measurement itself may be taxonomically patterned. A natural next analysis would combine a dated phylogeny with an explicit latent measurement model for species-level polymorphism, but adding such a model after seeing the present results would change the scope of this paper. Here, genus clustering is best treated as a replicated but ambiguity-sensitive secondary observation.

### From Camellia to a broader hierarchy of flower-colour variation

This paper also clarifies the relationship between species-specific mechanistic work and global comparative patterns. Work in *Camellia* asks how similar visible flower-colour outcomes can be generated repeatedly through partially reusable molecular modules. The present analysis asks a later-level question: once colour variation exists within species, how much of it is carried and how is it arranged geographically? The two problems need not share the same answer. Repeated generation of visible colour, persistence of alternative states and spatial organization are distinct levels of flower-colour evolution. Keeping them separate prevents a molecular recurrence result from being mistaken for a global ecological mechanism, and prevents a spatial pattern from being overinterpreted as evidence about how the colour states originated.

### Limitations and scope

Five limitations define the present claim. First, the discovery and reserve frames are fixed-photo analytical samples, not random samples of angiosperms; the observed fractions of polymorphic species are therefore descriptive of admitted frames. Second, `D` summarizes four broad visible colour states and does not represent continuous reflectance spectra, pigment chemistry or pollinator-perceived colour. Third, ambiguous-palette observations remain biologically unresolved. `D_min4` and `D_max4` are exact only under the stated assumption that ambiguity can be completed into the four existing states, and applying one endpoint uniformly across all species is not an adversarial optimization over arbitrary species-specific latent allocations. Fourth, geographic organization is observational and does not identify selection, migration, drift or environmental causation. Fifth, although the common frozen acquisition frame required iNaturalist Research Grade, georeferencing, flowering annotation and the other stated metadata gates, it did not impose a native-range filter or explicit `captive=false`/`wild=true` parameter. Horticultural, introduced or otherwise non-native-range records may therefore remain, so the spatial result is an estimand of organization in the observed community-photograph sample rather than a native-range-only cline.

Within those limits, the core result is unusually stable. Greater flower-colour polymorphism was associated with stronger species-specific geographic organization in discovery and independent reserve data, and the reserve relationship survived explicit controls for sampling geometry, observer/calendar structure, image background, clear technical failure and two extreme ambiguity completions. The simplest synthesis is therefore not a universal map of flower-colour evolution. It is a hierarchy: **species differ in how much colour variation they carry; greater variation is associated with stronger internal geographic organization; and the specific direction and location of that organization remain species dependent.**

---

## Data availability and reproducibility

All analytical contracts, frozen input ledgers, result receipts, null summaries, synchronized paper-number ledgers, figure-data tables and figure-generation code are versioned in the project repository. Acquisition provenance is independently reconstructed by `scripts/analysis/audit_polymorphism_acquisition_provenance_20260911.py` and frozen in `results/polymorphism_acquisition_provenance_20260911/`; that audit opens neither candidate image pixels nor flower-colour outcomes. The paper-facing numerical authority is `paper/polymorphism_v0_1/paper_numbers.json`; main figures are generated from the synchronized `paper/polymorphism_v0_1/figure_data/` tables. Repository/archive DOI and final public-data statement will be inserted at submission freeze.

## Author contributions

To be completed at submission freeze using the final author list and CRediT roles.

## Acknowledgements

To be completed at submission freeze.

## Conflict of interest

To be completed at submission freeze.

---

## References

Dalrymple, R. L. et al. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist*, 228, 1972–1985. https://doi.org/10.1111/nph.16737

Delph, L. F. & Kelly, J. K. (2014). On the importance of balancing selection in plants. *New Phytologist*, 201, 45–56. https://doi.org/10.1111/nph.12441

Farquhar, J. E., Pili, A. & Russell, W. (2023). Using crowdsourced photographic records to explore geographical variation in colour polymorphism of an Australian varanid. *Journal of Biogeography*, 50, 1409–1421. https://doi.org/10.1111/jbi.14500

Jansen, N., Pruijn, N. & Mayer, M. (2025). Citizen Observations Shed New Light on Geographic Variation in Colour Polymorphism of a Widespread Reptile. *Journal of Biogeography*, 52, 629–640. https://doi.org/10.1111/jbi.15062

Narbona, E., Wang, H., Ortiz, P. L., Arista, M. & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology*, 20(S1), 8–20. https://doi.org/10.1111/plb.12575

Schemske, D. W. & Bierzychudek, P. (2001). Evolution of flower color in the desert annual *Linanthus parryae*: Wright revisited. *Evolution*, 55, 1269–1282. https://doi.org/10.1111/j.0014-3820.2001.tb00650.x

Schemske, D. W. & Bierzychudek, P. (2007). Spatial differentiation for flower color in the desert annual *Linanthus parryae*: was Wright right? *Evolution*, 61, 2528–2543. https://doi.org/10.1111/j.1558-5646.2007.00219.x

Trunschke, J., Lunau, K., Pyke, G. H., Ren, Z.-X. & Wang, H. (2021). Flower Color Evolution and the Evidence of Pollinator-Mediated Selection. *Frontiers in Plant Science*, 12, 617851. https://doi.org/10.3389/fpls.2021.617851

Vaidya, P. et al. (2018). Ecological causes and consequences of flower color polymorphism in a self-pollinating plant (*Boechera stricta*). *New Phytologist*, 218, 380–392. https://doi.org/10.1111/nph.14998

Westerband, A. C., Funk, J. L. & Barton, K. E. (2021). Intraspecific trait variation in plants: a renewed focus on its role in ecological processes. *Annals of Botany*, 127, 397–410. https://doi.org/10.1093/aob/mcab011
