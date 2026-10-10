# FCP original 85,337 species×cell photos: nearest local same-species four-colour environmental orientation null (2026-10-10)

## Current source-supported bottleneck

PR #145 found that the original four-colour flower-photo phenotype is weakly predicted by thermal, moisture, solar, wind, vapor-pressure, elevation and SoilGrids block data on **20,546 heldout species photo cells** even beyond training-only species colour means and a low-rank real-source photographic great-circle spatial kernel. However local residual spatial autocorrelation remained, so it cannot be called an evolutionarily adaptive abiotic effect.

PR #146 directly measured whether source photographed labels could be exchanged among TRUE geographically close photos of the SAME species without violating original species identity or geographic labels. Among 38,968 climate-complete source photo-cell phenotypes, only **44 source photos in 22 species within 50km**, **144 photos/72 species within 100km**, **382 photos/187 species within 250km**, **908 photos/437 species within 500km** had different source photographed colour states in the same species' nearby complete-link group. This fails the predefined 10% group-based matched permutation power feasibility gate at every distance, and is a sampling identifiability limitation rather than negative biological adaptation evidence.

## Further exploratory fixed ONE nearest original source pair/species

To leverage the source without choosing favourable photographs based on their flower colour outcome, this task selects **exactly ONE original photographed site pair per taxon** from all already CLASSIFIABLE source species×cell photos with all climatic data, minimizing actual great-circle site distance; if ties, choose the lowest original photo ID order. Do NOT select pairs based on whether their colours differ or whether environment contrasts are extreme. Repeat exactly this pre-outcome fixed-pair set at 50km, 100km, 250km and 500km original photographed-site distance caps and report every cap, including HOLD.

No later source photograph shall be substituted after seeing its four-colour photo label. The 500km source pair set includes closer subsets and does not provide an independent replication of 100/250km.

### Four original hue classes, not presumed ordered pigment values

The original four coarse photographed classes are white, yellow-orange, red-pink and blue-purple. Do not assume yellow, red and blue represent greater anthocyanin or shared biochemical pigment intensity. For each actual original photo pair A/B ordered by photo ID, the response is the vector difference in the **four one-hot colour indicators**. Pairs with the same photographed colour have zero vector; they remain in the complete source denominator, not reassigned to white or called biologically monomorphic.

The paired environmental value difference is B - A for all WorldClim BIO1–BIO19, 12-month climatic radiation/wind/vapour annual mean and monthly CV, and ten modelled SoilGrids topsoil properties where source complete.

### Physical geography conditional test

All geographic and environmental differences derive only from the original publicly available photographed sites, not region-cell centroids. Before examining four-colour orientation, residualize every signed environmental difference against signed latitude, signed absolute-latitude, signed longitude sine/cosine and signed source photograph elevation differences **over ALL already geographically eligible nearest source pairs**. Environmental residualization is entirely phenotype-blind.

For each environment, test a scale-normalized four-dimensional original colour-gradient contrast statistic: the norm of the summed product of environment residual with photo-colour one-hot difference, standardized by the pairwise RMS null scale. Compare with **999 independent within-pair photo-colour orientation swaps** between original A and B; each swap preserves the species' original two four-colour photo-counts and photographed sites. Never shuffle labels across species, swap original photo coordinates, or let a specimen appear in more than one selected pair.

Use the same sign swap across ALL environmental variables, so each radius/environment complete-case sample yields familywise **maxT adjustment over ALL 25 climate and/or 35 full soil+climate named environmental variables**, not one-at-a-time uncorrected winner p-values. Additionally report fixed conservative Bonferroni ×8 over 4 nested radii ×2 correlated climate/soil cohorts. Report every named physical block and individual variable, including null/negative effects. At any radius, skip fitting a null with fewer than 60 original paired taxa, 30 photographic hue-discordant pairs, or 20 distinct genera.

### What this estimates and does NOT estimate

A positive result would mean a **source-observational four-colour gradient orientation associated with geographic/elevation-residualized environment among nearest photographed localities within named species** under a conditional within-pair label-exchange assumption. It is more directly geographically paired than the earlier source-global 5fold predictions, but it does NOT completely remove spatial autocorrelation. It is neither a genomic within-population flower-colour polymorphism, an experimentally randomized exposure, a pollinator/payoff test, proof of local adaptation, nor an independent new source taxon set.

A negative/maxT-nonsignificant result would indicate that these original paired photos do not independently support a repeatable local directional environmental contrast; it would not show a biological absence of pigment response because source flower photographs are sparse and single observations per region. Note that using only a pair per taxon sacrifices within-species replication to avoid dyadic pseudoreplication. Soil missingness changes which original PHOTO PAIRS are fully testable, not how those pairs were originally selected.

**Full source historical denominators retained**: 42,111 nominal species, 85,337 original species×cell photos, 39,075 classified source photos, 38,968 complete-climate source photos (20,546 with original test-heldout species intercept), and 27,003 completely selected-soil source photo cells (11,136 with heldout species intercept). Original 1,499 high-depth H1/H2, main and the independent unexposed 2,000+730 prospective taxa are unchanged.

This is an additional retrospective exploration prompted by earlier spatial autocorrelation tests, not preregistered confirmatory inference. The scientific gate for proving global environment beyond BOTH true spatial and independent phylogenetic residual covariance remains unfulfilled.
