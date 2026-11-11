# FCP white-pigment benefit–cost: necessary geographic footprint versus generic hue turnover
## Post-outcome, equal-species macroecological test — 2026-10-08

### Biological question, not the assumed answer

Petal pigment may have multiple benefits (pollinator display, vegetative/petal antiherbivore or abiotic functions) and costs (resource allocation, predator attraction, pleiotropy, mating and reproductive investment). These opposing selective components can change by site and thereby maintain colour variants across a species' range. But the FCP image data measure photographed floral visible colours, NOT biochemical pigment concentration, heritable morphs or plant fitness.

The already validated broad macro finding is nearby geographic clustering of within-species colour variation in 3 species-disjoint cohorts, and H2 independently confirms a photographic achromatic/chromatic phenotype-space direction conditional on colour-class construction. White/chromatic photo categories are frequent but observed white is strongly exposure coupled. Previously, a paired local-white/nonwhite versus nonwhite-hue 50-km mixture comparison across species with both opportunities was NOT statistically supported (41/36/49 species; signflip p .5667/.8123/.2553).

The *next necessary but not sufficient* test of a **white-axis-specific benefit–cost mosaic** is:

> Is white–chromatic colour sorting systematically larger than the sorting of two different nonwhite visible hues within the **same plant species and geographic opportunity**, after preserving each morph-pair's photo counts and month-specific composition?

This differs from the unsupported *signed* global chromatic-high-latitude/elevation cline. The analysis measures **unsigned geographic heterogeneity**. Species may enrich white at either end of their range; a common sign is NOT required.

### What is actually tested

For each original source-SHA-frozen species-disjoint photo cohort, require >=40 classifiable photos from the existing measured tables, with valid date, and at least **1 degree** of absolute-latitude sampled span. Define the geographic axis as three photo-count quantile bins of absolute latitude. **Identical latitudes always remain in the same bin**, even if that makes tertiles uneven and reduces testable opportunity. No pseudo-spatial split of photos from exactly the same locality is allowed.

For each species, form three fixed composition comparisons BEFORE checking geographic effects:

1. **W_all:** white versus all three other coarse visible colour classes pooled, requiring >=5 per side.
2. **W_lead:** white versus the most frequent individually observed nonwhite hue, requiring >=5 per side.
3. **hue_lead_next:** most frequent versus second-most-frequent nonwhite hue, requiring >=5 per individual hue.

The leading hue is selected by observed global species photo frequency with a deterministic label-name tie break; NEVER by the largest geographic contrast. W_lead versus hue_lead_next is the matched candidate versus coloured-hue control: the same species, the same fixed geographic bins, and a shared leading-hue category. Both must occupy >=2 bins with >=3 photos in each supported bin.

For a binary comparison with global target proportion p, photo-tertile group sizes n_b, and within-tertile target fraction p_b, calculate

    V2 = [sum_b n_b (p_b - p)^2] / [N p(1-p)].

This is a 0-to-1 Cramer's-V-squared-like normalized colour-label / sampled-latitude-bin association, **without a directional slope**. V2 cannot be called a selection coefficient.

The primary null shuffles the binary labels **only inside original species-pair-calendar-month photo blocks**, preserving the global two-state composition, each month-specific two-state composition, all fixed photo locations, exact geographic bins and sample size. Use 199 predetermined random permutations per species and morph contrast. A no-exchangeability null contributes exactly zero **identifiable excess** in the global W_all denominator and is clearly marked nonidentified, not biologically zero.

Compute observed-minus-mean-null V2 and its permutation-null SD standardized residual per morph pair, enabling comparison of the two quite different morph pair sample sizes. In species where W_lead and hue_lead_next are BOTH eligible, month-exchangeable and null-identifiable, form within-species

    matched difference = z(W_lead) - z(hue_lead_next).

Summarize species equally, three cohorts separately, 1999 species bootstrap intervals and 9999 paired signflips. Genus-balanced means are secondary nonphylogenetic robustness controls. A universal white-axis-specific result requires >=25 matched species in **each** cohort, a positive mean and bootstrap CI above zero and p<0.05 from a two-sided paired signflip in all three; no after-the-fact promotion of whichever cohort or pigmentation hue appears favourable. This is a stringent *same-photo-system* exploratory generality rule, not an independent new dataset.

### Interpretation and independent empirical alternatives

- **Stronger W-lead than hue-lead-next**: would be consistent with a white-associated geographic axis more strongly sorted than a generic nonwhite colour contrast, *which could also be caused by neutral population history, spatial distribution, distinct flower colour detection probabilities, or developmental plasticity*. It cannot establish opposing selective costs or evolutionary-maintenance benefits.
- **Both similarly spatially assorted**: no demonstrated white-specific universal benefit–cost footprint; generic geographically allocated floral-colour ITV remains viable.
- **Coverage limited**: inability to compare both pairs across species cannot be used to infer a zero relative effect or absence of trade-offs.
- **If only whole-species white state incidence is high**: category pooling and global colour-label marginal prevalence may trivially favour a white+any-nonwhite count. The counts 516 white+one individually nonwhite hue versus 166 two nonwhite hues are strictly *sample-descriptive*, with overlapping categories and known exposure issues; they DO NOT estimate genetic white morph maintenance.

Independent biochemistry/fitness evidence is stronger for specific organism systems than for the photo map: in **Boechera stricta** pigmented variants experience less herbivory while white variants can have higher fecundity under well-watered treatment (Vaidya et al. 2018 DOI 10.1111/nph.14998), and in **Silene littorea** petal-specific pigment loss (PAL) occurs at substantially higher natural frequencies than whole-plant pigment loss (WAL), while photosynthetic-tissue pigments are maintained (Del Valle et al. 2019 DOI 10.1186/s12870-019-2082-6). The latter suggests an **organ-specific pleiotropic-cost constraint**, not an experimentally established universal pigment-synthesis energy cost.

### STOP conditions

Do not recast this or any photographic result as confirmed fitness balancing, morph-genotype equilibrium, universal anthocyanin pigment chemistry or local adaptation. Do not change frozen New Phytologist H1/H2/older results or claim unbiased world polymorphism incidence from the 1,499 opportunity-selected species. If the matched comparison is not supported, report it explicitly rather than another story about white selection.
