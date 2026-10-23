# FCP: seven-climate axes explain only a small portion of photographed continuous hue sorting

**24 October 2026. Status: original-data computation and independent source-artifact common-species audit both PASS. Post hoc observational comparison, NOT inherited flower-morph maintenance or local adaptation.**

## Original source and separate analysis denominators

- Original 447,500 photographed records from 4,475 named plant species; original 4-state-classifiable nine-hue ROI: 227,107 photographs (50.75%).
- Original seven source-pinned WorldClim 2.1 climate normals: BIO1 annual mean temperature, BIO4 temperature seasonality, BIO5 warmest-month maximum temperature, BIO12 annual precipitation, BIO14 driest-month precipitation, BIO15 precipitation seasonality, and twelve-month average solar radiation.
- At least 40 original classifiable photos with complete climate, five 10-arcminute geographic raster cells and 20 independent observers: **3,296** original species; **3,227** after one photo per observer; all 4,475 remain explicitly listed in the species-level CSV.
- At least two original photographic colour labels, each supported by five original photographs: **1,592** species; one-observer sensitivity **1,509**. The full analysis does not silently classify all observed source colour bins as inherited morphs.
- Geographical/month residualisation retains only species whose climate residual retains at least 5% of original species-level climate variance. For the complete-source signed white BIO5 on same-multi-colour support, **1,578 species (710 genera)**, one-observer **1,494 species**.

Run source: https://github.com/zuizui0223/fcp/actions/runs/40029760904 (PASSED original photo/WorldClim SHA checks, exact 9D covariance synthetic control, 199 source/month-permutation worlds, 999 genus bootstraps).

Independent source-artifact join: https://github.com/zuizui0223/fcp/actions/runs/40030561641 (PASSED; full original taxon CSVs pinned by SHA and common-support 1,999-genus bootstrap).

Reproduction artifacts committed in this PR:
- scripts/analysis/fcp4475_nine_colour_seven_climate_covariance_20261024.py
- scripts/analysis/audit_fcp4475_nine_hue_climate_common_support_20261024.py
- results/fcp4475_nine_hue_climate_common_support_verified_20261024.json
- .github/workflows/fcp-4475-nine-hue-seven-climate-covariance-20261024.yml
- .github/workflows/fcp-4475-nine-hue-climate-common-support-20261024.yml

## Finding 1. Climate differences covary with all nine continuous hue fractions in photographed plants

Species-conditional covariance of seven standardized environmental variables with the complete nine original hue fractions gave small mean source-space hue-vector covariance norm excesses over fixed-photo-month permutation worlds. All seven mean vector norms were above their null expectation in raw photographic data; after projecting calendar-month and broad 3D geographic-coordinate effects out of the climate, they remained positive but much smaller (original taxon equal, source counts per axis). For example, BIO5 hue-vector excess decreased from **+0.00601 to +0.00221**; one-photo-per-observer sensitivity from **+0.00565 to +0.00213**. Each 199-world upper permutation p=0.005; Holm over the seven climate axes within the relevant norm endpoint gives 0.035.

**Critical limitation:** This is association of entire ROI photo colour with a long-term climate normal. It cannot establish that actual petals respond biologically, or that the categorical and photographic organ confounds have been removed. The 12 distinct exploratory endpoint families are not jointly corrected.

## Finding 2. Signed white-palette covariance has a recurrent climate direction mainly for temperature variables

Rather than 9D norm, take the **direction** of observed continuous white-palette share, not genetic white-morph occupancy. Compare exact original nine-hue covariances on the *same*, >=2-supported-photo-class species and calendar-month/XYZ-residualised climate:

| Climate predictor | Same-source species | Mean white-component source-month-null excess signed covariance | Simultaneous 7-climate genus-bootstrap 95% CI | Family7 excludes zero? |
|---|---:|---:|---:|---|
| BIO1 annual mean temperature | 1,467 | +0.003844 | [+0.002115,+0.005573] | Yes |
| BIO4 temperature seasonality | 1,426 | +0.002094 | [+0.000364,+0.003824] | Yes |
| BIO5 warmest-month maximum temperature | 1,578 | **+0.005278** | **[+0.003474,+0.007082]** | Yes |
| BIO12 annual precipitation | 1,576 | −0.001380 | [−0.003240,+0.000480] | No |
| BIO14 driest-month precipitation | 1,567 | −0.000538 | [−0.002243,+0.001168] | No |
| BIO15 precipitation seasonality | 1,539 | −0.000464 | [−0.002164,+0.001236] | No |
| SRADmean solar radiation | 1,133 | +0.000528 | [−0.001071,+0.002127] | No |

The signs and seven-family CI decisions for all seven variables are unchanged with one photo per observer. **BIO5 is the largest positive white-component signed covariance in these normalized observational estimands**, not a universal per-degree climate slope or proof of pigment adaptation.

Note that geo/month residual covariance is measured in **original species climate-SD-scaled units**, but is not a slope per one climate SD: the residual climate no longer has variance exactly one. It is *not* valid to interpret +0.00528 as “a 1°C increase causes white hue fraction to rise 0.53 percentage points”. Species, photo and climate sampling shape this covariance.

## Finding 3. Temperature–white association is largely between the four original photograph label groups, but some remains within labels

The equality

    total signed white-photo climate covariance =
       between original coarse photographic-label component
       + within original coarse photographic-label component

was verified **per species** and after synchronised source-month-matched null correction. Only the identical >=2-well-supported-label taxon sets were used to compute component contrasts.

| BIO5, geo/month residual | 1,578 original taxa | 1,494 observer-thinned taxa |
|---|---:|---:|
| Total signed original-white covariance excess | **+0.005278** | +0.005086 |
| Between existing photographed-label component | **+0.004002** | +0.003856 |
| Within existing photographed-label component | **+0.001276** | +0.001230 |
| Between-component share of **signed mean covariance** | **75.8%** | 75.8% |
| Between-component share genus bootstrap 95% | [67.2%,84.4%] | [65.8%,84.9%] |

Within this original source, the between-label excess is larger than within-label excess by **+0.002726** (1,999 genus-bootstrap CI [+0.001569,+0.003895]). The within-label component is also positive: **+0.001276** (marginal genus CI [+0.000756,+0.001816]).

**Meaning:** coarse *photographic colour-label frequency and associated group means* contribute more to the signed climate/white-component association than continuous palette changes *inside the existing coarse bins*. This is a **statistical identity conditional on colour classifications and observed photography**. The 75.8% is neither a proportion of heritable variance nor evidence that stable genetic morph frequencies respond to BIO5. Unmeasured smaller modes within bins, petal/blotch/flower-age and camera differences can enter either component.

## Which biological claims are *not* identified?

1. **Not climate-driven adaptation.** No reciprocal transplant, per-morph pollen transfer, physiological stress performance, matched neutral genetics, anthocyanin synthesis or fitness/cost measurements were collected in this world photo source.
2. **Not independence from spatial sampling.** We remove a low-dimensional linear spatial trend using 3D Earth coordinates plus photo month, not a full geographically structured latent field. BIO predictors correlate with one another and spatial history; the strict previously unsupported same-year/month <=5km independent observer tests have not been replaced.
3. **Not flower-organ error correction.** Source continuous nine-palette ROI segmentation and coarse labels have not all been verified by two botanical reviewers. The 589 original images from 27 suspected two-peak taxa already show plausible buds, fruits, post-floral sepals, pigment blotches and camera effects. A separate frozen three-cohort flower-minus-background BIO5 robustness test was unsupported (rho −0.00328, p=0.82); sparse same-observer tests were also unsupported.
4. **Not categorical species-wide genetic polymorphism.** Original 227,107 colour-classifiable photos are 50.75% of the original atlas, with strong classification missingness and no individual genetic matching. Within-existing-label changes are not necessarily continuous inherited phenotypes; between-label changes are not necessarily discrete genetically maintained morphs.

## Conservative manuscript-style inference

> Across the observed within-species photographic nine-colour trait composition, climatic variation, particularly the thermal BIO1/BIO4/BIO5 axes, showed small but repeatable covariation with the photographed white-palette fraction after simple geographic and photo-month controls. In taxa with adequately supported original photo-colour classes, approximately three quarters of the BIO5-associated signed white covariance arose between those photographed categories and the remainder arose within them. These terms are descriptive decompositions of photographed variation, not experimentally verified changes in genetic morph frequencies versus continuous pigment expression. Given prior negative flower-minus-background controls and unvalidated source flower organs, a flower-specific adaptive mechanism remains unestablished.

This is a supported observation bounded by the above limitations. Do **not** rewrite the historical New Phytologist main paper to claim a molecular pigment trade-off, universal climatic cline, temperature causation, or general genetic morph replacement.
