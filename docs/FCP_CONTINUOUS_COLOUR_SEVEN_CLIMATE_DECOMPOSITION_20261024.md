# FCP original 447,500-photo nine-colour × seven-climate decomposition

**October 24, 2026.** Source-bound exploratory analysis. The source image/photo measurement and original New Phytologist P1–P4 results are immutable. The new question is **how the photographed continuous nine-colour distribution covaries with local long-term climate normals within a named species**. It does **not** ask whether camera photos prove heritable colour morphs or that climate adaptation maintains them.

## Exact source and climate availability

The actual 4,475-species × 100 original photographs = 447,500 rows are reloaded using the pinned 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6 original compressed photo CSV SHA256. Exactly 227,107 originally classified photos contain normalized nine-colour flower-ROI fractions; 220,393 nonclassified photos remain missing and are **not assigned a hue**.

Source climate is the original WorldClim 2.1 10-arc-minute BIO archive (SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf) and monthly solar-radiation archive (SHA256 c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026). The seven historical axis labels are:
- BIO1 annual mean temperature;
- BIO4 temperature seasonality;
- BIO5 warmest-month maximum temperature;
- BIO12 annual precipitation;
- BIO14 driest-month precipitation;
- BIO15 precipitation seasonality;
- SRADmean mean of the original 12 solar-radiation normals.

All seven axes are used as **a family of separate source-exact univariate photo-climate relationships** with Holm multiplicity; they are not jointly interpreted as seven statistically independent causes. Sampling is conditional on the original 10-arc-minute photographed locations, not measured fine-scale weather at flowering time. Meteorology, drought, nectar, pollinators, soil or genotype measurements are **not** in this source.

## The distinction between a continuously shifting hue and colour-state frequency changes

The actual photo source has a 9-element normalized colour-palette composition Y and an originally frozen photo-derived coarse label C: white, yellow/orange, red/pink or blue/purple. Climate normal Z is centered and standardized within each species. For any climate axis a, the exact finite-sample covariance identity is:

    Cov(Z_a,Y) = E_C[ Cov(Z_a,Y | C) ] + Cov_C(E[Z_a|C], E[Y|C]).

The first term is the **within-existing-photographic-label hue-shift term** and the second is the **between-existing-photographic-label mean/frequency difference term**. This avoids manufacturing species-wide *genetic* 2/3/4 morphs by Gaussian mixtures. It is also more precise than claiming all colour changes are either a smooth curve or an admixture of stable modes: both terms can coexist or cancel, and original C is simply a derived label from the same photo colour.

Record both (a) the **signed white-palette fraction covariance per within-species climate SD**, and (b) the vector norm of covariance in all nine hue fractions. The latter has a nonzero finite-sample null expectation: subtract matched-null expectation before reporting the mean excess norm. **Do not divide norms of the two components and report them as additive percentages of genetic variance**: the components are signed 9D vectors and can cancel.

This is an identity conditional on original photo labels, not an identifiable biological split between (A) continuous within-genotype colour change and (B) independent inherited-morph frequency changes. The latter requires real petal/organ and genotype confirmation.

## Exact observation opportunity, matching and spatial sensitivity

The entire 4,475-species original roster remains in the per-species outputs. A species is estimable for source-image climate covariance only if there are at least 40 *original classified, valid georeferenced, complete seven-climate-axis* photographs, at least 5 different occupied 10-arc-minute WorldClim climate cells and at least 20 independent observer IDs. A species-specific climate axis must have nontrivial variation. The two-state frequency description is shown only for taxa with **at least two original coarse labels supported by at least five photographs each**.

Two analysis scales per species:
1. **Raw photographed climate covariance:** fitted on all eligible photos; geography and camera source may be responsible for the association, so this is explicitly descriptive.
2. **Calendar-month and broad-geography residual climate covariance:** regress each climate variable on calendar month and Earth-centred geographic XYZ coordinates, using a species-specific linear projection and requiring at least 5% source climate variance remains unexplained. Photo hue is not itself regressed, but covariance with orthogonalized climate equals covariance of orthogonalized hue with orthogonalized climate for the fitted nuisance basis. This is **not a proof of climate independent of nonlinear geography or population genetic history**; strong climate–space collinearity is reported as an eligibility limitation.

In each species and each analysis scale, **199 synchronized month-conditioned permutation worlds** shuffle the original complete seven-climate vectors among photographs **within the same species and calendar month**, while fixing the original photo hue vectors, labels, photos, coordinates, climate sample, and number of individuals. This is a month-matched comparison against chance assignment; it **does not match observers, geography or colour/life stage**, and is not a randomization of climate exposure in natural plants.

Recompute after **colour-blind one-photo-per-observer thinning**. Report the actual species/cell/observer opportunity for raw and residual climate analysis; no photo is imputed to strengthen statistical power.

Apply Holm over the seven climate axes separately for each combination of (source vs spatial-residual), (whole 9D vs signed white), and (total/within/among original coarse labels). These are **twelve exploratory endpoint families**; the Holm over seven does not silently control the entire hypothesis search. Add 999 resamples of genera for two-sided uncertainty intervals of equal-species mean excess effects.

## Predeclared interpretation gates

- A nonzero total effect does not imply both within and between effects. A positive magnitude is not a signed shift toward any named colour.
- A between-label effect does not establish discrete inherited modes or climate-driven morph frequency change; C was constructed from the same original photographs.
- A within-label effect indicates photographic continuous hue variation **inside a coarse bin** associated with climate but may still be lighting, flower organ, developmental stage, camera, cultivation, geography or compositional mixing.
- If the raw effect is positive but the geography/month residual effect vanishes or the eligible set collapses, **no independent climatic association has been established**.
- No fitness, pollinator-mediated selection, germplasm, anthocyanin investment cost or local adaptation is identified.
- Earlier three species-disjoint cohort nine-hue JSD×BIO5 partial correlations (+0.00925,+0.00826,+0.01019) were statistically positive, but flower-minus-background (rho −0.00328,p=0.82) and same-observer tests did not support a flower-specific observer-independent climatic cause. This new output cannot overrule those negative challenges.
- The strongest prior finding remains species-conditional broad-scale photograph-colour spatial structure, not a genetic morph prevalence claim.

## Reproducibility

- Code: scripts/analysis/fcp4475_nine_colour_seven_climate_covariance_20261024.py
- CI: .github/workflows/fcp4475-nine-hue-seven-climate-covariance-20261024.yml
- Output: two all-4,475-species tables (original and one-per-observer, **no exact coordinates/photos/observer IDs**) plus full seven-axis machine summary.
- All source photos and BIO/SRAD checksum verification, exact covariance synthetic fixtures, no imputed missing hues, source-denominator tests and timestamped GitHub Actions source work must succeed before any numerical result is asserted.
