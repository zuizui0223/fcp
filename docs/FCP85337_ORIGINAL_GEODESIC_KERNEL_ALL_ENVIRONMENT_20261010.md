# FCP 85,337 source photo-cells: within-species geodesic spatial-kernel control

The original source comprises 85,337 species-by-region flower photographs, of which 39,075 have a valid four-state photo-colour label and 46,262 remain unclassified. Source attributes for all 100,543 distinct original photo identities were attached only at genuine publicly available photographed coordinates, including WorldClim BIO1-19, elevation, 12-month solar/wind/vapour pressure summaries and SoilGrids ten properties.

## Question and design

Does measured flower-colour prediction across sites of the SAME named species improve when environmental variables are added to an explicit geographic spatial covariance baseline? Earlier source species-intercept analyses had 20,546 climate-complete test photos and 11,136 soil-complete test photos. Spherical third-degree polynomial spatial adjustment still left small positive environmental predictive information, but it was not a covariance field.

For this sensitivity, group-heldout fivefold splits exclude entire original equal-area geographic cells. Within every training fold, derive original-species colour intercepts from OTHER photographed cells only. Map original site latitude and longitude to unit-sphere coordinates and estimate a low-rank spatial feature basis using the exponential of negative great-circle distance divided by 100km or 500km. The 128 Nyström landmarks are sampled from TRAIN sites using an outcome-blind fixed seed. All models evaluate exactly the same heldout photo identities and folds.

Compare species intercept + latitude/longitude, then plus actual photographed-site spatial-kernel features, then plus six climate/elevation/weather blocks, and when source soil is complete, all seven physical blocks. Report each predeclared block's conditional drop-one predictive contribution and original species and geographic-cell conditional bootstrap intervals, never just the highest gain.

## Limits

This is a 128-landmark approximation to a geodesic spatial covariance, not a dense Gaussian spatial random field. Even with positive environmental prediction, residual spatial dependence and observation bias are not proven absent. The fixed species intercept absorbs species-invariant phylogenetic main effects, but does not estimate genetic flower-colour variation, population-level frequencies or adaptation. The directly-tipped LCVP species subset remains a separate phylogenetic study and the full original 42,111-taxon phylogenetic inference is HOLD. No change to original 1,499 high-depth paper or reserved 2,000+730 species.

The result is valid only after the exact 85,337/39,075 source checks and synthetic spatial-kernel tests pass in the dedicated GitHub Actions workflow.
