# Climate landscape macroecology: results, 2026-11-27

## Completed source analysis
Run 45971197518 (job 140411780000) passed all steps, including frozen photograph and WorldClim BIO archive authentication. Source 447,500 photos, 4,475 taxa; 447,446 climate-complete source locations and 227,107 original classifiable photographs. No observer splitting. Full numerical receipt: results/fcp4475_climate_landscape_verified_20261127.json. Artifact 15184882647; provider ZIP digest SHA256 a3398a827f9d8b5a3478c3096f94d15752d858cfd95c62be2fa4e48c0b45adf4.

## Macroecological findings
Range-wide four-label diversity: n=4,333 species, additional within-species sampled climate breadth (thermal/moisture) and BIO4/BIO15 mean seasonality increase R2 from .32519 to .32628 (delta .00109). All four marginal genus 95% intervals include zero. This is no strong evidence for these predictors of the *amount* of photographic colour variation.

Spatial colour allocation: species with at least two repeated photo categories and enough near/far pair opportunities. Regression on sampled climate partition (1 - near-pair climate squared distance / range-wide climate squared distance), conditional on sampled range size, climatic means and breadth, near-pair fractions, total colour diversity, photo counts, seasonality, source cohort, and latitude. Coefficients per predictor SD; genus percentile bootstrap 1,999.

| Panel | n | Thermal coefficient [95%] | Moisture coefficient [95%] | incremental R2 |
|---|---:|---|---|---:|
| 50 km four-category | 1315 | .009234 [.001912,.016526] | -.011363 [-.019421,-.003582] | .00581 |
| 250 km four-category | 2110 | .007261 [.002547,.012147] | -.001608 [-.006431,.003239] | .00458 |
| 50 km exact common species | 1258 | .009453 [.002106,.017258] | -.010602 [-.019183,-.003067] | .00572 |
| 250 km exact common species | 1258 | .009540 [.004892,.014212] | -.004537 [-.009857,.000498] | .01253 |
| 50 km continuous nine-palette Hellinger | 1315 | .004568 [-.000219,.009205] | -.005002 [-.008983,-.001169] | .00508 |
| 250 km continuous nine-palette Hellinger | 2110 | .003326 [.000961,.005664] | -.001360 [-.003889,.001108] | .00314 |

The thermal association is positive at both scales in the matched categorical sample and at 250 km in the continuous measure. The negative conditional moisture association is supported at 50 km in categorical and continuous outcomes but not at 250 km. No formal difference between distances or thermal vs moisture coefficients has been established by comparing marginal intervals; possible covariance between climatic partition predictors remains. These are small, within-sample exploratory effects, not causal mechanisms or validated predictive accuracy. The 50-km Hellinger thermal interval includes zero. Genus clustering is not full phylogenetic control.

## Next biological question
Does the climate-organization association preferentially involve white-versus-chromatic contrasts, rather than redistribution among chromatic hues? Test exact pairwise additive components, compare component-normalized effects on the same taxa and evaluate absolute contributions on those taxa. Primary relative component baseline >= .01, both scales 50/250 km, all original eligible classified photographs. Reuse the fixed climate artifact, not new source selection. Do not interpret white as molecularly verified pigment loss. This extension does not establish pigment production cost, plasticity, or inherited morphs. Existing manuscript and previous analyses remain unchanged.
