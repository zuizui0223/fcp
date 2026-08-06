#!/usr/bin/env python3
from pathlib import Path

path = Path("docs/jbi_manuscript_editorial_revision.md")
text = path.read_text(encoding="utf-8")

old = """Models were fitted in Python 3.12 using `statsmodels` 0.14.6 (Seabold & Perktold, 2010). Wald standard errors were estimated with family-clustered sandwich covariance. Reported 95% confidence intervals were calculated on the log-odds scale as the coefficient ±1.96 clustered standard errors and exponentiated for odds-ratio intervals. All 20 broader-set models converged in four iterations. The complete model matrix, including sample composition, coefficients, intervals, p-values, effort terms, convergence and fitted-probability ranges, is provided in Table S1.
"""
new = old + """
Because standardisation does not remove collinearity, we explicitly diagnosed the design matrices used for the focal baseline model and the coarse occurrence-cloud sensitivity models. We calculated variance inflation factors for each non-intercept predictor and the condition number of the intercept-inclusive design matrix. We did not use VIF-based stepwise selection: the climatic metric and effort covariates were retained a priori, and range or connectivity summaries were added one at a time except in the explicitly labelled integrated 100-km connectivity model. Diagnostic values and the exact predictor sets are archived in `docs/supporting/jbi_model_condition_diagnostics.csv`.
"""
if old not in text:
    raise SystemExit("Spatial-model insertion point not found")
text = text.replace(old, new, 1)

old = """We evaluated two sets of coarse alternatives to the focal relationship. First, models added sampled-range extent and fragmentation summaries derived from GBIF coordinates, including median nearest-neighbour distance, 95% spatial extent, numbers of components at 50- and 100-km thresholds, the fraction in the largest 100-km component, occupied 1° grid cells and an integrated 100-km connectivity model. These models included climate-cell effort and GBIF-record effort.
"""
new = """We evaluated two sets of coarse alternatives to the focal relationship. First, models added sampled-range extent and fragmentation summaries derived from GBIF coordinates, including median nearest-neighbour distance, 95% spatial extent, numbers of components at 50- and 100-km thresholds, the fraction in the largest 100-km component and occupied 1° grid cells. Each summary was added in a separate model containing moisture breadth, climate-cell effort and GBIF-record effort; the range summaries were not entered together. An additional, explicitly labelled integrated 100-km connectivity sensitivity model jointly included the 100-km component count and largest-component fraction. This one-at-a-time structure was used to limit redundancy and overfitting in the small comparative sample.
"""
if old not in text:
    raise SystemExit("Occurrence-cloud paragraph not found")
text = text.replace(old, new, 1)

old = """The baseline-unambiguous moisture model included 34 species from 25 families: 20 within-population and 14 among-population cases (Table 2). The standardised moisture-breadth coefficient was −0.854, corresponding to an odds ratio of 0.426. The family-clustered 95% confidence interval was 0.184–0.985, and the clustered Wald p-value was 0.0460. The model converged in four iterations, with fitted probabilities from 0.050 to 0.701.
"""
new = old + """
The two-predictor focal design showed little collinearity: both moisture breadth and climate-cell effort had VIF = 1.10, and the design-matrix condition number was 1.36. Thus, the focal coefficient was not made unstable by linear redundancy between the climatic metric and the effort covariate.
"""
if old not in text:
    raise SystemExit("Baseline-results insertion point not found")
text = text.replace(old, new, 1)

old = """The fragmentation dataset contained 55 model-complete species from 31 families. Added fragmentation and connectivity terms had odds ratios from 0.574 to 1.061; every 95% confidence interval included one and the smallest p-value among these added terms was 0.280 (Table S5). The moisture-breadth odds ratio remained below one across these models, ranging from 0.480 to 0.597, although its intervals generally included one.
"""
new = old + """
Collinearity was modest in most one-at-a-time occurrence-cloud models (maximum VIF 1.51–4.15; condition numbers 2.06–4.39). The occupied-1°-grid model exceeded the prespecified diagnostic threshold of five (maximum VIF 5.47), as did the integrated 100-km model marginally (maximum VIF 5.02). We therefore treat those two estimates as descriptive sensitivity checks rather than separable independent effects; neither is used to select or redefine the focal model.
"""
if old not in text:
    raise SystemExit("Fragmentation-results insertion point not found")
text = text.replace(old, new, 1)

path.write_text(text, encoding="utf-8")
print(path)
