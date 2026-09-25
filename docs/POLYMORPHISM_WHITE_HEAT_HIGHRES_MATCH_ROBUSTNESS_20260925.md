# White–heat high-resolution local matching robustness — 2026-09-25

Status: **post-result robustness analysis**.

The first local matching analysis used WorldClim 2.1 at 10 arc-minute resolution. That grid is too coarse to strongly discriminate many pairs separated by 25–50 km. This analysis repeats the already frozen geographic matching algorithm with **WorldClim 2.1 BIO5 at 2.5 arc-minute resolution**.

Nothing in the prospective BIO5 mechanism gate is changed. No matching threshold, colour eligibility rule, high-clip definition or model is retuned.

## Fixed analysis

Reuse `scripts/analysis/run_white_heat_local_match_20260925.py` unchanged.

Inputs:
- same third-cohort biological artifact;
- same response-blind high-clip exclusion and near-clip metric;
- WorldClim 2.1 `wc2.1_2.5m_bio_5.tif`.

Report the same fixed thresholds:
- 25 km;
- 50 km;
- 100 km;
- 250 km.

## Interpretation

Persistence of a positive white-minus-nonwhite BIO5 difference at 25–50 km would strengthen a local thermal-sorting interpretation.

Failure at 2.5 arc-minute resolution would strengthen the conclusion that the prospective 10 arc-minute BIO5 association is predominantly broad-scale geographic sorting rather than a local temperature contrast.

This remains post-result robustness and cannot establish causal heat selection or evolutionary transition direction.
