# Stable local bee versus heat status — 2026-09-25

## Completed discriminator

The pre-specified stable local bee versus heat workflow completed successfully in GitHub Actions run **36105413871** (artifact **10850813548**).

The analysis used the versioned Noori et al. (2026) curated GloBI bee–plant dataset and WorldClim BIO5. Local bee predictors were constructed outcome-blind before morph labels were opened.

## Outcome-blind coverage gate

Among 499 third-cohort flower species, 332 had exact species-level coordinated bee records after source filtering.

Qualified species by fixed candidate radius:
- 50 km: 20
- 100 km: 39
- 250 km: 92
- 500 km: 154

Per protocol, **500 km** was therefore selected as the smallest radius with >=100 coverage-qualified species. This is a broad regional bee-community scale, not a local flower-visitor measurement.

## Joint model B1 — heat + local bee richness

Informative panel:
- 167 species
- 5,836 flower records

BIO5:
- OR = **1.138** per within-species SD
- 95% CI = **1.061–1.220**
- p = **0.000300**

Local bee richness:
- OR = **0.946**
- 95% CI = **0.878–1.020**
- raw p = **0.149**
- Holm pollinator p = **0.297**

Thus the previously observed positive heat association persists after adjustment for stable local bee richness, whereas bee richness is not supported.

## Joint model B2 — heat + Bombus fraction + bee richness

Informative panel:
- 136 species
- 4,828 flower records

BIO5:
- OR = **1.155**
- 95% CI = **1.069–1.247**
- p = **0.000266**

Local Bombus fraction:
- OR = **1.019**
- 95% CI = **0.941–1.103**
- p = **0.647**

Local bee richness:
- OR = **0.928**
- 95% CI = **0.857–1.005**
- p = **0.0663**

Species-level audits likewise do not support recurrent white/non-white differences in local bee richness (p = 0.239) or Bombus fraction (p = 0.783).

## Interpretation

The stable joint test does **not** support the hypothesis that the broad BIO5–white association is an artefact of the measured local bee community.

Instead:
- BIO5 remains positive and highly resolved after bee adjustment;
- local bee richness is weak-to-null;
- Bombus composition is null.

The best current comparative interpretation is therefore:

> **Measured white flower states show a repeatable broad-scale association with warmer warmest-month climates that is not explained by the available stable bee-community predictors.**

This strengthens temperature relative to bee turnover as the current common environmental clue, but it still does not establish a causal heat-driven whitening mechanism.

The earlier same-observer and <=25–250 km white/non-white matching sensitivities remain critical: those tests did not retain the BIO5 association. Hence the evidence supports **broad geographic temperature sorting**, not universal local heat selection.

## Mechanistic implication

Together with experimental and physiological literature showing that high temperature can reduce anthocyanin accumulation or shift floral pigment composition, the FCP result motivates a specific next hypothesis:

**temperature-sensitive suppression or loss of pigmentation is a recurrent proximal route to white/pale floral states, while geography and lineage determine whether that route is realized and maintained.**

That hypothesis still requires an independent direction/genetic test. The present data do not identify whether white is derived, whether the effect is plastic or genetic, or whether anthocyanin loss is the relevant pathway in each species.

The frozen New Phytologist manuscript remains unchanged.
