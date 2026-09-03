# JBI ecological result synthesis — spatial organization of flower-colour variation

Status: ecological interpretation ledger, 2026-09-03.  
Purpose: make the paper result-led rather than method-led while preserving every frozen inferential boundary.

## 1. Ecological question

The central response is not flower colour itself. It is **how intraspecific flower-colour variation is organized in space**:

- local coexistence within populations (C);
- persistent geographic differentiation among places (S);
- both signals (C+S), where separately retained;
- unresolved documentation is never biological absence.

The ecological problem is therefore: **what kinds of ecological setting are associated with local maintenance of multiple colour variants versus geographic sorting of colour variation?**

This framing unifies the comparative and image-based results without numerically pooling their distinct samples.

## 2. Primary comparative hypothesis test — broad climatic sorting versus broad-niche coexistence

The frozen 34-species comparison provides the strongest current ecological test because the two competing expectations were explicit and the five climatic summaries were evaluated symmetrically.

### H1 — broad-gradient environmental sorting

If geographic differentiation commonly reflects colour variants sorting along broad macroclimatic gradients, S species should occupy climatic niches at least as broad as, and potentially broader than, species in which colour variants coexist locally.

**Prediction:** `S >= C` in species-level climatic breadth.

### H2 — broad-niche coexistence / niche-variation expectation

If broad ecological occupancy is associated with the maintenance of intraspecific diversity, while geographic differentiation can arise within narrower species-level climatic settings through fine-scale selection, historical isolation, restricted dispersal or drift, C species should occupy broader climatic niches than S species.

**Prediction:** `S < C` in species-level climatic breadth.

### Result

All five frozen climatic-niche estimates point in the H2 direction:

| Climatic summary | S vs C odds ratio | Ecological direction |
|---|---:|---|
| Temperature breadth | 0.817 | S narrower |
| Moisture breadth | 0.412 | S narrower |
| Climatic heterogeneity | 0.681 | S narrower |
| PCA dispersion | 0.712 | S narrower |
| PCA hull area | 0.577 | S narrower |

Moisture breadth is the largest observed contrast: OR = **0.412**, 95% CI **0.180–0.947**, family-clustered Wald `p = 0.0368`, permutation `p = 0.0423`. After Holm adjustment, the threshold evidence weakens (`p = 0.184` Wald; `p = 0.212` permutation), so moisture is not promoted as a unique mechanism.

Every leave-one-family-out estimate remains below one. Open Tree and time-scaled V.PhyloMaker2 models also retain negative point estimates for all five metrics, although their confidence intervals include one. Moisture-breadth CR2/Satterthwaite inference gives OR = **0.407**, `p = 0.063`.

### Ecological decision

**The observed comparative pattern is more consistent with H2 than with the simple H1 broad-gradient sorting expectation.**

The strongest result is the repeated direction across distinct climatic summaries, family deletion and phylogenetic treatments. The data do not establish that broad climatic niches cause local coexistence, nor that moisture is the causal agent.

A biologically important consequence is that **geographic flower-colour differentiation does not appear to require broad species-level macroclimatic breadth**. Geographic differentiation may instead arise inside relatively narrow macroclimatic niches through environmental differences at finer scales, spatially varying biotic interactions, demographic history, isolation or drift.

## 3. Mechanistic discrimination on the expanded display-core-v6

A later literature-driven screen used 74 core flower-colour-variation species, of which 66 were climate-eligible and 32 had positive C/S organization evidence. These tests are exploratory relative to the frozen 34-species comparison, but they directly ask which ecological explanation remains plausible.

### 3.1 Static temperature seasonality — candidate, not established

Mean occupied BIO4 initially distinguished S from C:

- S vs C OR = **0.349**, family-bootstrap 95% CI **0.085–0.909**;
- C+S vs C OR = **0.401**, CI **0.140–1.002**.

This is the only v1 screen whose S-vs-C bootstrap interval excluded one. However, the contrast is fragile to broad geographic position:

- adding mean BIO1: OR = 0.412, CI 0.100–1.227;
- adding median absolute latitude: OR = 0.425, CI 0.100–1.120;
- adding both: OR = 0.458, CI 0.112–1.342.

**Decision:** temperature seasonality is a candidate clue, not an established temperature mechanism. Its apparent signal can partly encode latitudinal/thermal geography.

### 3.2 Precipitation seasonality — prediction not supported

For S vs C, OR = **1.50**, CI **0.41–5.15**. This does not support the prediction that stronger precipitation seasonality favours local coexistence.

### 3.3 Geographic fragmentation — direction plausible, evidence unresolved

Fragmentation estimates point toward greater geographic differentiation but remain uncertain:

- S vs C OR = **1.80**;
- C+S vs C OR = **1.75**;
- confidence intervals are broad and include one.

**Decision:** fragmentation remains biologically plausible but is not supported as an explanatory axis in the present data.

### 3.4 Large-scale spatial climatic partitioning — not supported

Three independent diagnostics fail to show that S species have stronger large-scale climatic partitioning than C species:

- pairwise spatial environmental turnover: S vs C OR = **1.37**, CI **0.49–7.51**;
- normalized two-sector climate-centroid separation: OR = **1.62**, CI **0.56–3.48**;
- two-sector Gaussian climate overlap: OR = **1.00**, CI **0.48–3.36**.

**Decision:** the simple mechanism `persistent large-scale climatic partition -> S` is not supported by the current WorldClim/GBIF representation.

This negative discrimination is ecologically important. The primary 34-species pattern cannot simply be redescribed as geographically differentiated species occupying strongly separated macroclimates across their ranges.

### 3.5 Total climatic hypervolume size — not a useful discriminator

Rarefied PC1–PC3 hull volume gives:

- S vs C OR = **2.08**, CI **0.46–6.79**;
- C+S vs C OR = **0.47**, CI **0.26–1.16**.

**Decision:** total niche volume alone does not distinguish the spatial organization states. This reinforces the ecological shift from `how large is the niche?` toward `how is environmental variation distributed through space and time?`.

## 4. Dynamic temporal heterogeneity — decisive climatic follow-up remains not evaluable

Before dynamic climate outcomes were inspected, the TerraClimate analysis fixed 1,320 locations (66 species × 20 geographically balanced occupied cells), the 1958–2025 source interval, seven environmental metrics, directions, model family, family-bootstrap design, Holm family and seed.

The direct hypothesis was:

> stronger **temporal** environmental heterogeneity favours local coexistence C over persistent spatial segregation S.

Five direct temporal metrics were fixed, with a pre-result supported/unsupported/not-evaluable decision rule. The exact run `33599249288` did **not** produce a biological result. The fixed-point key passed, but the 6-hour job was cancelled during TerraClimate extraction after 1,000 / 6,545 point-variable requests; the prospective model step never ran and no dynamic result artifact exists.

**Decision: H1-dynamic = NOT_EVALUABLE because of incomplete data extraction, not biologically unsupported.**

No metric value may be inferred from that run, and no source, point set, direction, threshold or response state may be substituted post-result.

## 5. Biotic and life-history hypotheses — current evidence ceiling

The prospective mechanism registry also considered mating system, dispersal, phenology and biotic turnover. These hypotheses are biologically meaningful but the current independent comparative coverage is insufficient for promotion.

| Hypothesis | Prediction | Current decision |
|---|---|---|
| Pollinator geographic mosaic | persistent interaction turnover favours S; temporal turnover may favour C | **not evaluable** — no independent range-wide interaction-turnover dataset |
| Multi-agent balancing selection | varying pollinator/antagonist balance favours local coexistence or mosaics | **not evaluable** — inadequate comparative fitness-component data |
| Self-fertilization / effective gene flow | greater selfing shifts organization toward S | **not evaluable** — GIFT informative pure C/S coverage is far below the frozen gate |
| Dispersal × fragmentation | low-connectivity dispersal amplifies fragmentation -> S | **not evaluable** — dispersal coverage below gate; cannot replace with a post-hoc main effect |
| Flowering duration × temporal variability | temporal exposure across flowering season modifies C vs S | **not evaluable** — flowering coverage below gate and dynamic climate itself is NE |
| Drift / founder history | weak connectivity without environmental correspondence favours S | **not evaluable** — no comparative neutral population-genetic dataset |

These NE outcomes are useful boundaries, not ecological nulls.

## 6. Integrated ecological interpretation

The current evidence therefore supports a narrower and more interesting ecological story than a generic `flower colour tracks climate` claim:

1. **Intraspecific flower-colour variation has a real spatial organization.** In the separate six-species held-out image analysis, geographically neighbouring observations were more similar in continuous colour than expected under species-conditioned random labelling (`p = 0.0113`).
2. **That organization is not captured by one universal global transition geography.** The six-species shared-transition primary test was not confirmed (`p = 0.0906`).
3. **Across species, geographic differentiation tends to occur in narrower sampled climatic niches than local coexistence.** All five 34-species climatic breadth/heterogeneity estimates are below one, contrary to the simplest broad-gradient sorting expectation.
4. **Simple large-scale spatial climatic partitioning does not explain that contrast.** Environmental turnover, regional centroid displacement and regional climatic overlap are all unresolved/null in the expanded core.
5. **Total niche size is also insufficient.** The ecological object is more likely the allocation of environmental variation through space and time than total hypervolume alone.
6. **Temporal heterogeneity remains the clearest unresolved climatic mechanism**, but the exact prospective dynamic test is currently NE because extraction failed before modelling.
7. **Pollinator, mating-system, dispersal and neutral-history mechanisms remain open**, not rejected, because the present cross-species datasets do not meet the predeclared coverage requirements.

## 7. Manuscript-level claim

Recommended result-led statement:

> Across independent evidence streams, flower-colour variation was spatially organized within species, but geographic differentiation was neither concentrated along one universal global boundary nor associated with stronger large-scale climatic partitioning. In the frozen comparative sample, geographically differentiated cases instead occurred toward the narrower end of sampled species-level climatic niches. This pattern is inconsistent with a simple model in which broad macroclimatic gradients generally generate geographic colour sorting and points toward finer-scale, temporal, biotic or historical processes as the next mechanistic level.

The final clause is a hypothesis implication, not a demonstrated causal mechanism.

## 8. What should be foregrounded in the JBI paper

Main ecological results, in order:

1. competing H1/H2 test and the five concordant ORs;
2. moisture as the strongest observed axis but not unique mechanism;
3. rejection/non-support of the simple large-scale spatial climatic partition explanation;
4. six-species image result as independent evidence that colour is genuinely spatially organized within species but not along one common boundary;
5. dynamic temporal heterogeneity as the single highest-priority unresolved ecological test;
6. terminal 200-species atlas NE only as a measurement-transfer limitation, not as the paper's scientific centerpiece.

Methods, provenance, frozen gates and exact hashes remain essential for credibility but should support this ecological argument rather than dominate the Results narrative.
