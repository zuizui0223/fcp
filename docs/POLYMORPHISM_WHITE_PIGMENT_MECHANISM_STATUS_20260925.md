# White/pigmented proximal-mechanism audit status — 2026-09-25

## Current conclusion

The strongest commonality is now **proximal pathway convergence, not a universal ecological cause and not a universal evolutionary direction**.

Across the natural-FCP literature screen, at least **five independent natural systems** have direct molecular/biochemical evidence that white or strongly depigmented floral morphs map to reduced, blocked or rerouted anthocyanin/flavonoid pathway output:

1. *Gymnadenia rhellicani* — R2R3-MYB stop variant -> altered ANS expression -> reduced cyanidin;
2. *Ipomoea purpurea* — naturally occurring CHS loss-of-function and bHLH regulatory variants;
3. *Parrya nudicaulis* — petal-specific CHS downregulation / cis-regulatory control;
4. *Pleione limprichtii* — FLS / ANS / UFGT and candidate MBW-network expression differences;
5. *Silene littorea* — >42-fold F3h expression difference with candidate Myb1a regulation and matching flavonoid-profile evidence.

The five systems do **not** implicate one universal gene. They instead converge on the same biochemical network at different structural and regulatory levels.

## Evidence strength is not identical across the five systems

- *Gymnadenia* and *Ipomoea* provide especially direct genetic evidence.
- *Parrya* provides strong regulatory evidence at the pathway threshold.
- *Pleione* and *Silene* provide transcriptomic/metabolomic convergence with candidate regulators; these should not be described as if one causal nucleotide change had been proven.

Machine-readable manual review is split into:
- `results/polymorphism_natural_fcp_pigment_audit_20260925/manual_adjudication.csv` — initial A/B review;
- `results/polymorphism_natural_fcp_pigment_audit_20260925/candidate_c_rescue_adjudication.csv` — manual recovery of strong candidate-C cases.

## Important biological bounds

### Betalain counterexample

*Abronia fragrans* is a natural white/pink system in a betalain-producing lineage. Therefore **anthocyanin loss cannot be promoted to a universal angiosperm white-flower mechanism**.

The defensible higher-level statement is:

> Multiple molecular routes repeatedly reduce or reroute floral pigment output, often through the anthocyanin/flavonoid network in anthocyanic taxa, producing a recurrent achromatic endpoint.

### Heat can act directly in particular systems, but the comparative heat signal did not replicate

*Moricandia arvensis* provides a natural abiotic-plasticity example: spring flowers are lilac whereas summer flowers are white, and the colour/pigment shift is heat-linked.

The FCP comparative programme, however, does **not** support BIO5 as a common environmental cause:

- third cohort: broad long-term BIO5–white association passed its original gate;
- seasonal-heat follow-up: species-level and species-stratified signals were positive, but the prespecified same-cell local criterion failed narrowly (p=0.054995);
- original discovery cohort replication: 271 species, median within-species BIO5 contrast = -0.0068 SD, Wilcoxon p=0.743, conditional OR=1.033, p=0.131;
- original reserve cohort replication: 260 species, median contrast = +0.0662 SD, Wilcoxon p=0.0541, conditional OR=1.048, p=0.0319;
- the frozen legacy replication verdict is `LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST`.

Thus heat is a **context-dependent modifier/plasticity mechanism in some systems**, not a replicated global explanation of the recurrent white/nonwhite axis.

## Pollinator discrimination

Independent pollinator tests do not identify a common global pollinator cause:

- Sphingidae and Lepidoptera predictions were unsupported;
- effort-adjusted bee breadth was unsupported;
- stable local bee richness/Bombus predictors were unsupported as common white-state predictors.

The earlier third-cohort joint models retained BIO5 when bee predictors were added, but the subsequent discovery/reserve BIO5 replication failed. Therefore those joint models should not be read as evidence that heat is the universal alternative to pollinators.

The *Gymnadenia* case shows that pollinators can maintain a colour polymorphism in a particular system, but that is **system-specific mechanism evidence**, not a global white-flower rule.

## Evolutionary direction is separate and negative

The independent OpenTree endpoint test retained 196 species (116 nonwhite-dominant, 80 white-dominant). Across all 100 topology resolutions:

- median q(nonwhite -> white) / q(white -> nonwhite) = **0.701**;
- fraction with ratio > 1 = **0**;
- the prespecified nonwhite->white asymmetry gate failed.

Therefore the recurrent white/nonwhite axis must **not** be rewritten as a confirmed repeated pigmented->white evolutionary transition.

## Working two-level model

The current evidence supports:

1. **phenotypic accessibility / proximal convergence** — many molecular or regulatory routes can alter pigment-pathway output and reach an achromatic endpoint;
2. **context-dependent ecological sorting/maintenance** — temperature, herbivory, pollinators and other agents can favor or induce colour states in particular systems, but no single global ecological driver is confirmed.

This is stronger than the original “pollinator versus environment” dichotomy because it explains why the same white/nonwhite phenotype-space axis can recur even when the selective context differs among species.

The next discriminating test is therefore molecular rather than another broad abiotic/pollinator screen: ask whether **white-involving natural polymorphisms are disproportionately generated by pathway-loss/downregulation mechanisms relative to nonwhite hue–hue polymorphisms** among systems with direct molecular evidence.

The frozen New Phytologist manuscript remains unchanged.
