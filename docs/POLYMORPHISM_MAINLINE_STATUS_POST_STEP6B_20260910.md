# FCP polymorphism mainline after Step 6b

Date: 2026-09-10 JST

## Mainline decision

The paper should now be organized around three empirical facts:

1. flower-colour polymorphism is a measurable species-level property;
2. more polymorphic species tend to show stronger within-species geographic organization of flower colour;
3. polymorphism is clustered at a shallow taxonomic scale (genus), and that genus-level signal recurs in a species-disjoint reserve set.

The discovery-only positive sampled-span association is no longer a main species-attribute claim because it does not recur in the reserve set.

The old globally shared-boundary and privileged recurrent-colour-direction hypotheses remain rejected/unsupported and should be used only as contrasts that sharpen the estimand.

## Claim 1 — polymorphism is a measurable species-level property

Retain unchanged.

Discovery contains 369 eligible species with four-state Simpson diversity D spanning approximately 0 to 0.708; 46.61% have a second morph >=10% and 26.56% >=20% in the admitted/classifiable frame. These are properties of the sampled frame, not global prevalence estimates.

`mixed_uncertain` is structural/technical missingness plus a smaller chromatically ambiguous component. It cannot be numerically reclassified as hidden polymorphism after the frozen Step-2 gate failed. D remains a conservative four-state summary of admitted rows.

## Claim 2 — greater polymorphism is associated with stronger within-species geographic colour organization

Retain as a major claim, but with the Step-6/6b robustness boundary stated explicitly.

### Original discovery and reserve evidence

Discovery Step 5:

- >=10% second morph: 172 species, mean spatial rho 0.034565108, exact-randomization p=0.001;
- <10% complement: 197 species, mean rho 0.020434835;
- dilution contrast = 0.014130273, p=0.014;
- continuous rho(D, species spatial rho)=0.089213, p=0.034.

Reserve Step 5b:

- >=10% second morph: 151 species, primary mean rho 0.032299123, p=0.001;
- observer-pair exclusion p=0.001;
- quarter-stratified p=0.001;
- matched-background differential mean rho 0.013621726, p=0.005;
- >=10% minus <10% dilution contrast = 0.011671678, p=0.041;
- continuous rho(D, species primary spatial rho)=0.101601, p=0.025.

Discovery and reserve species overlap is **0**, so the recurrence is species-disjoint as well as photo-tranche-disjoint.

### Sampled-span robustness

Step 6 joint span+n diagnostic failed its deliberately strict two-tranche Holm gate because the reserve partial effect weakened when `n_classifiable` was conditioned on.

Step 6b separated the controls.

Span-only partial rank effects are:

- discovery: 0.103241, raw permutation p=0.0489476;
- reserve: 0.101872, raw permutation p=0.0500975.

Their Holm p-values are both 0.0978951, so the prospectively frozen two-tranche span-only gate is formally **not passed**.

However, sampled-span adjustment does not attenuate the effect sizes: discovery moves from raw 0.0892 to 0.1032 and reserve remains essentially unchanged from raw 0.1016 to 0.1019. Therefore sampled geographic span is not empirically behaving like an effect-size explanation of the D-spatial association, although the multiplicity-corrected robustness family is borderline and cannot be labeled confirmatory support.

The finite-sample-unbiased Simpson sensitivity gives essentially the same span-only partial effects (~0.103 discovery, ~0.102 reserve) and the same Holm outcome.

### `n_classifiable` is a different kind of variable

Raw photo opportunity is fixed at 100 photos/species in both discovery and reserve. `n_classifiable` is therefore post-classification yield, not raw sampling effort.

It is strongly related to polymorphism:

- discovery rho(D,n) = -0.537139;
- reserve rho(D,n) = -0.555901;
- the negative association persists for D_unbiased.

Conditioning on n strengthens the discovery D-spatial partial effect (0.1564) but weakens the reserve effect (0.0779, p~0.14). This sensitivity must be reported, but it should not be treated as equivalent to adjustment for pre-colour sampling opportunity: n may itself reflect colour-classification difficulty or complexity associated with polymorphism.

### Admissible wording

Use wording such as:

> Greater flower-colour polymorphism was associated with stronger within-species geographic organization in two species-disjoint global photo frames. The effect-size pattern was not attenuated by sampled geographic span, although a conservative multiplicity-corrected span-adjusted diagnostic was borderline; results were more sensitive to post-classification yield.

Do not claim that all measurement-opportunity bias is excluded.

## Claim 3 — polymorphism shows shallow taxonomic clustering

Promote this to the main species-attribute result.

### Discovery

Step 4 genus clustering:

- 61 repeated genera / 169 species;
- clustering gain ~20.3%;
- raw p=0.00229989;
- six-family Holm p=0.0114994.

Family-level clustering was not supported (gain ~3.6%, raw p~0.30), so the signal is shallow rather than a generic family-level taxonomic effect.

### Species-disjoint reserve recurrence

Step 6 reserve genus test:

- 54 repeated genera / 146 species;
- clustering gain ~15.3%;
- raw p=0.0227989;
- Holm p=0.0455977 in the two-attribute reserve recurrence family.

The reserve species set overlaps discovery by 0 species. Thus genus-level clustering recurs in an independent species set.

### Span decomposition

Genus clustering is not explained by sampled span.

Discovery:

- raw D gain 0.2031, p=0.00290;
- span-only residual-D gain 0.2113, p=0.00150;
- D_unbiased span-only gain 0.2123, p=0.00155.

Reserve:

- raw D gain 0.1545, p=0.0219;
- span-only residual-D gain 0.1521, p=0.0235;
- D_unbiased span-only gain 0.1541, p=0.0207.

Conditioning on post-classification n reduces the reserve clustering signal (gain ~0.085, p~0.11), so the same measurement-yield caveat applies. The span result itself is stable.

### Admissible wording

> Flower-colour polymorphism is taxonomically clustered at the genus level in two species-disjoint global frames, while family-level clustering is not supported.

Call this **taxonomic clustering**, not phylogenetic signal, because no explicit phylogeny is being modeled.

## Discovery-only species attributes

### Sampled geographic span

Discovery Step 4 showed rho(D, sampled span)=0.181031, Holm p=0.00479976, robust to D_unbiased and alternate span proxy.

But the reserve-photo recurrence is null:

- reserve rho(D, sampled span)=-0.002567;
- raw p=0.961002;
- Holm p=0.961002.

Therefore span must be downgraded from a general species-attribute claim to a discovery-specific association / sampling-frame diagnostic.

### Latitude

Not supported in discovery (rho approximately -0.034, Holm p=1). Do not emphasize.

### Pollination and life form

No biological null can be claimed because source gates failed:

- pollination: 144 source-backed species but overwhelmingly `mixed`, prospectively unusable for category comparison;
- direct GIFT Life_form_1: 2/369 FCP species, prospectively unusable.

These remain data-coverage limitations, not negative biological results.

## Rejected former claims

### Recurrent M1-M2 directions as privileged polymorphism directions

Reject. D is strongly related to total continuous colour variance, but M1-M2 variance share is not a privileged direction after total spread is controlled and does not outperform random two-dimensional orientations.

### Globally shared flower-colour boundary

Do not restore. Sharedness-v2 failed prospective identification power and prior cross-species shared-boundary analyses did not support a universal boundary.

This is now a useful contrast: **within-species geographic organization and genus-level taxonomic structure are recoverable without requiring species to share the same geographic boundary.**

## Current paper spine

### Claim A — species-level polymorphism

Flower-colour polymorphism is measurable as a continuous species attribute in a large global photo frame.

### Claim B — geographic organization

Species with greater polymorphism tend to have stronger within-species geographic colour organization, recurring in a species-disjoint reserve frame. Sampled geographic span does not attenuate the effect size, but strict multiplicity-corrected span robustness is borderline and post-classification-yield sensitivity must be disclosed.

### Claim C — shallow taxonomic organization

Polymorphism is clustered within genera, this genus-level signal recurs in a species-disjoint reserve set, and sampled-span adjustment leaves it intact; family-level clustering is not supported.

### Estimand contrast

The repeatable structure is **species-specific and genus-associated**, not a universal cross-species geographic boundary or a privileged global colour-space direction.

## Recommended one-sentence paper thesis

> Global flower-colour polymorphism is a genus-clustered species property whose magnitude predicts species-specific geographic colour organization, even though species do not share a detectable universal colour boundary.

The word `predicts` here is statistical/associational only; causal language should be avoided in the manuscript.
