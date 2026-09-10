# Polymorphism paper spine v0.1

Date: 2026-09-10 JST

## Working title

**Flower-colour polymorphism is genus-clustered and geographically organized within species at the global scale**

Alternative, stronger conceptual title:

**Species-specific geography without a universal boundary: global structure of flower-colour polymorphism**

## One-sentence thesis

Global flower-colour polymorphism is a genus-clustered species property whose magnitude is associated with species-specific geographic colour organization, even though species do not share a detectable universal colour boundary.

## Paper logic

### Question 1 — Is flower-colour polymorphism a measurable species property?

Yes. Use four-state Simpson diversity D over the fixed classifiable flower-colour states. In the 369-species discovery frame D is continuous from ~0 to 0.708; 46.61% have a second morph >=10% and 26.56% >=20% within the admitted frame.

The representation caveat is explicit: mixed/uncertain rows are not reclassified post hoc, so D is conservative.

### Question 2 — Is that species property taxonomically organized?

Yes, at the genus level rather than the family level.

Discovery:

- genus clustering gain ~20.3%, six-family Holm p=0.0115;
- family clustering not supported.

Species-disjoint reserve:

- 0 species overlap with discovery;
- genus clustering gain ~15.3%;
- reserve recurrence Holm p=0.0456.

Span-only adjustment leaves genus clustering essentially unchanged in both tranches.

This is the cleanest answer to **which kinds of species carry similar levels of polymorphism?**: closely related species at the genus scale resemble one another more than expected, while the signal does not extend clearly to families.

### Question 3 — What does greater polymorphism look like geographically?

More polymorphic species tend to show stronger within-species geographic organization of flower colour.

Discovery and reserve both show positive continuous D-spatial relationships and stronger spatial organization in the fixed >=10% polymorphic subset. The reserve reproduces the pattern under observer-pair exclusion, calendar-quarter stratification, and matched-background differential controls.

The species sets are disjoint.

### Opportunity robustness

Sampled geographic span does not attenuate the continuous effect size:

- discovery raw rho 0.0892 -> span-only partial rho 0.1032;
- reserve raw rho 0.1016 -> span-only partial rho 0.1019.

But the two span-only permutation tests are individually borderline (~0.049 and ~0.050) and their two-test Holm p is ~0.098, so the multiplicity-corrected robustness family does not pass.

Post-classification yield is a separate limitation. `n_classifiable` is strongly negatively associated with D (~-0.54 to -0.56); conditioning on it weakens the reserve effect. Because raw photo opportunity is fixed at 100/species and n is measured after the colour-classification gate, this should be reported as measurement-yield sensitivity rather than treated as a clean pre-colour confounder.

### Question 4 — Do species share the same global colour boundary?

No supported universal boundary is recovered, and the prospective sharedness-v2 qualification failed identification power. This is not an awkward negative result; it defines the scale of the positive result.

The recoverable regularity is:

> **how strongly each species is internally structured**

not:

> **where all species change colour together**.

That scale distinction is the conceptual payoff of the paper.

## Recommended Results order

### Result 1. Polymorphism forms a continuous species-level gradient

Show the 369-species D distribution, second-morph frequencies, and reserve/discovery measurement design. Keep M1-M3 out of the main result except for a short representation diagnostic.

### Result 2. Polymorphism is clustered within genera but not families

Lead with discovery, then species-disjoint reserve recurrence. Put span-only and D_unbiased robustness immediately after the primary result.

### Result 3. Greater polymorphism is associated with stronger within-species geographic organization

Lead with the discovery and reserve threshold contrasts, then the continuous D-spatial gradients. Follow with observer/time/background controls.

### Result 4. The pattern is species-specific rather than a shared global boundary

Use the failed sharedness qualification and earlier common-boundary result as an estimand contrast, not as another large analysis section.

### Result 5. Other candidate species attributes are not generalizable or not testable

Keep brief:

- sampled span positive in discovery but null in the species-disjoint reserve -> discovery-specific;
- latitude unsupported;
- pollination and life form fail source-coverage gates and cannot be interpreted biologically.

## Recommended figure architecture

### Figure 1 — A global species gradient of flower-colour polymorphism

A. schematic of four admitted colour states and D calculation;
B. ranked species distribution of D;
C. histogram/density of second-morph fraction;
D. discovery/reserve design and species-disjoint split.

Main message: polymorphism is measurable continuously rather than as an anecdotal list of polymorphic species.

### Figure 2 — Genus-level taxonomic clustering

A. observed equal-genus-weighted within-genus D difference against permutation null in discovery;
B. same in reserve;
C. effect-size comparison: discovery raw, span-only, D_unbiased; reserve raw, span-only, D_unbiased;
D. family-level null as a compact negative control.

Main message: species resemble congeners in polymorphism magnitude, and this recurs in a disjoint species set.

### Figure 3 — More polymorphic species are more geographically organized

A. discovery species D versus spatial rho;
B. reserve species D versus spatial rho;
C. >=10% polymorphic versus complement spatial effect in discovery and reserve;
D. reserve robustness estimates: primary, observer-pair exclusion, quarter stratification, matched-background differential.

Main message: polymorphism magnitude covaries with species-specific geographic structuring.

### Figure 4 — The scale of regularity

Conceptual/result synthesis:

- left: no identifiable universal cross-species colour boundary;
- center: each species can have its own geographic colour organization;
- right: strength of that organization increases with species-level polymorphism and polymorphism itself clusters by genus.

This figure should be a synthesis, not another hypothesis test.

## Discussion structure

### 1. The unit of regularity is the species, not a universal map boundary

The failed common-boundary route is reframed as a scale-identification result. Cross-species geography need not align for within-species organization to recur.

### 2. Genus clustering implies inherited or conserved constraints, without claiming phylogenetic signal

Possible mechanisms can be discussed as hypotheses: shared pigment architecture, developmental constraints, mating-system context, pollinator interactions, or common demographic histories. None is identified by the present analysis.

### 3. Geographic organization links polymorphism magnitude to spatial differentiation

Discuss mechanisms maintaining morphs differently across local environments, gene flow/selection balance, and spatially varying biotic interactions as candidate explanations. Keep them clearly mechanistic hypotheses rather than inferred causes.

### 4. Sampling and measurement limits define the current ceiling

Explicitly state:

- raw photo opportunity fixed at 100/species;
- sampled-span correction leaves effect sizes positive but multiplicity-corrected diagnostic is borderline;
- classifiable yield is strongly associated with D and produces reserve sensitivity;
- pollination/life-form global trait coverage is insufficient;
- no explicit phylogeny or population genetics is used.

## Claims not to make

Do not claim:

- a universal flower-colour transition zone;
- privileged M1-M2 polymorphism directions;
- causal effects of geographic range size;
- phylogenetic signal from genus/family labels alone;
- adaptation or environmental selection from spatial rho;
- pollination or life-form effects from failed trait-coverage gates;
- global prevalence of polymorphism from the selected 369-species frame.

## Current strongest novelty statement

The novelty is not merely that flower-colour polymorphism exists globally. It is the joint scale result:

> **Polymorphism behaves as a repeatable species attribute, resembles itself among congeners, and predicts the strength of species-specific geographic colour organization, while no common global boundary is required.**
