# Frequency dependence as a candidate mechanism of geographic flower-colour variation

Date: 2026-11-26. This separate branch develops a new biological hypothesis without changing the main FCP manuscript or reinterpreting photo-region labels as validated petals.

## Model
Let p be local frequency of a biologically verified white phenotype. Let the white-minus-coloured relative fitness difference be Δw = a(E) + b(1−2p), with b>0 for negative frequency dependence. In a simple haploid/phenotypic frequency dynamics model, dp/dt = p(1−p)Δw. The equilibrium p* = (1+a(E)/b)/2 is interior only when |a(E)|<b; when directional selection outweighs rare-morph advantage, fixation becomes possible. The equilibrium need not be 0.5. This is an illustrative hypothesis, not parameter estimates or a novel claim that opposing selection alone is new. For b<0 an interior crossing is unstable in this model. For b=0, the model is directional or neutral.

Predictions: environmental shifts in a(E) may move the equilibrium; negative frequency dependence can buffer loss of alternative phenotypes; positive frequency dependence can reinforce regional dominance. Pollinator switching in rewardless species versus colour constancy in rewarding species is a candidate source of variation in b. Fitness rather than visits alone is required.

## First direct data test
The frozen Table S8 has 104 revisited species-cell-calendar-month windows (73 species), 275 initial and 283 later observer-thinned photos, with no shared observers between two chosen years within each window. Of these, 66 windows (50 species) have both white and nonwhite labels across the pooled years. This previously restricted panel was defined by the overall repeated-colour criterion; the new primary test reconstructs equivalent windows from all 4,475 species without that colour-based inclusion gate.

For n0>=2 initial independent observer records, let p0=k0/n0, p1=k1/n1. The naive product (1−2p0)(p1−p0) is positively biased under constant true frequency because the same noisy p0 occurs twice. Averaging (1−2A)(p1−C) over distinct ordered initial-observer labels A,C yields T=(1−2p0)(p1−p0)−2p0(1−p0)/(n0−1). Under constant exchangeable local photographic frequencies, the conditional allocation of the pooled white counts is hypergeometric and E(T)=0. Use equal species weights and 19,999 independent conditional window randomizations with fixed sample sizes; the alternative is positive average convergence toward 0.5. This is one direct temporal test, not a new round of speculative sensitivity analyses.

Observed from the actual frozen Table S8 (Git blob 450b0df91eb8ba343c48f8d48732223c9191f08a):
- Naive score +0.1273795 (dimensionless; not a selection coefficient).
- Corrected primary score +0.0069920.
- Null central 95% [−0.076094,+0.077802]; conditional one-sided p=0.4218.
- No evidence of symmetric rare-photo-colour recovery in this small panel. The direction toward 0.5 is not a prediction of all possible negative-frequency-dependent equilibria. No biological fitness effect or absence of selection has been established.

## Next decisive measurement
After establishing genuine heritable floral variants, manipulate their local frequencies (e.g. 10:90, 50:50, 90:10) at fixed density across environments, measure maternal seed/fruit production and paternal siring or pollen export, and estimate the slope of relative fitness against starting frequency. A sign reversal in relative fitness and invasion when rare support a protected polymorphism. Photo geography can prioritize species and regions, but cannot estimate a(E), b or establish protected coexistence by itself.

Related evidence:
- Gigord, Macnair & Smithson (2001), PNAS, doi:10.1073/pnas.111162598 — manipulative rare-colour reproductive advantage in a rewardless orchid.
- Tenhumberg, Dellinger & Smith (2023), Journal of Ecology, doi:10.1111/1365-2745.14057 — theory combining pollinator and nonpollinator selection on pigmentation.
