# Sisyrinchium: does other white species co-occurrence predict M-purple sorting within recorded field lines?

**2026-10-13 JST. PROTOCOL LOCK before computing the new line-conditioned outcomes and null probabilities.** This is a retrospective, source-exposed extension of PR #182, NOT an external preregistration. Existing overall site-year and quadrat results (PR182) have been seen; no new field observation is gathered here.

## Biological target
The original Takahashi, Takakura & Kawata (2016) paper (PLOS ONE 11:e0164381; DOI 10.1371/journal.pone.0164381) reported that the purple M-species morph is more frequent where the two other white-flowered L- and S-species co-occur. The paper ALSO reported that density of those congeners and distance to the water edge were strongly correlated; it could not fit the two as independent causal explanatory factors. The **new question** is more restrictive:

1. **PRIMARY**: Across multiple named site-years, does the recorded purple M-species morph occur at quadrats having a higher **co-occurring L+S fraction WITHIN the same recorded line** than the white M-species morph, after conditioning on M-species colour inventory and counts in EACH (year,population,line) stratum?
2. **SECONDARY**: Does the strong PR182 within-site M-colour quadrat heterogeneity remain **within recorded lines** after fixing all line-specific M-colour totals, or is its apparent signal entirely from line-to-line shifts?

We will not retest an unconstrained GLM and call it an original finding. Positive association is not proof of interspecific reproductive interference; other species may covary with wetness, edge distance, microhabitat, density and species introductions. Negative/unsupported does not refute the original broader between-line association.

## Fixed original source and units
Takahashi et al. original Dryad DOI **10.5061/dryad.j8q34**, legitimately public Zenodo mirror **record 5003124**, exact `data.zip` size **127502 bytes**, MD5 **008ec196aa5acb769dbfd36e835e32ba**. Extract ONLY `data/frequency.txt` as a table (not `__MACOSX/._*` AppleDouble metadata). Original internal SHA256 **a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e**; exact 512 field quadrats in 16 (year,population) groups, named years **Y11/Y13**. Fields `v,w` are purple/white **within the polymorphic M-species**. Fields `l,s` are *individuals of other white-flowered L/S species*; **never append them to M-white counts**. `line` is the original row's recorded line identifier, which is only a coarse spatial grouping, NOT a verified independently replicated population. A line with one quadrat is noninformative for within-line tests. `distance` is the original distance from water's edge (units follow source; do not invent GPS metre point coordinates). 
- Preserve 512 rows, siteyear, line, quadrat IDs including M-absent quadrats (n=v+w=0) as `NO_TARGET_PLANTS`.
- For M-present quadrats set n=v+w, a=(l+s)/(n+l+s) as the **fraction of recorded plants in that quadrat belonging to the two other species**. This is the original paper's relative-abundance exposure; never condition on or analyze flower-hour counts as independent plants. This compositional exposure shares denominator n (a risk), but the conditional null FIXES n and the other-species counts exactly. No missing humidity values are converted to 0; this analysis does not model humidity.
- Plot keys `(year,population,line,quad)` must be unique, all v,w,l,s nonnegative integers. Sum original M colour counts and L/S counts unchanged.

## Source-oblivious eligibility
First apply the **same PR182** site-year eligibility `N_M >=20`, >=2 occupied M quadrats, >=5 each M morph. A site-year is **PRIMARY-INFORMATIVE** if at least TWO recorded lines each have >=2 occupied M quadrats, >=1 purple and >=1 white M plants in the line, >=10 M plants in the line, and **at least two different** `a` exposures among occupied quadrats; across these lines >=30 M plants and >=5 of each morph. Any failure is `HOLD_LINE_OR_EXPOSURE_SUPPORT`, NOT zero exposure/zero association. Excluded lines and site-years appear in ledgers, not silently omitted. Singletons (such as every `line` in Otto Y11) cannot test within-line contrast. Unlike the PR182 study-wide B, this analysis is explicitly conditional on the original line definitions.

## PRIMARY statistic and conditional null (association)
Within every informative line `h`, denote `N_h = Σ n_j`, `V_h = Σ v_j`, `p_h=V_h/N_h`; define `D_h=N_h p_h(1-p_h)>0`. Let other-species share `a_j` be fixed (site and quadrat size are fixed). Define the **within-line violet-vs-white difference in encountered L+S share**

`t_h = Σ_j (v_j-n_j p_h)*a_j / D_h`.

This identity is the mean a exposure among violet M plants minus mean a among white M plants inside that line. Within each site-year, average informative line scores using weights D_h, so `T_site = Σ_h D_h t_h / Σ_h D_h`. For study-level main summary, **EQUAL SITE-YEAR WEIGHT** across all eligible siteyears. Units are **percentage-point difference in L/S relative abundance**, not odds ratios, causal environmental effects, or flower colour phenotype diversity fraction.

Fixed-inventory NULL: independently allocate exactly V_h violet labels without replacement among the N_h M-plant slots in line h while preserving every quadrat's n_j; calculate t_h and the site-equal average. Use independent `numpy.random.Generator.multivariate_hypergeometric` per line with **R=9999**, seed **2026101309**. **One-sided positive hypothesis** predicted by earlier publication, plus-one `p=(1+#(null >= observed))/(R+1)`. p tests exchangeability *WITHIN lines conditional on all observed siteyear-line inventories*, not causal associations with an independent randomized competitor.

Potential confounding remains by continuous within-line distance. A **prelisted sensitivity** residualizes `a_j` by **n-weighted least-squares intercept+linear distance within each informative line**, keeping siteyear, line, n,V_h fixed. Recompute T with `a_residual`, same denominator D_h, and independent R=9999 line-conditional colour-label shuffles, seed **2026101310**. Lines with <3 M-occupied quadrats or no residual variation get excluded from distance sensitivity; require >=2 informative lines and >=30 M plants per included siteyear. This adjustment addresses **linear** water-distance trend, not nonlinear distance, habitat matching or unmeasured ecology. No posthoc selection of degree/order based on p.

## SECONDARY statistic and null (within-line clumping)
For **each** PR182-eligible siteyear, additionally require >=2 lines with >=2 occupied M quadrats, N_h>=2 and both M colours in line. Original `B_site=Σ_j n_j(p_j-p)^2/[N p(1-p)]` from PR182.
Given fixed line inventories `(N_h,V_h)`, the exact reference `B_linebetween=Σ_h N_h(p_h-p)^2/[N p(1-p)]` plus
`E[within-line B] = Σ_h N_h p_h(1-p_h)*(m_h-1)/(N_h-1)/[N p(1-p)]` (the latter 0 for N_h<=1 or m_h<=1).
Define `Δ_line=B_observed-[B_linebetween+E(within-line B)]`. Compute equal-site-year mean Δ_line, plus-one **one-sided** R=9999 site-line conditional multivariate-hypergeom null, seed **2026101311**. Report all line singleton/morph-monochrome holds. Expected E is analytic, not Monte Carlo plug-in.

For both primary/secondary, apply **Holm correction for the two planned primary-family comparisons** across the two predeclared endpoints. Distance-residual sensitivity is explicitly exploratory and not another independently discovered confirmatory test.

## Robustness and output ledger
- An **a priori sensitivity** restricts to **latest year per named population**, to avoid giving repeated places two votes. If <3 eligible distinct population names, mark HOLD rather than computing inferential p.
- Preserve all 16 source siteyears and every 512 quadrat, no multiplication of 13 siteyears as evolutionary replicates. Present sign count among siteyears (not as significance test), highest-leverage siteyear, and a **posthoc leave-one-siteyear-out audit**. If the primary result is driven by 1 site, downgrade species-wide generality.
- Explicitly compare whole-site PR182 B-null and new line-conditioned B-null; do not describe them as directly commensurate ecological distances. Lines have heterogeneous numbers and different sampling designs across populations. Small quadrats can produce chance mono-morphism; null corrects this at line scale.
- Publish source SHA/MD5 and no-label-manipulation original data provenance, siteyear/line diagnostic CSVs, original and line-conditional null results, valid/held counts, independent verifier and tests. Include zero M target plot ledger. Do not make claims about floral genotypes, pollinator behavior, mating success, geographic adaptation, or a novel original result already in Takahashi et al. 2016.

## Audit
Build code and synthetic tests in an isolated package **after this protocol commit**. Independently verify observed numerator, denominator, siteyear support and exact expected conditional B via `fractions.Fraction` using original TXT (not importing the primary code). Verify independent Monte Carlo reference fixed seeds and a hand-enumerated tiny case. Save source and endpoint SHA. Run tests locally, then source-recover and rerun remote CI on a **separate Draft PR stacked on PR182**. Leave FCP main H1/H2 and all preexisting PR170-182 outputs unchanged.
