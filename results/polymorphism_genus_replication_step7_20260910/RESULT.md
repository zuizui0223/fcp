# Polymorphism genus clustering — Step 7 result

**Raw genus clustering replicated in the species-disjoint reserve: `True`.**
**Opportunity-robust genus clustering replicated: `False`.**
**Discovery genus clustering survives span+n adjustment: `True`.**

## Discovery

- raw Step-4 reproduction: repeated genera **61**, species **169**, gain **0.202963**, p_lower **0.00229989**
- adjusted for sampled span + n: gain **0.154051**, p_lower **0.00754962**
- D_unbiased adjusted: gain **0.151979**, p_lower **0.00819959**

## Species-disjoint reserve

- repeated genera **54**, species **146**
- raw D: gain **0.153654**, p_lower **0.0216989**
- adjusted for sampled span + n: gain **0.033604**, p_lower **0.304385**
- D_unbiased adjusted: gain **0.035231**, p_lower **0.297685**

## Cross-tranche repeated-genus means

- shared genera with >=2 species in both tranches: **23**
- test run: **True**
- Spearman rho = **0.530632**, p_two_sided = **0.00989951**

## Measurement/opportunity diagnostics

- discovery rho(D, n_classifiable) = **-0.537139**
- reserve rho(D, n_classifiable) = **-0.555901**
- discovery rho(D, sampled span) = **0.181031**
- reserve rho(D, sampled span) = **-0.002567**

## Claim boundary

A replicated result is genus-level taxonomic clustering of photo-derived species polymorphism diversity. It is not formal phylogenetic signal, genetic determination, adaptation, or a causal effect of genus membership.
