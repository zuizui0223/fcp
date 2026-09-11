# RGFCA 42,111 estimand ladder

Date: 2026-09-11 JST
Status: frozen while the authorized 42,111-species breadth measurement is running and before the complete species-colour join is opened.

The expanded atlas is not one analysis with progressively smaller denominators. It is a ladder of distinct estimands whose sampling requirements differ.

| Layer | Frozen opportunity | Primary estimand | Explicit non-claim |
|---|---:|---|---|
| Species breadth | 42,111 species; one frozen anchor/species | species-equal distribution of one observed coarse flower-colour state | not modal species colour; not polymorphism prevalence |
| Repeated colour breadth | repeated random species subsets of the complete 42,111 result | minimum species N reproducing the complete one-anchor colour composition | not minimum N for polymorphism or maps |
| Cross-cell geographic breadth | 85,337 frozen taxon×cell anchors; 42,111 species; 128 cells | cell-wise species-equal observed colour composition | not pooled with the species-breadth frequency estimator |
| Observer-disjoint cross-cell pairs | subset of 13,564 species with >=2 frozen cells that pass observer recovery | probability that two frozen anchors from distinct cells and observers have different coarse colour states | not within-population coexistence; not automatically species-wide Simpson D |
| Multi-photo depth | metadata opportunity: >=2: 33,944; >=5: 24,612; >=10: 18,301; >=20: 12,985 species | future within-species diversity estimates at increasing precision after separately frozen acquisition | `after_observer_cap` alone does not guarantee observer-disjoint samples |
| Strong spatial depth | 2,636 species with >=20-photo opportunity, >=5 discovery cells and >=100 km span | future high-support within-species spatial organization | not currently measured in the 42,111 expansion |

## Weighting boundaries

- Species-breadth global composition: each discovered species receives equal weight through its single frozen anchor.
- Geographic map: within each cell, taxon×cell anchors receive equal species weight. Do not pool the 85,337 rows as a global frequency because widespread species then receive more weight.
- Cross-cell pair analysis: one frozen observer-disjoint pair per eligible species, giving equal species weight.
- Missing/unevaluable anchors stay in their original denominator and are never replaced after colour opening.

## Interpretation order

1. First report what is observable across 42,111 species.
2. Then report whether that species-equal composition saturates at a smaller N under repeated subsampling.
3. Separately map geographic composition using the taxon×cell estimator.
4. Separately quantify cross-cell within-species discordance using observer-disjoint pairs.
5. Only multi-photo depth analyses may support species-level polymorphism/diversity estimates; their sampling frame must be frozen before additional pixels are opened.

The ladder prevents the very large breadth denominator from being used to overstate what one photograph per species can identify.
