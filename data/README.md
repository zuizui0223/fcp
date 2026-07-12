# Empirical seed dataset

`seed_outcomes.csv` is the first literature-verified empirical layer for the FCP framework. It is intentionally small and conservative.

## Outcome classes

- `monomorphic`: only one floral colour morph is documented at the target biological scale.
- `geographic_mosaic`: multiple colour morphs occur across the species, but sampled local populations are effectively single-morph or strongly segregated.
- `mixed`: both among-population differentiation and at least some within-population coexistence are documented.
- `local_coexistence`: multiple colour morphs coexist within local populations; strong geographic segregation is not required for this code.

These classes describe **spatial organization of observed colour variation**, not the evolutionary mechanism maintaining it.

## Coding principles

1. `protected_polymorphism` from the mathematical model is not automatically equivalent to `local_coexistence` in field data.
2. Range size is not coded as a balancing mechanism. It will be analysed as an opportunity axis that can affect environmental heterogeneity, mutation/origin opportunity, persistence, and detection effort.
3. A species is not coded `monomorphic` merely because one source reports one colour. Negative states require adequate sampling evidence.
4. The seed table contains only positive polymorphism cases whose spatial organization can be inferred from primary studies.
5. `not_coded` means the source used for the current row was not sufficient to make that claim; it does not mean absence.

## Next empirical expansion

The next table should add literature-matched controls and explicit covariates:

- range size / occupied area;
- number of sampled populations and individuals;
- environmental heterogeneity across the occupied range;
- connectivity or dispersal proxies;
- life-history buffering / seed-bank traits;
- pollination system;
- publication and observation effort.

Any control species should be matched phylogenetically and by study effort before it is treated as a meaningful negative comparison.
