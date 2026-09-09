# RGFCA M3 genus-structure observer-disjoint stress test — result

Date: 2026-09-09 JST

Protocol frozen before outcome opening at `49332d087811713a1f33afa290146db392e281c1`.

## Fixed target

Only the previously identified flower-minus-background D, mode M3 genus-activity structure was tested. Exactly 20 previously frozen observer SHA256 partitions were used; no mode, salt, support threshold or taxonomic subset was selected after outcome opening.

## Result

- observer partitions tested: **20/20**
- salts where both observer-disjoint halves were evaluable and positive: **3/20**
- pass fraction: **0.15**
- prespecified requirement: **>=16/20**
- `observer_disjoint_genus_structure_robust`: **false**

Across halves, median observed genus R² remained moderate (A **0.499**, B **0.502**), but permutation calibration was unstable: median p-like A **0.178**, B **0.154**. Different observer partitions often concentrated the apparent genus structure in only one half.

## Reproduction

A complete second execution produced byte-identical outputs for all three retained CSV files:

- `m3_genus_observer_by_salt.csv`: `8be2fe0a599a3cf66bc004fc0c8365a88b867caea8ec68a0a8bf4ac9eed1d1aa`
- `m3_genus_observer_global_summary.csv`: `6f2883d21388c62f00b40223b049c58470052ba8f710f83637a13cc3691bfe45`
- `m3_genus_observer_half_results.csv`: `e931d1979cd11597a0d4c769f2923d4d1e7756cd3ac2fc55f63b622ac5e914a1`

## Interpretation boundary

The exploratory pooled-reserve D-M3 genus structure is not robust to complete observer-composition separation. It should therefore not be promoted as a lineage-level biological pattern from the current reserve. The strongest retained result in this exploratory lane remains the recurrence of the global M1–M3 signal axes themselves under observer-disjoint partitions; assignment of those axes to specific geography or taxonomic groups is substantially less stable.
