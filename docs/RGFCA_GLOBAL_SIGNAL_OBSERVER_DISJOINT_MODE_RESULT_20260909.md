# RGFCA global recurrent signal modes — observer-disjoint validation result

Date: 2026-09-09 JST
Protocol freeze: `bfc5f07e02b0ee658568df799a041aa794e2ef81`

## Result

The global flower-minus-background mode basis survives complete observer partitioning even though the sole local geographic candidate does not.

| mode | salts passing both halves | salt-pass fraction | median half cosine | min–max half median cosine | observer-disjoint robust |
|---|---:|---:|---:|---:|---|
| M1 | 20/20 | 1.00 | 0.9882 | 0.9604–0.9934 | yes |
| M2 | 20/20 | 1.00 | 0.9825 | 0.9448–0.9912 | yes |
| M3 | 20/20 | 1.00 | 0.9732 | 0.9332–0.9875 | yes |

Each half-data covariance was rebuilt from that observer half alone across the first 20 frozen grid shifts, using >=3 photos and >=2 observers per half-state, >=500 km within-species state separation, equal-observer state means, and species-equal covariance. The top five half-data axes were optimally aligned to the frozen full-data reference; M1–M3 were the validation targets.

Half-data support remained substantial: median retained species per salt-half was about 188 (range of salt-half minima 164–195), with median around 706 admitted states.

## Reproduction

A complete second execution was byte-identical for all retained outputs:

- `mode_observer_disjoint_summary.csv`: `f048bc65417a29edd83ce08dd0a01d35755561db37cdb4e36d228dc325766773`
- `mode_observer_disjoint_by_half.csv`: `35a985648d466b7b075465b21cede697a960b0934deddc910cc6964ae829908e`
- `mode_observer_disjoint_realizations.csv`: `2badc1482dc500c7dba49b88f8c95042dfd6ecedf58f04c886dbdd040b5ba7b8`

## Interpretation

This separates two levels of the exploratory atlas:

1. **Mode structure is robust:** the main recurrent colour-signal axes M1–M3 are recovered from disjoint observer sets.
2. **Local geography is not yet robust:** the strongest flower-side local candidate failed the 20-way observer-partition stress test (2/20 salts passed).

The current reserve therefore supports treating M1–M3 as stable descriptive axes of within-species flower/background signal variation, but it does not yet support interpreting any particular geographic hotspot or transition as biological. The next target is an observer-consensus geography: retain only cells whose mode direction recurs across disjoint observer partitions rather than using the pooled map alone.
