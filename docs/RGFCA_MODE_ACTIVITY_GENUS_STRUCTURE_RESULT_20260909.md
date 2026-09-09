# RGFCA recurrent-mode activity genus-structure exploration — result

Date: 2026-09-09 JST

Protocol frozen before outcome opening at `a2c6d65fd788dbe0c050ba2838a885197b9f08ff`.

## Result

Using species-equal median RMS activity across the first 20 fixed grid shifts, 122 retained species in 50 genera with >=2 species contributed to each component/mode analysis. The >=3-species sensitivity retained 50 species in 14 genera.

For the primary flower-minus-background signal D:

- M1: genus R² = **0.429**, permutation p-like = **0.349**
- M2: genus R² = **0.437**, permutation p-like = **0.296**
- M3: genus R² = **0.560**, permutation p-like = **0.0088**

The >=3-species sensitivity for D-M3 gave R² = **0.462**, p-like = **0.0209**.

Thus only M3 shows clear genus structure in the primary D analysis.

Component diagnostics indicate that this is not cleanly flower-only. In the >=3-species sensitivity, flower M3 had R² = 0.475 (p-like 0.0140) and background M3 had R² = 0.472 (p-like 0.0170). This suggests genus structure in the broader flower/background signal configuration rather than a demonstrated flower-only lineage effect.

## Reproduction

A complete second execution produced byte-identical outputs for all four retained CSV files:

- `genus_activity_summary.csv`: `63f08c36dbc52cda376646fa44e288c0925763b299893dacf3cda7aa358fcb68`
- `genus_structure_summary.csv`: `fd46562336722e1e36c209c2ef0dd31468f37a6339cd21fd92e29c2e5dde5022`
- `mode_activity_by_species_shift.csv`: `bdfff7de264b652096a84bea29220ce55af6c10718dda025f51ae4d6655ee22e`
- `mode_activity_species_summary.csv`: `6517842a7c82c66718ccff60791053fcb9ca27aaee2bb505ed9189381ce415cd`

## Interpretation boundary

The current reserve supports exploratory genus structure in the amount of geographic M3 signal activity. This is not a phylogenetic heritability estimate and does not establish shared mechanism, adaptation, or lineage-level causality. Observer-disjoint replication is the next required stress test.
