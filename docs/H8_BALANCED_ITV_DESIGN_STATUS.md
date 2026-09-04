# H8 balanced species-by-cell replication — design audit terminal state

## Status

H8 is **not opened as an experiment** because its metadata-only design audit found no species-by-cell design that met the frozen minimum of 35 species while retaining at least five cells per species.

No H8 iNaturalist query and no H8 image pixel was opened.

## Audit result

The H8 audit used only H7 fresh metadata support. H7 itself had stopped before pixels, so this audit did not use flower colour.

For at least five qualifying cells per species:

- ≥10 retained fresh metadata records/cell: **14 species**;
- ≥15/cell: **10 species**;
- ≥20/cell: **7 species**;
- ≥25/cell: **6 species**;
- ≥30/cell: **5 species**.

For six qualifying cells the counts were still lower. The frozen requirement was at least **35 species**.

Binding design decision:

`h8_metadata_design_feasible = false`

## Consequence

The species-by-cell replication architecture is not relaxed further. In particular, the program does not retreat to three or four cells per species, because the completed H6/H6b diagnostics showed that low-cell, low-information species were precisely where nominal spatial signal was strongest and least information-weight robust.

The next prospective design therefore changes the sampling unit rather than weakening the cell threshold: fixed-n, range-spanning individual photographs per species with continuous geographic distance and a photo-level permutation null.
