# JBI Chapter 1 — continuous spatial colour-organization protocol

## Central question

The spatial arm does not begin by asking which environment predicts each flower colour. It separates two ordered questions:

1. **Stage A — local organization:** after conditioning on each species' observed distribution, are geographically neighbouring flower-colour observations more similar than expected under random labelling?
2. **Stage B — shared concentration:** if Stage A is supported, do independent species concentrate their relatively strongest continuous colour transitions in the same globally evaluable cells more than expected under the complete species-conditioned null pipeline?

The global display may suppress species identity, but species identity is never removed from the inferential null model.

## Frozen photograph design

The development sample contains six species:

- *Antirrhinum majus*;
- *Dactylorhiza sambucina*;
- *Gentiana lutea*;
- *Ipomoea purpurea*;
- *Lysimachia arvensis*;
- *Raphanus sativus*.

For each species, 200 photographs were acquired and assigned deterministically to:

- 80 calibration photographs;
- 120 held-out evaluation photographs.

Total sample:

- calibration: 480;
- evaluation: 720;
- total: 1,200.

The split is outcome-blind, hash-frozen and joined to coordinates by `photo_id`.

## Frozen continuous colour representation

Calibration geometry did not justify forcing one universal discrete-morph model across all six species. The primary representation was therefore frozen before evaluation as a **species-specific continuous colour vector**.

For every photograph:

1. Florence localizes the operational flower region using the frozen prompt/model path;
2. the frozen species-specific feature subset is extracted;
3. every feature is standardized with calibration-only means and population standard deviations;
4. the complete standardized vector is retained as one indivisible observational label.

The representation is neither a biological morph label nor an inferred causal mechanism. It is an operational continuous description of visible flower-colour structure.

No evaluation-derived feature selection, scaling, ROI rule or species-specific dimension change is allowed.

## Evaluation opening rule

After the representation and analysis contracts were frozen, the same extraction path was applied unchanged to all 720 evaluation photographs.

A post-opening rule change prompted by the evaluation values would invalidate this confirmatory evaluation cycle. No such rule change was made.

## Stage A — species-conditioned local organization

### Colour-blind graph

Within each species, construct a symmetrized spherical k-nearest-neighbour graph using only observation latitude and longitude.

Primary graph:

- `k = 5`.

Predeclared sensitivities:

- `k = 3`;
- `k = 8`.

Colour values do not choose neighbours, remove graph edges or set the graph scale.

### Edge discontinuity

For standardized colour vectors `z_i` and `z_j` connected by edge `(i,j)`, define continuous colour discontinuity as the root-mean-square component difference:

`d_ij = sqrt(mean((z_i - z_j)^2))`.

For each species, `Q_i` is the arithmetic mean of its edge discontinuities. The global statistic is the arithmetic mean of the six species-specific `Q_i` values, giving every species one vote regardless of its retained edge count.

### Stage-A null

For every permutation replicate:

1. keep observation locations and the graph fixed;
2. permute complete colour vectors strictly within species;
3. recompute every edge discontinuity;
4. recompute every `Q_i`;
5. recompute the equal-species global mean.

Primary and sensitivity analyses each use 9,999 permutations with the frozen random seed.

The one-sided confirmatory alternative is lower observed discontinuity than the null: neighbouring observations are more similar than expected.

### Stage-A gate

Stage B is run only if the frozen primary Stage-A lower-tail Monte Carlo p-value is at most 0.05.

That gate was passed (`p = 0.0113`).

## Stage B — shared continuous transition concentration

Stage B does not convert the continuous vectors into post-hoc discrete colour states. It derives a relative transition-intensity surface for each species while preserving the full Stage-A random-labelling logic.

### Base graph

Use the frozen Stage-A primary graph (`k = 5`) for every species.

### Geometry-only primary support selection

Before any observed colour edge score is computed, evaluate a fixed ordered set of spatial supports using only:

- species identity;
- observation coordinates;
- base-graph edges;
- great-circle edge distances.

Candidate edge caps, in priority order:

1. 500 km;
2. 1,000 km;
3. 2,000 km.

Candidate equal-area longitude–sin(latitude) grids, in priority order within each edge cap:

1. 36×18;
2. 24×12;
3. 18×9.

For each species/configuration:

1. retain base-graph edges no longer than the candidate cap;
2. assign each retained edge to the fixed equal-area cell containing its great-circle midpoint;
3. define a species/cell as detectable only when at least two retained geometry edges fall in that cell.

Let `D_i(x)=1` when species `i` is detectable in cell `x`, and zero otherwise. Then:

`A(x) = sum_i D_i(x)`

is the opportunity denominator: the number of species for which a transition could have been measured in cell `x` under the frozen geometry rule.

A configuration passes only if all frozen support criteria hold:

- at least 30 retained edges for every species;
- at least 8 cells with `A(x) ≥ 2`;
- at least 2 cells with `A(x) ≥ 3`;
- at least 4 species contribute to at least one cell with shared opportunity.

The first passing configuration is the primary. If no configuration passes, Stage B is reported as not estimable without inspecting colour-derived concentration statistics.

All nine configurations passed. The frozen primary is therefore:

**500-km edge cap + 36×18 equal-area grid**.

### Species transition intensity

For one species and one fixed configuration:

1. compute RMS continuous colour discontinuity on retained edges;
2. average-rank transform the edge scores within species to `[0,1]`;
3. average those ranked edge intensities within each detectable cell.

This produces species-cell transition intensity `T_i(x)`.

Ranking makes Stage B compare the geographic placement of relatively strong versus weak transitions within each species rather than treating raw vector distances from species-specific feature spaces as directly commensurate.

### Shared transition surface

For cells with at least two detectable species:

`S(x) = sum_i D_i(x) T_i(x) / A(x)`.

Cells with `A(x) < 2` are `not evaluable` and stored as `NaN`, never as zero.

### Confirmatory concentration statistic

The primary statistic is the opportunity-weighted variance of `S(x)` across cells with `A(x) ≥ 2`, using `A(x)` as the cell weight.

Larger values indicate stronger geographic concentration of shared transition intensity.

Cell-wise uncorrected p-values are not the primary inference.

### Complete Stage-B null

For each of 9,999 primary replicates:

1. keep locations, base graphs, local edge filters, grids and detectability masks fixed;
2. permute complete colour vectors strictly within species;
3. recompute edge RMS discontinuities;
4. recompute within-species ranks;
5. recompute species-cell intensities;
6. recompute `S(x)`;
7. recompute opportunity-weighted concentration.

The confirmatory p-value is the plus-one Monte Carlo upper-tail probability.

Every non-primary configuration that independently passed the same geometry-only support criteria was evaluated as a predeclared descriptive sensitivity using 1,999 permutations. Sensitivity configurations do not replace the selected primary after results are observed.

## Ordered interpretation

### If Stage A is not supported

Stop. Do not estimate shared boundaries or search for geographic causes.

### If Stage A is supported but Stage B is not

Conclude that flower colour shows repeated within-species spatial organization, but that the present six-species sample does not establish one common global transition geography at the frozen supports.

This is the realized outcome:

- Stage A primary: supported (`p = 0.0113`);
- Stage B primary: not supported (`p = 0.0906`).

### If both Stage A and Stage B are supported

Only then freeze a geographic reference library and test post-discovery correspondence without modifying the discovered transition surface.

Because the primary Stage-B gate was not passed, geographic correspondence is not part of the confirmatory main line. Any later environmental or historical overlay must be labelled exploratory and cannot rescue the Stage-B hypothesis.

## What the completed analysis establishes

It can establish:

- non-random local organization of operational continuous flower-colour vectors within species;
- whether relatively strong transitions are globally concentrated across independent species under frozen supports.

It cannot by itself establish:

- climatic causation;
- historical causation;
- identical mechanisms among species;
- morph-specific adaptation or fitness;
- physiological colour measurements from uncalibrated community photographs.

## Governing rule

> Species may disappear from the map display, but species must never disappear from the null model.

An unconditional shuffle of all photographs would test global species composition and range geography, not the spatial organization of flower colour.
