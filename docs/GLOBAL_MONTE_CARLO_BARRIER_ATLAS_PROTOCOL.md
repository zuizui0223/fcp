# Global Monte Carlo flower-colour barrier atlas — prospective protocol v1

Status: **prospective design frozen on 2026-09-04 before any biological H9 result is available**. The branch was forked from `362cdcc949f1421a9a5bb0532453914a23b4be83` while the H9 location-blind measurement workflow was still queued. This protocol does not alter, rescue or reinterpret the frozen H1–H6 outcomes, and it does not use a future H9 colour result to choose its design.

## 1. Why a new design is necessary

The first random photo-first atlas was globally broad but extremely sparse within species.

Its measured frame contained:

- 20,845 measured photographs;
- 10,103 classifiable colour photographs;
- 8,989 candidate species;
- 5,144 species with at least one classifiable photograph;
- only 1,480 species with at least two classifiable photographs;
- 752 with at least three;
- 356 with at least five;
- 111 with at least ten;
- 35 with at least twenty.

Thus the total candidate average was only 2.32 photographs per species and the classifiable average was 1.12 per candidate species. Among the 5,144 species with any classifiable colour, 3,664 (71.2%) were singletons.

This matters because the frozen H1 null permutes classifiable colour labels only within species. A species represented by one classifiable photograph cannot change under that null, and a species represented by two or three photographs supplies only weak within-species geographic information. Therefore the H1 result (`p = 0.413`) is retained exactly as a failure to detect excess synchronized boundary concentration **in that globally broad but within-species sparse sampling frame**. It is not promoted into evidence that terrestrial geographic barriers do not organize flower colour.

The new experiment separates two goals that the first atlas conflated:

1. **global coverage:** many lineages and regions must have non-zero inclusion probability;
2. **within-species information:** each inferential species must contribute enough spatially dispersed photographs to estimate colour geography.

The computational solution is not to load every photograph of every species into one analysis. It is a fixed, reproducible **Monte Carlo coverage design** in which globally eligible species and their photographs are repeatedly subsampled under balanced inclusion rules.

## 2. Questions

### G1 — recurrent barrier field

After requiring adequate within-species replication, do independent plant species place stronger flower-colour discontinuities in the same broad geographic regions more often than expected under species-conditioned colour permutation?

This is deliberately weaker and more realistic than requiring the exact same grid edge to be the strongest transition in many species.

### G2 — stability of the inferred field

If a recurrent field exists, is its geography stable to independent species/photo resampling, taxonomic and geographic deletion, and reasonable predeclared spatial supports?

### G3 — prevalence and heterogeneity of species-specific spatial organization

Regardless of whether G1 passes, how common and how heterogeneous is within-species spatial colour structure across the globally sampled inferential species set?

G3 prevents a null G1 from being misread as a null statement about all species-specific geography. The primary outputs are the mean, median, positive fraction and between-species heterogeneity of species-level spatial effects.

### G4 — barrier interpretation, hierarchically gated

Only if G1 is supported may an independent, pre-frozen barrier-overlay extension ask whether the recurrent field aligns with broad biogeographic, topographic or environmental discontinuities. Barrier overlays cannot rescue a null G1.

## 3. Two-layer global sampling frame

### 3.1 Census layer

The census layer is metadata-only and outcome-blind. It is intended to be as globally inclusive as practical and can be processed in hash-partitioned chunks. No flower-colour pixel is opened at this stage.

The census records species identity, coordinates, observation/photo identifiers, observer identity, licence and metadata quality. It is allowed to estimate sampling capacity and geographic geometry but not flower colour.

Species that cannot support within-species inference remain part of the descriptive global census. They are not recoded as biologically unstructured.

### 3.2 Inferential layer

The inferential layer contains only species with enough metadata support for fixed-n, range-spanning measurement. Eligibility is determined without colour.

Before any fresh pixels are opened, a metadata-only feasibility audit evaluates the following predeclared raw-photo targets per species:

- 60;
- 80;
- 100.

Within each candidate species, prior experiment IDs are excluded, observer contribution is capped at two retained photographs, and retained raw photographs are selected by a deterministic geographic maximin rule with deterministic hash tie-breaking.

The automatic primary raw-photo target is the **largest** of 100, 80 or 60 that yields at least 300 metadata-eligible species globally. If no target yields at least 300 species, the experiment stops as `not_evaluable_global_inferential_species_capacity_before_pixels`; the threshold is not relaxed after colour opening.

The metadata audit must also report representation by major terrestrial realm, latitude band and family. These are coverage diagnostics, not outcome-driven inclusion criteria and are never used to select a favourable colour result.

## 4. Location-blind colour measurement

Fresh candidate images are measured without species name, coordinates, climate or barrier context using the already-qualified ROI-v4 flower-mask path plus the fixed generic CIELAB palette.

No failed image is replaced after measurement.

A species enters Monte Carlo inference only if it retains at least **40 classifiable photographs**. This number is fixed because each analysis draw uses 20 photographs and the design needs genuine within-species photo resampling rather than repeatedly analyzing the same fixed set.

The experiment requires at least **250 inferential species** after measurement. If fewer remain, G1–G3 are `not_evaluable`; no classifiable-photo threshold is relaxed.

## 5. Monte Carlo coverage schedule

The compute budget is bounded independently of the total census size.

Primary schedule:

- outer observed resamples: **200**;
- species per resample: **250**;
- photographs per included species: **20**;
- species inclusion: seeded balanced random blocks so eligible species differ in inclusion count by at most one whenever mathematically possible;
- photo inclusion: seeded balanced random blocks within species so all classifiable photographs have repeated non-zero opportunity while no single resample uses more than 20 photographs per species;
- the complete species and photo inclusion-count audit is saved.

If exactly 250 species survive, every species is included in every outer resample but its 20-photo subset still varies. If more than 250 survive, the balanced random schedule spreads computation across the full eligible pool.

This design uses the full eligible dataset **in probability and over repeated coverage**, rather than attempting an all-at-once Cartesian analysis.

## 6. Species-level geometry

For every included species in every observed resample:

1. construct a colour-blind spherical `k = 3` nearest-neighbour graph on the 20 coordinates;
2. compute soft four-group Jensen–Shannon colour divergence on each undirected graph edge;
3. rank-standardize edge divergence within species so species with intrinsically larger colour variance do not dominate the global map;
4. preserve the complete edge geometry as the opportunity denominator.

A companion continuous-distance statistic is computed for G3:

`rho_i = Spearman(pairwise great-circle distance, pairwise colour JSD)`

using all 190 photograph pairs from the fixed 20-photo draw.

G3 reports the Monte Carlo distribution of each species' `rho_i`, the equal-species global mean and median, the fraction of species with positive median `rho_i`, and a hierarchical between-species variance estimate. A near-zero global mean with large heterogeneity is interpreted as heterogeneous species-specific organization, not universal absence.

## 7. Global recurrent-barrier field

Each species graph edge contributes to a global equal-area field at its great-circle midpoint, with spatial kernel smoothing.

Primary support:

- global grid: 36 × 18 equal-area longitude × sin(latitude);
- kernel support: 500 km.

Predeclared sensitivities:

- 24 × 12 with 1,000-km kernel;
- 72 × 36 with 250-km kernel.

For every field cell:

- **opportunity** = species-normalized geometric edge support near the cell;
- **colour numerator** = the same support weighted by within-species rank-standardized colour divergence;
- **field intensity** = numerator divided by opportunity where opportunity passes the frozen minimum support.

Every species is normalized to equal total contribution within a resample. Abundant species cannot dominate merely by having more raw photographs.

The primary G1 statistic is opportunity-weighted spatial concentration of the consensus field. A recurrent barrier is therefore a region with stronger-than-null colour discontinuity **conditional on the fact that sampled species had geographic edge opportunity there**.

## 8. Species-conditioned null

The null keeps fixed:

- species identities;
- coordinates;
- observer-filtered photograph pools;
- all Monte Carlo species/photo sampling schedules;
- structural measurement missingness;
- graph geometry and global opportunity fields.

Only classifiable colour vectors are permuted within species.

Primary null size: **999 permutations** with a pre-frozen seed. Each permutation is run through the exact same 200-resample schedule as the observed data.

Support requires:

- postmeasurement inferential gate passed;
- observed G1 concentration greater than the null expectation;
- upper-tail Monte Carlo `p < 0.05`.

No spatial support, kernel, species weight, photo threshold or subset may be substituted after seeing G1.

## 9. Stability and robustness are first-class results

A positive G1 is not sufficient by itself. The following are reported regardless of significance where estimable:

1. **Monte Carlo convergence:** running consensus-field distance and concentration versus 25, 50, 100, 150 and 200 resamples;
2. **split-resample reproducibility:** correlation of odd- versus even-resample consensus fields;
3. **photo-resampling stability:** within-species spread of `rho_i` and boundary contribution across the 200 draws;
4. **leave-one-realm-out:** recompute the consensus after deleting each represented terrestrial realm;
5. **leave-one-family-out for major families:** delete families contributing at least 5% of inferential species;
6. **equal-species versus inclusion-probability weighting:** weighting sensitivity only; equal-species remains primary;
7. **predeclared spatial-support sensitivity:** 24×12/1,000 km and 72×36/250 km cannot replace the primary 36×18/500 km result.

A stability claim requires that the sign of the G1 excess is unchanged, the odd/even field correlation is positive and substantial, and no single realm or major family creates the result. Exact numerical stability thresholds, if used as pass/fail gates, must be frozen in a separate synthetic-validation amendment before fresh colour pixels open.

## 10. Computational interpretation

With 250 species × 20 photographs, the complete pairwise G3 calculation contains 47,500 photograph pairs per outer resample. Across 200 observed resamples this is 9.5 million pair evaluations, before null reuse/vectorization. Graph-based G1 is smaller still.

Thus the scientific sample can eventually contain far more than 250 species and far more than 20 photographs per species while each compute unit remains bounded. Additional species/photos improve coverage and resampling stability rather than exploding one monolithic analysis.

This is the intended distinction:

- **more Monte Carlo resamples** reduce computational/sampling Monte Carlo error conditional on the measured pool;
- **more adequately replicated species and photographs** add biological information.

Repeated resampling cannot manufacture within-species information from singleton species, which is why the new metadata and postmeasurement gates are binding.

## 11. Claim boundaries

A supported G1 would establish recurrent broad geographic concentration of flower-colour discontinuity across globally sampled, adequately replicated species under the frozen opportunity-conditioned design.

It would not by itself establish:

- one exact universal boundary line;
- climate causation;
- topographic causation;
- local adaptation;
- pollinator mediation;
- pigment physiology;
- absence of species-specific structure outside recurrent regions.

A null G1 would mean no recurrent field was detected at the frozen global support among the inferential species. G3 remains independently interpretable and can show substantial species-specific heterogeneity even when G1 is null.

## 12. Relationship to prior experiments

- Six-species Chapter 1 remains evidence that within-species spatial colour organization can occur under dense focal sampling.
- Random-atlas H1 remains the valid result for its original globally broad/sparse frame and is not rerun with relaxed rules.
- H6/H6b motivated the need to distinguish low-information equal weighting from reliable within-species replication.
- H7/H8 established, before pixels, that balanced species-by-cell sampling was not feasible from that fresh metadata frame.
- H9 changed the unit to fixed-n individual photographs with continuous distance. The present protocol generalizes that logic to a globally broader Monte Carlo coverage design and was frozen before H9 supplied a biological outcome.
