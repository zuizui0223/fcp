# Global Monte Carlo flower-colour barrier atlas — prospective protocol v1

Status: **prospective design frozen on 2026-09-04 before any global Monte Carlo flower-colour field was available**. The branch was forked from `362cdcc949f1421a9a5bb0532453914a23b4be83` while the H9 location-blind measurement workflow was still queued. The recurrent-zone extraction rule, ecological/geographic overlay family and external source families were subsequently frozen while the new global species-discovery workflow was still metadata-only. This protocol does not alter, rescue or reinterpret the frozen H1–H6 outcomes.

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

### G1 — recurrent flower-colour barrier field

After requiring adequate within-species replication, do independent plant species place stronger flower-colour discontinuities in the same broad geographic regions more often than expected under species-conditioned colour permutation?

This is deliberately weaker and more realistic than requiring the exact same grid edge to be the strongest transition in many species.

### G2 — stability and zone reproducibility

If a recurrent field exists, is its geography stable to independent species/photo resampling, taxonomic and geographic deletion, and reasonable predeclared spatial supports? If so, can persistent high-discontinuity areas be extracted reproducibly as neutral flower-colour zones (`Z01`, `Z02`, ...)?

### G3 — prevalence, heterogeneity and edge-level mechanism alignment

Regardless of whether G1 passes, how common and heterogeneous is within-species spatial colour structure? On the same fixed photo-graph edges, are stronger colour discontinuities associated with stronger independently frozen pollinator, climate, terrain or geographic-barrier contrasts?

This edge-level mechanism analysis is independent of the existence of a shared global zone and **cannot rescue a null G1**.

### G4 — interpretation of shared flower-colour zones, hierarchically gated

Only if G1 is supported and G2 stability passes may the neutral recurrent zones be compared with independently frozen pollinator-turnover, climate-turnover, terrain, marine-gap, river, mountain and established biogeographic-boundary surfaces. These overlays annotate a supported colour zone; they cannot create one.

## 3. Global sampling frame

### 3.1 Repeated metadata-only species-discovery census

The original 8,989-species pool came from only one random metadata page in each 18 × 9 equal-area cell. It is therefore retained as **round 0**, not treated as the complete global species frame.

Before any new colour pixel is opened, species discovery itself is repeated prospectively:

- the same 18 × 9 equal-area grid;
- 20 fresh metadata-only discovery rounds;
- one independently random 200-record iNaturalist page per cell per round;
- 3,240 fixed fresh request attempts in total;
- the same research-grade, flowering-annotation, georeferenced, positional-accuracy and licence filters as the original atlas;
- no colour, climate, barrier or H9 outcome enters discovery;
- no early stopping when species accumulation appears to plateau;
- failed cell-round requests are recorded and not replaced or rerun.

The census saves a compact observation-ID/species index, the cumulative species frame, per-round and per-cell accumulation audits, odd-versus-even fresh-round species-set overlap, baseline-versus-fresh overlap and the species-accumulation curve. Thus the geographic/taxonomic sampling frame itself has a reproducible stability diagnostic rather than depending on one random global draw.

All species discovered in round 0 or any of the 20 fresh rounds enter the next **metadata capacity scan**. Discovery does not imply inclusion in colour inference.

### 3.2 Capacity census and inferential layer

The inferential layer contains only species with enough metadata support for fixed-n, range-spanning measurement. Eligibility is determined without colour.

Before any fresh pixels are opened, a metadata-only capacity audit evaluates the following predeclared raw-photo targets per species:

- 100;
- 80;
- 60.

Within each candidate species, prior experiment IDs are excluded, observer contribution is capped at two retained photographs, and retained raw photographs are selected by a deterministic geographic maximin rule with deterministic hash tie-breaking.

The automatic primary raw-photo target is the **largest** of 100, 80 or 60 that yields at least 300 metadata-eligible species globally. If no target yields at least 300 species, the experiment stops as `not_evaluable_global_inferential_species_capacity_before_pixels`; the threshold is not relaxed after colour opening.

The capacity scan and the eventual candidate-image acquisition are distinct outcome-blind stages. Capacity may be estimated from a metadata-only draw, but actual image acquisition uses one separately frozen fresh draw at the automatically selected target and has its own premeasurement gate. A favourable second target may not be substituted if the acquisition draw underperforms.

The metadata audit must also report representation by major terrestrial realm, latitude band and family. These are coverage diagnostics, not outcome-driven inclusion criteria and are never used to select a favourable colour result.

Species that cannot support within-species inference remain part of the descriptive global census. They are not recoded as biologically unstructured.

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
- kernel bandwidth: 500 km;
- evaluated support extends to 3 bandwidths (1,500 km).

Predeclared sensitivities:

- 24 × 12 with 1,000-km bandwidth;
- 72 × 36 with 250-km bandwidth.

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
- graph geometry and global opportunity fields;
- all external pollinator, climate, terrain and geographic predictor values.

Only classifiable colour vectors are permuted within species.

Primary null size: **999 permutations** with a pre-frozen seed. Each permutation is run through the exact same 200-resample schedule as the observed data. The same null therefore supplies both the G1 field test and the downstream edge/overlay alignment tests without introducing a different favourable spatial null.

G1 support requires:

- postmeasurement inferential gate passed;
- observed G1 concentration greater than the null expectation;
- upper-tail Monte Carlo `p < 0.05`.

No spatial support, kernel, species weight, photo threshold or subset may be substituted after seeing G1.

## 9. Stability, persistent zone extraction and robustness

A positive G1 is not sufficient by itself. The following are reported regardless of significance where estimable:

1. **species-discovery stability:** accumulation across the fixed 20 metadata rounds and odd/even fresh-round overlap;
2. **Monte Carlo convergence:** running consensus-field distance and concentration versus 25, 50, 100, 150 and 200 resamples;
3. **split-resample reproducibility:** correlation of odd- versus even-resample consensus fields;
4. **photo-resampling stability:** within-species spread of `rho_i` and boundary contribution across the 200 draws;
5. **leave-one-realm-out:** recompute the consensus after deleting each represented terrestrial realm;
6. **leave-one-family-out for major families:** delete families contributing at least 5% of inferential species;
7. **equal-species versus inclusion-probability weighting:** weighting sensitivity only; equal-species remains primary;
8. **predeclared spatial-support sensitivity:** 24×12/1,000 km and 72×36/250 km cannot replace the primary 36×18/500 km result.

If G1 and the G2 stability gate are supported, persistent flower-colour zones are extracted without consulting named geography or external ecological surfaces:

- within each of the 200 observed resamples, mark the top **10%** of evaluable field cells;
- a primary zone seed must be in that per-resample top decile in at least **60%** of evaluable observed resamples;
- require at least **100** evaluable resamples for a seed cell;
- join seed cells by 8-neighbour connectivity with longitude wrap across the international date line;
- discard components smaller than **3 cells**;
- assign neutral IDs `Z01`, `Z02`, ... in descending opportunity-weighted integrated colour-field intensity;
- inspect geographic names or candidate explanations only after these neutral components are frozen.

Top-5% and top-15% hotspot definitions are sensitivity analyses only and cannot replace the primary top-10% rule.

## 10. External ecological and geographic concordance

The external source family is fixed before the global colour outcome: CHELSA v2.1 climate, EarthEnv topography, GMBA Mountain Inventory v2, RESOLVE Ecoregions 2017, HydroRIVERS/HydroSHEDS, Natural Earth land-water geometry, a stable versioned GloBI interaction dataset, and an independently frozen GBIF occurrence download. Exact external download identifiers/checksums must be recorded before they enter colour inference. Missing external coverage is `not_evaluable`, never biological zero.

### 10.1 G3 edge-level mechanism alignment — does not require a shared zone

For every externally evaluable predictor, calculate an external contrast on the same fixed within-species photo-graph edges used for colour discontinuity. Each species contributes one within-species Spearman correlation between colour-edge score and external-edge score. The primary statistic is the **equal-species mean** of those correlations.

Primary requirements:

- at least 5 evaluable graph edges per species;
- at least 30 evaluable species for a predictor;
- the same 999 within-species colour permutations as the main analysis;
- one-sided upper-tail test;
- Holm correction across the fixed primary predictor family.

The fixed primary edge family is:

1. direct pollinator-interaction community turnover from GloBI where coverage permits;
2. GBIF-based pollinator-guild occurrence turnover (bees, Syrphidae, Lepidoptera proxy, and a flower-visiting bird proxy);
3. multivariate CHELSA climate turnover using standardized BIO1, BIO4, BIO12 and BIO15;
4. EarthEnv terrain-barrier intensity from elevation, slope, terrain roughness index and vector ruggedness at the frozen 50-km grain;
5. marine-gap fraction along the great-circle edge;
6. major-river crossing intensity;
7. crossing intensity of RESOLVE realm/biome/ecoregion boundaries;
8. GMBA mountain-system boundary crossing intensity.

A positive edge result means that species-specific flower-colour geography is concordant with that external contrast more strongly than expected after species-conditioned colour permutation. It **does not** establish a shared global flower-colour zone and cannot rescue a null G1.

### 10.2 G4 shared-zone overlay — requires G1 + G2

For a supported recurrent colour field, each independently frozen external surface is compared with the colour field using opportunity-weighted Spearman correlation. The exact same 999 species-conditioned null colour fields define the spatial null. At least 50 jointly evaluable cells are required, and Holm correction is applied across the fixed primary overlay family.

Only globally supported predictors may annotate a neutral `Z` component. A zone receives a descriptor such as `climate-associated`, `pollinator-associated`, `terrain-associated`, `marine-associated`, `river-associated`, `biogeographic-boundary-associated`, or `mountain-boundary-associated` only when:

1. that predictor passes the global overlay test; and
2. the opportunity-weighted predictor mean inside the neutral zone is at or above the 75th percentile of that predictor over all evaluable cells.

Multiple descriptors are allowed. Causal labels such as “pollinator-driven” or “mountain-caused” are forbidden from these spatial concordance tests alone.

Direct GloBI interaction turnover and broad GBIF guild turnover are retained as complementary pollinator layers. Agreement strengthens ecological coherence; disagreement is reported as interaction-specificity/sampling uncertainty rather than selecting whichever layer fits better.

## 11. Computational interpretation

With 250 species × 20 photographs, the complete pairwise G3 calculation contains 47,500 photograph pairs per outer resample. Across 200 observed resamples this is 9.5 million pair evaluations, before null reuse/vectorization. Graph-based G1 is smaller still.

Thus the scientific sample can eventually contain far more than 250 species and far more than 20 photographs per species while each compute unit remains bounded. Additional species/photos improve coverage and resampling stability rather than exploding one monolithic analysis.

This is the intended distinction:

- **more metadata discovery rounds** broaden and stabilize which species have an opportunity to enter the global programme;
- **more Monte Carlo analysis resamples** reduce computational/sampling Monte Carlo error conditional on the measured pool;
- **more adequately replicated species and photographs** add biological information.

Repeated resampling cannot manufacture within-species information from singleton species, which is why the new metadata and postmeasurement gates are binding.

## 12. Claim boundaries and outcome logic

A supported G1+G2 would establish recurrent broad geographic concentration of flower-colour discontinuity across globally sampled, adequately replicated species under the frozen opportunity-conditioned design. Stable neutral components could then reasonably be described as **recurrent flower-colour biogeographic zones**.

The strongest possible result would be a stable neutral zone whose colour discontinuity exceeds the species-conditioned null and whose geography is independently concordant with one or more frozen ecological/geographic surfaces. It still would not by itself establish causal adaptation, pollinator-mediated selection or historical vicariance.

The interpretation matrix is fixed:

- **G1+G2 positive; ecological/geographic overlay positive:** recurrent flower-colour zone with independent ecological/geographic concordance;
- **G1+G2 positive; overlays null:** recurrent colour zone detected, mechanism/geographic interpretation unresolved;
- **G1 null; G3 edge alignment positive:** species-specific colour geography repeatedly covaries with external contrasts but is not globally synchronized;
- **G1 null; G3 heterogeneous but edge alignment null:** substantial species-specific geography without a common global mechanism at the tested scales;
- **G1 and G3 essentially null:** no general flower-colour geographic organization detected under the frozen design and scales.

No result by itself establishes:

- one exact universal boundary line;
- climate causation;
- topographic causation;
- local adaptation;
- pollinator mediation;
- pigment physiology;
- absence of species-specific structure outside recurrent regions.

## 13. Relationship to prior experiments

- Six-species Chapter 1 remains evidence that within-species spatial colour organization can occur under dense focal sampling.
- Random-atlas H1 remains the valid result for its original globally broad/sparse frame and is not rerun with relaxed rules.
- H6/H6b motivated the need to distinguish low-information equal weighting from reliable within-species replication.
- H7/H8 established, before pixels, that balanced species-by-cell sampling was not feasible from that fresh metadata frame.
- H9 changed the unit to fixed-n individual photographs with continuous distance. The present protocol generalizes that logic to a globally broader Monte Carlo coverage design while preserving outcome blindness at each new gate.
