# Global Monte Carlo flower-colour barrier atlas — prospective protocol v1

Status: **prospective design frozen on 2026-09-04 before any global Monte Carlo flower-colour field was available**. The branch was forked from `362cdcc949f1421a9a5bb0532453914a23b4be83` while the H9 location-blind measurement workflow was still queued. The recurrent-zone extraction rule, ecological/geographic overlay family and external source families were subsequently frozen while the new global species-discovery workflow was still metadata-only. On 2026-09-05, after the metadata-only capacity scan showed that thousands of species can satisfy the largest raw-photo target but before fresh candidate acquisition or any new global flower-colour pixels were opened, the compute envelope was prospectively bounded to at most 1,000 fresh-candidate species and at most 500 pixel-measured species. This preserves the repeated-world-map estimand rather than allowing metadata abundance to turn the analysis into an all-at-once mega-analysis. This protocol does not alter, rescue or reinterpret the frozen H1–H6 outcomes.

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

### 3.2 Capacity census and bounded inferential layer

The inferential layer contains only species with enough metadata support for fixed-n, range-spanning measurement. Eligibility is determined without colour.

Before any fresh pixels are opened, a metadata-only capacity audit evaluates the following predeclared raw-photo targets per species:

- 100;
- 80;
- 60.

Within each candidate species, prior experiment IDs are excluded and observer contribution is capped at two retained photographs. The automatic primary raw-photo target is the **largest** of 100, 80 or 60 that yields at least 300 metadata-eligible species globally. If no target yields at least 300 species, the experiment stops as `not_evaluable_global_inferential_species_capacity_before_pixels`; the threshold is not relaxed after colour opening.

The capacity scan and the eventual candidate-image acquisition are distinct outcome-blind stages. Capacity may be estimated from a metadata-only draw, but actual image acquisition uses one separately frozen fresh draw at the automatically selected target and has its own premeasurement gate. A favourable second target may not be substituted if the acquisition draw underperforms.

The metadata-only capacity census eventually showed that thousands of species can satisfy the largest target. That information changed only the **computational envelope**, not the biological thresholds or colour model. Before fresh candidate acquisition and before any global candidate pixels were opened, the following bounded hierarchy was frozen:

1. **global capacity universe:** all discovered species are retained for metadata provenance and target selection;
2. **fresh candidate species:** if more than 1,000 species pass the chosen capacity target, select exactly 1,000 by deterministic SHA256 rank of `inat_taxon_id` under seed `20260916`; otherwise query all passing species;
3. **pixel-measured species:** after the fixed fresh acquisition and its ≥300-species gate, if more than 500 full-target species remain, select exactly 500 by an independent deterministic SHA256 rank under seed `20260918`; otherwise measure all full-target species;
4. **one world-map realization:** draw 250 measured species and 20 classifiable photographs per included species under balanced inclusion schedules.

No species is replaced if the fresh candidate draw or measurement underperforms. No extra species is added after seeing flower colour. The 1,000- and 500-species bounds use taxon identity only and cannot read colour, climate, pollinator, barrier, measurement-success or G1/G3 outcomes.

This means the global census defines the universe and non-zero inclusion frame, while bounded prospectively random stages define the computational sample. The biological sample size is not multiplied by the number of Monte Carlo repetitions.

The metadata audit must also report global/equal-area coverage diagnostics. These are coverage diagnostics, not outcome-driven inclusion criteria and are never used to select a favourable colour result.

Species that cannot support within-species inference remain part of the descriptive global census. They are not recoded as biologically unstructured.

## 4. Location-blind colour measurement

Fresh candidate images are measured without species name, coordinates, climate or barrier context using the already-qualified ROI-v4 flower-mask path plus the fixed generic CIELAB palette.

At most **500 species** enter pixel measurement. At the largest 100-photo target this caps the location-blind image workload at **50,000 photographs**. Every frozen candidate photograph of a selected measurement species is attempted; failed image downloads or measurements are terminal missingness and are not replaced.

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

If exactly 250 species survive, every species is included in every outer resample but its 20-photo subset still varies. If more than 250 survive, the balanced random schedule spreads computation across the full measured eligible pool. With the 500-species measurement ceiling, a fully populated frame gives each species essentially symmetric long-run inclusion opportunity across the 200 realizations.

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

Only complete classifiable four-group colour vectors are permuted within species. One null assignment is reused across all 200 outer realizations for that permutation, preserving the observed geometry and the dependence induced by repeated appearances of the same measured photographs.

Primary null size: **999 permutations** with a pre-frozen seed. Each permutation is run through the exact same 200-resample schedule as the observed data. The same null therefore supplies both the G1 field test and the downstream edge/overlay alignment tests without introducing a different favourable spatial null.

G1 support requires:

1. the postmeasurement species/photo gate to pass;
2. observed recurrent-field concentration above the null mean;
3. primary one-sided Monte Carlo `p < 0.05`.

If G1 is null, external overlays do not rescue it.

## 9. Recurrent-zone extraction and stability

Persistent zones are extracted only if G1 is supported.

Within each evaluable observed realization, hotspot cells are the upper 10% of that realization's field. A cell may seed a persistent zone only if it is evaluable in at least 100 of 200 realizations and is a hotspot in at least 60% of those evaluable realizations. Adjacency includes longitude wrap. Connected components with at least three seed cells are assigned neutral IDs `Z01`, `Z02`, ... before geographic names or external predictors are opened.

Predeclared stability diagnostics include running consensus at 25, 50, 100, 150 and 200 realizations, odd/even consensus correlation, leave-one-realm-out, major-family deletion, support-scale sensitivity, high-observation-effort-cell deletion and observer-unique sensitivity. Repetition estimates sampling and map stability; it does not multiply biological sample size.

## 10. Ecological interpretation hierarchy

G3 edge-level mechanism alignment is reported regardless of G1 but cannot rescue a null G1. Primary within-species colour/external alignment controls great-circle distance so that ordinary distance-decay in both variables is not mistaken for a mechanism.

Only after G1 support and G2 stability may neutral zones be annotated using the independently frozen pollinator, climate, soil, terrain and geographic-barrier surfaces. Candidate explanations include pollinator-community turnover, CHELSA climate turnover, SoilGrids edaphic turnover, EarthEnv terrain, marine gaps, major HydroRIVERS crossings, RESOLVE realm/biome/ecoregion boundaries and GMBA mountain-system boundaries.

G5 separately compares sympatric species with matched allopatric controls. Negative `Delta = D_colour(sympatric) - mean D_colour(matched allopatric)` indicates convergence in sympatry; positive Delta indicates divergence/character displacement. These species-pair analyses are separate from G1 and do not create a global zone.

## 11. Claim ceiling

A positive G1 supports **broad recurrent geographic concentration of within-species flower-colour discontinuity within the public-photo sampling frame under the frozen observation-bias controls**. It does not demonstrate one exact universal line, does not imply that every species changes there, and does not by itself establish a causal mechanism.

A null G1 with heterogeneous G3 effects supports the narrower conclusion that spatial colour organization is species- or context-specific rather than globally synchronized. A null G1 and near-zero heterogeneous/prevalence signal supports weak general geographic organization at the tested scales.

No result is described as an unbiased census of all world flowers. Regions or taxa absent from the public-photo frame cannot be reconstructed by Monte Carlo resampling.
