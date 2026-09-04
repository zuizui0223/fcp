# Repeated Global Flower-Colour Atlas (RGFCA)

Status: **methodological estimand frozen before any RGFCA global flower-colour field is opened**. This note makes explicit what the repeated global Monte Carlo design estimates, what repetitions do not mean biologically, and how observation availability is handled.

## 1. Core idea: repeated world-map realizations

The global analysis is not one enormous matrix containing every available photograph from every species. Instead it constructs bounded, balanced **world-map realizations**.

For realization `r = 1, ..., R`:

1. select a balanced random subset `S_r` of globally eligible species;
2. for every included species `i`, select a balanced random subset `P_ir` of classifiable photographs;
3. place those photographs at their actual geographic coordinates on the world map;
4. calculate within-species colour-discontinuity geometry without pooling species identities;
5. project equal-species discontinuity support and geographic opportunity onto the same global equal-area field;
6. save the field and the exact species/photo inclusion ledger;
7. repeat under the same frozen schedule.

Primary values remain:

- `R = 200` observed realizations;
- 250 species per realization;
- 20 photographs per included species;
- seeded balanced species inclusion;
- seeded balanced within-species photograph inclusion.

The bounded analysis unit is therefore **one realization of a global flower-photograph map**, not the complete global photograph archive.

## 2. Balanced repeated coverage

Let `c_i(R)` be the number of the first `R` realizations in which eligible species `i` appears. The balanced schedule is constructed so that, whenever mathematically possible,

`max_i c_i(R) - min_i c_i(R) <= 1`.

The same principle is used for photographs within species. High-volume public-photo species therefore do not gain larger inferential weight merely because more photographs exist for them.

If more than 250 species are eligible, no single realization contains the entire inferential species pool, but all eligible species have repeated non-zero inclusion opportunity. If exactly 250 survive, species composition is fixed while the 20-photo realization changes within species.

This is **coverage through repetition**, not an all-at-once analysis.

## 3. World-field estimands

For every evaluable global cell `x` and realization `r`, let

- `O_r(x)` = species-normalized geographic edge opportunity;
- `N_r(x)` = the same support weighted by within-species rank-standardized flower-colour discontinuity;
- `F_r(x) = N_r(x) / O_r(x)` where the frozen opportunity gate is passed.

The repeated atlas reports:

### 3.1 Consensus intensity

`F_bar(x) = mean_r F_r(x)` over evaluable realizations.

### 3.2 Realization uncertainty

For every cell, report across-realization SD and MAD of `F_r(x)` together with the number of evaluable realizations.

### 3.3 Hotspot recurrence

Within each realization, mark the frozen top 10% of evaluable field cells and define

`P_hot(x) = (# evaluable realizations in which x is a top-decile hotspot) / (# realizations in which x is evaluable)`.

The persistent-zone rule remains `P_hot >= 0.60`, at least 100 evaluable realizations, followed by the frozen connected-component rule. Neutral IDs `Z01`, `Z02`, ... are assigned before named geography or mechanism overlays are inspected.

This recurrence probability is the direct translation of the intuitive procedure: **scatter a balanced set of flower photographs over the world map many times and ask which flower-colour transition regions keep coming back**.

### 3.4 Monte Carlo convergence

The atlas reports change in the consensus field and recurrent-zone support at `R = 25, 50, 100, 150, 200`.

## 4. Matched species-conditioned null

Repetition alone is not evidence of biological structure. Every null replicate uses the same species/photo sampling opportunity and the same fixed coordinates/graph geometry as the observed program. Only complete classifiable colour vectors are permuted within species.

The comparison is therefore recurrent observed flower-colour geography versus recurrent colour geography expected from exactly the same sampled geometry after within-species colour-location association is broken.

## 5. What repetition does and does not buy

Repetition provides bounded computation, equal-species sampling opportunity, broad coverage of a larger eligible species pool, direct stability estimates, and separation of instability caused by taxon composition, photograph composition or spatial support.

It does **not** create 200 independent biological worlds, multiply biological sample size by 200, or recover species/regions that have essentially zero public-photo inclusion probability.

## 6. Repeated global discovery

The method also treats the upstream species universe as a sampling problem. The original one-pass baseline contained 8,989 species. A first 20-round metadata-only census discovered 30,393 species in the combined union, but a technical audit showed that adjacent `order_by=random, page=1` rounds could be nearly duplicated by upstream random-page persistence. Those records remain valid discoveries, but the V1 odd/even overlap is not used as independent-repetition evidence.

A cache-resistant V2 was therefore frozen before any RGFCA colour field using stable ID ordering, deterministic cell-specific pages and explicit V1 ID exclusion. The capacity universe is the deduplicated V1+V2 species union after successful V2 completion.

## 7. Observation-bias architecture

The repeated-atlas framework reduces **dominance bias**, not all observation bias.

Primary design protections include:

- equal-area global discovery attempts;
- observer cap of two retained candidate photographs per observer/species;
- deterministic geographic maximin raw-photo selection;
- balanced species and photograph inclusion across realizations;
- equal-species contribution to the global field;
- division of colour support by geographic edge opportunity;
- leave-one-realm-out and major-family deletion;
- location-blind image measurement.

### 7.1 Empirical target-group sampling availability

Before any RGFCA global colour outcome, three count-only surfaces were frozen over the same 162 equal-area cells. All 486 requests succeeded with zero errors.

The public-photo frame is extremely concentrated:

- all target-group research photo records: 53,270,601; equal-area-cell Gini **0.885**;
- reusable-licence records: 40,554,393; Gini **0.882**;
- flowering-annotated reusable records: 4,049,721; Gini **0.885**;
- the 17 cells in the pre-frozen top observation-effort decile contain **84.8%**, **84.2%** and **83.3%** of these three record pools, respectively.

Median reusable-licence retention across defined cells is 0.790. Median flowering-annotation retention conditional on reusable licence is only 0.095. Geographic ranks remain highly similar across the filters (`rho=0.999` all-vs-licence; `rho=0.983` all-vs-flowering), showing that filtering reduces volume far more than it removes the geographic observation concentration.

These are target-group **sampling-availability proxies**, not pure observer-effort probabilities, because biological plant availability also varies geographically.

### 7.2 Platform-wide activity control

A separate all-taxa research-photo count surface is prospectively frozen before RGFCA colour to provide a broader platform-activity proxy. This allows target-group availability to be distinguished, imperfectly, from areas that are simply intensively used by the observation platform. The all-taxa layer is a negative control, not a correction that replaces the primary RGFCA field.

### 7.3 Frozen robustness checks

A positive primary G1 is not automatically called observation-robust. The atlas also reports:

1. an observer-unique sensitivity with at most one classifiable photograph per observer/species/realization when at least 20 distinct observers are available;
2. complete reanalysis after deleting the 17 pre-frozen highest-observation cells;
3. opportunity-weighted field correlations with target-group record density, reusable-licence fraction, flowering-annotation fraction and platform-wide activity, evaluated against the same 999 species-conditioned null fields;
4. observer concentration diagnostics and primary-versus-observer-unique field similarity.

A colour zone that disappears after high-effort deletion, or that simply follows a sampling-availability layer without exceptional separation from the matched null, is classified as observation-sensitive rather than robust biogeography.

### 7.4 Residual bias that cannot be claimed away

The method cannot identify unrecorded taxa/regions, fully remove preferential uploading of unusual colour morphs, remove colour-dependent identification/annotation errors, or perfectly separate biological abundance from human observer effort. Consequently the claim ceiling is recurring flower-colour geography **within the measurable public-photo sampling frame**, robust to the explicitly measured observation-availability structure—not an unbiased census of all flowers on Earth.

## 8. Biological hierarchy

RGFCA separates:

1. **global recurrence:** do independent species repeatedly place strong colour discontinuity in the same broad geographic regions? (`G1/G2`)
2. **species-specific organization:** how heterogeneous are within-species spatial colour effects? (`G3`)
3. **mechanistic concordance:** are colour discontinuities aligned, beyond geographic distance alone, with pollinator turnover, climate, soil, terrain, marine gaps, rivers, mountain systems or established biogeographic boundaries? (`G3/G4`)
4. **sympatry/allopatry:** are co-occurring species more colour-similar or more colour-divergent than matched non-co-occurring controls? (`G5`)

A null shared global field does not imply that flower colour has no biogeography. Species-specific effects or ecological concordance cannot retroactively rescue a null recurrent-zone test.

## 9. Methodological claim ceiling

The defensible methodological claim is not that Monte Carlo, bagging or spatial resampling is new. It is:

> We operationalize global flower-colour biogeography as a balanced repeated-atlas problem: bounded multispecies realizations of public flower photographs are repeatedly mapped at their true coordinates, converted to opportunity-corrected colour-discontinuity fields, and summarized by geographic hotspot recurrence against a species-conditioned matched null, with observation-availability surfaces frozen independently as negative controls.
