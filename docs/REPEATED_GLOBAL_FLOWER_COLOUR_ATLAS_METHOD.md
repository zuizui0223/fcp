# Repeated Global Flower-Colour Atlas (RGFCA)

Status: **methodological estimand frozen after metadata-only species discovery and before any RGFCA flower-colour field is opened**. This note does not change the numerical gates in `GLOBAL_MONTE_CARLO_BARRIER_ATLAS_PROTOCOL.md`; it makes explicit what the repeated global Monte Carlo design estimates and what the repetitions do not mean biologically.

## 1. Core idea: repeated world-map realizations

The global analysis is not defined as one enormous matrix containing every available photograph from every species. Instead, it constructs a sequence of bounded, balanced **world-map realizations**.

For realization `r = 1, ..., R`:

1. select a balanced random subset `S_r` of globally eligible species;
2. for every included species `i`, select a balanced random subset `P_ir` of classifiable photographs;
3. place those photographs at their actual geographic coordinates on the world map;
4. calculate within-species colour-discontinuity geometry without pooling species identities;
5. project equal-species discontinuity support and geographic opportunity onto the same global equal-area field;
6. save the field and the exact species/photo inclusion ledger;
7. repeat under the same frozen schedule.

Primary values remain those already frozen in the prospective protocol:

- `R = 200` observed realizations;
- 250 species per realization;
- 20 photographs per included species;
- seeded balanced species inclusion;
- seeded balanced within-species photograph inclusion.

The bounded analysis unit is therefore **one realization of a global flower-photograph map**, not the complete global photograph archive.

## 2. Balanced repeated coverage

Let `c_i(R)` be the number of the first `R` realizations in which eligible species `i` appears. The balanced schedule is constructed so that, whenever mathematically possible,

`max_i c_i(R) - min_i c_i(R) <= 1`.

The same principle is used for photographs within species. Thus high-volume iNaturalist species do not gain larger inferential weight merely because more photographs exist for them.

If more than 250 species are eligible, no single realization contains the entire inferential species pool, but all eligible species have repeated non-zero inclusion opportunity. If exactly 250 survive, species composition is fixed while the 20-photo realization changes within species.

This is **coverage through repetition**, not an all-at-once analysis.

## 3. World-field estimands

For every evaluable global cell `x` and realization `r`, let

- `O_r(x)` = species-normalized geographic edge opportunity;
- `N_r(x)` = the same support weighted by within-species rank-standardized flower-colour discontinuity;
- `F_r(x) = N_r(x) / O_r(x)` where the frozen opportunity gate is passed.

The repeated atlas reports at least four distinct quantities:

### 3.1 Consensus intensity

`F_bar(x) = mean_r F_r(x)` over evaluable realizations.

This asks where colour discontinuity is consistently strong conditional on sampled geographic opportunity.

### 3.2 Realization uncertainty

For every cell, report the across-realization SD and MAD of `F_r(x)` together with the number of evaluable realizations. A visually strong mean field with high realization variance is not a stable biogeographic result.

### 3.3 Hotspot recurrence

Within each realization, mark the frozen top 10% of evaluable field cells. Define

`P_hot(x) = (# evaluable realizations in which x is a top-decile hotspot) / (# realizations in which x is evaluable)`.

The prospective persistent-zone rule remains `P_hot >= 0.60`, at least 100 evaluable realizations, followed by the frozen connected-component rule. Neutral IDs `Z01`, `Z02`, ... are assigned before named geography or mechanism overlays are inspected.

This recurrence probability is the most direct statistical translation of the intuitive procedure: **scatter a balanced set of flower photographs over the world map many times and ask which flower-colour transition regions keep coming back**.

### 3.4 Monte Carlo convergence

The atlas reports the change in the consensus field and recurrent-zone support at fixed checkpoints (`R = 25, 50, 100, 150, 200`). A proposed global band is not considered robust simply because it appears in the final realization.

## 4. Matched species-conditioned null

Repetition alone is not evidence of biological structure. Every null replicate uses the same species/photo sampling opportunity and the same fixed coordinates/graph geometry as the observed program. Only complete classifiable colour vectors are permuted within species.

Therefore the comparison is not

`observed global map` versus `spatially random points`.

It is

`repeated observed flower-colour maps on the sampled geometry`

versus

`repeated within-species colour-permuted maps on exactly the same sampled geometry`.

This preserves observer/geographic sampling geometry, species range geometry, measured missingness and the balanced inclusion schedule. A recurrent observed region is interesting only if comparable recurrence is not produced by this matched null.

## 5. What repetition does and does not buy

### Repetition does provide

1. **bounded computation** — the cost of one realization stays fixed even as the measured global pool grows;
2. **sampling fairness** — no data-rich species is allowed to dominate one monolithic global fit;
3. **coverage of a large species pool** — species not present in one realization can enter later realizations;
4. **direct stability estimates** — field intensity, zones and species effects have observable resampling distributions;
5. **failure localization** — instability caused by taxon composition, photograph composition, realm deletion or spatial support can be separated.

### Repetition does not provide

1. 200 biological replicates from one species;
2. 200 independent worlds;
3. additional biological information for a species represented by one photograph;
4. permission to multiply nominal sample size by the number of realizations;
5. a substitute for enough independently sampled species and enough within-species photographs.

Biological replication is carried by the number and coverage of independently sampled species/photographs. The outer repetitions estimate **sampling/field stability conditional on that biological pool**.

## 6. Why repeated global discovery is part of the method

The method repeats not only colour analysis but the upstream global species census before pixels are opened.

The completed metadata-only discovery provides a strong empirical reason for doing so:

- the original one-pass baseline contained 8,989 species;
- 20 fresh equal-area discovery rounds found 29,304 fresh species;
- the combined census contains 30,393 species;
- 21,404 species were new beyond the baseline;
- odd-versus-even fresh-round species-set Jaccard = 0.9844;
- baseline-versus-fresh Jaccard = 0.2599;
- all 3,240 fixed requests completed without error;
- no candidate image pixel or flower colour was opened during discovery.

Thus a single random global pass was not an adequate description of the available taxonomic sampling frame, whereas independent halves of the repeated fresh census converged on nearly the same species set. This is a methodological result in its own right: **global photographic biodiversity coverage should be estimated through repeated, auditable realizations rather than assumed from one random draw**.

## 7. Relationship to biological hypotheses

RGFCA separates three levels that must not be conflated:

1. **global recurrence:** do independent species repeatedly place strong colour discontinuity in the same broad geographic regions? (`G1/G2`)
2. **species-specific organization:** even if no shared zone exists, which species show spatial colour structure and how heterogeneous are those effects? (`G3`)
3. **mechanistic concordance:** are colour discontinuities aligned, beyond geographic distance alone, with pollinator turnover, climate, soil, terrain, marine gaps, rivers, mountain systems or established biogeographic boundaries? (`G3/G4`)

The sympatry/allopatry extension (`G5`) is separate again: it asks whether co-occurring species are more colour-similar or more colour-divergent than matched non-co-occurring controls.

A null shared global field therefore does not imply that flower colour has no biogeography. Conversely, species-specific effects or ecological concordance cannot be used to retroactively declare a null recurrent-zone test positive.

## 8. Methodological claim ceiling

Until a literature comparison is completed, the defensible methodological claim is not “the first method ever to resample a global biodiversity map.” The claim is narrower:

> We operationalize global flower-colour biogeography as a balanced repeated-atlas problem: rather than fitting one all-at-once model to the complete opportunistic photograph archive, we generate bounded random world-map realizations, preserve equal-species contribution and sampling opportunity, and estimate the recurrence and stability of geographic colour structure against a species-conditioned matched null.

The scientific novelty must ultimately be evaluated against existing spatial bootstrap, bagging, biodiversity-atlas and Monte Carlo mapping literature; the implementation above is frozen independently of that novelty assessment.
