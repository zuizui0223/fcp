# RGFCA independent global signal-axis replication — prospective protocol

Date frozen: 2026-09-09 JST

This protocol replaces the earlier Canadian-Rockies local-signal replication as the primary ecological replication target. The current-reserve exploration has converged: the recurrent global M1–M3 flower-minus-background signal axes survived repeated grid/observer/photo perturbation and complete observer-disjoint partitions, whereas assignment of those modes to specific geography, the Canadian Rockies and South Florida local candidates, and pooled M3 genus structure did not survive the corresponding observer/species stress tests. No local hotspot, island hypothesis, latitude effect or lineage effect is confirmatory in this tranche.

The current 50,000-photo discovery/reserve frame is outcome-opened and cannot replicate itself.

## Prerequisite measurement-validity gate

Ecological replication images and colour outcomes remain sealed until the independently frozen human-adjudicated focal-flower measurement gate at commit `1ef5c3a0983dee37b57e94d14ab4a7233eafe57a` has a retained result with `measurement_gate_passed=true` on the complete 300-image frame.

That gate requires all three prospective floors on the exact frozen ROI-v4 estimator:

1. pooled focal-target precision >= 0.70;
2. pooled focal-target recall >= 0.35;
3. median annotated-image focal-target precision >= 0.70.

No passing result for that gate is currently admitted by this protocol. Metadata-only ecological feasibility work may proceed, but ecological image pixels, segmentation outputs, palette values and M1–M3 outcomes may not be opened before the retained measurement-gate result passes. If the measurement gate fails, this ecological replication remains blocked; do not lower either gate or substitute a successful subset/model.

## Superseded local-replication target

`RGFCA_INDEPENDENT_LOCAL_SIGNAL_REPLICATION_DESIGN_20260909.md`, frozen before the later observer-multipart audits, fixed Canadian-Rockies cell `(-23,11), M2+` as the confirmatory target. Subsequent current-reserve observer-composition stress tests did not reproduce that local pattern robustly. It is therefore retained as historical design evidence only and must not be executed as the primary replication after this protocol.

The same applies to South Florida M2+, observer-consensus hotspot cells, island/mainland effects, latitudinal flower-side effects and M3 genus structure: all are secondary exploratory targets at most after the global-axis primary analysis is irreversibly closed.

## Primary scientific target

The confirmatory target is the existence of the same three-dimensional within-species flower-minus-background colour-signal subspace in a genuinely fresh, species-disjoint global tranche.

Let each local state have the 12-dimensional palette-difference vector

`D = flower_fraction_12 - matched_background_fraction_12`

with palette order:

`white, yellow, orange, red, pink, magenta, purple, blue, bronze, green, brown, black`.

Within each species, centre local-state `D` vectors by that species' mean. Form the species-equal covariance by giving each admitted species equal total weight, exactly as in the frozen discovery recurrence analysis.

The frozen discovery reference vectors are immutable and are copied here so fresh data cannot alter them:

| colour | M1 | M2 | M3 |
|---|---:|---:|---:|
| white | -0.64122299159389695 | 0.66632447750791313 | -0.014038124608976001 |
| yellow | 0.13482211192318391 | -0.26168571050850942 | -0.40184906212408278 |
| orange | 0.0098330606136441002 | -0.0044890970509536996 | -0.0067753323268837 |
| red | 0.0229963472672206 | 0.0019883085711544001 | 0.0064525068998573 |
| pink | 0.031752156249569501 | -0.089318865681069295 | -0.018878795110177399 |
| magenta | 0.019778592536036602 | -0.043515115615847602 | -0.0025633611642378001 |
| purple | 0.053498486545234702 | -0.052827773403501997 | -0.028055157245154099 |
| blue | 0.0072704874687645003 | -0.064885360637869496 | 0.0075572587118776001 |
| bronze | 0.0033793834867418002 | -0.0127002060820953 | 0.18885476761721459 |
| green | 0.70828978337774184 | 0.49478713731348112 | 0.157150536416549 |
| brown | -0.2115387641304306 | -0.42838818105298831 | 0.67670709703927578 |
| black | -0.1388586537438114 | -0.20528961335971219 | -0.56456233410526158 |

Retained discovery reference-loading CSV SHA-256: `d9dfddc4cc18052921c2b2a4262ab09da903277757f9ab251222627282c172ad`.

Discovery explained fractions for context only were M1 0.320310, M2 0.219613 and M3 0.133120. They are not replication thresholds.

## Independent species frame

The ecological replication species are fully species-disjoint from the 500 discovery/reserve species. This is deliberate: the primary claim is that the signal system generalizes beyond the taxa used to discover it.

Before any ecological pixels are opened, construct and hash a metadata-only candidate manifest using species identity, taxon ID, photo/observation/observer IDs, coordinates, date, licence and source metadata only.

Every ecological replication `photo_id`, `observation_id` and exact `observer_id` must also be disjoint from all outcome-opened or model-development/validation frames, including:

- the 50,000-photo discovery/reserve atlas;
- the 110-image Monarda validation frame;
- the 300-image independent human-reference frame;
- JRC detector development/qualification frames;
- previously opened FlowerMask examples and any later measurement-development frame.

Known overlap is disqualifying. Species selection cannot use image pixels, masks, palette values, M1–M3 scores, island status, climate, hotspot geography or any ecological outcome.

## Metadata-only ecological qualification and sampling

Target exactly **300 species** selected outcome-blind from the species-disjoint candidate pool.

A species is metadata-qualified only if, on the fixed 300-km equal-area lattice framework, it has at least **3 spatial support groups** whose centroids include a pair separated by >=500 km. Each support group must contain at least **6 distinct photographs from >=3 distinct observers**, with distinct observation IDs.

If more than 300 species qualify, rank by the following metadata-only tuple and take the first 300:

1. descending number of qualifying spatial support groups;
2. descending number of distinct observers;
3. descending maximum support-centroid separation;
4. ascending stable taxon ID as the deterministic tie break.

For each selected species retain 3–5 support groups, chosen to maximize geographic dispersion using coordinates only. Within each retained support group sample exactly 6 photos using deterministic round-robin observer balancing, with no observer contributing more than 2 photos to that group. Thus the planned ecological frame contains 18–30 photos per species (5,400–9,000 photos total).

If the metadata pool cannot supply the exact 300-species frame under these fixed rules, record `ecological_metadata_frame_qualified=false` and stop. Do not reduce the species, spatial-group, photo or observer floors after opening outcomes.

Metadata acquisition/qualification may happen before the measurement gate passes; image-pixel opening may not.

## Measurement and analytical support after the gate passes

Only after the independent measurement gate passes may the exact hashed ecological manifest be acquired/measured with the unchanged validated measurement pipeline.

The continuous 12-palette analysis uses the same admission semantics as the discovery signal analysis. A measured local state is analytically admitted only with >=3 admitted photographs from >=2 observers. A species is analytically admitted only if it retains >=2 admitted local states with centroid separation >=500 km.

The primary replication is support-evaluable only if at least **200 species** remain analytically admitted. If fewer remain, report `ecological_replication_support_adequate=false`; do not relax QC or state/species floors.

## Fresh covariance and reference comparison

The fresh tranche is allowed to estimate its own top-three eigenvectors only for comparison with the already frozen reference; fresh eigenvectors do not redefine M1–M3.

Let `U_old` be the fixed 12x3 matrix above and `U_fresh` the top-three orthonormal eigenvectors from the species-equal fresh covariance.

Primary subspace similarity is the three singular values of `U_old.T @ U_fresh`.

For mode-specific diagnostics, optimally permute the three fresh axes to maximize total absolute cosine with M1–M3, then orient signs to positive cosine. This matching rule is fixed before ecological outcome opening.

## Repeated fresh-world analysis

Run exactly **500 realizations**, seed `20260920`.

For every realization:

1. shift the 300-km lattice independently in x and y by Uniform(0,300 km);
2. rebuild local states from coordinates;
3. require >=3 admitted photos and >=2 observers per state;
4. within each state bootstrap observers with replacement and then draw one photo from each selected observer;
5. require >=2 states separated >=500 km per species;
6. construct the species-equal within-species covariance;
7. estimate the fresh top-three subspace;
8. compare it to the immutable old reference.

A realization is support-evaluable only if >=200 species remain admitted.

Report support count, species/state support distribution, all three singular values, mean singular value, minimum singular value, optimally matched M1–M3 cosines and fresh explained fractions.

## Primary replication gate

The **global three-dimensional signal subspace** is independently replicated only if all conditions hold:

1. independent measurement gate passed before ecological outcome opening;
2. exact 300-species metadata frame qualified;
3. >=400/500 fresh realizations are support-evaluable with >=200 species;
4. median minimum singular value between old and fresh three-dimensional subspaces is >=0.80;
5. at least 80% of support-evaluable realizations have minimum singular value >=0.70;
6. the observer-disjoint robustness gate below passes.

This is the primary result. A failure is retained as non-replication; no other dimensionality, palette subset, species subset or geographic subset may replace it.

## Mode-specific replication labels

M1, M2 and M3 are each labelled individually replicated only if, after the fixed optimal assignment:

- median absolute cosine with its frozen reference axis is >=0.80; and
- at least 70% of support-evaluable realizations have cosine >=0.80.

Individual-axis failure does not retroactively rotate or rename the reference. It is possible for the three-dimensional subspace to replicate while one or more individual PCA directions rotate within that subspace; report that distinction explicitly.

## Observer-disjoint fresh robustness

Because local geography and lineage structure in the discovery reserve were observer-sensitive, complete observer separation is mandatory.

Use exactly 20 fixed SHA-256 partitions of fresh exact observer IDs with salts:

`rgfca-fresh-axis-observer-00` through `rgfca-fresh-axis-observer-19`.

For each salt, assign each observer wholly to A/B by digest parity. Recompute the fresh three-dimensional covariance separately in A and B using the same fixed 20 lattice shifts (the first 20 shifts generated by seed `20260920`), with no observer/photo bootstrap inside this stress test. A half is evaluable only if >=100 species retain analytical support.

A salt passes only if both halves are evaluable and both halves have:

- minimum singular value vs `U_old` >=0.70; and
- each optimally matched M1–M3 cosine >=0.75.

The observer-disjoint robustness gate passes only if >=16/20 salts pass. No salt may be dropped or replaced after outcome opening.

## Component diagnostics

Using exactly the same fresh states/weights, retain companion projections for:

- flower-only `F`;
- matched background `B`;
- the identity `D = F - B`.

These diagnose whether a replicated D mode is predominantly flower-side, background-side or mixed. They do not replace the D-based primary replication gate and no F/B result may be used to redefine the reference axes.

## Calibration diagnostics

For context only, compare the observed old-vs-fresh subspace similarity against 10,000 deterministic random 3-dimensional orthonormal subspaces in 12 dimensions, seed `20260921`. Report the fraction whose mean and minimum singular values equal or exceed the observed fresh values. This diagnostic is not a substitute endpoint and cannot rescue failure of the fixed replication gates.

## Opening order

1. freeze this protocol and immutable old reference values;
2. metadata-only discover/qualify the exact 300-species ecological frame and hash it;
3. verify complete species/photo/observation/observer disjointness;
4. wait for retained independent human-reference result with `measurement_gate_passed=true`;
5. only then acquire/open ecological image pixels and run the unchanged validated measurement pipeline;
6. freeze measurement outputs and QC ledger;
7. execute the 500-realization primary replication exactly once;
8. execute the fixed observer-disjoint stress test;
9. independently recompute all support and similarity gates from retained outputs;
10. only after the primary replication is irreversibly reported may old local geography, islands, climate or lineage hypotheses be examined as secondary exploration.

## Stop rules

Do not:

- run the old Canadian-Rockies local replication as the primary target;
- reopen South Florida or another hotspot because global replication fails;
- add old discovery species to reach the 300-species frame;
- lower state/species/observer thresholds;
- change the 12-colour palette or frozen M1–M3 vectors;
- rotate the old axes using fresh outcomes;
- open ecological image outcomes before measurement validation passes;
- tune the measurement model on ecological replication images.

## Claim boundary

A successful result would show that a three-dimensional photo-derived flower-minus-background colour-signal system discovered in the original 500 species recurs in a new observer-, photo- and **species-disjoint** global taxon set under a separately validated focal-flower measurement target. It would materially strengthen the existence of a general recurrent signal representation.

It would still not establish pollinator-perceived colour distance, adaptation, causal environmental selection, a universal spatial boundary, an island syndrome, or fixed hotspot geography. Those require separately designed evidence after this independent axis replication is closed.