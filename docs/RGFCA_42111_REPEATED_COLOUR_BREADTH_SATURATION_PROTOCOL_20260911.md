# RGFCA 42,111 repeated colour-breadth saturation protocol

Date: 2026-09-11 JST
Status: frozen while the authorized blind Step-8F partitions were still running, before the complete 42,111-species species-colour join or global colour composition was available.

## Question

What is the smallest number of species needed to reproduce the species-equal global distribution of one frozen observed coarse flower-colour state per discovered species, conditional on the same fixed measurement/classification pipeline?

This is a measurement-saturation question, not a claim that one photograph identifies a species' modal colour or polymorphism.

## Full reference

The reference is the complete terminal 42,111-species breadth experiment. All 42,111 species remain in the denominator. The one anchor that failed before pixel opening and all later acquisition/ROI/colour-admission failures remain missing; none may be replaced.

The reference four-state composition is calculated among classifiable anchors, with the classifiable fraction against all 42,111 species reported alongside it.

## Repeated subsampling

Species counts:

`N = 100, 250, 500, 1,000, 2,000, 5,000, 10,000, 20,000, 30,000, 42,111`

For every N < 42,111, draw 200 samples without replacement from the full 42,111 species, using NumPy `default_rng(20260911 + replicate)` and a fresh permutation for each replicate. The sample is defined before classifiability is consulted. N=42,111 is the full reference.

For each sample report:

- classifiable fraction;
- four-state composition among its classifiable anchors;
- L1 distance from the full four-state composition;
- cosine similarity to the full four-state composition;
- maximum absolute deviation across the four colour-group proportions;
- absolute deviation in classifiable fraction from the full experiment.

## Prespecified saturation criterion

A species count N is `composition_saturated` only when at least 95% of its 200 repeated samples simultaneously satisfy:

- four-state L1 distance <= 0.05;
- maximum absolute colour-group proportion deviation <= 0.025;
- absolute classifiable-fraction deviation <= 0.03.

Cosine similarity is descriptive and is not part of the pass gate because it can be insensitive to some compositional shifts.

The **minimum repeated composition dataset** is the smallest N satisfying this criterion. No intermediate N may be added after seeing the outcome to manufacture a lower minimum. If no tested N below 42,111 passes, the result is `no_subsampled_minimum_passed` rather than relaxing the gate.

## Interpretation ceiling

This saturation result identifies how many sampled species are sufficient to reproduce this 42,111-species one-anchor observed-state composition. It does not establish how many species are sufficient for within-species D, polymorphism prevalence, C*/S*, taxonomic mechanisms, environmental mechanisms, or the global geographic map. Those are separate estimands with separate depth and taxon×cell requirements.
