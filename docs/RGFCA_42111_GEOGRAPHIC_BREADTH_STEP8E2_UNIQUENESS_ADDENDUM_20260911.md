# RGFCA Step 8E2 addendum — unique image identity for taxon×cell anchors

Date: 2026-09-11 JST
Status: frozen after metadata-only conflict diagnosis and before any Step-8E image pixel or flower-colour outcome is opened.

## Trigger

The first taxon×cell anchor builder applied the already frozen within-pair minimum-hash rule independently to all 85,337 taxon×cell pairs, then correctly failed because selected photo IDs were not globally unique.

A metadata-only diagnostic established:

- taxon×cell pairs: 85,337;
- duplicated selected observation IDs: 0;
- duplicated selected photo IDs: 1;
- selected pair rows affected by any duplicate ID: 2;
- source observation IDs spanning more than one pair: 0;
- source photo IDs spanning more than one pair: 2;
- deterministic constrained greedy matching can assign globally unique observation and photo IDs to all 85,337 pairs;
- unresolved pairs under that matching: 0;
- image pixels opened: false;
- flower colour used: false.

Thus the failure is an identity-allocation conflict in frozen metadata, not a colour or biological result.

## Frozen constrained selection rule

Retain the original candidate universe and candidate hash:

`candidate_hash = SHA256('20260911|cell_id|inat_taxon_id|observation_id|photo_id')`.

Define a deterministic pair processing order:

`pair_hash = SHA256('20260911|pair|cell_id|inat_taxon_id')`.

Process the 85,337 taxon×cell pairs in ascending `pair_hash`. Within each pair, process candidate rows in ascending `candidate_hash`, then `observation_id`, then `photo_id`. Select the first candidate whose observation ID and photo ID have not already been assigned to an earlier pair.

If no unused candidate exists for a pair, retain that pair as unresolved; do not relax identity constraints and do not substitute after colour opening.

The completed metadata diagnostic showed that this rule can match all 85,337 pairs without unresolved pairs. That diagnostic is a feasibility check only; the frame is frozen by a separate execution of this exact rule.

## Why this is not outcome tuning

The rule is introduced because the original frame violated an intended technical requirement: one durable image identity must not be counted as two separate taxon×cell observations. The correction uses only pre-existing species/taxon/cell/observation/photo metadata and deterministic hashes. No image byte, flower colour, ecological effect, literature label, climate variable, or downstream inference has been inspected.

## Estimator boundary unchanged

- the 42,111 one-anchor-per-species layer remains the species-equal global composition estimator;
- the 85,337 unique taxon×cell layer remains the cell-wise geographic estimator;
- the two layers are not pooled into one frequency estimate;
- neither layer establishes within-species polymorphism, C*, S*, modal species colour, or focal-species petal truth.
