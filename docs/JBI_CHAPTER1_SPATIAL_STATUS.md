# JBI Chapter 1 spatial execution status

## Active branch and PR

- branch: `analysis/jbi-global-colour-boundaries`
- stacked base: `audit/upstream-spatial-reaudit` (PR #17)
- active PR: #20

## Primary question

The primary Chapter 1 inference is now:

**global flower-colour map → species-conditioned random-labelling test → species-level transition boundaries → shared-boundary strength → post-discovery geographic correspondence**.

The previous C/S-to-climate comparison is retained as secondary/provenance context and is not the primary spatial discovery analysis.

## Completed on this branch

- species-conditioned permutation core: implemented;
- no between-species colour-label exchange: guarded by tests;
- label-independent boundary detectability denominator `A(x)`: implemented;
- `A(x)=0` / insufficient support: treated as not evaluable rather than biological zero;
- detectability invariance under colour-label permutation: guarded;
- Chapter 1 protocol: frozen;
- iNaturalist measurement-gate protocol: frozen;
- deterministic 6 x 200 photo split contract: implemented;
- calibration/evaluation allocation: 80/120 per species, 480/720 total;
- split assignment basis: species + stable photo ID only;
- row-order/geography/observer/date invariance: tested;
- post-outcome source manifests: rejected for recognized colour/visibility/segmentation fields;
- source/split/freeze-manifest all-or-none CI contract: implemented;
- interrupted-worktree recovery scanner: implemented;
- recovery scanner accepts only an exact 6 x 200 pre-measurement manifest with globally unique photo IDs;
- equivalent duplicate copies are tolerated when their species/photo-ID hash is identical;
- conflicting eligible ID sets fail closed instead of choosing one silently.

## Empirical state

- acquired development photographs: 1,200 (6 species x 200);
- spatial/observer/season overlap acquisition controls: passed before this branch;
- flower-colour classifications: 0;
- calibration photographs measured: 0/480;
- held-out photographs measured: 0/720;
- random versus non-random spatial placement: `not_evaluated`;
- species-level colour boundaries: `not_evaluated`;
- shared-boundary surface: `not_evaluated`;
- geographic correspondence: `not_evaluated`.

## Current blocking input

The acquired 1,200-photo source manifest is not present in the remote repository or in the available GitHub Actions artifacts inspected for PR #17. The current execution environment also does not mount the original Windows worktree at `C:/Users/zuizui/.codex/worktrees/fcp-spacetime-recovery`.

The canonical import target remains:

`data/frozen/jbi_ch1_photo_source_manifest.csv`

The local worktree no longer requires the original filename to be known. From the root of `fcp-spacetime-recovery`, run:

```bash
python scripts/data/recover_jbi_ch1_photo_source_manifest.py --root .
```

The scanner recursively inspects CSV/TSV acquisition tables and writes the canonical source only when the candidate is unambiguous under the frozen contract. It also writes:

`docs/supporting/jbi_ch1_photo_source_recovery_report.json`

The report records every inspected table, rejection reason, eligible candidate, canonical species/photo-ID hash, and any equivalent duplicate paths.

After source recovery, materialize the frozen split with:

```bash
python scripts/data/freeze_jbi_ch1_photo_split.py \
  data/frozen/jbi_ch1_photo_source_manifest.csv \
  data/frozen/jbi_ch1_photo_split_v1.csv \
  data/frozen/jbi_ch1_photo_split_v1_manifest.json
```

This deterministically materializes:

- `data/frozen/jbi_ch1_photo_split_v1.csv`;
- `data/frozen/jbi_ch1_photo_split_v1_manifest.json`.

No colour labels, geographic predictors, observer information, or dates participate in assignment.

## Next gate

1. recover/import the original 1,200-photo acquisition manifest without adding measurement outcomes;
2. materialize and commit the deterministic 480/720 split plus hashes;
3. run the 480-photo blinded visibility/segmentation/colour-code calibration;
4. freeze the measurement codebooks and failure rules;
5. only then open the 720-photo evaluation set.
