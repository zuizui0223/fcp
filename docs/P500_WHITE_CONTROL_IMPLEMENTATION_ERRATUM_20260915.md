# P500 white-control implementation correction (2026-09-15)

## Outcome first

This change qualifies **offline implementation safeguards**, not P500 image
opening. The original protocol and freeze receipt are preserved byte-for-byte.
The receipt's asserted historical flags are not independently verified chronology.
No pixels, biological response, new coupling model or H2 test are opened here.

## Input-schema reconciliation

The path called `preflight_result_path` in the original receipt actually contains
a **candidate metadata acquisition** result. It has `decision` and
`outcome_firewall`, not the `gate_results` and `blindness` schema required by
`validate_preopening_firewall`. The earlier test incorrectly expected that real
file to be accepted, and another test attempted to mutate a nonexistent key.

The strict runtime validator remains fail-closed on that file. We have not
invented missing flags, rewritten the acquisition record, or relaxed the original
opening contract to turn a failed check into an authorization.

`inspect_candidate_metadata` now recognizes the actual acquisition schema and
reports its documented metadata support separately. It **always** reports
`execution_chronology_verified=false` and `opening_authorized=false`.
Its scope is the supplied result JSON, not a complete metadata-census or
prior-ID exclusion audit. An explicit, evidence-backed execution handoff and
chronology contract remain necessary before image acquisition.

The corrected tests assert rejection of the real acquisition file as an opening
record, acceptance of the unchanged receipt with an explicitly synthetic complete
schema, and rejection of corrupted inputs. Synthetic declarations are never saved
as research provenance. A passing test suite means this separation is enforced,
not that the missing historical evidence has been supplied.

## Fixed-condition and numerical guards

- Reject receipt changes to the primary predictor, diagnostic list, equivalence
  interval, retention floor, decision states, counts, status, and opened outcomes.
- Reject impossible sensitivity counts above the primary count, non-integral
  counts, truthy strings/numbers instead of booleans, and invalid OR intervals.
- Preserve the existing valid-input CLEAR / FLAGGED / INDETERMINATE rule and
  its 0.80–1.25 and 90% thresholds. No precedence rule is reinterpreted.
- Require integral decoded RGB samples and integer decode bit depth 1–32;
  unsupported depths raise an error rather than truncate silently.

## CI scope

The new `p500-white-control-offline.yml` runs only unit tests and static JSON
inspection. Repository permission is read-only and checkout credentials are not
persisted. The artifact explicitly preserves `opening_authorized=false`.

The inherited recovery workflow's duplicate step-level `env` was consolidated
without changing either variable, the command, trigger, or scientific rule.
That workflow is **not dispatched** by this correction; closed acquisitions are
not restarted. CI success does not validate coupling-model implementation,
image segmentation, raw-data chronology, or prospective H2 results.
