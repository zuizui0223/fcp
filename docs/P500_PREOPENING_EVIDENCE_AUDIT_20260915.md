# P500 pre-opening evidence: unresolved, not an ecological result

**Scope correction:** this is a comparison of the older PR 33 control records,
not a complete current execution inventory. After the user clarified that all
execution was in Actions, a later authorized measurement run with successful
pixel-processing jobs was found. The
[Actions correction](FCP_H123_P500_ACTIONS_CORRECTION_20260915.md) supersedes
any statement here that P500 globally remains unopened or has only metadata.
The missing fields remain a real schema difference, not proof of absent
execution evidence. No missing declaration has been fabricated or accepted.

Audit date: 2026-09-15. Current PR 33 head was checked as
`9e87458185bba8b6df901280b7a7078ce3544de8`. No acquisition, image access,
response classification or prospective test was dispatched in this audit.

## Decision and scope

The older control validator still grants no opening authorization. Its
pre-opening evidence state is INDETERMINATE;
this is **not a computed measurement-control verdict** and is not H2
NON-SUPPORT. Both scientific result fields in the original receipt remain null.
There is no evidence here that an undisclosed opening actually occurred;
absence of adequate non-access evidence is not proof of a chronology violation.

## Field-level comparison of the retained records

The strict validator requires all nine entries below in `preflight.blindness`.
The candidate record instead provides `outcome_firewall`, while the receipt
declares all nine false in `required_blind_flags`.

| Required field | Same-name candidate declaration | Receipt declaration | Independent chronology established |
| --- | --- | --- | --- |
| image_results_requested | missing | false | no |
| filtered_specimens_opened | missing | false | no |
| image_urls_opened | missing | false | no |
| image_bytes_opened | missing | false | no |
| flower_colour_opened | false | false | no |
| palette_opened | false | false | no |
| D_opened | false | false | no |
| H1_opened | missing | false | no |
| H2_W_opened | false | false | no |

The candidate's `image_pixels_opened=false` is related evidence, but is not
silently renamed to `image_bytes_opened`; the scope of each event requires an
explicit evidence-backed interpretation. The terms `image_results_requested`,
`filtered_specimens_opened` and `image_urls_opened` are not operationally defined
in enough detail by the control protocol to infer them from API metadata
acquisition. Metadata containing a photo URL is not by itself proof that its
image bytes were accessed. Conversely, a false declaration is not an exhaustive
record of all access paths. H1 here must refer to the P500 cohort, not the
already completed H1 analyses in different cohorts.

There is a separate schema mismatch: the candidate decision is expressed in
`decision.verdict`, whereas the opening validator expects
`gate_results.overall_pass`. This audit does not create a converted preflight
or change the validator to accept these records.

## Sources actually inspected

- [Candidate result](../results/polymorphism_h2_p500_candidate_metadata_20260913/result.json):
  500 selected species, 499 full-100 species, 49,998 metadata rows, zero request
  errors. Its own next-step wording requires separate pixel authorization.
- [Original receipt](../results/polymorphism_h2_p500_white_measurement_control_freeze_20260913/receipt.json)
  and [original control protocol](POLYMORPHISM_H2_P500_WHITE_MEASUREMENT_CONTROL_PROTOCOL_20260913.md).
  Git records the protocol commit as `8604e7a9250b0dfe74e1fbac50e8f8104f45be16`
  at 2026-09-14 23:36:31 +09:00 and receipt commit as
  `681fb8381f540f80a51aa17c036ea5f3517cb269` at 23:37:22.
  The date in a filename is not substituted for these commit timestamps;
  neither is authenticated proof of all historical non-access.
- [Acquisition runner](../scripts/acquisition/run_polymorphism_h2_p500_candidate_metadata_20260913.py),
  [acquisition workflow](../.github/workflows/polymorphism-h2-p500-candidate-metadata.yml)
  and [identity preflight](../scripts/analysis/preflight_polymorphism_h2_p500_acquisition_20260913.py).
  Their inspected declarations/checks use the earlier metadata schema, not all
  nine later required fields. This is a source inspection, not proof that every
  dependency or external access channel was exhaustively audited.
- [Current strict validator](../fcp_pipeline/p500_white_measurement_control.py)
  still rejects the real candidate as an opening record; offline metadata
  inspection always returns opening_authorized=false.
- [PR 33](https://github.com/zuizui0223/fcp/pull/33) explicitly retains the
  chronology blocker. Publication-branch implementation work does not supersede it.

## What can and cannot close this gate

Technical recording, response-join verification and an artificial-data-qualified
conditional fitter now exist on the publication branch. These improve future
execution integrity but cannot supply earlier events or missing historical
evidence. The original protocol, receipt and candidate result are unchanged.

To evaluate an execution handoff, obtain contemporaneous records tied to the
exact cohort/artifacts and revisions, explain each field's event meaning, and
identify which access paths and periods the records cover. Any gaps remain
unknown. Do not author an all-false replacement, backdate a receipt, rerun
acquisition to overwrite provenance, or reinterpret missingness as false.
Even complete provenance would not by itself validate focal-petal measurement
or justify calling a digital-highlight control physical colour calibration.

This comparison alone supplies neither prospective white-axis confirmation nor
a biological negative. The later Actions measurement must be evaluated from
its actual authorization and terminal records, without equating location
blindness to the older highlight-control gate. Do not restart it or retrofit
this branch's control onto already opened pixels. Existing H1–H3 reporting can
continue while complete, immutable prospective outputs remain pending.
