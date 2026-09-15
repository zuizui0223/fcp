# Correction: P500 measurement is already running in GitHub Actions

Checked at 2026-09-15 05:29:32 UTC (14:29:32 JST). This supersedes the
publication branch's statements that P500 pixels remain unopened. The user
clarified that execution was entirely in GitHub Actions. The earlier audit
looked at PR 33 and an older checkout, not the full current Actions inventory.
That was insufficient to establish the current state.

## Authoritative live evidence

[Run 34919485994](https://github.com/zuizui0223/fcp/actions/runs/34919485994),
`Prospective H2 P500 location-blind measurement`, uses commit
`f403606b5611d895988a6a4f15d9ca1e281ce848` on
`analysis/polymorphism-42111-h1-h2-gates-20260912`.
The top-level run was reported queued, but its job inventory contained:

- authorization/firewall job: completed successfully;
- measurement partitions: 100 completed successfully, 8 in progress, 148 queued;
- total specified partitions: 256; planned population: 499 species / 49,900 rows.

These are job counts, **not verified measured-image or classifiable-species
counts**. No partial colour/palette output was downloaded or joined in this
audit. The snapshot is time-specific; poll this same run rather than restart it.

## Observed chronology, not reconstructed false flags

1. Metadata run [34728801295](https://github.com/zuizui0223/fcp/actions/runs/34728801295)
   completed its draw, metadata checks and artifact upload on September 13.
   Its overall failure was at the Git commit/rebase step (`unstaged changes`),
   not image measurement or acquisition-capacity failure. Artifact 10308723712
   has digest `b4918cd9d597edc5a9a744344ee7efa25eed66dbbbc30dad20ace788a9a073ef`.
2. Recovery [34730873300](https://github.com/zuizui0223/fcp/actions/runs/34730873300)
   verified original artifact hashes and committed it without reacquisition.
3. The later measurement protocol was committed as
   `b90dc53e90cc26d7c797bdc26d4e07848ad16016` at September 15 10:57:03 JST.
   Its one-shot authorization was committed as the run head at 11:00:18 JST.
4. The authorization/firewall job ran 11:00:25–11:01:17 JST. It validated the
   exact protocol revision, 499/49,900 population and frozen artifact hashes.
5. Partition (0,0,0), job 104224376226, acquired its images during
   11:02:23–11:02:56 JST and ran ROI/palette measurement during
   11:02:56–11:15:34 JST. It sealed and uploaded its terminal output successfully.
   This alone contradicts a current assertion of globally unopened P500 pixels.

Artifact API metadata confirms a firewall artifact (10376569553, digest
`c650c55ae8cd804fdacd36035ac9aab25b5201638a31804a0ce96e94c84e3668`) and
partition (0,0,0) artifact (10378280183, digest
`d07b863be0df21d361a395c463713ad76e1161bde6820202d7317445a43c4b04`).
Their contents and all terminal row counts have not yet been independently
verified. The individual partition is a timing example, not a claim that it
was the earliest acquisition across all workers.

## Which protocol is executing

The run's [exact workflow](https://github.com/zuizui0223/fcp/blob/f403606b5611d895988a6a4f15d9ca1e281ce848/.github/workflows/polymorphism-h2-p500-prospective-measurement.yml)
and [measurement protocol](https://github.com/zuizui0223/fcp/blob/f403606b5611d895988a6a4f15d9ca1e281ce848/docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260915.md)
specify a separate one-shot location-blind measurement route. The workflow
withholds source/species/location fields from workers, requires all 256 terminal
receipts and 49,900 unique terminal results before metadata-colour joining,
then runs the fixed q_white/W structured-null test. These are inspected design
requirements; final execution compliance is not yet established.

The inspected workflow and its H2 runner do not invoke the publication branch's
new near-clip coupling/snapshot helpers or the older PR 33 measurement-control
gate. Location blindness is not equivalent to response-blind highlight control.
Do not claim this run satisfies every condition of the earlier control contract,
or interpret a future H2 support label as a CLEAR artifact-control verdict.
The [source-level compatibility audit](P500_PROTOCOL_COMPATIBILITY_AUDIT_20260915.md)
documents the missing correspondence; terminal artifact verification remains required.

The publication branch's synthetic coupling work was authored **after actual
P500 pixel opening elsewhere in the repository**. It remains synthetic-only
development, but must not be described as a global pre-opening freeze. No
retroactive addition, fabricated historical flag, remeasurement or threshold
change is authorized by this correction.

## Next step

Leave the active run undisturbed. Check terminal execution and complete receipts
before reading/joining scientific output; preserve every failure and missing row.
If the run completes, audit the immutable aggregate and its claim ceiling,
including the unresolved highlight-control compatibility and Monarda failure.
Current H1–H3 frozen outcomes remain unchanged. P500's scientific result is still
pending, not confirmed, negative, or cleared of measurement artifacts.
