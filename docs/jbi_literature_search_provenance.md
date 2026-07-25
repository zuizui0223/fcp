# Journal of Biogeography literature-search provenance

This record reconstructs the automated literature-discovery chronology from repository commits. Dates and times are commit timestamps converted to Japan Standard Time (UTC+09:00). A commit documents that the corresponding code or generated output entered the repository; it does not by itself identify a human screener or prove that all manual review occurred at that exact time.

## Verified chronology

| Date and time (JST) | Stage | Commit | Verified repository evidence |
|---|---|---|---|
| 16 July 2026, 11:38:52 | Global discovery workflow introduced | `5147a4f4f3a1cbf8a3be4566fad6baff5fefaafe` | Added the OpenAlex workflow using eight broad English queries, two cursor pages per query and up to 200 works per page. |
| 16 July 2026, 13:24:47 | First global discovery output and QC committed | `9e8f5322e8deb67cf70e27528eb1a03950a72e45` | Committed 631 mapped works, 430 candidate species and QC covering 96 families. |
| 16 July 2026, 14:29:18 | Additional global discovery/QC output | `bb81b0e251a81892ac3901898a5b03ca64e4cb7d` | Recorded another generated global-discovery/QC state. |
| 16 July 2026, 14:30:04 | Additional global discovery/QC output | `da0db3217585c26c2889d07d3b375a6e5290e1c5` | Recorded a subsequent generated global-discovery/QC state. |
| 16 July 2026, 19:45:37 | Expanded evidence dataset committed | `84eeeed5e4c5b7e32d5478c9cc688cfc94905271` | Expanded the global flower-colour evidence dataset following vocabulary and evidence-pipeline revisions. |
| 16 July 2026, 21:08:00 | Targeted follow-up for strong deferred candidates | `1131b67134c8cb7dc0c6a83af3fe8b654187314d` | Committed targeted follow-up results for high-priority deferred candidates. |
| 17 July 2026, 11:43:30 | Complementary pigment-term acquisition added to PR workflow | `f958a582b7dab31ace532ce98817d9513c5c65ad` | Added a second targeted OpenAlex pass using pigment terminology and retained its works and QC as workflow artifacts. |
| 17 July 2026, 12:26:26 | Ecological follow-up search code committed | `fbfafbc3e9dd3832190a565fd39f9b571fa270f8` | Added a third-pass OpenAlex search using morph-frequency, geographic-variation and pollinator-selection terminology. This commit verifies the search implementation, not a separately committed result table. |
| 17 July 2026, 12:52:29 | Classic flower-colour variation follow-up code committed | `6e27d33e1c06249055ad17ce482a6d926a598ec1` | Added another targeted follow-up vocabulary pass. This commit verifies the implementation, not a separately committed result table. |
| 17 July 2026, 13:06:48 | Follow-up evidence merged | `f2ff13874c0954561eb540b8341416ee9f5b10e2` | Merged retained follow-up evidence into the deferred-species ranking. |
| 17 July 2026, 13:07:43 | Full deferred pool resolved from merged evidence | `3040253841f7dd694d1cca331f50530b7a4de778` | Reclassified the deferred pool using the aggregated follow-up evidence. |
| 19 July 2026, 12:40:38 | Ambiguous confirmed-species enrichment committed | `c237e3f64fac130acdc952f8716852f26a6db44f` | Committed targeted evidence enrichment for confirmed species whose spatial-scale evidence remained ambiguous. |

## Search period used in the manuscript

The repository therefore supports the statement that the automated global discovery and recorded targeted follow-up/enrichment activity used for the evidence pipeline occurred from **16 to 19 July 2026**. The initial global corpus was generated on 16 July; deferred-candidate follow-up and aggregation occurred on 16–17 July; and targeted enrichment of ambiguous confirmed cases was committed on 19 July.

## Review and adjudication record

The repository documents operational evidence rules rather than a conventional two-reviewer systematic-review process:

- a baseline case required an explicit source passage supporting the assigned category;
- the passage had to be linked to the species and classification in the source-level audit;
- conflicting evidence required a `mixed` classification rather than forced binary assignment;
- climatic model results could not be used to assign or retain the classification;
- the baseline manifest was frozen before final manuscript models were rerun; and
- unsupported cases remained flagged for manual review and were excluded from confirmed classifications.

**Not verified:** the repository does not record the names or number of human screeners, independent duplicate screening, inter-reviewer agreement, or a formal disagreement-resolution procedure. These details must not be reconstructed retrospectively. If classification was performed by one named author with code-assisted evidence extraction, the title page and Methods should state that only after author confirmation.

## Interpretation boundary

The dates above are provenance dates for the implemented evidence pipeline. They do not turn the literature-derived candidate sample into a prospectively registered systematic review or a random sample of angiosperms.
