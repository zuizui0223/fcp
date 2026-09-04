# Random photo-first inference execution handoff

This branch was frozen before the complete fresh measurement result existed. It does not contain an H1 or H2 biological outcome at freeze time.

The post-measurement handoff is deliberately mechanical:

1. the measurement branch must first complete and seal all 128 partitions;
2. reassembly must return exactly the frozen 20,845 candidate rows and commit the measured table plus measurement-result manifest;
3. the sync workflow imports those exact Git blobs into this branch and verifies their hashes and denominator;
4. that two-file import triggers the frozen H1 then hierarchical H2 workflow;
5. H1 uses 10,000 colour-blind sampled photographs per replicate, 200 replicates, cap 2 per species×cell, and 999 within-species null permutations with seeds 20260903 / 20260904;
6. H2 reuses the matched 999 H1 edge-persistence null maps and cannot rescue a non-significant or not-evaluable H1.

The JSON-safe launcher is a technical serializer only. It maps non-finite not-evaluable/QC scalar placeholders to JSON null and does not alter sampling, morph states, transition edges, null permutations, climate variables, statistics, or decisions.
