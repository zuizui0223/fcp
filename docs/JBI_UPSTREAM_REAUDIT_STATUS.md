# JBI upstream re-audit status

## Current stage

The historical 34-species freeze remains unchanged. A separate upstream re-audit has started from the archived full systematic-search corpus.

## Established facts

- archived systematic search: 105,249 raw records and 79,242 deduplicated records;
- the archived search manifest reports 19 truncated query/database shards;
- final natural eligibility and spatial classification were never human-adjudicated in that corpus;
- the historical exploratory candidate builder dropped species with both within- and among-population evidence before ecological modelling;
- mixed evidence is therefore retained explicitly in the new review representation;
- taxon validation now precedes species-level state aggregation;
- climatic metrics and historical model results are not inputs to the classification workflow.

## Active output target

The first deliverable is a GBIF-validated species review queue with separate evidence axes for:

1. local coexistence of discrete natural floral colour variants;
2. geographic or among-population colour structuring.

No new manuscript result should be inferred until source-level duplicate review and adjudication are complete.

## Remaining gates

1. run the mixed-preserving review-queue workflow and inspect its diagnostics;
2. repair or explicitly bound the 19 truncated search shards;
3. complete independent source-level review and adjudication;
4. only then create a replacement statistical freeze and rerun ecological analyses.
