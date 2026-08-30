# JBI Chapter 1 resume marker — 2026-08-30

The frozen 720-photograph evaluation set was explicitly opened after the calibration-stage feature representation and spatial null model were frozen.

Active evaluation workflow: `JBI Chapter 1 Florence evaluation features`, run `33281907575`.

The workflow contains one preflight job and 36 deterministic 20-photograph extraction shards (6 species × 6 shards), with no review-dependent gate. Each shard must retain `evaluation_row=true`, `calibration_only=false`, `final_label=false`, and `evaluation_feature_measurement=true`.

After extraction, the required next gate is an automated harvest that verifies exactly 720 unique evaluation photograph IDs, exactly 120 per species, zero calibration-ID overlap, and exact equality with the frozen evaluation subset before Stage A or Stage B is interpreted.

Four provisional QA figures (C1–C4) have been rendered for Stage A global/species results and Stage B surface/sensitivity. They remain provisional until their numerical source manifest and the exact 720-row evaluation harvest are verified together.
