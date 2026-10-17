#!/usr/bin/env bash
# Restore historical 149,900-photo classifier EXACTLY, not new model code.
set -euo pipefail
SOURCE="9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
THIRD="7e538e5c51c05a7cc47b2fcf53eea92634c8a863"
git fetch --no-tags --depth=1 origin "$SOURCE" "$THIRD"
mkdir -p scripts/analysis fcp_pipeline
git show "$SOURCE:fcp_pipeline/photo_first_measurement.py" > fcp_pipeline/photo_first_measurement.py
git show "$SOURCE:scripts/analysis/materialize_random_photo_first_h9_measurement_infrastructure.sh" > scripts/analysis/materialize_random_photo_first_h9_measurement_infrastructure.sh
bash scripts/analysis/materialize_random_photo_first_h9_measurement_infrastructure.sh
for f in \
    fcp_pipeline/third_cohort_prospective_execution_gate.py \
    scripts/analysis/run_polymorphism_h1_observer_disjoint_reliability_20260913.py \
    scripts/analysis/run_polymorphism_delta_geometry_validation_20260912.py \
    scripts/analysis/run_polymorphism_delta_geometry_structured_null_20260912.py \
    scripts/analysis/run_polymorphism_h2_third_cohort_prospective_white_axis_20260917.py
do
    git show "$THIRD:$f" > "$f"
done
test "$(git hash-object fcp_pipeline/photo_first_measurement.py)" = "b2761c7fd2f8af615c1d96a615942e7af67d492a"
echo FCP_WAVE4_ORIGINAL_IMAGE_AND_STATISTICAL_METHODS_RESTORED
