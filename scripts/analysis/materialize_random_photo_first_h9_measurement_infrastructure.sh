#!/usr/bin/env bash
set -euo pipefail

SOURCE_REPO="https://raw.githubusercontent.com/zuizui0223/fcp"
SOURCE_COMMIT="9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
DETECTOR_SHA256="f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5"

fetch() {
  local path="$1"
  mkdir -p "$(dirname "$path")"
  curl --fail --location --retry 4 --retry-delay 2 \
    "${SOURCE_REPO}/${SOURCE_COMMIT}/${path}" \
    --output "$path"
}

fetch fcp_pipeline/flower_roi_v4.py
fetch fcp_pipeline/flower_roi_v4_runtime.py
fetch fcp_pipeline/photo_first_measurement_execution.py
fetch docs/supporting/random_photo_first_measurement_contract_v1.json
fetch docs/supporting/random_photo_first_measurement_execution_v1.json
fetch docs/supporting/jbi_atlas_roi_estimator_contract_v4.json
fetch data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json
fetch requirements-atlas-roi-v4.txt
fetch scripts/acquisition/acquire_random_photo_first_measurement_partition.py
fetch scripts/analysis/measure_random_photo_first_blind_partition.py
fetch scripts/analysis/seal_random_photo_first_measurement_partition.py
fetch data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt

# Git-blob identities were independently frozen by the completed measurement pipeline.
test "$(git hash-object fcp_pipeline/flower_roi_v4.py)" = "9d1d2a1848f0c798b53c4df51543dc2682342377"
test "$(git hash-object fcp_pipeline/flower_roi_v4_runtime.py)" = "6f3e71f7216b1635a2e06628336f5ea9ce0105d6"
test "$(git hash-object fcp_pipeline/photo_first_measurement.py)" = "b2761c7fd2f8af615c1d96a615942e7af67d492a"
test "$(git hash-object docs/supporting/jbi_atlas_roi_estimator_contract_v4.json)" = "b043bc83d0d69489f38771e8f4bbe524963b70ea"
test "$(git hash-object data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json)" = "b25425f1f738bcea75cef5f60980332a9857610b"
test "$(git hash-object docs/supporting/random_photo_first_measurement_contract_v1.json)" = "27b5933025d33d5e91e957fdf892b7409c605ef2"
test "$(git hash-object docs/supporting/random_photo_first_measurement_execution_v1.json)" = "9ac42ddd6175d971a076961ff18a96bcf1355a5b"
echo "${DETECTOR_SHA256}  data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt" | sha256sum -c -

echo "h9_measurement_infrastructure_materialized=${SOURCE_COMMIT}"
