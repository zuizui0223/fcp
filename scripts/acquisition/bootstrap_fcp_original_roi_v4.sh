#!/usr/bin/env bash
# Restore the *method*, never historical measured-photo rows, into the FCP checkout.
set -euo pipefail
usage() { echo "Usage: $0 --code-only | --full [--sam-dir PATH]" >&2; exit 2; }
[[ $# -ge 1 ]] || usage
mode="$1"; shift
[[ "$mode" == '--code-only' || "$mode" == '--full' ]] || usage
sam_dir="${RUNNER_TEMP:-/tmp}/fcp-fresh-original-efficient-sam"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sam-dir) [[ $# -ge 2 ]] || usage; sam_dir="$2"; shift 2 ;;
    *) usage ;;
  esac
done
root="$(git rev-parse --show-toplevel)"
cd "$root"
source_commit='9fae6ccdf684a46026f72ba12e98de2c5c54bf2a'
# Exact historical blobs used by the 1,499-species analysis.
files=(
 'fcp_pipeline/flower_roi_v4.py|9d1d2a1848f0c798b53c4df51543dc2682342377'
 'fcp_pipeline/flower_roi_v4_runtime.py|6f3e71f7216b1635a2e06628336f5ea9ce0105d6'
 'fcp_pipeline/photo_first_measurement.py|b2761c7fd2f8af615c1d96a615942e7af67d492a'
 'fcp_pipeline/photo_first_atlas_v2.py|d2984b5a462bf746ea5646332b74be55411b78f8'
 'fcp_pipeline/photo_first_atlas.py|6dced732d315f166565a0dfaf3507c90f967255e'
 'fcp_pipeline/shared_transition_surface.py|3d69648c6649270b196d8aa46f6420cb5b9d9e4a'
 'fcp_pipeline/photo_first_measurement_execution.py|b6e76efb4db8e7529b5ce5a91c60b6ae68637fb5'
 'scripts/analysis/measure_random_photo_first_blind_partition.py|b1d7091db44682559bbe4b48a9f3b8db648501a6'
 'docs/supporting/random_photo_first_measurement_contract_v1.json|27b5933025d33d5e91e957fdf892b7409c605ef2'
 'docs/supporting/random_photo_first_measurement_execution_v1.json|9ac42ddd6175d971a076961ff18a96bcf1355a5b'
 'docs/supporting/jbi_atlas_roi_estimator_contract_v4.json|b043bc83d0d69489f38771e8f4bbe524963b70ea'
 'data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json|b25425f1f738bcea75cef5f60980332a9857610b'
 'requirements-atlas-roi-v4.txt|197e9227bbaee2f3b83edc8c857b636b0a1ca35d'
)
# Read-only historical git history; source-composition and fingerprints are non-negotiable.
git fetch --no-tags --depth=1 origin "$source_commit" >&2
for element in "${files[@]}"; do
  path="${element%%|*}"; digest="${element##*|}"
  mkdir -p "$(dirname "$path")"
  if [[ -e "$path" ]]; then
    [[ "$(git hash-object "$path")" == "$digest" ]] || {
      echo "Refusing to overwrite nonidentical existing source $path" >&2; exit 1;
    }
  else
    git show "$source_commit:$path" > "$path"
  fi
  [[ "$(git hash-object "$path")" == "$digest" ]] || {
    echo "Frozen ROI-v4 method Git blob mismatch: $path" >&2; exit 1;
  }
done
if [[ "$mode" == '--full' ]]; then
  detector='data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt'
  mkdir -p "$(dirname "$detector")" "$sam_dir"
  if [[ ! -f "$detector" ]]; then
    git show "$source_commit:$detector" > "$detector"
  fi
  echo 'f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5  '"$detector" | sha256sum -c -
  revision='d525f622e6f640acf5a0fc37c7ca1f243da5bde0'
  encoder="$sam_dir/efficient_sam_vitt_encoder.onnx"
  decoder="$sam_dir/efficient_sam_vitt_decoder.onnx"
  if [[ ! -f "$encoder" ]]; then
    curl --fail --location --retry 3 --retry-delay 3 \
      "https://raw.githubusercontent.com/yformer/EfficientSAM/$revision/weights/efficient_sam_vitt_encoder.onnx" \
      --output "$encoder"
  fi
  if [[ ! -f "$decoder" ]]; then
    curl --fail --location --retry 3 --retry-delay 3 \
      "https://raw.githubusercontent.com/yformer/EfficientSAM/$revision/weights/efficient_sam_vitt_decoder.onnx" \
      --output "$decoder"
  fi
  echo '84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951  '"$encoder" | sha256sum -c -
  echo 'a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11  '"$decoder" | sha256sum -c -
fi
python - <<'PY'
import json,hashlib
from pathlib import Path
p=Path('scripts/analysis/measure_random_photo_first_blind_partition.py')
b=p.read_bytes()
assert hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()=='b1d7091db44682559bbe4b48a9f3b8db648501a6'
print(json.dumps({'status':'FROZEN_ROI_V4_CODE_VERIFIED','method_replaced':False,'historical_measurements_read':False}))
PY
[[ "$mode" != '--code-only' ]] || echo 'Code-only verification passed; original detector/EfficientSAM model weights are NOT staged.'
