from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "finalize_fcp_v2_biological_partition_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_bio_finalizer", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _condition_rows(measurement_id: str, *, heavy: bool) -> pd.DataFrame:
    rows = [
        {
            "measurement_id": measurement_id,
            "condition_id": "baseline",
            "condition_family": "baseline",
            "morph": "white",
            "measurement_status": "classified_four_state_morph",
            "counterfactual_status": "available",
            "image_sha256": "a" * 64,
        }
    ]
    for ev in ("m1_0", "m0_5", "p0_0", "p0_5", "p1_0"):
        rows.append(
            {
                "measurement_id": measurement_id,
                "condition_id": f"fixed_ev_{ev}",
                "condition_family": "fixed_mask_exposure",
                "morph": "white",
                "measurement_status": "classified_four_state_morph",
                "counterfactual_status": "available",
                "image_sha256": "a" * 64,
            }
        )
    if heavy:
        for ev in ("m1_0", "m0_5", "p0_5", "p1_0"):
            rows.append(
                {
                    "measurement_id": measurement_id,
                    "condition_id": f"full_ev_{ev}",
                    "condition_family": "full_pipeline_exposure",
                    "morph": "white",
                    "measurement_status": "classified_four_state_morph",
                    "counterfactual_status": "available",
                    "image_sha256": "a" * 64,
                }
            )
        rows.append(
            {
                "measurement_id": measurement_id,
                "condition_id": "neutral_bg",
                "condition_family": "background_neutralization",
                "morph": "white",
                "measurement_status": "classified_four_state_morph",
                "counterfactual_status": "available",
                "image_sha256": "a" * 64,
            }
        )
        for label in ("base", "x_m5", "x_p5", "y_m5", "y_p5", "scale_m10", "scale_p10"):
            rows.append(
                {
                    "measurement_id": measurement_id,
                    "condition_id": f"jitter_{label}",
                    "condition_family": "roi_prompt_jitter",
                    "morph": "white",
                    "measurement_status": "classified_four_state_morph",
                    "counterfactual_status": "available",
                    "image_sha256": "a" * 64,
                }
            )
    return pd.DataFrame(rows)


def test_expected_condition_count_contract() -> None:
    m = load_module()
    assert m.BASE_CONDITIONS == 6
    assert m.HEAVY_EXTRA_CONDITIONS == 12


def test_source_drift_must_not_have_biological_rows(tmp_path: Path) -> None:
    m = load_module()
    expected = pd.DataFrame(
        {
            "measurement_id": ["x"],
            "heavy_counterfactual": [False],
            "expected_source_sha256": ["a" * 64],
        }
    )
    receipt = pd.DataFrame(
        {
            "measurement_id": ["x"],
            "biological_acquisition_status": ["source_byte_drift"],
            "expected_source_sha256": ["a" * 64],
            "observed_source_sha256": ["b" * 64],
            "heavy_counterfactual": [False],
            "failure_reason": ["drift"],
        }
    )
    bio = _condition_rows("x", heavy=False)
    e=tmp_path/"e.csv"; r=tmp_path/"r.csv"; b=tmp_path/"b.csv"; o=tmp_path/"o.csv"
    expected.to_csv(e,index=False); receipt.to_csv(r,index=False); bio.to_csv(b,index=False)
    import sys
    old=sys.argv
    sys.argv=["prog","--expected-worker-manifest",str(e),"--acquisition-receipt",str(r),"--biological-long",str(b),"--output-terminal-receipt",str(o)]
    try:
        with pytest.raises(RuntimeError, match="failed/drifted"):
            m.main()
    finally:
        sys.argv=old


def test_complete_heavy_row_has_18_conditions(tmp_path: Path) -> None:
    m = load_module()
    expected = pd.DataFrame(
        {
            "measurement_id": ["x"],
            "heavy_counterfactual": [True],
            "expected_source_sha256": ["a" * 64],
        }
    )
    receipt = pd.DataFrame(
        {
            "measurement_id": ["x"],
            "biological_acquisition_status": ["acquired_sha_match"],
            "expected_source_sha256": ["a" * 64],
            "observed_source_sha256": ["a" * 64],
            "heavy_counterfactual": [True],
            "failure_reason": [""],
        }
    )
    bio = _condition_rows("x", heavy=True)
    assert len(bio)==18
    e=tmp_path/"e.csv"; r=tmp_path/"r.csv"; b=tmp_path/"b.csv"; o=tmp_path/"o.csv"
    expected.to_csv(e,index=False); receipt.to_csv(r,index=False); bio.to_csv(b,index=False)
    import sys
    old=sys.argv
    sys.argv=["prog","--expected-worker-manifest",str(e),"--acquisition-receipt",str(r),"--biological-long",str(b),"--output-terminal-receipt",str(o)]
    try:
        m.main()
    finally:
        sys.argv=old
    out=pd.read_csv(o)
    assert len(out)==1
    assert out.loc[0,"pass_b_terminal_status"]=="biological_measurement_complete"
    assert int(out.loc[0,"condition_rows"])==18
