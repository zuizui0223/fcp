import hashlib
import json

from scripts.data.audit_inaturalist_automated_colour_cache import finite_number


def test_finite_number_rejects_null_and_nonfinite():
    assert finite_number(1.5) is True
    assert finite_number(None) is False
    assert finite_number(float("nan")) is False


def test_cache_json_contract_can_be_strict_standard_json(tmp_path):
    record = {"background_L_mean": None, "background_features_available": False}
    path = tmp_path / "record.json"
    path.write_text(json.dumps(record, allow_nan=False), encoding="utf-8")
    assert hashlib.sha256(path.read_bytes()).hexdigest()
    assert json.loads(path.read_text(encoding="utf-8"))["background_L_mean"] is None
