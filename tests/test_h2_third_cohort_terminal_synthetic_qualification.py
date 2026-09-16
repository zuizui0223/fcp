import hashlib
import importlib.util
import json
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analysis/qualify_polymorphism_h2_terminal_transport_20260916.py"
PRODUCTION_WRITER = (
    "scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py::write_result"
)


def _load_module():
    assert SCRIPT.exists(), "synthetic terminal qualifier has not been implemented"
    spec = importlib.util.spec_from_file_location("h2_terminal_qualifier", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def test_synthetic_h2_complete_round_trips_through_production_writer_and_zip(tmp_path):
    module = _load_module()
    package_receipt = module.qualify(tmp_path)

    result_path = tmp_path / "result.json"
    qualification_path = tmp_path / "qualification_receipt.json"
    zip_path = tmp_path / "terminal_package.zip"
    package_path = tmp_path / "package_receipt.json"
    for path in (result_path, qualification_path, zip_path, package_path):
        assert path.is_file(), path

    result = json.loads(result_path.read_text(encoding="utf-8"))
    qualification = json.loads(qualification_path.read_text(encoding="utf-8"))
    persisted_package = json.loads(package_path.read_text(encoding="utf-8"))

    assert result["stage"] == "H2_COMPLETE"
    assert result["synthetic_only"] is True
    assert result["biological_inputs_opened"] is False
    assert result["stage_receipt"]["stage"] == "H2_COMPLETE"
    assert qualification["verdict"] == "PASS"
    assert qualification["production_writer"] == PRODUCTION_WRITER
    assert qualification["result_sha256"] == _sha256(result_path)
    assert persisted_package == package_receipt
    assert package_receipt["package_sha256"] == _sha256(zip_path)

    with ZipFile(zip_path, "r") as archive:
        assert sorted(archive.namelist()) == ["qualification_receipt.json", "result.json"]
        archived_result = json.loads(archive.read("result.json").decode("utf-8"))
        archived_qualification = json.loads(
            archive.read("qualification_receipt.json").decode("utf-8")
        )
    assert archived_result == result
    assert archived_qualification == qualification
