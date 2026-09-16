#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

from fcp_pipeline.p500_prospective_execution_gate import (
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    validate_stage_transition,
)


ROOT = Path(__file__).resolve().parents[2]
PRODUCTION_SCRIPT = (
    ROOT / "scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py"
)
PRODUCTION_WRITER = (
    "scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py::write_result"
)
OUTPUT_NAMES = (
    "result.json",
    "qualification_receipt.json",
    "terminal_package.zip",
    "package_receipt.json",
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_production_writer():
    spec = importlib.util.spec_from_file_location("p500_production_h2_writer", PRODUCTION_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load production H2 script: {PRODUCTION_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    writer = getattr(module, "write_result", None)
    if not callable(writer):
        raise RuntimeError("production H2 script does not expose write_result")
    return writer


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_deterministic_zip(zip_path: Path, members: list[tuple[str, Path]]) -> None:
    with ZipFile(zip_path, "w") as archive:
        for arcname, source in sorted(members):
            info = ZipInfo(arcname, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, source.read_bytes())


def qualify(output_dir: Path) -> dict:
    """Exercise only the already-existing terminal serializer/transport path.

    No candidate table, measured photo table, colour vector, D/W value, or actual
    third-cohort outcome is read by this function. The P500 stage constants are
    used only as the existing validated H2_COMPLETE schema profile; this is a
    transport qualification, not a biological re-analysis or cohort-size claim.
    """

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name in OUTPUT_NAMES:
        path = out / name
        if path.exists():
            path.unlink()

    support_stage = {
        "stage": "SUPPORT_GATE_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "measurement_evaluable_species": 300,
        "support_decision": "PASS",
    }
    h2_stage = {
        "stage": "H2_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "primary_h2_vector_species": 40,
        "primary_structured_null_p": 0.5,
        "h2_decision": "NOT_CONFIRMED",
    }
    validate_stage_transition(support_stage, h2_stage)

    result = {
        "analysis": "polymorphism_h2_terminal_transport_synthetic_qualification_20260916",
        "status": "synthetic_only_terminal_transport_qualification",
        "stage": "H2_COMPLETE",
        "synthetic_only": True,
        "biological_inputs_opened": False,
        "third_cohort_biological_outcome_opened": False,
        "scientific_claims_authorized": False,
        "schema_profile": "existing_p500_h2_complete_terminal_schema_only",
        "qualification_scope": [
            "production_write_result_serialization",
            "filesystem_persistence",
            "independent_json_readback",
            "terminal_stage_validation",
            "deterministic_zip_packaging",
        ],
        "synthetic_support_stage": support_stage,
        "stage_receipt": h2_stage,
    }

    writer = _load_production_writer()
    writer(out, result)

    result_path = out / "result.json"
    persisted = json.loads(result_path.read_text(encoding="utf-8"))
    if persisted != result:
        raise RuntimeError("production writer disk read-back differs from in-memory synthetic result")
    if persisted.get("synthetic_only") is not True or persisted.get("biological_inputs_opened") is not False:
        raise RuntimeError("synthetic outcome firewall changed during serialization")
    validate_stage_transition(persisted["synthetic_support_stage"], persisted["stage_receipt"])

    qualification_receipt = {
        "qualification": "h2_third_cohort_terminal_transport_preopening_gate",
        "verdict": "PASS",
        "synthetic_only": True,
        "biological_inputs_opened": False,
        "third_cohort_biological_outcome_opened": False,
        "scientific_claims_authorized": False,
        "production_writer": PRODUCTION_WRITER,
        "schema_profile": "existing_p500_h2_complete_terminal_schema_only",
        "schema_profile_note": (
            "P500 counts are synthetic schema fixtures only; this gate makes no claim "
            "about the third-cohort denominator or biological result."
        ),
        "result_sha256": _sha256(result_path),
    }
    qualification_path = out / "qualification_receipt.json"
    _write_json(qualification_path, qualification_receipt)

    zip_path = out / "terminal_package.zip"
    _write_deterministic_zip(
        zip_path,
        [
            ("result.json", result_path),
            ("qualification_receipt.json", qualification_path),
        ],
    )

    with ZipFile(zip_path, "r") as archive:
        names = sorted(archive.namelist())
        if names != ["qualification_receipt.json", "result.json"]:
            raise RuntimeError(f"unexpected terminal package members: {names}")
        archived_result = json.loads(archive.read("result.json").decode("utf-8"))
        archived_receipt = json.loads(archive.read("qualification_receipt.json").decode("utf-8"))
    if archived_result != persisted:
        raise RuntimeError("ZIP read-back result differs from persisted result")
    if archived_receipt != qualification_receipt:
        raise RuntimeError("ZIP read-back qualification receipt differs from persisted receipt")
    validate_stage_transition(archived_result["synthetic_support_stage"], archived_result["stage_receipt"])

    package_receipt = {
        "qualification": "h2_third_cohort_terminal_transport_package",
        "verdict": "PASS",
        "synthetic_only": True,
        "biological_inputs_opened": False,
        "scientific_claims_authorized": False,
        "members": names,
        "package_sha256": _sha256(zip_path),
        "package_bytes": zip_path.stat().st_size,
        "result_sha256": qualification_receipt["result_sha256"],
    }
    _write_json(out / "package_receipt.json", package_receipt)
    return package_receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    receipt = qualify(args.output_dir)
    print(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
