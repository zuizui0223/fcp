"""Pre-pixel identity gate and durable single-execution ledger. No image API."""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import types
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[2]
AMENDMENT = "docs/supporting/rgfca_monarda_execution_amendment_v1.json"
AUTHORIZATION = "docs/supporting/rgfca_monarda_execution_authorization_v1.json"
PARENT = "docs/supporting/rgfca_monarda_region_agreement_contract_v1.json"
RUNTIME_COMMIT = "9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
RUNTIME_BLOBS = {
    "fcp_pipeline/flower_roi_v4.py": "9d1d2a1848f0c798b53c4df51543dc2682342377",
    "fcp_pipeline/flower_roi_v4_runtime.py": "6f3e71f7216b1635a2e06628336f5ea9ce0105d6",
    "fcp_pipeline/photo_first_measurement.py": "b2761c7fd2f8af615c1d96a615942e7af67d492a",
    "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json": "b043bc83d0d69489f38771e8f4bbe524963b70ea",
    "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json": "b25425f1f738bcea75cef5f60980332a9857610b",
}
QUALIFIED_PATHS = (
    "scripts/analysis/run_rgfca_monarda_region_agreement.py",
    "scripts/analysis/monarda_execution_guard.py",
    "scripts/analysis/audit_rgfca_monarda_archive.py",
    "tests/test_rgfca_monarda_region_agreement.py",
    "tests/test_monarda_execution_guard.py",
    ".github/workflows/rgfca-monarda-region-agreement.yml",
    AMENDMENT, PARENT,
    "docs/supporting/rgfca_monarda_archive_intake_v1.json",
    "docs/supporting/rgfca_monarda_source_mapping_v1.json",
)


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def blob_id(raw: bytes, *, text: bool = False) -> str:
    if text:
        raw = raw.replace(b"\r\n", b"\n")
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def git_bytes(revision: str, relative: str, root: Path = ROOT) -> bytes:
    return subprocess.check_output(["git", "show", f"{revision}:{relative}"], cwd=root)


def environment_snapshot() -> dict:
    return {
        "python": platform.python_version(), "system": platform.system(),
        "machine": platform.machine(), "platform": platform.platform(),
        "python_build": sys.version,
        "packages": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()
                     if d.metadata.get("Name")},
    }


def check_environment(expected: dict, observed: dict) -> None:
    for key in ("python", "system", "machine"):
        if observed[key] != expected[key]:
            raise RuntimeError(f"execution environment mismatch: {key}")
    packages = {k.lower().replace("_", "-"): v for k, v in observed["packages"].items()}
    for name, version in expected["packages"].items():
        if packages.get(name.lower().replace("_", "-")) != version:
            raise RuntimeError(f"execution dependency mismatch: {name}")


def atomic_json(path: Path, value: dict | list) -> None:
    temporary = path.with_suffix(path.suffix + ".pending")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def verify_authorization(root: Path = ROOT) -> tuple[dict, dict, dict]:
    """Bind a committed authorization to committed qualified source, before decode."""
    path = root / AUTHORIZATION
    if not path.is_file():
        raise RuntimeError("No committed Monarda execution authorization; pixels remain closed")
    raw = path.read_bytes()
    if raw != git_bytes("HEAD", AUTHORIZATION, root):
        raise RuntimeError("execution authorization is not the exact committed file")
    auth = json.loads(raw)
    if auth.get("status") != "authorized_once_after_pre_pixel_qualification":
        raise RuntimeError("execution is not authorized")
    qualification = auth["qualification"]
    if qualification.get("conclusion") != "success" or not qualification.get("run_id"):
        raise RuntimeError("missing successful qualification receipt")
    revision = qualification["head_sha"]
    subprocess.run(["git", "merge-base", "--is-ancestor", revision, "HEAD"],
                   cwd=root, check=True, capture_output=True)
    if set(auth["qualified_source_blobs"]) != set(QUALIFIED_PATHS):
        raise RuntimeError("qualified source manifest is incomplete")
    identities = {}
    for relative, expected in auth["qualified_source_blobs"].items():
        frozen = git_bytes(revision, relative, root)
        observed = (root / relative).read_bytes()
        if blob_id(frozen) != expected or blob_id(observed, text=True) != expected:
            raise RuntimeError(f"qualified source identity mismatch: {relative}")
        identities[relative] = {"git_blob": expected, "canonical_sha256": sha_bytes(frozen)}
    amendment = json.loads((root / AMENDMENT).read_text(encoding="utf-8"))
    if blob_id((root / PARENT).read_bytes(), text=True) != amendment["parent_contract_blob"]:
        raise RuntimeError("original contract was altered")
    observed_environment = environment_snapshot()
    check_environment(amendment["execution_environment"], observed_environment)
    return auth, amendment, {
        "authorization_sha256": sha_bytes(raw), "qualified_sources": identities,
        "environment": observed_environment, "qualification": qualification,
        "monarda_pixels_decoded": False, "model_loaded": False,
    }


def materialize_runtime(root: Path = ROOT) -> tuple[Path, dict]:
    """Copy only verified historical blobs into an isolated cache, never repo code."""
    cache = root / ".artifacts/monarda-runtime-v1/immutable-source"
    receipts = {}
    for relative, expected in RUNTIME_BLOBS.items():
        raw = git_bytes(RUNTIME_COMMIT, relative, root)
        if blob_id(raw) != expected:
            raise RuntimeError(f"immutable runtime source mismatch: {relative}")
        destination = cache / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if destination.read_bytes() != raw:
                raise RuntimeError(f"existing immutable cache differs: {relative}")
        else:
            with destination.open("xb") as stream:
                stream.write(raw)
        receipts[relative] = {"git_blob": expected, "sha256": sha_bytes(raw)}
    return cache, receipts


def import_runtime(cache: Path):
    """Import relative runtime dependencies from the verified isolated namespace."""
    name = "_monarda_frozen_roi_v4"
    if any(k == name or k.startswith(name + ".") for k in sys.modules):
        raise RuntimeError("isolated runtime already imported in this process")
    package = types.ModuleType(name)
    package.__path__ = [str(cache / "fcp_pipeline")]
    sys.modules[name] = package
    for module in ("flower_roi_v4", "flower_roi_v4_runtime"):
        source = cache / f"fcp_pipeline/{module}.py"
        expected = RUNTIME_BLOBS[f"fcp_pipeline/{module}.py"]
        if blob_id(source.read_bytes()) != expected:
            raise RuntimeError("runtime cache changed before import")
        spec = importlib.util.spec_from_file_location(f"{name}.{module}", source)
        loaded = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = loaded
        spec.loader.exec_module(loaded)
    return loaded


class ExecutionLedger:
    """One canonical run; no automatic resume or overwrite entry point."""

    def __init__(self, output: Path, census: list[dict], receipt: dict):
        keys = [(r["split"], r["coco_image_id"]) for r in census]
        if len(keys) != 110 or len(set(keys)) != 110:
            raise RuntimeError("execution start requires 110 unique census rows")
        # mkdir is the exclusive run claim; even an empty directory is a hold.
        output.mkdir(parents=True, exist_ok=False)
        self.output = output
        self.rows = {key: {
            "split": row["split"], "coco_image_id": row["coco_image_id"],
            "member_path": row["member_path"], "annotation_count": row["annotation_count"],
            "alignment_status": "not_started", "scoring_status": "not_started",
        } for key, row in zip(keys, census)}
        atomic_json(output / "execution_start.json", {
            **receipt, "started_utc": datetime.now(timezone.utc).isoformat(),
            "all_images": 110, "status": "execution_claimed_no_pixels_yet",
        })
        self.checkpoint()

    def checkpoint(self):
        atomic_json(self.output / "all_110_checkpoint.json", list(self.rows.values()))

    def event(self, phase: str, key=None):
        record = {"utc": datetime.now(timezone.utc).isoformat(), "phase": phase,
                  "image_key": list(key) if key is not None else None}
        with (self.output / "events.jsonl").open("a", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(record) + "\n")
            stream.flush()
            os.fsync(stream.fileno())

    def update(self, phase: str, row: dict):
        key = (row["split"], int(row["coco_image_id"]))
        if key not in self.rows:
            raise RuntimeError("image outside fixed ledger")
        if phase == "score" and self.rows[key]["scoring_status"] != "not_started":
            raise RuntimeError("duplicate image score")
        self.rows[key].update(row)
        if phase == "score":
            self.rows[key]["scoring_status"] = "terminal"
        self.checkpoint()
        self.event(phase + "_terminal", key)

    def fail(self, reason: str):
        self.checkpoint()
        atomic_json(self.output / "execution_terminal.json", {
            "status": "not_evaluable_execution_incomplete", "reason": reason,
            "all_110_images_accounted_for": len(self.rows) == 110,
            "terminal_scoring_rows": sum(r["scoring_status"] == "terminal" for r in self.rows.values()),
            "ecological_results_changed": False, "automatic_restart_permitted": False,
        })

    def complete(self, result_path: Path):
        atomic_json(self.output / "execution_terminal.json", {
            "status": "terminal_result_retained",
            "result_sha256": sha_bytes(result_path.read_bytes()),
            "checkpoint_sha256": sha_bytes((self.output / "all_110_checkpoint.json").read_bytes()),
            "events_sha256": sha_bytes((self.output / "events.jsonl").read_bytes()),
            "automatic_restart_permitted": False,
        })
