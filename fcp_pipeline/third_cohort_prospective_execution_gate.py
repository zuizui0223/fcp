"""Forward-only execution contract for the third-cohort prospective H2 test."""
from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
from typing import Any

import pandas as pd

METADATA_FREEZE_COMMIT = "ca69930986e39e3ea9b2d2f97ab247045ba9db0b"
AUTHORIZED_METADATA_SHA256 = "13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290"
AUTHORIZED_SPECIES_SHA256 = "36a866b040b835e20539d318b0533bcb905cbe760d23df29e7da6587a054e593"
MEASUREMENT_SOURCE_COMMIT = "9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
BLIND_SALT = "FCP_H2_THIRD_COHORT_PROSPECTIVE_20260917_V1"
EXPECTED_SPECIES = 499
EXPECTED_ROWS = 49_900
EXPECTED_TERMINAL_PARTITIONS = 256
MIN_MEASUREMENT_EVALUABLE_SPECIES = 250
MIN_H2_VECTOR_SPECIES = 20

STAGES = (
    "PREOPENING",
    "FIREWALL_FROZEN",
    "PARTITION_MEASUREMENT_COMPLETE",
    "REASSEMBLY_COMPLETE",
    "SUPPORT_GATE_COMPLETE",
    "H2_COMPLETE",
)


def _require_bool(value: Any, label: str) -> None:
    if type(value) is not bool:
        raise RuntimeError(f"{label} must be boolean")


def measurement_id(photo_id: object) -> str:
    payload = f"{BLIND_SALT}\x1fphoto\x1f{photo_id}".encode("utf-8")
    return "FCPH2T3-" + hashlib.sha256(payload).hexdigest().upper()[:24]


def measurement_batch(mid: str) -> int:
    digest = hashlib.sha256(f"fcp-h2-third-batch\x1f{mid}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % 2


def validate_authorized_denominator(
    metadata: pd.DataFrame,
    authorized_species: pd.DataFrame,
    *,
    expected_species: int = EXPECTED_SPECIES,
    target_per_species: int = 100,
) -> pd.DataFrame:
    if expected_species < 1 or target_per_species < 1:
        raise RuntimeError("expected denominator must be positive")
    expected_rows = int(expected_species) * int(target_per_species)
    required = {
        "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
        "photo_license", "latitude", "longitude",
    }
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise RuntimeError(f"authorized metadata lacks measurement fields: {missing}")
    sp_col = "query_species" if "query_species" in metadata.columns else "species"
    if sp_col not in metadata.columns:
        raise RuntimeError("authorized metadata lacks species identity")
    if "species" not in authorized_species.columns or "inat_taxon_id" not in authorized_species.columns:
        raise RuntimeError("authorized species table lacks identity fields")

    x = metadata.copy()
    x["species"] = x[sp_col].astype(str)
    x["inat_taxon_id"] = pd.to_numeric(x["inat_taxon_id"], errors="raise").astype("int64")
    sp = authorized_species[["species", "inat_taxon_id"]].copy()
    sp["species"] = sp["species"].astype(str)
    sp["inat_taxon_id"] = pd.to_numeric(sp["inat_taxon_id"], errors="raise").astype("int64")

    if len(x) != expected_rows:
        raise RuntimeError(f"authorized row denominator mismatch: {len(x)} != {expected_rows}")
    if len(sp) != expected_species or sp["species"].nunique() != expected_species or sp["inat_taxon_id"].nunique() != expected_species:
        raise RuntimeError("authorized species denominator or uniqueness mismatch")
    if x["species"].nunique() != expected_species or x["inat_taxon_id"].nunique() != expected_species:
        raise RuntimeError("metadata species denominator mismatch")
    metadata_pairs = set(zip(x["species"], x["inat_taxon_id"], strict=False))
    species_pairs = set(zip(sp["species"], sp["inat_taxon_id"], strict=False))
    if metadata_pairs != species_pairs:
        raise RuntimeError("metadata species set differs from authorized species set")
    counts = x.groupby(["species", "inat_taxon_id"], observed=True).size()
    if not (counts.astype(int) == int(target_per_species)).all():
        raise RuntimeError("not every authorized species has the exact frozen row target")
    if x["observation_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate observation IDs in authorized denominator")
    if x["photo_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate photo IDs in authorized denominator")
    return x.reset_index(drop=True)


def validate_candidate_metadata(candidate: Mapping[str, Any]) -> None:
    if candidate.get("analysis") != "polymorphism_h2_third_cohort_candidate_metadata_acquisition":
        raise RuntimeError("unexpected third-cohort candidate metadata analysis")
    if candidate.get("status") != "complete_fresh_metadata_draw_before_pixels":
        raise RuntimeError("third-cohort candidate metadata was not frozen before pixels")
    decision = candidate.get("decision", {})
    if decision.get("verdict") != "THIRD_COHORT_METADATA_GATE_PASS":
        raise RuntimeError("third-cohort candidate metadata gate did not pass")
    if decision.get("pixels_may_be_authorized_in_separate_next_step") is not True:
        raise RuntimeError("metadata result does not permit a separate later pixel step")
    if decision.get("biological_pixel_opening_authorized_by_this_gate") is not False:
        raise RuntimeError("metadata gate improperly authorized biological pixel opening")
    if decision.get("no_species_replacement") is not True or decision.get("no_target_relaxation") is not True:
        raise RuntimeError("metadata replacement/relaxation prohibition missing")
    capacity = candidate.get("capacity", {})
    if capacity.get("authorized_species") != EXPECTED_SPECIES:
        raise RuntimeError("authorized species denominator mismatch")
    if capacity.get("authorized_metadata_rows") != EXPECTED_ROWS:
        raise RuntimeError("authorized row denominator mismatch")
    if candidate.get("transport", {}).get("request_errors") != 0:
        raise RuntimeError("metadata receipt contains request errors")
    files = candidate.get("files", {})
    if files.get("authorized_metadata_sha256") != AUTHORIZED_METADATA_SHA256:
        raise RuntimeError("authorized metadata hash mismatch")
    if files.get("authorized_species_sha256") != AUTHORIZED_SPECIES_SHA256:
        raise RuntimeError("authorized species hash mismatch")
    firewall = candidate.get("outcome_firewall", {})
    for key in (
        "image_pixels_opened", "flower_colour_opened", "morph_opened", "palette_opened",
        "D_opened", "H2_vectors_opened", "H2_W_opened", "structured_null_opened",
        "P500_recovered_H2_read",
    ):
        if firewall.get(key) is not False:
            raise RuntimeError(f"metadata outcome firewall not closed: {key}")


def validate_authorization(auth: Mapping[str, Any]) -> None:
    if auth.get("status") != "authorize_exactly_one_third_cohort_location_blind_h2_execution":
        raise RuntimeError("unexpected third-cohort biological authorization status")
    expected = {
        "metadata_freeze_commit": METADATA_FREEZE_COMMIT,
        "authorized_metadata_sha256": AUTHORIZED_METADATA_SHA256,
        "authorized_species_sha256": AUTHORIZED_SPECIES_SHA256,
        "frozen_species": EXPECTED_SPECIES,
        "frozen_rows": EXPECTED_ROWS,
        "target_rows_per_species": 100,
        "measurement_machine_source_commit": MEASUREMENT_SOURCE_COMMIT,
        "H2_target": "fixed_q_white_W_structured_null",
        "primary_threshold": 0.10,
        "strict_sensitivity_threshold": 0.20,
        "structured_null_replicates": 999,
        "primary_seed": 20260915,
        "strict_seed": 20261015,
        "minimum_classifiable_photos_per_species": 40,
        "minimum_measurement_evaluable_species": MIN_MEASUREMENT_EVALUABLE_SPECIES,
        "minimum_primary_H2_vector_species": MIN_H2_VECTOR_SPECIES,
    }
    for key, value in expected.items():
        if auth.get(key) != value:
            raise RuntimeError(f"authorization field mismatch: {key}")
    for key in ("species_replacement_allowed", "row_replacement_allowed", "target_relaxation_allowed", "axis_refit_allowed"):
        if auth.get(key) is not False:
            raise RuntimeError(f"authorization must forbid {key}")


def build_start_record(*, branch: str, head_sha: str, known_output_paths_absent: bool,
                       working_tree_inputs_verified: bool) -> dict[str, Any]:
    _require_bool(known_output_paths_absent, "known_output_paths_absent")
    _require_bool(working_tree_inputs_verified, "working_tree_inputs_verified")
    if not head_sha or len(head_sha) != 40:
        raise RuntimeError("head_sha must be a 40-character commit SHA")
    now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    authorized = known_output_paths_absent and working_tree_inputs_verified
    return {
        "schema": "third_cohort_prospective_execution_start_v1",
        "stage": "PREOPENING",
        "created_at_utc": now,
        "branch": branch,
        "head_sha": head_sha,
        "metadata_freeze_commit": METADATA_FREEZE_COMMIT,
        "authorized_metadata_sha256": AUTHORIZED_METADATA_SHA256,
        "authorized_species_sha256": AUTHORIZED_SPECIES_SHA256,
        "measurement_source_commit": MEASUREMENT_SOURCE_COMMIT,
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "image_pixels_opened_before_start": False,
        "flower_colour_opened_before_start": False,
        "H2_opened_before_start": False,
        "known_output_paths_absent": known_output_paths_absent,
        "working_tree_inputs_verified": working_tree_inputs_verified,
        "opening_authorized": authorized,
        "claim_scope": "forward chronology only; no retrospective certification of unrecorded external access",
    }


def validate_start_record(record: Mapping[str, Any], *, expected_head_sha: str) -> None:
    if record.get("schema") != "third_cohort_prospective_execution_start_v1":
        raise RuntimeError("unexpected third-cohort start-record schema")
    if record.get("stage") != "PREOPENING" or record.get("head_sha") != expected_head_sha:
        raise RuntimeError("third-cohort start record identity mismatch")
    if record.get("metadata_freeze_commit") != METADATA_FREEZE_COMMIT:
        raise RuntimeError("metadata-freeze commit mismatch")
    if record.get("authorized_metadata_sha256") != AUTHORIZED_METADATA_SHA256:
        raise RuntimeError("start-record authorized metadata mismatch")
    if record.get("authorized_species_sha256") != AUTHORIZED_SPECIES_SHA256:
        raise RuntimeError("start-record authorized species mismatch")
    if record.get("measurement_source_commit") != MEASUREMENT_SOURCE_COMMIT:
        raise RuntimeError("start-record measurement source mismatch")
    for key in ("image_pixels_opened_before_start", "flower_colour_opened_before_start", "H2_opened_before_start"):
        if record.get(key) is not False:
            raise RuntimeError(f"preopening outcome flag not false: {key}")
    if record.get("known_output_paths_absent") is not True or record.get("working_tree_inputs_verified") is not True:
        raise RuntimeError("preopening input/output verification failed")
    if record.get("opening_authorized") is not True:
        raise RuntimeError("start record does not authorize opening")
    raw_time = record.get("created_at_utc")
    if not isinstance(raw_time, str) or not raw_time.endswith("Z"):
        raise RuntimeError("invalid start-record timestamp")
    datetime.fromisoformat(raw_time[:-1] + "+00:00")


def validate_stage_transition(previous: Mapping[str, Any], current: Mapping[str, Any]) -> None:
    prev = previous.get("stage")
    cur = current.get("stage")
    if prev not in STAGES or cur not in STAGES:
        raise RuntimeError("unknown third-cohort execution stage")
    if STAGES.index(cur) != STAGES.index(prev) + 1:
        raise RuntimeError(f"illegal third-cohort stage transition: {prev} -> {cur}")
    if current.get("species") != EXPECTED_SPECIES or current.get("rows") != EXPECTED_ROWS:
        raise RuntimeError("third-cohort denominator changed across stages")
    if current.get("replacement_rows", 0) != 0 or current.get("replacement_species", 0) != 0:
        raise RuntimeError("replacement is forbidden")

    if cur == "FIREWALL_FROZEN":
        if current.get("measurement_ids") != EXPECTED_ROWS:
            raise RuntimeError("firewall measurement-ID census mismatch")
        if current.get("terminal_partitions") != EXPECTED_TERMINAL_PARTITIONS:
            raise RuntimeError("firewall partition census mismatch")
        if current.get("candidate_pixels_opened") is not False:
            raise RuntimeError("firewall was not frozen before pixels")
    elif cur == "PARTITION_MEASUREMENT_COMPLETE":
        if current.get("terminal_acquisition_rows") != EXPECTED_ROWS or current.get("terminal_measurement_rows") != EXPECTED_ROWS:
            raise RuntimeError("partition measurement coverage is incomplete")
        if current.get("terminal_partitions") != EXPECTED_TERMINAL_PARTITIONS:
            raise RuntimeError("terminal partition census mismatch")
        if current.get("persisted_image_pixels") is not False:
            raise RuntimeError("partition execution persisted image pixels")
    elif cur == "REASSEMBLY_COMPLETE":
        if current.get("unique_measurement_ids") != EXPECTED_ROWS or current.get("duplicate_measurement_ids") != 0:
            raise RuntimeError("reassembly measurement-ID census mismatch")
    elif cur == "SUPPORT_GATE_COMPLETE":
        n = current.get("measurement_evaluable_species")
        if type(n) is not int or not 0 <= n <= EXPECTED_SPECIES:
            raise RuntimeError("invalid measurement-evaluable species count")
        expected = "PASS" if n >= MIN_MEASUREMENT_EVALUABLE_SPECIES else "NOT_EVALUABLE"
        if current.get("support_decision") != expected:
            raise RuntimeError("measurement-support decision mismatch")
    elif cur == "H2_COMPLETE":
        if previous.get("support_decision") != "PASS":
            raise RuntimeError("H2 cannot open after a failed support gate")
        nvec = current.get("primary_h2_vector_species")
        if type(nvec) is not int or nvec < 0:
            raise RuntimeError("invalid H2 vector count")
        if nvec < MIN_H2_VECTOR_SPECIES:
            if current.get("h2_decision") != "NOT_EVALUABLE_VECTOR_SUPPORT":
                raise RuntimeError("H2 vector-support decision mismatch")
        else:
            p = current.get("primary_structured_null_p")
            if not isinstance(p, (int, float)) or not 0 <= float(p) <= 1:
                raise RuntimeError("invalid H2 p-value")
            expected = "CONFIRMED" if float(p) < 0.05 else "NOT_CONFIRMED"
            if current.get("h2_decision") != expected:
                raise RuntimeError("H2 decision mismatch")


__all__ = [
    "AUTHORIZED_METADATA_SHA256", "AUTHORIZED_SPECIES_SHA256", "BLIND_SALT", "EXPECTED_ROWS", "EXPECTED_SPECIES",
    "EXPECTED_TERMINAL_PARTITIONS", "MEASUREMENT_SOURCE_COMMIT", "METADATA_FREEZE_COMMIT",
    "MIN_H2_VECTOR_SPECIES", "MIN_MEASUREMENT_EVALUABLE_SPECIES", "STAGES", "build_start_record",
    "measurement_batch", "measurement_id", "validate_authorization", "validate_authorized_denominator",
    "validate_candidate_metadata", "validate_stage_transition", "validate_start_record",
]
