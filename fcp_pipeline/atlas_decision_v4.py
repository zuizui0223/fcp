"""Strict pre-image decision state machine for the FCP image-first atlas v4.

This module is deliberately outcome-only: it has no image loading, ROI, coordinate,
environment, or pollinator interface.  It converts already frozen branch outcomes
into the confirmatory decision tree and preserves any downstream computations as
diagnostic-only when an upstream branch did not authorize the fallback.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping, Sequence


PROTOCOL = "jbi-image-first-global-flower-colour-atlas-inference-v4"
BRANCHES = (
    "species_conditioned_shared_transition",
    "environmental_boundary_concordance",
    "pollinator_biogeographic_concordance",
)
OUTCOMES = ("supported", "unsupported", "not_evaluable")
ENVIRONMENTAL_FAMILIES = ("macroclimate", "terrain", "land_cover", "ecoregion")


def git_blob_sha1(path: Path) -> str:
    """Return the Git blob SHA-1 of a working-tree file without invoking Git."""

    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def validate_v4_contract(contract: Mapping[str, Any], *, root: Path = Path(".")) -> None:
    """Fail closed if the post-v3, pre-image freeze contract drifts."""

    if contract.get("protocol") != PROTOCOL:
        raise ValueError("unexpected atlas v4 inference protocol")
    if contract.get("status") != "prospectively_frozen_before_any_scaleout_candidate_pixel":
        raise ValueError("v4 inference contract must be frozen before scale-out pixels")
    if contract.get("candidate_pixels_opened_at_freeze") is not False:
        raise ValueError("candidate pixels were already opened at v4 freeze")

    parent = contract.get("immutable_parent_v3", {})
    parent_path = root / str(parent.get("path", ""))
    if not parent_path.is_file() or git_blob_sha1(parent_path) != parent.get("git_blob_sha1"):
        raise ValueError("immutable v3 parent identity mismatch")

    invariants = contract.get("immutable_legacy_invariants", [])
    required_roles = {"stage_a_six_species", "stage_b_six_species", "literature_34_species", "negative_validation_three_species"}
    observed_roles = {str(row.get("role")) for row in invariants}
    if observed_roles != required_roles:
        raise ValueError("legacy 6/34/3 invariant registry changed")
    for row in invariants:
        path = root / str(row.get("path", ""))
        expected = str(row.get("git_blob_sha1", ""))
        if not path.is_file() or git_blob_sha1(path) != expected:
            raise ValueError(f"immutable legacy artifact changed: {row.get('role')}")

    sampling = contract.get("sampling_repetitions", {})
    expected_sampling = {
        "cohort_count": 8,
        "species_per_cohort": 25,
        "observations_per_species": 300,
        "total_species": 200,
        "total_observations": 60000,
    }
    for key, expected in expected_sampling.items():
        if int(sampling.get(key, -1)) != expected:
            raise ValueError(f"sampling repetition dimension changed: {key}")
    if sampling.get("species_overlap_between_frozen_cohorts") is not False:
        raise ValueError("the eight frozen cohorts must remain species-disjoint")
    if sampling.get("resample_after_pixels") is not False:
        raise ValueError("post-pixel sampling repetitions are prohibited")
    if sampling.get("all_cohorts_required") is not True:
        raise ValueError("all eight sampling repetitions must finish")

    roi = contract.get("roi_validity", {})
    if roi.get("estimator_version") != "v4" or roi.get("retuning_after_scaleout_pixels") is not False:
        raise ValueError("ROI v4 must stay frozen after scale-out pixels open")
    if roi.get("independent_development_passed") is not True or roi.get("locked_test_passed") is not True:
        raise ValueError("ROI v4 requires both development and locked-test qualification")

    cascade = contract.get("strict_confirmatory_cascade", {})
    if tuple(cascade.get("branches", ())) != BRANCHES:
        raise ValueError("strict confirmatory branch order changed")
    if tuple(cascade.get("outcomes", ())) != OUTCOMES:
        raise ValueError("strict confirmatory outcome vocabulary changed")
    if cascade.get("advance_only_on") != "unsupported":
        raise ValueError("fallback must advance only after an evaluable unsupported result")
    if cascade.get("not_evaluable_advances") is not False:
        raise ValueError("not_evaluable must stop the confirmatory cascade")
    if cascade.get("supported_advances") is not False:
        raise ValueError("supported must stop the confirmatory cascade")
    if tuple(cascade.get("environmental_family_order", ())) != ENVIRONMENTAL_FAMILIES:
        raise ValueError("environmental family order changed")
    if cascade.get("same_frozen_colour_field_required") is not True:
        raise ValueError("all branches must reuse the same frozen colour field")
    if cascade.get("blocked_downstream_results_are") != "diagnostic_only":
        raise ValueError("blocked downstream computations must remain diagnostic-only")

    current = contract.get("current_preimage_state", {})
    if current.get("shared_transition_qualification") != "not_evaluable":
        raise ValueError("failed shared-transition signal recovery must remain not_evaluable")
    if current.get("confirmatory_cascade_can_advance") is not False:
        raise ValueError("current not_evaluable spatial branch cannot authorize environment")

    ledger = contract.get("decision_ledger_schema", {})
    required_fields = {
        "sampling_repetition_id",
        "cohort_id",
        "species_blind_id",
        "branch",
        "outcome",
        "reason_code",
        "effective_n",
        "data_version",
        "config_sha256",
        "manifest_sha256",
        "confirmatory_status",
    }
    if set(ledger.get("required_fields", ())) != required_fields:
        raise ValueError("decision ledger schema changed")


def strict_cascade_decisions(
    branch_results: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Apply the v4 cascade without treating unknown evidence as a negative result.

    `branch_results` may contain precomputed results for all branches.  Results from
    a branch that was not authorized by the confirmatory cascade are retained in
    `diagnostic_results`, but can never become the promoted conclusion.
    """

    unknown = set(branch_results) - set(BRANCHES)
    if unknown:
        raise ValueError(f"unknown atlas inference branches: {sorted(unknown)}")

    confirmatory: dict[str, dict[str, Any]] = {}
    diagnostic: dict[str, dict[str, Any]] = {}
    authorized = True
    blocker: str | None = None
    promoted: str | None = None

    for branch in BRANCHES:
        raw = dict(branch_results.get(branch, {}))
        raw_outcome = raw.get("outcome")
        if raw_outcome is not None and raw_outcome not in OUTCOMES:
            raise ValueError(f"invalid outcome for {branch}: {raw_outcome!r}")

        if not authorized:
            if raw_outcome is not None:
                diagnostic[branch] = {
                    **raw,
                    "confirmatory_status": "diagnostic_only",
                    "blocked_by": blocker,
                }
            confirmatory[branch] = {
                "outcome": "not_evaluable",
                "reason_code": f"blocked_by_upstream_{blocker}",
                "confirmatory_status": "blocked",
            }
            continue

        if raw_outcome is None:
            outcome = "not_evaluable"
            reason = "missing_frozen_branch_result"
        else:
            outcome = str(raw_outcome)
            reason = str(raw.get("reason_code") or f"branch_{outcome}")

        confirmatory[branch] = {
            **raw,
            "outcome": outcome,
            "reason_code": reason,
            "confirmatory_status": "evaluated" if raw_outcome is not None else "not_evaluable",
        }

        if outcome == "unsupported":
            continue
        if outcome == "supported":
            promoted = branch
            blocker = branch
            authorized = False
        elif outcome == "not_evaluable":
            blocker = branch
            authorized = False

    if promoted is not None:
        terminal = "supported"
    elif all(confirmatory[branch]["outcome"] == "unsupported" for branch in BRANCHES):
        terminal = "unsupported"
    else:
        terminal = "not_evaluable"

    return {
        "branch_order": list(BRANCHES),
        "advance_only_on": "unsupported",
        "confirmatory_results": confirmatory,
        "diagnostic_results": diagnostic,
        "promoted_conclusion": promoted,
        "terminal_outcome": terminal,
    }


def validate_decision_ledger_rows(rows: Sequence[Mapping[str, Any]]) -> None:
    """Validate the future species × repetition × branch audit ledger."""

    required = {
        "sampling_repetition_id",
        "cohort_id",
        "species_blind_id",
        "branch",
        "outcome",
        "reason_code",
        "effective_n",
        "data_version",
        "config_sha256",
        "manifest_sha256",
        "confirmatory_status",
    }
    seen: set[tuple[str, str, str, str]] = set()
    for row in rows:
        missing = required - set(row)
        if missing:
            raise ValueError(f"decision ledger row missing fields: {sorted(missing)}")
        if row["branch"] not in BRANCHES:
            raise ValueError("decision ledger contains an unknown branch")
        if row["outcome"] not in OUTCOMES:
            raise ValueError("decision ledger contains an invalid outcome")
        if row["confirmatory_status"] not in {"evaluated", "blocked", "diagnostic_only", "not_evaluable"}:
            raise ValueError("decision ledger contains an invalid confirmatory status")
        key = (
            str(row["sampling_repetition_id"]),
            str(row["cohort_id"]),
            str(row["species_blind_id"]),
            str(row["branch"]),
        )
        if key in seen:
            raise ValueError("decision ledger contains duplicate species-repetition-branch rows")
        seen.add(key)
