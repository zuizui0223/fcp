#!/usr/bin/env python3
"""Finalize the strict terminal v5 confirmatory cascade.

This command does not open coordinates or measure colour.  It only verifies and
assembles outcomes that were legally reached by the prospectively frozen v5
state machine.  A scientifically valid measurement-completeness
``not_evaluable`` outcome closes the cascade without requiring any coordinate or
colour-inference artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_inference_cascade import (
    STOP_CONFIRMATORY,
    TERMINAL,
    walk_confirmatory_cascade,
    validate_contract,
)


PROTOCOL = "jbi-atlas-terminal-confirmatory-cascade-v5-result-v1"
MEASUREMENT_PASS = "pass_scaleout_measurement_completeness"
MEASUREMENT_NE = "not_evaluable_scaleout_measurement_completeness"
PRIMARY_PROTOCOL = "jbi-atlas-real-colour-inference-amendment-v5-v1"
ENV_PROTOCOL = "jbi-atlas-environmental-concordance-v5-v1"
BOMBUS_NE = "pollinator_biogeographic_concordance_not_evaluable_precolour_source_gate"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _branch_outcomes(primary: Mapping[str, Any]) -> dict[str, str]:
    if (
        primary.get("protocol") != PRIMARY_PROTOCOL
        or primary.get("status") != "complete_reached_primary_spatial_branches_v5"
        or primary.get("coordinate_colour_join_performed") is not True
        or primary.get("same_standardized_colour_field_required_for_all_downstream_branches") is not True
    ):
        raise ValueError("primary real-colour result is not an authorized v5 result")
    outcomes: dict[str, str] = {}
    for row in primary.get("branches", []):
        if not isinstance(row, Mapping):
            raise ValueError("primary branch ledger contains a non-object row")
        branch = str(row.get("branch") or "")
        outcome = str(row.get("outcome") or "")
        if branch in outcomes:
            raise ValueError(f"duplicate primary branch: {branch}")
        outcomes[branch] = outcome
    if "species_conditioned_spatial_organization" not in outcomes:
        raise ValueError("primary result is missing the spatial branch")
    return outcomes


def build_final_decision(
    *,
    measurement_gate: Mapping[str, Any],
    inference: Mapping[str, Any],
    primary: Mapping[str, Any] | None = None,
    environmental: Mapping[str, Any] | None = None,
    bombus_gate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return one legal frozen-cascade decision without reading protected keys."""

    validate_contract(inference)
    measurement_status = str(measurement_gate.get("status") or "")
    if (
        measurement_gate.get("frozen_measurements") != 60000
        or measurement_gate.get("coordinates_opened") is not False
    ):
        raise ValueError("measurement gate changed the frozen denominator or opened coordinates")

    outcomes: dict[str, str]
    coordinates_opened = False
    evidence_state = "measurement_gate_only"

    if measurement_status == MEASUREMENT_NE:
        if measurement_gate.get("coordinate_join_permitted") is not False:
            raise ValueError("not-evaluable measurement gate cannot permit coordinate join")
        if primary is not None or environmental is not None or bombus_gate is not None:
            raise ValueError("downstream evidence supplied after measurement not_evaluable")
        outcomes = {"species_conditioned_spatial_organization": "not_evaluable"}
    elif measurement_status == MEASUREMENT_PASS:
        if measurement_gate.get("coordinate_join_permitted") is not True:
            raise ValueError("passed measurement gate must permit the protected join")
        if primary is None:
            raise ValueError("passed measurement gate requires the primary v5 inference result")
        outcomes = _branch_outcomes(primary)
        coordinates_opened = bool(primary.get("coordinate_colour_join_performed"))
        evidence_state = "primary_real_colour"

        primary_next = str(primary.get("next_confirmatory_branch") or "")
        if primary_next == "environmental_concordance":
            if environmental is None:
                raise ValueError("environmental result is required because that branch was reached")
            if (
                environmental.get("protocol") != ENV_PROTOCOL
                or environmental.get("status") != "complete_environmental_concordance_v5"
                or environmental.get("branch") != "environmental_concordance"
                or environmental.get("standardized_flower_field_remeasured") is not False
                or environmental.get("standardized_flower_field_restandardized") is not False
                or environmental.get("standardized_flower_colour_field_sha256")
                != primary.get("standardized_colour_field_sha256")
            ):
                raise ValueError("environmental result does not preserve the frozen v5 colour field")
            env_outcome = str(environmental.get("outcome") or "")
            outcomes["environmental_concordance"] = env_outcome
            evidence_state = "environmental_concordance"
            if env_outcome == "unsupported":
                if bombus_gate is None or bombus_gate.get("status") != BOMBUS_NE:
                    raise ValueError("reached pollinator branch must use the binding pre-colour Bombus gate")
                outcomes["pollinator_biogeographic_concordance"] = "not_evaluable"
                evidence_state = "pollinator_precolour_source_gate"
            elif bombus_gate is not None:
                raise ValueError("pollinator evidence supplied when pollinator branch was not reached")
        else:
            if environmental is not None or bombus_gate is not None:
                raise ValueError("downstream result supplied for a branch not reached by primary inference")
    else:
        raise ValueError(f"unexpected measurement gate status: {measurement_status!r}")

    decisions = walk_confirmatory_cascade(outcomes, inference)
    if not decisions:
        raise RuntimeError("frozen cascade produced no decision")
    final_state = decisions[-1].next_step
    if final_state not in (STOP_CONFIRMATORY, TERMINAL):
        raise RuntimeError("finalized cascade did not terminate")
    branch_ledger = [
        {"branch": item.branch, "outcome": item.outcome, "next_step": item.next_step}
        for item in decisions
    ]
    supported = [row["branch"] for row in branch_ledger if row["outcome"] == "supported"]
    return {
        "protocol": PROTOCOL,
        "inference_version": inference["version"],
        "status": "complete_terminal_confirmatory_cascade_v5",
        "measurement_completeness_status": measurement_status,
        "coordinate_colour_join_performed": coordinates_opened,
        "branch_ledger": branch_ledger,
        "supported_branches": supported,
        "final_confirmatory_state": final_state,
        "evidence_state": evidence_state,
        "not_evaluable_advanced_as_unsupported": False,
        "post_colour_pollinator_source_substitution_used": False,
        "claim_ceiling": (
            "Only branches explicitly listed as supported in the frozen v5 cascade may be "
            "reported as confirmatory terminal-atlas support. not_evaluable never advances."
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-gate", type=Path, required=True)
    parser.add_argument("--primary-result", type=Path)
    parser.add_argument("--environmental-result", type=Path)
    parser.add_argument("--bombus-source-gate", type=Path)
    parser.add_argument(
        "--inference-v5",
        type=Path,
        default=ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    measurement_gate = load_json(args.measurement_gate)
    inference = load_json(args.inference_v5)
    primary = load_json(args.primary_result) if args.primary_result else None
    environmental = load_json(args.environmental_result) if args.environmental_result else None
    bombus = load_json(args.bombus_source_gate) if args.bombus_source_gate else None
    result = build_final_decision(
        measurement_gate=measurement_gate,
        inference=inference,
        primary=primary,
        environmental=environmental,
        bombus_gate=bombus,
    )
    parents = {
        "measurement_gate": sha256(args.measurement_gate),
        "inference_v5": sha256(args.inference_v5),
    }
    if args.primary_result:
        parents["primary_result"] = sha256(args.primary_result)
    if args.environmental_result:
        parents["environmental_result"] = sha256(args.environmental_result)
    if args.bombus_source_gate:
        parents["bombus_source_gate"] = sha256(args.bombus_source_gate)
    result["parents_sha256"] = parents
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
