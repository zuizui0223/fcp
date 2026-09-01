"""Frozen terminal-scaleout inference cascade for the FCP image-first atlas.

This module is intentionally small and policy-like.  It does not measure colour or
choose spatial geometry.  It enforces the prospectively frozen order in
``jbi_image_first_atlas_inference_contract_v5.json`` and keeps ``not_evaluable``
separate from ``unsupported``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import json
import re


BRANCHES = (
    "species_conditioned_spatial_organization",
    "shared_transition",
    "environmental_concordance",
    "pollinator_biogeographic_concordance",
)
OUTCOMES = ("supported", "unsupported", "not_evaluable")
EXECUTION_STATES = (
    "run_confirmatory",
    "not_reached",
    "stopped_not_evaluable",
    "stopped_supported",
    "terminal_complete",
    "diagnostic_only",
)
STOP_CONFIRMATORY = "STOP_CONFIRMATORY"
TERMINAL = "TERMINAL"
_GIT_BLOB_RE = re.compile(r"^[0-9a-f]{40}$")


@dataclass(frozen=True)
class BranchDecision:
    branch: str
    outcome: str
    next_step: str


class CascadeContractError(ValueError):
    """Raised when the frozen inference contract is internally inconsistent."""


def load_contract(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_contract(value)
    return value


def validate_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("version") != "jbi_image_first_atlas_inference_v5":
        raise CascadeContractError("unexpected inference contract version")
    if contract.get("pixel_status_at_freeze") != "not_revealed":
        raise CascadeContractError("v5 must have been frozen before terminal pixels were revealed")
    if contract.get("legacy_analyses_are_immutable") is not True:
        raise CascadeContractError("legacy analyses must remain immutable")

    if tuple(contract.get("confirmatory_sequence", ())) != BRANCHES:
        raise CascadeContractError("confirmatory branch order is not frozen as expected")
    if tuple(contract.get("state_vocabulary", ())) != OUTCOMES:
        raise CascadeContractError("biological outcome vocabulary changed")

    execution = tuple(contract.get("execution_state_vocabulary", ()))
    if execution != EXECUTION_STATES:
        raise CascadeContractError("execution-state vocabulary changed")

    expected_transitions = {
        "species_conditioned_spatial_organization": {
            "supported": "shared_transition",
            "unsupported": "environmental_concordance",
            "not_evaluable": STOP_CONFIRMATORY,
        },
        "shared_transition": {
            "supported": STOP_CONFIRMATORY,
            "unsupported": "environmental_concordance",
            "not_evaluable": STOP_CONFIRMATORY,
        },
        "environmental_concordance": {
            "supported": STOP_CONFIRMATORY,
            "unsupported": "pollinator_biogeographic_concordance",
            "not_evaluable": STOP_CONFIRMATORY,
        },
        "pollinator_biogeographic_concordance": {
            "supported": TERMINAL,
            "unsupported": TERMINAL,
            "not_evaluable": TERMINAL,
        },
    }
    if contract.get("transition_table") != expected_transitions:
        raise CascadeContractError("transition table differs from the frozen strict cascade")

    sampling = contract.get("terminal_scaleout_sampling", {})
    species_total = int(sampling.get("species_total", -1))
    cohort_count = int(sampling.get("cohort_count", -1))
    species_per_cohort = int(sampling.get("species_per_cohort", -1))
    if (species_total, cohort_count, species_per_cohort) != (200, 8, 25):
        raise CascadeContractError("terminal sampling dimensions changed")
    if cohort_count * species_per_cohort != species_total:
        raise CascadeContractError("frozen cohorts do not exhaust the species universe")
    for key in (
        "cohorts_are_disjoint",
        "cohorts_exhaust_species_universe",
        "all_cohorts_must_be_reported",
        "favourable_cohort_selection_forbidden",
    ):
        if sampling.get(key) is not True:
            raise CascadeContractError(f"sampling safeguard is not frozen: {key}")

    for branch in ("species_conditioned_spatial_organization", "shared_transition"):
        rule = contract.get(branch, {})
        if float(rule.get("alpha", -1.0)) != 0.05:
            raise CascadeContractError(f"{branch} alpha changed")
        if int(rule.get("minimum_evaluable_cohorts", -1)) != 8:
            raise CascadeContractError(f"{branch} must retain all eight cohorts")
        if int(rule.get("minimum_directionally_concordant_cohorts", -1)) != 6:
            raise CascadeContractError(f"{branch} replication threshold changed")

    if contract.get("colour_field_freeze", {}).get("downstream_remeasurement_or_restandardization_forbidden") is not True:
        raise CascadeContractError("downstream branches may not redefine the colour field")
    if contract.get("diagnostic_firewall", {}).get("diagnostic_results_cannot_promote_confirmatory_claims") is not True:
        raise CascadeContractError("diagnostic firewall is not enforced")

    blobs = contract.get("immutable_git_blobs", {})
    if not isinstance(blobs, Mapping) or not blobs:
        raise CascadeContractError("immutable Git-blob registry is empty")
    for path, blob in blobs.items():
        if not str(path).strip() or not _GIT_BLOB_RE.fullmatch(str(blob)):
            raise CascadeContractError(f"invalid immutable blob pin: {path}={blob}")


def next_step(branch: str, outcome: str, contract: Mapping[str, Any]) -> str:
    validate_contract(contract)
    if branch not in BRANCHES:
        raise ValueError(f"unknown branch: {branch}")
    if outcome not in OUTCOMES:
        raise ValueError(f"unknown outcome: {outcome}")
    return str(contract["transition_table"][branch][outcome])


def branch_decision(branch: str, outcome: str, contract: Mapping[str, Any]) -> BranchDecision:
    return BranchDecision(branch=branch, outcome=outcome, next_step=next_step(branch, outcome, contract))


def classify_spatial_organization(
    *,
    pooled_p_value: float | None,
    cohort_directions: Sequence[float | None],
    contract: Mapping[str, Any],
) -> str:
    """Classify the frozen pooled spatial-organization test.

    ``cohort_directions`` are observed statistic minus null center for C01-C08.
    Missing/non-finite values make the confirmatory branch not evaluable.
    """
    validate_contract(contract)
    rule = contract["species_conditioned_spatial_organization"]
    return _classify_pooled_with_replication(
        pooled_p_value=pooled_p_value,
        cohort_directions=cohort_directions,
        alpha=float(rule["alpha"]),
        minimum_evaluable=int(rule["minimum_evaluable_cohorts"]),
        minimum_positive=int(rule["minimum_directionally_concordant_cohorts"]),
    )


def classify_shared_transition(
    *,
    pooled_p_value: float | None,
    cohort_directions: Sequence[float | None],
    qualification_passed: bool,
    contract: Mapping[str, Any],
) -> str:
    """Classify shared transition only when its exact preimage qualification passed."""
    validate_contract(contract)
    if qualification_passed is not True:
        return "not_evaluable"
    rule = contract["shared_transition"]
    return _classify_pooled_with_replication(
        pooled_p_value=pooled_p_value,
        cohort_directions=cohort_directions,
        alpha=float(rule["alpha"]),
        minimum_evaluable=int(rule["minimum_evaluable_cohorts"]),
        minimum_positive=int(rule["minimum_directionally_concordant_cohorts"]),
    )


def _classify_pooled_with_replication(
    *,
    pooled_p_value: float | None,
    cohort_directions: Sequence[float | None],
    alpha: float,
    minimum_evaluable: int,
    minimum_positive: int,
) -> str:
    if pooled_p_value is None:
        return "not_evaluable"
    try:
        p_value = float(pooled_p_value)
    except (TypeError, ValueError):
        return "not_evaluable"
    if not (0.0 <= p_value <= 1.0):
        return "not_evaluable"
    if len(cohort_directions) != minimum_evaluable:
        return "not_evaluable"

    directions: list[float] = []
    for value in cohort_directions:
        if value is None:
            return "not_evaluable"
        try:
            number = float(value)
        except (TypeError, ValueError):
            return "not_evaluable"
        if number != number or number in (float("inf"), float("-inf")):
            return "not_evaluable"
        directions.append(number)

    positive = sum(value > 0.0 for value in directions)
    if p_value <= alpha and positive >= minimum_positive:
        return "supported"
    return "unsupported"


def walk_confirmatory_cascade(
    outcomes: Mapping[str, str], contract: Mapping[str, Any]
) -> list[BranchDecision]:
    """Walk only branches that are legally reached by the frozen cascade.

    The caller must provide an outcome for every branch that is actually reached.
    Outcomes supplied for an unreached branch are rejected, preventing post-hoc
    selection of a favourable downstream overlay.
    """
    validate_contract(contract)
    decisions: list[BranchDecision] = []
    current = BRANCHES[0]
    visited: set[str] = set()
    while current not in (STOP_CONFIRMATORY, TERMINAL):
        if current in visited:
            raise RuntimeError("cascade contains a cycle")
        visited.add(current)
        if current not in outcomes:
            raise ValueError(f"missing outcome for reached branch: {current}")
        outcome = outcomes[current]
        decision = branch_decision(current, outcome, contract)
        decisions.append(decision)
        current = decision.next_step

    reached = {item.branch for item in decisions}
    illegal = set(outcomes) - reached
    if illegal:
        raise ValueError(
            "outcome supplied for branch that was not confirmatorily reached: "
            + ", ".join(sorted(illegal))
        )
    return decisions


def validate_decision_ledger_rows(
    rows: Iterable[Mapping[str, Any]], contract: Mapping[str, Any]
) -> None:
    """Validate branch-order semantics for a persisted decision ledger."""
    validate_contract(contract)
    required = tuple(contract["required_decision_ledger_columns"])
    rows_list = list(rows)
    if not rows_list:
        raise ValueError("decision ledger is empty")
    for index, row in enumerate(rows_list, start=1):
        missing = [column for column in required if column not in row]
        if missing:
            raise ValueError(f"decision ledger row {index} missing: {missing}")
        branch = str(row["branch"])
        if branch not in BRANCHES:
            raise ValueError(f"decision ledger row {index} has unknown branch: {branch}")
        state = str(row["execution_state"])
        if state not in EXECUTION_STATES:
            raise ValueError(f"decision ledger row {index} has unknown execution state: {state}")
        outcome = str(row["outcome"])
        if state in {"run_confirmatory", "stopped_not_evaluable", "stopped_supported", "terminal_complete"}:
            if outcome not in OUTCOMES:
                raise ValueError(f"decision ledger row {index} lacks a valid biological outcome")
        elif state == "diagnostic_only" and outcome not in OUTCOMES:
            raise ValueError(f"diagnostic row {index} lacks a valid descriptive outcome")

    confirmatory = [
        row
        for row in rows_list
        if str(row["execution_state"])
        in {"run_confirmatory", "stopped_not_evaluable", "stopped_supported", "terminal_complete"}
    ]
    if not confirmatory:
        raise ValueError("decision ledger contains no confirmatory branch")
    outcome_map: dict[str, str] = {}
    for row in confirmatory:
        branch = str(row["branch"])
        outcome = str(row["outcome"])
        prior = outcome_map.get(branch)
        if prior is not None and prior != outcome:
            raise ValueError(f"conflicting confirmatory outcomes for {branch}")
        outcome_map[branch] = outcome
    walk_confirmatory_cascade(outcome_map, contract)

    if any(
        str(row["execution_state"]) == "diagnostic_only"
        and str(row.get("advance_reason", "")).lower().startswith("confirmatory")
        for row in rows_list
    ):
        raise ValueError("diagnostic-only results cannot advance confirmatory inference")
