#!/usr/bin/env python3
"""Validate the prospectively frozen terminal image-first atlas inference contract.

The validator deliberately operates without opening candidate images.  It checks
branch semantics and verifies that all legacy/reference artifacts pinned by the v5
contract still have their frozen Git-blob identities.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_inference_cascade import (  # noqa: E402
    BRANCHES,
    STOP_CONFIRMATORY,
    TERMINAL,
    load_contract,
    next_step,
)


def git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def validate_immutable_blobs(root: Path, contract: dict) -> list[dict[str, str]]:
    checked: list[dict[str, str]] = []
    for relative, expected in contract["immutable_git_blobs"].items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"missing frozen/reference artifact: {relative}")
        observed = git_blob_sha(path)
        if observed != expected:
            raise RuntimeError(
                f"immutable artifact changed: {relative}: expected {expected}, observed {observed}"
            )
        checked.append({"path": relative, "blob": observed})
    return checked


def validate_transition_semantics(contract: dict) -> None:
    if next_step(BRANCHES[0], "supported", contract) != BRANCHES[1]:
        raise RuntimeError("spatial support must proceed to shared-transition test")
    if next_step(BRANCHES[0], "unsupported", contract) != BRANCHES[2]:
        raise RuntimeError("spatial non-support must proceed to environment")
    if next_step(BRANCHES[1], "unsupported", contract) != BRANCHES[2]:
        raise RuntimeError("shared-transition non-support must proceed to environment")
    if next_step(BRANCHES[2], "unsupported", contract) != BRANCHES[3]:
        raise RuntimeError("environmental non-support must proceed to pollinator geography")

    for branch in BRANCHES[:-1]:
        if next_step(branch, "not_evaluable", contract) != STOP_CONFIRMATORY:
            raise RuntimeError(f"not_evaluable illegally advances from {branch}")
    if next_step(BRANCHES[-1], "not_evaluable", contract) != TERMINAL:
        raise RuntimeError("terminal pollinator not_evaluable state is malformed")


def validate_freeze_firewalls(contract: dict) -> None:
    if contract["pixel_status_at_freeze"] != "not_revealed":
        raise RuntimeError("contract is not a pre-pixel freeze")
    if contract["roi_freeze"]["candidate_image_retuning_forbidden"] is not True:
        raise RuntimeError("ROI may be retuned on terminal candidate images")
    if contract["colour_field_freeze"]["downstream_remeasurement_or_restandardization_forbidden"] is not True:
        raise RuntimeError("downstream branches can redefine the frozen colour field")
    if contract["terminal_scaleout_sampling"]["favourable_cohort_selection_forbidden"] is not True:
        raise RuntimeError("favourable cohort selection is not prohibited")
    if contract["diagnostic_firewall"]["diagnostic_results_cannot_promote_confirmatory_claims"] is not True:
        raise RuntimeError("diagnostic-only results can contaminate confirmatory inference")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        default="docs/supporting/jbi_image_first_atlas_inference_contract_v5.json",
    )
    parser.add_argument("--repo-root", default=str(ROOT))
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    contract_path = root / args.contract
    contract = load_contract(contract_path)
    validate_transition_semantics(contract)
    validate_freeze_firewalls(contract)
    blobs = validate_immutable_blobs(root, contract)

    result = {
        "status": "pass",
        "contract": str(Path(args.contract)),
        "version": contract["version"],
        "pixel_status_at_freeze": contract["pixel_status_at_freeze"],
        "immutable_blob_count": len(blobs),
        "confirmatory_sequence": list(BRANCHES),
        "not_evaluable_advances": False,
        "legacy_artifacts_unchanged": True,
        "terminal_scaleout": {
            "species": contract["terminal_scaleout_sampling"]["species_total"],
            "photos_target": contract["terminal_scaleout_sampling"]["photo_target_total"],
            "cohorts": contract["terminal_scaleout_sampling"]["cohort_count"],
        },
    }
    if args.output_json:
        output = root / args.output_json
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
