#!/usr/bin/env python3
"""Gate original FCP 4,730 U100 downstream statistics on BOTH real ROI terminals.

No image opening, photo classification, ecological estimator or threshold here.
The two immutable workflow IDs are independent; whichever completed second
is responsible for starting the single all-source statistical replay.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

EXPECTED = {
    "new_U100": {
        "id": 38997169400,
        "head_sha": "b79d57686394c1bcbb32fcc2142f0dd19db4ced9",
        "name": "FCP full 4730 original flower ROI and H1 D 50km W replay",
        "event": "push",
        "head_branch": "main",
    },
    "wave4": {
        "id": 38992543115,
        "head_sha": "ef448c9514e0c702e7887e00b288ee966728b0d5",
        "name": "FCP wave4 identical original photo ROI and D W spatial replay",
        "event": "pull_request",
        "head_branch": "analysis/fcp-wave4-roi-import-repair-20261017",
    },
}
SOURCE_RUN = 38996743273


def _verify_run(label: str, data: dict) -> None:
    expected = EXPECTED[label]
    for name, original in expected.items():
        if data.get(name) != original:
            raise RuntimeError(f"{label}: immutable original ROI {name} drift")
    if data.get("run_attempt") != 1:
        raise RuntimeError(f"{label}: original ROI first-attempt identity drift")
    if data.get("head_repository", {}).get("full_name") != "zuizui0223/fcp":
        raise RuntimeError(f"{label}: different source repository cannot authorise results")


def decide(
    new_run: dict,
    wave4_run: dict,
    trigger_id: int,
    *,
    manual: bool = False,
) -> dict:
    """Pure protocol gate; does NOT infer biological colour results."""
    both = {"new_U100": new_run, "wave4": wave4_run}
    for cohort, data in both.items():
        _verify_run(cohort, data)
    if trigger_id not in {r["id"] for r in both.values()}:
        raise RuntimeError("event did not originate from either pinned genuine original ROI run")
    status = {cohort: {"status": r.get("status"), "conclusion": r.get("conclusion")} for cohort,r in both.items()}
    complete = all(r.get("status") == "completed" and r.get("conclusion") == "success"
                   for r in both.values())
    if not complete:
        return {
            "ready": False,
            "reason": "BOTH_ACTUAL_PIXEL_JOBS_NOT_YET_SUCCESSFUL",
            "upstream": status,
            "source_run_id": SOURCE_RUN,
        }
    # A 'workflow_run' delivery may arrive after the other source also finished.
    # Only the source which actually completed LAST owns the automatic join.
    # Prevents double 447,800-photo statistical runs on delayed delivery.
    latest = max(
        both.values(),
        key=lambda r: (r.get("updated_at") or "", int(r["id"]))
    )["id"]
    if not manual and trigger_id != latest:
        return {
            "ready": False,
            "reason": "COMPLETED_BUT_OTHER_ROI_EVENT_IS_THE_LAST_FINISHER",
            "upstream": status,
            "last_finisher_id": latest,
            "source_run_id": SOURCE_RUN,
        }
    return {
        "ready": True,
        "reason": "BOTH_PINNED_ROI_RUNS_FINISHED_SUCCESSFULLY",
        "source_run_id": SOURCE_RUN,
        "full_U100_ROI_run_id": EXPECTED["new_U100"]["id"],
        "wave4_ROI_run_id": EXPECTED["wave4"]["id"],
        "upstream": status,
        "last_finisher_id": latest,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--new-run", required=True, type=Path)
    p.add_argument("--wave4-run", required=True, type=Path)
    p.add_argument("--trigger-id", type=int, required=True)
    p.add_argument("--manual", action="store_true")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--github-output", required=True, type=Path)
    a = p.parse_args()
    rec = decide(json.loads(a.new_run.read_text()), json.loads(a.wave4_run.read_text()),
                 a.trigger_id, manual=a.manual)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
    with a.github_output.open("a") as fh:
        fh.write("ready=" + str(rec["ready"]).lower() + "\n")
    print("FCP_4730_GENUINE_DUAL_ROI_TERMINAL_GATE", json.dumps(rec, sort_keys=True))


if __name__ == "__main__":
    main()
