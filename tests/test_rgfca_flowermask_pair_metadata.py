"""Offline artificial metadata only; never requests a provider or opens images."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from scripts.analysis import audit_rgfca_flowermask_pair_metadata as audit


def artificial_metadata(row, *, missing=False):
    files = []
    for index, (name, kind) in enumerate([("img.png", "image/png"), ("label.png", "image/png"),
                                          ("label_names.txt", "text/plain")]):
        if missing and name == "label.png":
            continue
        files.append({"filename": name, "id": f"00000000-0000-4000-8000-{index:012d}",
                      "folder_id": row["candidate_folder_id"], "status": "COMPLETED", "size": 3,
                      "content_details": {"size": 3, "sha256_hash": "a" * 64, "content_type": kind,
                                          "download_url": "https://invalid.example/do-not-request-image"}})
    return json.dumps(files).encode()


@pytest.fixture
def plan():
    return audit.build_plan()


def test_complete_plan_is_deterministic_and_reuses_exact_cached_source(plan):
    assert audit.build_plan() == plan
    assert len(plan["rows"]) == 300
    assert sum(row["prior_schema_opened"] for row in plan["rows"]) == 6
    assert sum(row["source"] == "request" for row in plan["rows"]) == 299
    cached = [row for row in plan["rows"] if row["source"] == "cached"]
    assert len(cached) == 1
    assert cached[0]["candidate_folder_id"] == audit.CACHED_FOLDER
    assert all(row["metadata_url"].startswith(audit.BASE) for row in plan["rows"])
    assert audit.validate_plan(plan) == plan["rows"]


@pytest.mark.parametrize("mutation", ["omit", "reorder", "duplicate", "url", "opened", "retries", "cap"])
def test_plan_mutations_fail_closed(plan, mutation):
    if mutation == "omit":
        plan["rows"].pop()
    elif mutation == "reorder":
        plan["rows"].reverse()
    elif mutation == "duplicate":
        plan["rows"][1] = deepcopy(plan["rows"][0])
    elif mutation == "url":
        plan["rows"][0]["metadata_url"] = "https://data.mendeley.com/public-files/image.png"
    elif mutation == "opened":
        plan["rows"][0]["prior_schema_opened"] = False
    elif mutation == "retries":
        plan["rules"]["retries"] = 1
    else:
        plan["rules"]["maximum_new_requests"] = 300
    with pytest.raises(ValueError):
        audit.validate_plan(plan)


def test_metadata_reports_advertised_pair_without_opening_payload(plan):
    out = audit.summarize_metadata(artificial_metadata(plan["rows"][0]), plan["rows"][0])
    assert out["required_pair_metadata_present"]
    assert not out["advertised_payload_hashes_verified"]
    assert not out["pagination_completeness_verified"]
    assert all("download_url" not in row for row in out["files"])


def test_missing_pair_is_retained_as_unresolved_not_source_error(plan):
    out = audit.summarize_metadata(artificial_metadata(plan["rows"][0], missing=True), plan["rows"][0])
    assert not out["required_pair_metadata_present"]
    assert out["file_count"] == 2


@pytest.mark.parametrize("mutation", ["duplicate_id", "duplicate_name", "casefold", "folder", "size",
                                        "float_size", "hash", "incomplete", "imageData", "detail_imageData"])
def test_bad_metadata_is_not_silently_repaired(plan, mutation):
    row = plan["rows"][0]
    entries = json.loads(artificial_metadata(row))
    if mutation == "duplicate_id":
        entries[1]["id"] = entries[0]["id"]
    elif mutation == "duplicate_name":
        entries[1]["filename"] = entries[0]["filename"]
    elif mutation == "casefold":
        entries[1]["filename"] = entries[0]["filename"].upper()
    elif mutation == "folder":
        entries[0]["folder_id"] = "other"
    elif mutation == "size":
        entries[0]["size"] = 8
    elif mutation == "float_size":
        entries[0]["size"] = 3.0
    elif mutation == "hash":
        entries[0]["content_details"]["sha256_hash"] = "bad"
    elif mutation == "incomplete":
        entries[0]["status"] = "PENDING"
    elif mutation == "imageData":
        entries[0]["imageData"] = "opaque"
    else:
        entries[0]["content_details"]["imageData"] = "opaque"
    with pytest.raises(ValueError):
        audit.summarize_metadata(json.dumps(entries).encode(), row)


@pytest.mark.parametrize("raw", [b"", b"{}", b"not JSON", b"[1]", b"null",
                                 pytest.param(b"x" * 250_001, id="oversized-metadata")])
def test_non_metadata_response_is_rejected(plan, raw):
    with pytest.raises(ValueError):
        audit.summarize_metadata(raw, plan["rows"][0])


def test_redirect_and_payload_url_blocked_before_request(plan):
    with pytest.raises(ValueError, match="redirects"):
        audit.NoRedirect().redirect_request(None, None, 302, "", {}, "https://example.org/img.png")
    row = {**plan["rows"][0], "metadata_url": "https://example.org/img.png"}
    with pytest.raises(ValueError, match="non-metadata"):
        audit.fetch_metadata(row)


def test_complete_artificial_inspection_checkpoints_all_300_and_never_fetches_cached(plan, tmp_path):
    output = tmp_path / "run"
    calls, pauses = [], []
    def fetch(row):
        checkpoint = json.loads((output / "result.json").read_bytes())
        assert len(checkpoint["terminal"]) == 300
        assert checkpoint["new_requests_started"] == len(calls) + 1
        assert row["candidate_folder_id"] != audit.CACHED_FOLDER
        calls.append(row["annotation_id"])
        return artificial_metadata(row)
    out = audit.run_inspection(plan, output, fetch=fetch, pause=pauses.append, clock=lambda: 0.0)
    assert len(calls) == 299
    assert len(pauses) == 298 and all(seconds == 1 for seconds in pauses)
    assert out["status_counts"] == {"metadata_inspected": 300}
    assert out["pair_metadata_present_count"] == 300
    assert audit.verify_saved_result(plan, output) == out
    with pytest.raises(FileExistsError):
        audit.run_inspection(plan, output, fetch=lambda row: pytest.fail("repeated request"))


def test_source_failure_stops_without_retry_and_retains_every_row(plan, tmp_path):
    calls = []
    def fetch(row):
        calls.append(row["annotation_id"])
        raise TimeoutError("artificial transport failure")
    out = audit.run_inspection(plan, tmp_path / "run", fetch=fetch, pause=lambda _: None)
    assert len(calls) == 1
    assert out["status"] == "stopped_source_failure"
    assert out["status_counts"] == {"source_inspection_failed": 1, "not_attempted_after_source_failure": 299}
    assert audit.verify_saved_result(plan, tmp_path / "run") == out
    assert not list((tmp_path / "run/responses").iterdir())


@pytest.mark.parametrize("mutation", ["authorize", "summary", "census", "requests", "status", "post_failure", "fake_summary"])
def test_saved_failure_cannot_be_promoted_or_relabelled(plan, tmp_path, mutation):
    def fail(row):
        raise TimeoutError()
    output = tmp_path / "run"
    out = audit.run_inspection(plan, output, fetch=fail, pause=lambda _: None)
    if mutation == "authorize":
        out["benchmark_execution_authorized"] = True
    elif mutation == "summary":
        out["pair_metadata_present_count"] = 1
    elif mutation == "census":
        out["terminal"].pop()
    elif mutation == "requests":
        out["new_requests_started"] += 1
    elif mutation == "status":
        out["status"] = "completed_metadata_only"
    elif mutation == "post_failure":
        out["terminal"][1]["status"] = "source_inspection_failed"
    else:
        out["terminal"][1]["summary"] = {"required_pair_metadata_present": True}
    audit.write_json(output / "result.json", out)
    with pytest.raises(ValueError):
        audit.verify_saved_result(plan, output)


def dispatch_history():
    return {"total_count": 1, "workflow_runs": [{"id": 123, "head_sha": "a" * 40,
            "event": "workflow_dispatch", "run_attempt": 1,
            "path": ".github/workflows/rgfca-flowermask-pair-metadata.yml"}]}


def test_one_first_attempt_dispatch_is_the_only_allowed_history():
    audit.validate_dispatch_history(dispatch_history(), 123, "a" * 40, 1)


@pytest.mark.parametrize("mutation", ["prior", "incomplete", "rerun", "wrong_id", "wrong_sha",
                                        "wrong_event", "wrong_workflow", "recorded_rerun"])
def test_prior_parallel_or_unverified_dispatch_cannot_acquire(mutation):
    history = dispatch_history()
    attempt = 1
    if mutation == "prior":
        history["workflow_runs"].append(deepcopy(history["workflow_runs"][0]))
        history["total_count"] = 2
    elif mutation == "incomplete":
        history["total_count"] = 101
    elif mutation == "rerun":
        attempt = 2
    else:
        key, value = {"wrong_id": ("id", 456), "wrong_sha": ("head_sha", "b" * 40),
                      "wrong_event": ("event", "push"), "wrong_workflow": ("path", "other.yml"),
                      "recorded_rerun": ("run_attempt", 2)}[mutation]
        history["workflow_runs"][0][key] = value
    with pytest.raises(ValueError):
        audit.validate_dispatch_history(history, 123, "a" * 40, attempt)


def test_committed_plan_bytes_match_workflow_and_protocol(plan):
    root = Path(__file__).resolve().parents[1]
    raw = (root / "docs/supporting/rgfca_flowermask_pair_metadata_plan_v1.json").read_bytes()
    assert json.loads(raw) == plan
    expected = "2ce2a578db6c470a9021189241e2e5c9fae561c2435896f64e61d50e0ccd4d83"
    assert audit.sha(raw) == expected
    workflow = (root / ".github/workflows/rgfca-flowermask-pair-metadata.yml").read_text(encoding="utf-8")
    protocol = (root / "docs/RGFCA_FLOWERMASK_PAIR_METADATA.md").read_text(encoding="utf-8")
    assert expected in workflow and expected in protocol
    assert "event=workflow_dispatch&per_page=100" in workflow
    assert "cancel-in-progress: false" in workflow
    assert "contents: read" in workflow and "actions: read" in workflow
    assert "test ! -f docs/supporting/rgfca_flowermask_pair_metadata_result_v1.json" in workflow
