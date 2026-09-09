"""Closed-result corruption checks; no images, models or acquisition."""
import json
import shutil

import pytest

from scripts.analysis import verify_rgfca_flowermask_pair_result as retained


BUNDLE = retained.ROOT / "data/validation/flowermask_pair_metadata_result_v1"
SUMMARY = retained.ROOT / "docs/supporting/rgfca_flowermask_pair_metadata_result_v1.json"


def test_actual_result_reconstructs_incomplete_pairing_with_full_denominator():
    result = retained.verify(BUNDLE)
    assert result == json.loads(SUMMARY.read_bytes())
    assert result["retained_files_verified"] == 303
    assert result["metadata_rows_inspected"] == 300
    assert result["pair_metadata_present_count"] == 299
    assert result["prior_schema_opened_rows_retained"] == 6
    assert len(result["source_group_row_counts"]) == 6
    assert set(result["source_group_row_counts"].values()) == {50}
    assert not result["all_300_pair_metadata_present"]
    assert not result["benchmark_execution_authorized"]
    assert result["unresolved"][0]["observed_filenames"] == ["label.png", "label_names.txt", "label_viz.png"]


@pytest.mark.parametrize("mutation", ["drop_response", "extra_payload", "response_bytes", "receipt_bytes",
                                     "history_bytes", "start_bytes", "drop_row", "promote_pair", "authorize"])
def test_corrupted_or_incomplete_closed_bundle_is_rejected(tmp_path, mutation):
    bundle = tmp_path / "saved"
    shutil.copytree(BUNDLE, bundle)
    response = next((bundle / retained.INNER / "responses").glob("*.json"))
    if mutation == "drop_response":
        response.unlink()
    elif mutation == "extra_payload":
        (bundle / "unexpected.png").write_bytes(b"artificial non-image fixture")
    elif mutation == "response_bytes":
        response.write_bytes(response.read_bytes() + b" ")
    elif mutation in {"receipt_bytes", "history_bytes", "start_bytes"}:
        relative = {"receipt_bytes": f"{retained.INNER}/result.json",
                    "history_bytes": "flowermask-pair-dispatch-history.json",
                    "start_bytes": f"{retained.INNER}/execution_start.json"}[mutation]
        path = bundle / relative
        path.write_bytes(path.read_bytes() + b" ")
    else:
        path = bundle / retained.INNER / "result.json"
        result = json.loads(path.read_bytes())
        if mutation == "drop_row":
            result["terminal"].pop()
        elif mutation == "promote_pair":
            result["all_300_pair_metadata_present"] = True
        else:
            result["benchmark_execution_authorized"] = True
        path.write_text(json.dumps(result), encoding="utf-8")
    with pytest.raises(ValueError):
        retained.verify(bundle)


def test_workflow_checks_saved_result_without_manual_or_network_acquisition():
    workflow = (retained.ROOT / ".github/workflows/rgfca-flowermask-pair-result.yml").read_text(encoding="utf-8")
    assert "--verify-summary" in workflow
    assert "workflow_dispatch:" not in workflow
    assert "contents: read" in workflow
    assert "audit_rgfca_flowermask_pair_metadata inspect" not in workflow
    assert "retention-days: 90" in workflow
