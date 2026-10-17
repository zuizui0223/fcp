"""Synthetic first-attempt identity gates, not biological flower-colour outputs."""
import copy
import pytest
from scripts.analysis import fcp_rescue_original_256_roi_terminals_20261017 as m


def fixture(label="new_U100", final="failure", missing=None, conclusion="failure"):
    run_id=next(i for i,r in m.PINNED.items() if r["label"]==label)
    expected=m.PINNED[run_id]
    run={
        "id":run_id,"name":expected["name"],"event":expected["event"],
        "head_branch":expected["head_branch"],"head_sha":expected["head_sha"],
        "head_repository":{"full_name":"zuizui0223/fcp"},
        "run_attempt":1,"status":"completed","conclusion":conclusion,
    }
    jobs=[{"id":1,"name":"freeze-new-original-pixel-IDs","conclusion":"success"}]
    if expected["preflight_name"]:
        jobs.append({"id":2,"name":expected["preflight_name"],"conclusion":"success"})
    jobs.append({"id":3,"name":m.FINAL,"conclusion":final})
    for b in range(2):
        for s in range(32):
            for p in range(4):
                name=f"unmodified-original-roi-workers ({b}, {s}, {p})"
                jobs.append({"id":100+b*128+s*4+p,
                             "name":name,"conclusion":None if name==missing else "success"})
    return run,[{"jobs":jobs[:98]},{"jobs":jobs[98:198]},{"jobs":jobs[198:]}]


@pytest.mark.parametrize("cohort",["new_U100","wave4"])
def test_rescue_exact_256_actual_success_and_only_failed_statistics(cohort):
    run,pages=fixture(label=cohort)
    record=m.decide(run,pages)
    assert record["ready"] is True
    assert record["terminal_worker_jobs"]==256
    assert record["original_worker_jobs_successful"]==256
    assert record["source_original_photo_IDs"]==(248100 if cohort=="new_U100" else 49800)
    assert record["original_H1_D_50km_H2_W_not_rerun"] is True
    assert record["original_pixel_measurement_method_replaced"] is False


def test_no_rescue_on_successful_original_completed_roi_science():
    run,pages=fixture(final="success",conclusion="success")
    record=m.decide(run,pages)
    assert record["ready"] is False
    assert record["reason"]=="FULL_ACTUAL_ORIGINAL_PIXEL_AND_SCIENCE_ALREADY_COMPLETE"


def test_no_rescue_if_any_original_pixel_partition_not_complete():
    run,pages=fixture(missing="unmodified-original-roi-workers (1, 31, 3)")
    rec=m.decide(run,pages)
    assert rec["ready"] is False
    assert rec["original_worker_jobs_successful"]==255


def test_reject_incomplete_original_matrix_not_forged_by_completion_flag():
    run,pages=fixture()
    pages[0]["jobs"].pop(7)
    with pytest.raises(RuntimeError,match="identities changed"):
        m.decide(run,pages)


def test_reject_repeated_original_partition_ids():
    run,pages=fixture()
    pages[0]["jobs"][4]["name"]=pages[0]["jobs"][5]["name"]
    with pytest.raises(RuntimeError,match="duplicate original photo ROI job records"):
        m.decide(run,pages)


def test_reject_forged_head_sha_repo_and_inference_run_attempt():
    run,pages=fixture()
    for key,bad in [("head_sha","changed"),("run_attempt",2)]:
        forged=copy.deepcopy(run)
        forged[key]=bad
        with pytest.raises(RuntimeError,match="drift"):
            m.decide(forged,pages)
    forged=copy.deepcopy(run)
    forged["head_repository"]={"full_name":"other/fcp"}
    with pytest.raises(RuntimeError,match="integrity drift"):
        m.decide(forged,pages)


def test_reject_invalid_source_run_without_pixel_processing():
    run,pages=fixture()
    run["id"]=38993836606
    with pytest.raises(RuntimeError,match="two exact frozen"):
        m.decide(run,pages)


def test_rescue_when_only_terminal_statistics_cancelled_after_all_source_pixels():
    run,pages=fixture(final="cancelled",conclusion="cancelled")
    assert m.decide(run,pages)["ready"] is True


def test_no_rescue_if_run_still_running():
    run,pages=fixture()
    run["status"]="in_progress"
    run["conclusion"]=None
    assert m.decide(run,pages)["ready"] is False
