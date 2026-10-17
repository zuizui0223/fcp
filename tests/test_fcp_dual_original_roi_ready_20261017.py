"""Synthetic protocol tests only: no flower images or biological results."""
import copy
import pytest
from scripts.analysis import fcp_dual_original_roi_ready_20261017 as m


def valid(label, *, timestamp="2026-10-17T05:00:00Z", status="completed", conclusion="success"):
    return {
        **m.EXPECTED[label],
        "run_attempt": 1,
        "head_repository": {"full_name": "zuizui0223/fcp"},
        "updated_at": timestamp,
        "status": status,
        "conclusion": conclusion,
    }


def test_wait_for_wave4_even_when_full_4730_new_pixel_measurement_succeeded():
    new=valid("new_U100")
    wave=valid("wave4", status="in_progress", conclusion=None)
    rec=m.decide(new,wave,m.EXPECTED["new_U100"]["id"])
    assert rec["ready"] is False
    assert rec["reason"]=="BOTH_ACTUAL_PIXEL_JOBS_NOT_YET_SUCCESSFUL"


def test_wait_for_full_4730_even_when_wave4_completed():
    new=valid("new_U100", status="queued",conclusion=None)
    wave=valid("wave4")
    assert not m.decide(new,wave,wave["id"])["ready"]


def test_trigger_on_later_completed_source_new_and_never_duplicate():
    new=valid("new_U100",timestamp="2026-10-17T10:00:00Z")
    wave=valid("wave4",timestamp="2026-10-17T09:00:00Z")
    assert m.decide(new,wave,new["id"])["ready"] is True
    old=m.decide(new,wave,wave["id"])
    assert old["ready"] is False
    assert old["last_finisher_id"]==new["id"]


def test_trigger_on_later_completed_source_wave4_and_never_duplicate():
    new=valid("new_U100",timestamp="2026-10-17T09:00:00Z")
    wave=valid("wave4",timestamp="2026-10-17T10:00:00Z")
    assert m.decide(new,wave,wave["id"])["ready"] is True
    assert m.decide(new,wave,new["id"])["ready"] is False


def test_manual_replay_allowed_only_after_both_authentic_source_runs_completed():
    new=valid("new_U100",timestamp="2026-10-17T10:00:00Z")
    wave=valid("wave4",timestamp="2026-10-17T09:00:00Z")
    assert m.decide(new,wave,wave["id"],manual=True)["ready"] is True
    wave["status"]="in_progress"
    wave["conclusion"]=None
    assert m.decide(new,wave,new["id"],manual=True)["ready"] is False


def test_reject_modified_sha_first_attempt_and_wrong_source_repo():
    new=valid("new_U100")
    wave=valid("wave4")
    for field,bad in [("head_sha","incorrect"),("run_attempt",2),("head_branch","wrong")]:
        other=copy.deepcopy(new)
        other[field]=bad
        with pytest.raises(RuntimeError,match="drift"):
            m.decide(other,wave,new["id"])
    other=copy.deepcopy(wave)
    other["head_repository"]={"full_name":"a-fork/fcp"}
    with pytest.raises(RuntimeError,match="different source repository"):
        m.decide(new,other,wave["id"])


def test_reject_wrong_trigger_or_run_name():
    new=valid("new_U100")
    wave=valid("wave4")
    with pytest.raises(RuntimeError,match="did not originate"):
        m.decide(new,wave,38993836606)
    wave["name"]="empty code-only preflight"
    with pytest.raises(RuntimeError,match="drift"):
        m.decide(new,wave,wave["id"])


def test_same_completion_timestamp_tiebreak_is_deterministic():
    new=valid("new_U100")
    wave=valid("wave4")
    bigger=max(new["id"],wave["id"])
    assert m.decide(new,wave,bigger)["ready"] is True
    assert m.decide(new,wave,min(new["id"],wave["id"]))["ready"] is False
