"""FCP AI source24: test no ground-truth / human inflation or case-data leakage."""
import hashlib
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

f=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_fcp_source24_blinded_ai_vs_original_20261010.py"
spec=importlib.util.spec_from_file_location("fcp_blind_source24_cmp",f)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def specimen():
    labels=["white"]*10+["yellow_orange"]*7+["red_pink"]*2+["blue_purple"]*5
    preds=["white"]*4+["yellow_orange"]*2+["red_pink"]+["blue_purple"]+[""]*2+[
        "yellow_orange"]*5+["white"]+["red_pink"]+["red_pink"]*2+["blue_purple"]*5
    assert len(preds)==len(labels)==24
    d=pd.DataFrame({
        "audit_case_id":[f"FCQ-{i:016x}" for i in range(24)],
        "cohort":["discovery"]*8+["validation"]*8+["third"]*8,
        "kind":["high_leverage"]*9+["same_species_same_original_colour_matched"]*9+
              ["original_source_photo_random_reference"]*6,
        "species":[f"Example{i}" for i in range(23)]+["Example22"],
        "original_algorithm_colour_label":labels,
        "assistant_preliminary_four_state":preds,
        "visual_support":["unclear" if p=="" else "clear" for p in preds],
    })
    return d


def test_source_pilot_labels_are_not_called_human_ground_truth():
    a=m.aggregate(specimen())
    assert (a["n_blind_AI_colour_assigned"],a["n_AI_source_colour_agreement"],
            a["n_AI_source_colour_disagreement"],a["n_AI_abstention"])==(22,16,6,2)
    assert a["n_review_trigger"]==8
    assert a["human_expert_reannotations"]==0
    assert a["measured_true_classifier_error_rate"] is None


def test_morph_specific_selected_source_counts_and_abstentions():
    x=m.aggregate(specimen())
    assert x["by_original_pipeline_colour_class"]["white"]=={
        "n":10,"agreement":4,"disagreement":4,"abstention":2,"review_triggers":6}
    assert x["by_original_pipeline_colour_class"]["blue_purple"]["agreement"]==5


def test_published_aggregate_does_not_include_photo_ids_species_or_machine_labels_by_photo():
    a=specimen()
    payload=str(m.aggregate(a))
    for case in a.audit_case_id:
        assert case not in payload
    for species in a.species:
        assert species not in payload
    assert m.AI_SHA in payload and m.KEY_SHA in payload


def test_sha_mismatch_is_fatal():
    with pytest.raises(RuntimeError,match="SHA256 changed"):
        m.check(b"missing/mutated original",m.AI_SHA,"prior blind annotation")


def test_new_ai_or_original_source_labels_cannot_silently_change_count():
    d=specimen()
    d.loc[0,"assistant_preliminary_four_state"]="red_pink"
    with pytest.raises(RuntimeError,match="frozen first-pass"):
        m.aggregate(d)
