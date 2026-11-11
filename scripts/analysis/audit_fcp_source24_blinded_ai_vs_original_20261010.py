#!/usr/bin/env python3
"""Post-outcome, selected-24 source photo AI-vs-algorithm colour concordance.

Only aggregate counts can be shared. Original source labels and AI case-level
colours NEVER go into reviewer queues or committed numerical summaries.
Neither classifier is expert-verified biological ground truth.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile
import pandas as pd

AI_SHA="f232f3e710793ae2ae60821638c73f74cd28af3f3d9851bf628046cb2f8189be"
PILOT_SHA="06191b36eca9862f6241b5110ae13aee0ff7e7414a5caaaa1025daffca33eadd"
KEY_SHA="ef0275fedf7ab4e3e0aa5309fd80c9a9aee97de733b79a93bc2137a04f6f4b5c"
BLIND_SHA="8cc44a71dfb64601e23f68b68c25e8cbc534a332e0a17ff51a51bd8737a35c27"
STATES=("white","yellow_orange","red_pink","blue_purple")
AI_COLUMNS=["audit_case_id","source_image_order","assistant_preliminary_four_state",
            "visual_support","image_only_comment","review_status"]
PILOT_COLUMNS=["audit_case_id","photo_id","frozen_original_image_sha256"]


def check(payload,sha,name):
    if sha256(payload).hexdigest()!=sha:
        raise RuntimeError(f"{name}: original source SHA256 changed")
    return pd.read_csv(BytesIO(payload),dtype=str).fillna("")


def load_original_and_frozen_blind_ai(ai_path,pilot_path,archive_path):
    ai=check(Path(ai_path).read_bytes(),AI_SHA,"previously blinded AI source")
    pilot=check(Path(pilot_path).read_bytes(),PILOT_SHA,"original 24 image panel")
    with ZipFile(archive_path) as archive:
        key=check(archive.read("photo_review/UNBLINDING_KEY_do_not_show_reviewers.csv"),
                  KEY_SHA,"archived source-sealed 405 original label key")
        blind=check(archive.read("photo_review/BLINDED_reannotation_photo_queue.csv"),
                    BLIND_SHA,"original 405 human-review panel")
    if list(ai.columns)!=AI_COLUMNS or list(pilot.columns)!=PILOT_COLUMNS:
        raise RuntimeError("prior to unblinding, source AI/pilot reviewer schema must be immutable")
    if len(ai)!=24 or len(pilot)!=24 or not ai.audit_case_id.equals(pilot.audit_case_id):
        raise RuntimeError("24 pilot image case identity/order drift")
    if ai.audit_case_id.duplicated().any() or pilot.photo_id.duplicated().any():
        raise RuntimeError("source case/photo ID reuse")
    if ai.source_image_order.tolist()!=[str(i) for i in range(1,25)]:
        raise RuntimeError("original blinded contact-sheet order cannot be reconstructed")
    if not ai.review_status.eq("AI_SINGLE_PASS_NOT_HUMAN_ADJUDICATION").all():
        raise RuntimeError("source labels cannot be silently upgraded to human evaluation")
    if not ai.assistant_preliminary_four_state.isin(("",*STATES)).all():
        raise RuntimeError("AI four-colour source labels changed")
    if not ai.visual_support.isin(["clear","moderate","ambiguous","unclear"]).all():
        raise RuntimeError("AI visual reliability classes changed")
    if not ai.loc[ai.assistant_preliminary_four_state.eq(""),"visual_support"].eq("unclear").all():
        raise RuntimeError("AI abstentions must be explicitly unclassifiable")
    if not {"audit_case_id","photo_id","cohort","kind","species",
            "original_algorithm_colour_label"}.issubset(key):
        raise RuntimeError("archived original source colour schema changed")
    if not key.original_algorithm_colour_label.isin(STATES).all():
        raise RuntimeError("original source colour categories changed")
    a=ai.merge(pilot[["audit_case_id","photo_id"]],on="audit_case_id",validate="one_to_one")
    a=a.merge(blind[["audit_case_id","photo_id"]],on=["audit_case_id","photo_id"],
              validate="one_to_one")
    a=a.merge(key[["audit_case_id","photo_id","cohort","kind","species",
                   "original_algorithm_colour_label"]],on=["audit_case_id","photo_id"],
              validate="one_to_one")
    if len(a)!=24 or a.audit_case_id.duplicated().any():
        raise RuntimeError("original 24 cases cannot all be matched to historical source")
    return a


def aggregate(a):
    if len(a)!=24:
        raise RuntimeError("selected AI pilot is not the exact 24 original cases")
    assigned=a.assistant_preliminary_four_state.isin(STATES)
    agree=assigned & a.assistant_preliminary_four_state.eq(a.original_algorithm_colour_label)
    disagree=assigned & ~agree
    abstain=~assigned
    review=disagree|abstain
    counts=(int(agree.sum()),int(disagree.sum()),int(abstain.sum()))
    if counts!=(16,6,2):
        raise RuntimeError("the frozen first-pass AI-screening receipt changed")
    by_source={}
    for state,g in a.groupby("original_algorithm_colour_label"):
        ass=g.assistant_preliminary_four_state.isin(STATES)
        same=ass & g.assistant_preliminary_four_state.eq(state)
        by_source[state]={"n":len(g),"agreement":int(same.sum()),
                         "disagreement":int((ass & ~same).sum()),
                         "abstention":int((~ass).sum()),
                         "review_triggers":int((~same).sum())}
    by_quality={}
    for q,g in a.groupby("visual_support"):
        g_ass=g.assistant_preliminary_four_state.isin(STATES)
        by_quality[q]={"n":len(g),"disagreement":int((g_ass &
            ~g.assistant_preliminary_four_state.eq(g.original_algorithm_colour_label)).sum()),
            "abstention":int((~g_ass).sum())}
    by_design={}
    for name,g in a.groupby("kind"):
        is_agree=g.assistant_preliminary_four_state.eq(g.original_algorithm_colour_label)
        by_design[name]={"n":len(g),"review_triggers":int((~is_agree).sum())}
    return {
        "schema":"fcp_source24_AI_blind_vs_source_pipeline_aggregate_v1",
        "n_selected_source_photographs":24,"n_unique_species":int(a.species.nunique()),
        "n_blind_AI_colour_assigned":int(assigned.sum()),
        "n_AI_source_colour_agreement":counts[0],
        "n_AI_source_colour_disagreement":counts[1],
        "n_AI_abstention":counts[2],
        "n_review_trigger":int(review.sum()),
        "by_original_pipeline_colour_class":by_source,
        "by_blind_AI_visual_support":by_quality,
        "by_sealed_original_selection_group":by_design,
        "checksum_previous_blind_AI_CSV":AI_SHA,
        "checksum_archived_original_405_source_key":KEY_SHA,
        "human_expert_reannotations":0,"measured_true_classifier_error_rate":None,
        "source_biological_photo_colour_codes_changed":False,
        "hard_nonclaims":[
          "Agreement with an AI photo screening is not biological photo-colour truth.",
          "The source high-leverage/controls are a targeted 24-photo sample, not a prevalence denominator.",
          "An unobserved focal plant/organ or background flowering species may cause disagreement.",
          "No species-wide or local 50km colour depletion has been recalculated.",
          "The original per-photo labels must not leak to independent human reviewers."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ai-blind",required=True,type=Path)
    p.add_argument("--pilot24",required=True,type=Path)
    p.add_argument("--original405-artifact",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    z=p.parse_args()
    a=load_original_and_frozen_blind_ai(z.ai_blind,z.pilot24,z.original405_artifact)
    payload=aggregate(a)
    z.out.parent.mkdir(parents=True,exist_ok=True)
    z.out.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2))


if __name__=="__main__":main()
