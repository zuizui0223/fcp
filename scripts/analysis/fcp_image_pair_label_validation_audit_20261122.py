#!/usr/bin/env python3
"""Audit the 128-pair exploratory visual-organ correction and produce blind review sheets.

This is a *source-table replay* of one AI-reader pilot, not independent
botanical or original-photo validation. It does not re-label the 447,500 photos.

Run:
 python scripts/analysis/fcp_image_pair_label_validation_audit_20261122.py \
   --out-dir /tmp/fcp-image-validation --bootstrap 5000

Source: manuscripts/climate_colour/Table_S11_pair_audit.csv
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "manuscripts/climate_colour/Table_S11_pair_audit.csv"
# Approximate full-pair-frame contrasts as disclosed in Supporting Information.
# Use only for an anchored difference-in-differences replay, not an estimate
# from the selected 128 pairs.
FULL_FRAME_WHITE_DIFF_PP = 2.63
FULL_FRAME_MISMATCH_DIFF_PP = -0.01
SEED = 20261122
COLOURS = {"white", "yellow_orange", "red_pink", "blue_purple", "unknown"}
READABILITY = {"yes", "no", "uncertain"}


def require(ok, msg):
    if not ok:
        raise AssertionError("FCP_PAIR_VISUAL_AUDIT: " + str(msg))


def mean(values, weights):
    total = sum(weights)
    require(total > 0, "empty or zero-weight denominator")
    return sum(x * w for x, w in zip(values, weights)) / total


def load():
    with SOURCE.open(encoding="utf-8", newline="") as f:
        src = list(csv.DictReader(f))
    require(len(src) == 128, "expected exactly 128 visual-pair audit draws")
    require([r["id"] for r in src] == [f"P{i:03d}" for i in range(1, 129)],
            "unexpected pair IDs or ordering")
    counts = Counter((int(r["tail"]), int(r["original_pair_mismatch"])) for r in src)
    require(counts == {(-1, 0): 32, (-1, 1): 32, (1, 0): 32, (1, 1): 32},
            "2x2 thermal-tail by original-agreement strata changed")
    seen = set()
    for r in src:
        require(r["valid"] in ("True", "False"), "reader status must be explicit")
        require(r["tail"] in ("-1", "1"), "unrecognized thermal tail")
        require(r["original_pair_mismatch"] in ("0", "1"), "bad original pair code")
        require(int(r["oa"]) != int(r["ob"]), "same original photographer")
        for key in ("left_photo_id", "right_photo_id"):
            require(int(r[key]) > 0, "invalid original source photo id")
            seen.add(int(r[key]))
        require(0 < float(r["stratum_weight"]) < 1, "invalid support weight")
        require(float(r["w"]) > 0, "invalid source frame pair weight")
        if r["valid"] == "True":
            require(r["visdiff"] in ("True", "False"), "readable pair has no visual mismatch")
            require(r["viswhite"] in ("0.0", "0.5", "1.0"),
                    "readable pair needs visual white share")
        else:
            require(r["visdiff"] == "" and r["viswhite"] == "",
                    "unreadable status must not be used as a corrected zero")
    require(len(seen) == 250, "unexpected photo ID duplication structure")
    require(len({r["inat_taxon_id"] for r in src}) == 116,
            "unexpected number of represented species")
    return src, counts


def status(r):
    return r["valid"] == "True"


def original_white(r):
    return float(r["white"])


def corrected_white(r):
    return float(r["viswhite"]) if status(r) else original_white(r)


def original_mismatch(r):
    return float(r["original_pair_mismatch"])


def corrected_mismatch(r):
    return float(r["visdiff"] == "True") if status(r) else original_mismatch(r)


def tail_contrast(rows, outcome, readable_only=False):
    out = {}
    for tail in (-1, 1):
        values = [r for r in rows if int(r["tail"]) == tail and
                  (not readable_only or status(r))]
        weights = [float(r["stratum_weight"]) for r in values]
        out[tail] = mean([outcome(r) for r in values], weights)
    return 100 * (out[1] - out[-1])


def contrasts(rows):
    orig_white = tail_contrast(rows, original_white)
    visual_white = tail_contrast(rows, corrected_white)
    orig_diff = tail_contrast(rows, original_mismatch)
    visual_diff = tail_contrast(rows, corrected_mismatch)
    change_white = visual_white - orig_white
    change_diff = visual_diff - orig_diff
    return {
        "audit_original_white_pp": orig_white,
        "audit_unreadable_imputed_white_pp": visual_white,
        "audit_original_mismatch_pp": orig_diff,
        "audit_unreadable_imputed_mismatch_pp": visual_diff,
        "visual_correction_white_pp": change_white,
        "visual_correction_mismatch_pp": change_diff,
        "full_frame_original_white_pp_rounded": FULL_FRAME_WHITE_DIFF_PP,
        "full_frame_original_mismatch_pp_rounded": FULL_FRAME_MISMATCH_DIFF_PP,
        "full_frame_corrected_white_pp_approx": FULL_FRAME_WHITE_DIFF_PP + change_white,
        "full_frame_corrected_mismatch_pp_approx": FULL_FRAME_MISMATCH_DIFF_PP + change_diff,
        "readable_only_white_pp": tail_contrast(rows, corrected_white, True),
        "readable_only_mismatch_pp": tail_contrast(rows, corrected_mismatch, True),
    }


def wilson(nsuccess, ntotal, z=1.959963984540054):
    p = nsuccess / ntotal
    z2 = z * z
    center = (p + z2 / (2 * ntotal)) / (1 + z2 / ntotal)
    h = z * math.sqrt(p*(1-p)/ntotal + z2/(4*ntotal**2))/(1+z2/ntotal)
    return [center-h, center+h]


def bootstrap(rows, n=5000):
    rng = random.Random(SEED)
    strata = {}
    for r in rows:
        strata.setdefault((r["tail"], r["original_pair_mismatch"]), []).append(r)
    require(len(strata) == 4 and all(len(v) == 32 for v in strata.values()),
            "strata disrupted")
    white, mismatch = [], []
    ordered = [strata[k] for k in sorted(strata)]
    for i in range(n):
        draw = [rng.choice(group) for group in ordered
                for _ in range(len(group))]
        c = contrasts(draw)
        white.append(c["full_frame_corrected_white_pp_approx"])
        mismatch.append(c["full_frame_corrected_mismatch_pp_approx"])
    def pct(v, p):
        v = sorted(v)
        i = (len(v)-1)*p
        lo = int(math.floor(i))
        hi = min(lo+1, len(v)-1)
        return v[lo]+(i-lo)*(v[hi]-v[lo])
    return {
        "reps": n, "seed": SEED,
        "stratification": "fixed 32 draws in each of four source audit strata",
        "ci_scope": "fixed full-frame rounded reference and exploratory AI-reader label; not photo-level fitting / independent validation",
        "corrected_white_pp_95": [pct(white,.025),pct(white,.975)],
        "corrected_mismatch_pp_95": [pct(mismatch,.025),pct(mismatch,.975)],
    }


def reviewer_packet(rows, out):
    """No thermal side, original colour, provisional visual score in reviewer sheets."""
    rng = random.Random(SEED)
    ids = list(range(len(rows)))
    rng.shuffle(ids)
    # Flip order once per underlying pair, shared by both reviewers,
    # so their scoring can be compared without mapping left/right.
    flip = {k: rng.choice([False, True]) for k in range(len(rows))}
    fields = ["reviewer", "blind_pair_id", "photo_A_url", "photo_B_url",
              "photo_A_open_corolla_visible", "photo_B_open_corolla_visible",
              "photo_A_primary_open_petal_colour", "photo_B_primary_open_petal_colour",
              "same_homologous_open_floral_organ", "same_flower_stage",
              "multiple_flowers_or_stages_mixed", "bud_fruit_sepal_dominance",
              "image_quality_sufficient", "confidence_low_medium_high",
              "reviewer_comment"]
    originals = []
    audit_id = {}
    for i in range(len(rows)):
        # Stable nonsequential review IDs. The key file is held separate.
        audit_id[i] = hashlib.sha256((str(SEED)+":"+rows[i]["id"]).encode()).hexdigest()[:10]
    for reviewer in ("A", "B"):
        with (out / f"reviewer_{reviewer}_blind_128_pairs.csv").open(
                "w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=fields)
            wr.writeheader()
            for i in ids if reviewer == "A" else list(reversed(ids)):
                r = rows[i]
                lhs, rhs = r["left_photo_id"], r["right_photo_id"]
                if flip[i]:
                    lhs, rhs = rhs, lhs
                wr.writerow({
                    "reviewer": reviewer, "blind_pair_id": audit_id[i],
                    "photo_A_url": f"https://www.inaturalist.org/photos/{lhs}",
                    "photo_B_url": f"https://www.inaturalist.org/photos/{rhs}",
                })
    # Keep this decoding key out of blinded reviewer packets.
    with (out / "UNBLIND_AFTER_TWO_REVIEWS_original_128_pair_key.csv").open(
            "w", newline="", encoding="utf-8") as f:
        fields_key = ["blind_pair_id", "original_source_pair_id",
                      "inat_taxon_id", "left_photo_id", "right_photo_id",
                      "pair_A_is_original_left", "thermal_tail",
                      "original_pair_mismatch", "original_colour_left",
                      "original_colour_right", "prior_AI_reader_code",
                      "prior_AI_readable", "source_sampling_weight",
                      "audit_stratum_weight"]
        wr=csv.DictWriter(f, fieldnames=fields_key)
        wr.writeheader()
        for i,r in enumerate(rows):
            wr.writerow({
                "blind_pair_id":audit_id[i], "original_source_pair_id":r["id"],
                "inat_taxon_id":r["inat_taxon_id"],
                "left_photo_id":r["left_photo_id"], "right_photo_id":r["right_photo_id"],
                "pair_A_is_original_left":not flip[i], "thermal_tail":r["tail"],
                "original_pair_mismatch":r["original_pair_mismatch"],
                "original_colour_left":r["ca"], "original_colour_right":r["cb"],
                "prior_AI_reader_code":r["blind_score"],
                "prior_AI_readable":r["valid"],
                "source_sampling_weight":r["w"],
                "audit_stratum_weight":r["stratum_weight"],
            })
    return {"n_queue_pairs":len(rows),
            "reviewers":"two separate human botanical adjudicators REQUIRED; no annotations provided",
            "source_pages_not_copied":"links only; verify current photo access and credits at review time",
            "scoring_protocol":{
                "readable_pair":"both open floral corollas visibly resolved, matched homologous organ and stage; otherwise no",
                "valid_categories":sorted(COLOURS),
                "same_homologous_open_floral_organ":"yes/no/uncertain",
                "same_flower_stage":"yes/no/uncertain",
                "confidence":"low/medium/high",
                "rule":"do not infer colours from source label or metadata; review images independently; collect uncertainty as missing, never assign false agreement",
                "adjudication":"blind independent scores first, reconcile discordance with a third botanist; retain original raw scores"
            }}


def main(out, reps):
    rows, strata = load()
    point = contrasts(rows)
    require(128 == len(rows), "wrong total rows")
    nvalid = sum(status(r) for r in rows)
    require(nvalid == 66 and 128 - nvalid == 62, "pilot reader completeness changed")
    same = [r for r in rows if original_mismatch(r) == 0 and status(r)]
    different = [r for r in rows if original_mismatch(r) == 1 and status(r)]
    require(len(same) == 34 and len(different) == 32, "readable stratum counts drifted")
    remain_same = sum(r["visdiff"] == "False" for r in same)
    remain_diff = sum(r["visdiff"] == "True" for r in different)
    require(remain_same == 32 and remain_diff == 3,
            "AI reader agreement data do not match S11")
    require(abs(point["full_frame_corrected_white_pp_approx"] - .37) < .015,
            "white delta does not reproduce manuscript rounded corrected point")
    require(abs(point["full_frame_corrected_mismatch_pp_approx"] - (-1.74)) < .02,
            "mismatch delta does not reproduce rounded corrected point")
    require(abs(point["readable_only_white_pp"] - (-2.64)) < .01,
            "readable white contrast mismatch")
    require(abs(point["readable_only_mismatch_pp"] - 9.44) < .01,
            "readable disagreement contrast mismatch")
    out.mkdir(parents=True, exist_ok=True)
    forms = reviewer_packet(rows, out)
    receipt = {
        "schema": "fcp_visual_organ_proxy_pilot_quality_gate_v1",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "status": "pilot_archive_replay_and_BLANK_two_human_review_forms",
        "n_original_photographed_flower_pairs":128,
        "n_distinct_species":116,
        "n_distinct_photo_ids":250,
        "strata": {f"tail_{k[0]}_original_mismatch_{k[1]}": n for k, n in strata.items()},
        "pilot_n_readable":nvalid,
        "pilot_n_unreadable":128-nvalid,
        "pilot_readable_rate":nvalid/128,
        "pilot_readable_rate_wilson95":wilson(nvalid,128),
        "pilot_different_label_readable":len(different),
        "pilot_different_label_still_different":remain_diff,
        "pilot_different_label_still_different_wilson95":wilson(remain_diff,len(different)),
        "pilot_same_label_readable":len(same),
        "pilot_same_label_still_same":remain_same,
        "pilot_same_label_still_same_wilson95":wilson(remain_same,len(same)),
        "point_replay":point,
        "bootstrap":bootstrap(rows,reps),
        "review_packet":forms,
        "limitations":[
            "The original pilot is a single AI visual reader, not botanically independent ground truth.",
            "Source 128 draws sampled with replacement, 116 taxa/250 distinct photo IDs, not a global error prevalence sample.",
            "The reader could not resolve 62/128 pairs; the unreadable-is-original assumption is just one sensitivity.",
            "Full-frame original 2.63 pp / -0.01 pp are rounded published anchors; replayed corrected points are approximate.",
            "The target is matched open floral organs, not genetically inherited colour morphs or pollinator fitness.",
            "This thermal-tail frame cannot calibrate the global 3,600-species spatial coefficient without an independent geographically sampled audit."
        ]
    }
    (out/"visual_organ_pilot_replay_and_validation_gate.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({k:receipt[k] for k in [
        "n_original_photographed_flower_pairs","n_distinct_species",
        "n_distinct_photo_ids","pilot_n_readable","pilot_n_unreadable",
        "pilot_different_label_still_different","pilot_same_label_still_same",
        "point_replay","bootstrap"]}, indent=2))


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=ROOT/"output/fcp_image_validation")
    p.add_argument("--bootstrap",type=int,default=5000)
    args=p.parse_args()
    main(args.out_dir,args.bootstrap)
