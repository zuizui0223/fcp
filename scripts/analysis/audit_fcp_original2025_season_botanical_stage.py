#!/usr/bin/env python3
"""Stage/flower-organ QC for ACTUAL original 2025 FCP seasonal image ROIs.

All independent reviewers receive only the provider photo/observation URL and
botanical identity, never predicted morph, season, weather, or observer.
No completed human review means UNKNOWN anthesis, never a negative flower.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from collections import Counter
from pathlib import Path

SOURCE="postblind/first_original_source_matched_season_64_photo_roi_QC.csv"
RECEIPT="postblind/first_original_source_matched_season_64_roi_receipt.json"
TERMINAL="measure/terminal_p000_b0000000.csv"
REVIEW=("measurement_id","photo_id","species","observation_url","photo_url_large",
        "photo_license","reviewer1_id","reviewer1_stage","reviewer1_tissue",
        "reviewer2_id","reviewer2_stage","reviewer2_tissue",
        "third_reviewer_id","adjudicated_stage","adjudicated_tissue","review_notes")
STAGE={"open_flower","buds_only","fruit_only","no_flower_visible","uncertain"}
TISSUE={"petal_tepal","other_open_flower_organ","spikelet_glume_bract",
        "fruit_seed","leaf_stem","background","uncertain"}

def sha(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()

def rows(data:bytes)->list[dict]:
    with io.StringIO(data.decode("utf-8")) as h:
        return list(csv.DictReader(h))

def original(zip_path:Path):
    with zipfile.ZipFile(zip_path) as f:
        receipt=json.loads(f.read(RECEIPT))
        payload=f.read(SOURCE)
        terminal=f.read(TERMINAL)
    if (receipt.get("status")!="FIRST_ORIGINAL_SEASONAL_PIXELS_MEASURED_BOTANICAL_ANTHESIS_NOT_VALIDATED"
        or receipt.get("genuine_original_source_date")!="2026-09-27"
        or receipt.get("postblind_join_sha256")!=sha(payload)
        or receipt.get("terminal_sha256")!=sha(terminal)
        or receipt.get("original_photos_measured")!=64
        or receipt.get("unique_original_observers")!=64
        or receipt.get("photo_id_join_complete") is not True
        or receipt.get("no_genetic_polymorphism_or_flowering_colour_ecological_effect_estimated") is not True):
        raise ValueError("Not sealed genuine original ROI-v4 2025 seasonal photo evidence")
    source=rows(payload)
    measured=rows(terminal)
    if len(source)!=64 or len(measured)!=64 or (
        len({r["measurement_id"] for r in source})!=64 or
        {r["measurement_id"] for r in measured}!={r["measurement_id"] for r in source}
    ):
        raise ValueError("Missing, duplicated or unmatched original ROI photo ID")
    if len({r["observer_id"] for r in source})!=64:
        raise ValueError("Observer-disjoint original source was changed")
    for season in ("summer","autumn"):
        sub=[r for r in source if r["local_season"]==season]
        if len(sub)!=32 or len({r["observation_id"] for r in sub})!=32:
            raise ValueError("Original seasonal source-observation support changed")
        n=sum(r["measurement_status"]=="classified_four_state_morph" for r in sub)
        if n!=receipt["first_work_unit_by_season"][season]["automated_four_state_classifiable"]:
            raise ValueError("Original ROI classifications were altered")
    return source,receipt

def create_blank(source_zip:Path,output:Path):
    if output.exists():raise ValueError("Cannot overwrite botanical review source form")
    source,receipt=original(source_zip)
    with output.open("x",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=REVIEW,lineterminator="\n")
        writer.writeheader()
        for row in sorted(source,key=lambda x:x["measurement_id"]):
            writer.writerow({
                "measurement_id":row["measurement_id"],"photo_id":row["photo_id"],
                "species":row["species"],
                "observation_url":"https://www.inaturalist.org/observations/"+row["observation_id"],
                "photo_url_large":row["photo_url_large"],
                "photo_license":row["photo_license"],
            })
    return {"status":"BLIND_BOTANICAL_REVIEW_FORM_UNFILLED",
            "original_photo_ids":64,"signed_source_postblind_sha256":receipt["postblind_join_sha256"],
            "blank_review_sha256":sha(output.read_bytes()),
            "model_colour_and_season_withheld_from_reviewers":True}

def botanical_state(row:dict)->str:
    s1=row["reviewer1_stage"].strip();s2=row["reviewer2_stage"].strip()
    t1=row["reviewer1_tissue"].strip();t2=row["reviewer2_tissue"].strip()
    u1=row["reviewer1_id"].strip();u2=row["reviewer2_id"].strip()
    third=row["third_reviewer_id"].strip()
    sa=row["adjudicated_stage"].strip();ta=row["adjudicated_tissue"].strip()
    if not any((s1,s2,t1,t2,u1,u2,third,sa,ta)):return "NOT_REVIEWED_UNKNOWN"
    if not all((s1,s2,t1,t2,u1,u2)):
        if any((third,sa,ta)):raise ValueError("Premature third-review evidence")
        return "REVIEW_INCOMPLETE_UNKNOWN"
    if u1==u2 or s1 not in STAGE or s2 not in STAGE or t1 not in TISSUE or t2 not in TISSUE:
        raise ValueError("Invalid independent reviewer ID, anthesis or flower-organ code")
    if (s1,t1)!=(s2,t2):
        if not any((third,sa,ta)):return "REVIEW_CONFLICT_UNRESOLVED"
        if not all((third,sa,ta)) or third in (u1,u2) or sa not in STAGE or ta not in TISSUE:
            raise ValueError("Third independent botanical adjudication required")
        s,t=sa,ta
    else:
        if any((third,sa,ta)):raise ValueError("Adjudication cannot rewrite agreement")
        s,t=s1,t1
    if s=="uncertain" or t=="uncertain":return "REVIEWED_BOTANICAL_UNCERTAINTY"
    if s!="open_flower":return "REVIEWED_NOT_OPEN_FLOWER"
    if t=="petal_tepal":return "REVIEWED_OPEN_PETAL_TEPAL"
    if t=="other_open_flower_organ":return "REVIEWED_OPEN_OTHER_ORGAN_SECONDARY_ONLY"
    return "REVIEWED_NONCOMPARABLE_OR_BACKGROUND_TISSUE"

def audit(source_zip:Path,review_csv:Path,outdir:Path):
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite original botanical review audit")
    source,receipt=original(source_zip)
    photos={r["measurement_id"]:r for r in source}
    with review_csv.open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=REVIEW:
            raise ValueError("Reviewer file schema changed or exposes floral-colour outcomes")
        reviews=list(reader)
    if len(reviews)!=64 or len({r["measurement_id"] for r in reviews})!=64 or {r["measurement_id"] for r in reviews}!=set(photos):
        raise ValueError("Review does not exactly cover 64 original photo IDs")
    summary=Counter();checks=[];verified=Counter()
    for rev in reviews:
        r=photos[rev["measurement_id"]]
        if (rev["photo_id"]!=r["photo_id"] or rev["species"]!=r["species"]
            or rev["photo_url_large"]!=r["photo_url_large"]
            or rev["photo_license"]!=r["photo_license"]
            or rev["observation_url"]!="https://www.inaturalist.org/observations/"+r["observation_id"]):
            raise ValueError("Photo identity/URL was altered after original ROI measurement")
        state=botanical_state(rev)
        eligible=(state=="REVIEWED_OPEN_PETAL_TEPAL"
            and r["measurement_status"]=="classified_four_state_morph"
            and r["same_roi_white_yellow_organ_mixed"].strip().lower()!="true")
        season=r["local_season"]
        verified[season]+=int(eligible)
        summary[(season,state)]+=1
        checks.append({
            "measurement_id":r["measurement_id"],"photo_id":r["photo_id"],
            "season":season,"botanical_review_status":state,
            "automated_colour_classifiable":r["measurement_status"]=="classified_four_state_morph",
            "same_roi_white_yellow_mixed":r["same_roi_white_yellow_organ_mixed"],
            "primary_reviewed_open_flower_photo_colour":eligible,
            "unreviewed_image_not_a_nonflower":True,
        })
    outdir.mkdir(parents=True,exist_ok=True)
    dest=outdir/"summer_autumn_2025_original64_botanical_review_qc.csv"
    with dest.open("x",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=list(checks[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(checks)
    pending=sum(v for (_,s),v in summary.items() if s in {
        "NOT_REVIEWED_UNKNOWN","REVIEW_INCOMPLETE_UNKNOWN",
        "REVIEW_CONFLICT_UNRESOLVED","REVIEWED_BOTANICAL_UNCERTAINTY"})
    report={
        "schema":"fcp_original2025_observer_disjoint_botanical_flower_stage_gate_v1",
        "status":("BOTANICAL_REVIEW_PENDING_NO_WHITE_RATE" if pending else "REVIEWED_BOTANICAL_STAGE_NO_CAUSAL_CLAIM"),
        "original_source_pixels":64,"pending_or_uncertain_review_photos":pending,
        "open_petaloid_classifiable_by_season":{"summer":verified["summer"],"autumn":verified["autumn"]},
        "per_season_stage_outcomes":{f"{a}__{b}":n for (a,b),n in sorted(summary.items())},
        "signed_original_source_sha256":receipt["postblind_join_sha256"],
        "review_sheet_sha256":sha(review_csv.read_bytes()),
        "review_qc_sha256":sha(dest.read_bytes()),
        "reviewers_blinded_to_original_morph_and_season":True,
        "not_reviewed_is_unknown_not_absent_flower":True,
        "model_only_white_label_not_biological_flower_white_rate":True,
    }
    (outdir/"summer_autumn_2025_original64_botanical_review_receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--signed-real-roi-zip",type=Path,required=True)
    p.add_argument("--blank-review",type=Path)
    p.add_argument("--review-sheet",type=Path)
    p.add_argument("--outdir",type=Path)
    a=p.parse_args()
    if a.blank_review and not a.review_sheet and not a.outdir:
        result=create_blank(a.signed_real_roi_zip,a.blank_review)
    elif a.review_sheet and a.outdir and not a.blank_review:
        result=audit(a.signed_real_roi_zip,a.review_sheet,a.outdir)
    else:
        p.error("Choose --blank-review alone OR --review-sheet together with --outdir")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
