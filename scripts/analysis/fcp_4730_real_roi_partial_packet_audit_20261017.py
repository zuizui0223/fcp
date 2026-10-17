#!/usr/bin/env python3
"""Audit COMPLETED genuine FCP original ROI-v4 partition terminals mid-run.

This is technical measurement quality control only: no new photo classification,
no source-photo replacement, no cohort combination and no ecological 50-km/W
significance test from incomplete data. The original sealed terminal packets
contain blinded IDs and colours, never taxon/location/source URLs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

TARGET = {
    "new_U100": {"run_id":38997169400,"original_photo_ids":248100,"expected_partitions":256},
    "wave4": {"run_id":38992543115,"original_photo_ids":49800,"expected_partitions":256},
}
ALLOWED = {"white", "yellow_orange", "red_pink", "blue_purple"}
PATTERN = re.compile(r"b([01])_partition_s(\d{2})_p(\d{2})\.json")
FORBIDDEN = {"photo_id","inat_taxon_id","species","observer_id","latitude",
             "longitude","photo_url_large","query_species","observation_id"}


def digest(file: Path) -> str:
    h=hashlib.sha256()
    with file.open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()


def inspect_group(root:Path, *, label:str, global_ids:set[str])->tuple[dict,list[dict]]:
    if label not in TARGET:
        raise RuntimeError("Unrecognised frozen original ROI cohort")
    receipts=list(root.rglob("b*_partition_s*_p*.json"))
    if not receipts:
        raise RuntimeError(f"{label}: no completed original sealed ROI receipts")
    summary={"cohort":label,"source_run_id":TARGET[label]["run_id"],
             "total_original_source_photo_ids":TARGET[label]["original_photo_ids"],
             "completed_partitions":0,"possible_partitions":256,
             "terminal_photo_rows":0,"classified_four_state_rows":0,
             "mixed_uncertain_rows":0,"acquisition_failure_rows":0}
    seen_parts=set()
    detail=[]
    for js in sorted(receipts):
        match=PATTERN.fullmatch(js.name)
        if match is None:raise RuntimeError("Invalid original partition basename "+js.name)
        b,s,p=(int(x) for x in match.groups())
        if not (0<=b<2 and 0<=s<32 and 0<=p<4):
            raise RuntimeError("Original frozen semantic partition outside 256")
        if (b,s,p) in seen_parts:
            raise RuntimeError("Duplicate sealed photo partition")
        seen_parts.add((b,s,p))
        csv_file=js.with_suffix(".csv")
        if not csv_file.is_file():
            raise RuntimeError("Missing original terminal CSV for "+js.name)
        receipt=json.loads(js.read_text())
        if (receipt.get("status")!="complete_random_photo_first_terminal_partition"
            or receipt.get("semantic_shard")!=s
            or receipt.get("compute_partition")!=p
            or receipt.get("source_urls_present") is not False
            or receipt.get("species_present") is not False
            or receipt.get("coordinates_present") is not False):
            raise RuntimeError("Historical blind original colour terminal contract violated")
        total=classified=uncertain=failed=0
        with csv_file.open(newline="",encoding="utf8") as f:
            reader=csv.DictReader(f)
            if not reader.fieldnames or not {"measurement_id","morph","measurement_status"} <= set(reader.fieldnames):
                raise RuntimeError("Original real ROI CSV does not have required fields")
            if FORBIDDEN & set(reader.fieldnames):
                raise RuntimeError("Original blind terminal leaked source metadata")
            for row in reader:
                measurement=row["measurement_id"].strip()
                if not measurement or measurement in global_ids:
                    raise RuntimeError("Missing/duplicated original blind measurement ID")
                global_ids.add(measurement)
                total+=1
                if row["morph"] in ALLOWED:
                    if row["measurement_status"]!="classified_four_state_morph":
                        raise RuntimeError("Colour label not emitted by original ROI classifier")
                    classified+=1
                elif row["morph"]=="mixed_uncertain":
                    uncertain+=1
                else:
                    raise RuntimeError("Unrecognised original four-state/mixed morph")
                failed+=(row["measurement_status"]=="image_acquisition_failed")
        if (total!=receipt.get("terminal_rows")
            or total!=receipt.get("selected_rows")
            or classified!=receipt.get("classified_rows")
            or uncertain!=receipt.get("mixed_uncertain_rows")
            or failed!=receipt.get("acquisition_failed_rows")
            or total!=classified+uncertain):
            raise RuntimeError("Original sealed partition terminal count disagrees with actual pixels")
        summary["completed_partitions"]+=1
        summary["terminal_photo_rows"]+=total
        summary["classified_four_state_rows"]+=classified
        summary["mixed_uncertain_rows"]+=uncertain
        summary["acquisition_failure_rows"]+=failed
        detail.append({"original_batch":b,"semantic_shard":s,"compute_partition":p,
                       "terminal_photo_rows":total,
                       "classified_four_state_rows":classified,
                       "mixed_uncertain_rows":uncertain,
                       "acquisition_failure_rows":failed,
                       "csv_sha256":digest(csv_file),"receipt_sha256":digest(js)})
    if summary["completed_partitions"]>256 or summary["terminal_photo_rows"]>summary["total_original_source_photo_ids"]:
        raise RuntimeError("Original photo source denominator exceeded")
    if summary["classified_four_state_rows"]+summary["mixed_uncertain_rows"]!=summary["terminal_photo_rows"]:
        raise RuntimeError("Mis-accounted biological and technical original photograph terminals")
    summary["classifiable_fraction_of_terminal_rows"]=(
        summary["classified_four_state_rows"]/summary["terminal_photo_rows"])
    summary["source_photo_IDs_not_yet_in_completed_partitions"]=(
        summary["total_original_source_photo_ids"]-summary["terminal_photo_rows"])
    summary["technical_assay_only_no_inference"]=True
    return summary,detail


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--new",type=Path,required=True)
    p.add_argument("--wave4",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--source-run-records",type=Path,required=True)
    a=p.parse_args()
    pinned=json.loads(a.source_run_records.read_text())
    for group in TARGET:
        record=pinned[group]
        if (record.get("id")!=TARGET[group]["run_id"]
            or record.get("run_attempt")!=1
            or record.get("conclusion")=="failure"):
            raise RuntimeError("Untrusted original photographic run preflight")
    seen=set()
    new,detail_new=inspect_group(a.new,label="new_U100",global_ids=seen)
    wave,detail_wave=inspect_group(a.wave4,label="wave4",global_ids=seen)
    result={
        "schema":"fcp_4730_genuine_original_roi_partial_terminal_quality_audit_v1",
        "snapshot_is_partial_until_both_256_partition_runs_succeed":True,
        "original_h9_source_frame_taxa":4730,
        "original_source_photo_IDs_total":447800,
        "original_149900_photo_measurements_not_opened_or_changed":True,
        "new_source_plus_wave4_original_IDs":297900,
        "technical_photo_terminal_audit_only":True,
        "original_species_conditional_50km_or_H1W_not_inferred":True,
        "source_run_ids":{k:v["run_id"] for k,v in TARGET.items()},
        "cohorts":{"new_U100":new,"wave4":wave},
        "both_new_cohorts_completed_unique_blind_measurement_IDs":len(seen),
        "new_plus_wave4_classified_four_state_rows":new["classified_four_state_rows"]+wave["classified_four_state_rows"],
        "new_plus_wave4_mixed_uncertain_rows":new["mixed_uncertain_rows"]+wave["mixed_uncertain_rows"],
        "new_plus_wave4_image_acquisition_failure_rows":new["acquisition_failure_rows"]+wave["acquisition_failure_rows"],
        "cross_cohort_duplicate_blind_measurement_ids":0,
    }
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"real_original_partial_roi_quality.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    with (a.output/"real_original_partial_roi_partition_ledger.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["cohort",*detail_new[0].keys()],lineterminator="\n")
        writer.writeheader()
        for label,details in (("new_U100",detail_new),("wave4",detail_wave)):
            for item in details:
                writer.writerow({"cohort":label,**item})
    print("FCP_ALL4730_REAL_ROI_ACTUAL_COMPLETED_TERMINAL_QUALITY",
          json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
