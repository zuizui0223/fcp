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


IMMUTABLE_CORRECTED_SOURCE_RECEIPT_SHA256 = (
    "96d5ac8c2650843f47ad28bf3f6f3bf25979c8ebe9ff4dc61da32327ea66fe0b"
)


def verify_corrected_five_cohort_analysis_eligibility(receipt:dict)->dict:
    """Authenticate actual independent-taxonomy eligibility, not partial flower colour."""
    expected={
        "schema":"fcp_4730_original_h9_source_vs_independent_447500_method_eligibility_v1",
        "full_original_U100_frame_taxa":4730,
        "raw_original_source_photo_rows":447800,
        "raw_distinct_photo_IDs_across_all_5_source_cohorts":447797,
        "actual_original_100photo_eligible_taxa_after_cross_cohort_audit":4475,
        "actual_original_independent_100photo_analysis_photo_rows":447500,
        "preserved_original_old_measured_photos":149900,
        "preserved_wave4_original_source_photo_rows":49800,
        "new_U100_source_100photo_taxa":2481,
        "new_U100_independent_analysis_100photo_taxa":2478,
        "under100_original_H9_taxa":252,
        "excluded_new_U100_source_photo_reuse_taxa":[62367,75633,1045161],
        "reused_photo_IDs":[45948761,311814882,517553374],
        "reused_observation_IDs":[287821814],
        "new_photo_pixels_measured_in_this_source_correction":False,
        "previous_H9_source_archive_byte_immutable":True,
        "original_50km_H1_D_W_statistical_method_changed":False,
        "successful_source_correction_run_id":39007962571,
    }
    for key,value in expected.items():
        if receipt.get(key)!=value:
            raise RuntimeError("real 5-cohort source correction receipt drift: "+key)
    return {
        "signed_original_five_cohort_correction_verified":True,
        "original_raw_source_rows":447800,
        "actual_distinct_source_photo_IDs":447797,
        "original_independent_100photo_analysis_taxa":4475,
        "original_independent_100photo_analysis_rows":447500,
        "taxa_rejected_whole_source_photo_reuse":3,
        "raw_photo_IDs_preserved_in_unchanged_running_ROI_workflow":True,
    }


def digest(file: Path) -> str:
    h=hashlib.sha256()
    with file.open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()


def authenticate_original_blind_firewall(root:Path, label:str)->tuple[dict,dict]:
    """Verify the original 256-partition *source assignments*, not only ROI receipts.

    The firewall is the pixel-stage original SHA-pinned manifest, created
    without opening photographs. No biological or geographical fields are
    copied into the audit outputs.
    """
    if label not in TARGET:
        raise RuntimeError("unknown immutable original photo cohort")
    manifest=root/"measurement_firewall_manifest.json"
    frozen=root/"partition_assignments.csv"
    sealed=root/"sealed_keys/metadata_join_key.csv"
    if not all(path.is_file() for path in (manifest,frozen,sealed)):
        raise RuntimeError("original pre-pixel photograph-ID firewall incomplete")
    receipt=json.loads(manifest.read_text())
    schema=("fcp_all4730_original_v4_blind_firewall_v1" if label=="new_U100"
            else "fcp_same_149900_method_wave4_firewall_v1")
    if (receipt.get("schema")!=schema
        or receipt.get("new_photo_rows")!=TARGET[label]["original_photo_ids"]
        or receipt.get("partitions")!=256
        or receipt.get("blind_worker_columns")!=["measurement_id","image_filename","photo_license"]
        or receipt.get("photo_pixels_opened") is not False):
        raise RuntimeError("original source ROI firewall source identity or method drift")
    if label=="new_U100":
        if (receipt.get("source_run_id")!=38996743273
            or receipt.get("new_100photo_species")!=2481
            or receipt.get("coordinates_or_taxa_in_pixel_worker") is not False):
            raise RuntimeError("original FCP full-U100 source identity was changed")
        sha_field="metadata_join_SHA256"
    else:
        if (receipt.get("source_metadata_run_id")!=38989666486
            or receipt.get("new_species")!=498
            or receipt.get("source_photo_urls_in_worker") is not False):
            raise RuntimeError("original wave4 source colour-blind firewall drift")
        sha_field="metadata_join_sha256"
    if digest(sealed)!=receipt.get(sha_field):
        raise RuntimeError("immutable original sealed metadata join SHA drift")
    assignments={}
    seen=set()
    expected_by_batch={0:set(),1:set()}
    with frozen.open(newline="",encoding="utf8") as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or not {"measurement_id","batch","semantic_shard","compute_partition"} <= set(reader.fieldnames):
            raise RuntimeError("original 256-worker pre-pixel assignment schema drift")
        for item in reader:
            mid=item["measurement_id"].strip()
            b,s,p=(int(item[k]) for k in ("batch","semantic_shard","compute_partition"))
            if not mid or mid in seen or not(0<=b<2 and 0<=s<32 and 0<=p<4):
                raise RuntimeError("duplicated original photo assignment/invalid 256-worker index")
            seen.add(mid)
            assignments.setdefault((b,s,p),set()).add(mid)
            expected_by_batch[b].add(mid)
    if len(seen)!=TARGET[label]["original_photo_ids"]:
        raise RuntimeError("original photo firewall assignment denominator changed")
    for b in (0,1):
        packet=root/f"batch_{b}/worker_packet/measurement_manifest.csv"
        if not packet.is_file():
            raise RuntimeError("missing original blind photo worker manifest")
        ids=set()
        with packet.open(newline="",encoding="utf8") as f:
            reader=csv.DictReader(f)
            if reader.fieldnames!=["measurement_id","image_filename","photo_license"]:
                raise RuntimeError("original pixel worker interface no longer blind")
            for item in reader:
                mid=item["measurement_id"].strip()
                if not mid or mid in ids or item["image_filename"]!=mid+".jpg":
                    raise RuntimeError("frozen original measurement manifest ID/filename drift")
                ids.add(mid)
        if ids!=expected_by_batch[b]:
            raise RuntimeError("original blinded worker photo-ID coverage differs from frozen assignment")
    sealed_ids=set()
    with sealed.open(newline="",encoding="utf8") as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or "measurement_id" not in reader.fieldnames:
            raise RuntimeError("original sealed join key invalid")
        for item in reader:
            mid=item["measurement_id"].strip()
            if not mid or mid in sealed_ids:
                raise RuntimeError("original sealed metadata join has duplicate photo IDs")
            sealed_ids.add(mid)
    if sealed_ids!=seen:
        raise RuntimeError("sealed original photo IDs do not match blind worker firewall")
    return assignments,{"frozen_firewall_photo_IDs":len(seen),
                        "frozen_firewall_sha256":digest(manifest),
                        "frozen_partition_assignments_sha256":digest(frozen),
                        "original_256_partition_assignments_validated":True,
                        "each_batch_worker_IDs_match_frozen_firewall":True,
                        "sealed_metadata_join_matches_original_blind_IDs":True}


def inspect_group(root:Path, *, label:str, global_ids:set[str],
                  assignment:dict)->tuple[dict,list[dict]]:
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
        this_part_ids=set()
        frozen_part=assignment.get((b,s,p),set())
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
                if measurement not in frozen_part or measurement in this_part_ids:
                    raise RuntimeError("original pixel terminal ID not assigned to this frozen partition")
                this_part_ids.add(measurement)
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
        if this_part_ids!=frozen_part:
            raise RuntimeError("original ROI completed packet omits a frozen pre-pixel photo ID")
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
    p.add_argument("--new-firewall",type=Path,required=True)
    p.add_argument("--wave4-firewall",type=Path,required=True)
    p.add_argument("--eligibility-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--source-run-records",type=Path,required=True)
    a=p.parse_args()
    if digest(a.eligibility_receipt)!=IMMUTABLE_CORRECTED_SOURCE_RECEIPT_SHA256:
        raise RuntimeError("actual corrected 447500 original eligibility Release SHA256 drift")
    source_correction=verify_corrected_five_cohort_analysis_eligibility(
        json.loads(a.eligibility_receipt.read_text()))
    pinned=json.loads(a.source_run_records.read_text())
    for group in TARGET:
        record=pinned[group]
        if (record.get("id")!=TARGET[group]["run_id"]
            or record.get("run_attempt")!=1
            or record.get("conclusion")=="failure"):
            raise RuntimeError("Untrusted original photographic run preflight")
    new_assign,new_firewall=authenticate_original_blind_firewall(a.new_firewall,"new_U100")
    wave_assign,wave_firewall=authenticate_original_blind_firewall(a.wave4_firewall,"wave4")
    new_all_ids=set().union(*new_assign.values())
    wave_all_ids=set().union(*wave_assign.values())
    if new_all_ids & wave_all_ids:
        raise RuntimeError("new-U100 and wave4 original frozen photo IDs overlap")
    seen=set()
    new,detail_new=inspect_group(a.new,label="new_U100",global_ids=seen,assignment=new_assign)
    wave,detail_wave=inspect_group(a.wave4,label="wave4",global_ids=seen,assignment=wave_assign)
    result={
        "schema":"fcp_4730_genuine_original_roi_partial_terminal_quality_audit_v1",
        "snapshot_is_partial_until_both_256_partition_runs_succeed":True,
        "original_h9_source_frame_taxa":4730,
        # Raw H9 photographic source occurrences are NOT all distinct across
        # the old discovery/new-U100 boundary: 3 photo IDs are repeated.
        # Report the corrected real five-cohort source audit explicitly;
        # actual ROI workers still process their immutable full source input.
        "original_source_photo_IDs_total":447800,
        "raw_five_cohort_source_photo_row_occurrences":447800,
        "original_source_globally_distinct_photo_IDs":447797,
        "original_cross_cohort_photo_ID_reuses":3,
        "original_cross_cohort_observation_ID_reuses":1,
        "independent_original_100photo_analysis_taxa":4475,
        "independent_original_100photo_analysis_photo_rows":447500,
        "new_U100_100photo_taxa_wholly_excluded_from_final_inference":[62367,75633,1045161],
        "excluded_new_U100_photo_rows_for_independence":300,
        "historical_source_photo_analytical_eligibility_ledger_release":
            "fcp-original-4730-independent-447500-analysis-eligibility-20261017",
        "this_quality_audit_does_not_remove_or_relabel_original_ROI_packets":True,
        "original_149900_photo_measurements_not_opened_or_changed":True,
        "new_source_plus_wave4_original_IDs":297900,
        "technical_photo_terminal_audit_only":True,
        "original_species_conditional_50km_or_H1W_not_inferred":True,
        "source_run_ids":{k:v["run_id"] for k,v in TARGET.items()},
        "independent_five_cohort_source_correction":source_correction,
        "independent_five_cohort_release_receipt_sha256":
            IMMUTABLE_CORRECTED_SOURCE_RECEIPT_SHA256,
        "cohorts":{"new_U100":new,"wave4":wave},
        "original_blind_firewall_audits":{"new_U100":new_firewall,"wave4":wave_firewall},
        "all_original_photo_source_blind_firewall_IDs_unique_and_disjoint":True,
        "every_completed_ROI_partition_has_exact_frozen_photo_IDs":True,
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
