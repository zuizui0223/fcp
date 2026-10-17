"""Protocol-only synthetic tests; no source photos or ecological estimators."""
from __future__ import annotations

import csv
import json

import pytest

from scripts.analysis import fcp_4730_real_roi_partial_packet_audit_20261017 as m


def csvout(path,rows,fieldnames):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames,lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fake_original_source_firewall(tmp_path,monkeypatch,label="new_U100"):
    monkeypatch.setitem(m.TARGET[label],"original_photo_ids",3)
    root=tmp_path/("src_"+label)
    ids=["FCPR-AAA","FCPR-BBB","FCPR-CCC"]
    assignments=[
        {"measurement_id":ids[0],"batch":0,"semantic_shard":0,"compute_partition":0},
        {"measurement_id":ids[1],"batch":0,"semantic_shard":0,"compute_partition":0},
        {"measurement_id":ids[2],"batch":1,"semantic_shard":1,"compute_partition":1},
    ]
    csvout(root/"partition_assignments.csv",assignments,
           ["measurement_id","batch","semantic_shard","compute_partition"])
    for b in [0,1]:
        csvout(root/f"batch_{b}/worker_packet/measurement_manifest.csv",
               [{"measurement_id":a["measurement_id"],"image_filename":a["measurement_id"]+".jpg",
                 "photo_license":"cc-by"} for a in assignments if a["batch"]==b],
               ["measurement_id","image_filename","photo_license"])
    join=root/"sealed_keys/metadata_join_key.csv"
    csvout(join,[{"measurement_id":i,"inat_taxon_id":f"PRIVATE_{i}"} for i in ids],
           ["measurement_id","inat_taxon_id"])
    record={
        "schema":("fcp_all4730_original_v4_blind_firewall_v1" if label=="new_U100" else
                  "fcp_same_149900_method_wave4_firewall_v1"),
        "new_photo_rows":3,"partitions":256,"blind_worker_columns":[
           "measurement_id","image_filename","photo_license"],
        "photo_pixels_opened":False,
    }
    if label=="new_U100":
        record.update({"source_run_id":38996743273,"new_100photo_species":2481,
                       "coordinates_or_taxa_in_pixel_worker":False,
                       "metadata_join_SHA256":m.digest(join)})
    else:
        record.update({"source_metadata_run_id":38989666486,"new_species":498,
                       "source_photo_urls_in_worker":False,
                       "metadata_join_sha256":m.digest(join)})
    (root/"measurement_firewall_manifest.json").write_text(json.dumps(record))
    return root,ids


def terminal(tmp_path,ids,*,bad=None):
    folder=tmp_path/"sealed"
    raw=[
        {"measurement_id":ids[0],"morph":"white",
         "measurement_status":"classified_four_state_morph"},
        {"measurement_id":ids[1] if bad is None else bad,"morph":"mixed_uncertain",
         "measurement_status":"not_evaluable_roi_or_flip_gate"},
    ]
    f=folder/"b0_partition_s00_p00.csv"
    csvout(f,raw,["measurement_id","morph","measurement_status"])
    (folder/"b0_partition_s00_p00.json").write_text(json.dumps({
       "status":"complete_random_photo_first_terminal_partition",
       "semantic_shard":0,"compute_partition":0,"selected_rows":2,
       "terminal_rows":2,"classified_rows":1,"mixed_uncertain_rows":1,
       "acquisition_failed_rows":0,
       "source_urls_present":False,"species_present":False,
       "coordinates_present":False,
    }))
    return folder


@pytest.mark.parametrize("label",["new_U100","wave4"])
def test_actual_terminal_photo_ids_match_source_partition_not_just_row_count(tmp_path,monkeypatch,label):
    root,ids=fake_original_source_firewall(tmp_path,monkeypatch,label)
    mapping,audit=m.authenticate_original_blind_firewall(root,label)
    assert sum(len(group) for group in mapping.values())==3
    assert audit["each_batch_worker_IDs_match_frozen_firewall"] is True
    folder=terminal(tmp_path,ids)
    summary,rows=m.inspect_group(folder,label=label,global_ids=set(),assignment=mapping)
    assert summary["completed_partitions"]==1
    assert summary["terminal_photo_rows"]==2
    assert summary["classified_four_state_rows"]==1
    assert summary["source_photo_IDs_not_yet_in_completed_partitions"]==1
    assert len(rows)==1


def test_forged_original_roi_blind_id_fails_even_when_csv_receipt_counts_match(tmp_path,monkeypatch):
    root,ids=fake_original_source_firewall(tmp_path,monkeypatch)
    mapping,_=m.authenticate_original_blind_firewall(root,"new_U100")
    folder=terminal(tmp_path,ids,bad="FCPR-UNAUTHENTIC")
    with pytest.raises(RuntimeError,match="not assigned to this frozen partition"):
        m.inspect_group(folder,label="new_U100",global_ids=set(),assignment=mapping)


def test_original_roi_packet_cannot_omit_one_frozen_photo_id(tmp_path,monkeypatch):
    root,ids=fake_original_source_firewall(tmp_path,monkeypatch)
    mapping,_=m.authenticate_original_blind_firewall(root,"new_U100")
    folder=terminal(tmp_path,ids)
    valid=folder/"b0_partition_s00_p00.csv"
    lines=valid.read_text().splitlines()
    valid.write_text("\n".join(lines[:2])+"\n")
    data=json.loads((folder/"b0_partition_s00_p00.json").read_text())
    data["terminal_rows"]=1
    data["selected_rows"]=1
    data["mixed_uncertain_rows"]=0
    (folder/"b0_partition_s00_p00.json").write_text(json.dumps(data))
    with pytest.raises(RuntimeError,match="omits a frozen"):
        m.inspect_group(folder,label="new_U100",global_ids=set(),assignment=mapping)


def test_original_photo_source_assigned_only_once_and_blind(tmp_path,monkeypatch):
    root,_=fake_original_source_firewall(tmp_path,monkeypatch)
    packed=root/"batch_0/worker_packet/measurement_manifest.csv"
    packed.write_text(packed.read_text().replace("FCPR-BBB", "FCPR-AAA"))
    with pytest.raises(RuntimeError,match="filename drift"):
        m.authenticate_original_blind_firewall(root,"new_U100")


def test_original_sealed_geo_metadata_must_match_same_pre_pixel_ids(tmp_path,monkeypatch):
    root,_=fake_original_source_firewall(tmp_path,monkeypatch)
    record=root/"measurement_firewall_manifest.json"
    data=json.loads(record.read_text())
    data["metadata_join_SHA256"]="bad"
    record.write_text(json.dumps(data))
    with pytest.raises(RuntimeError,match="SHA drift"):
        m.authenticate_original_blind_firewall(root,"new_U100")
