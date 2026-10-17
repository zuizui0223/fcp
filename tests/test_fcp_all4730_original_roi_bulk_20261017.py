"""No actual florals are analysed by these synthetic pipeline checks."""
import hashlib
import json

import pandas as pd
import pytest

from scripts.analysis import fcp_all4730_original_roi_bulk_20261017 as m


def make_source(tmp):
    source=tmp/"original_H9_source"
    source.mkdir()
    records=[]
    for i in range(100):
        records.append({
            "inat_taxon_id":12345678,
            "species":"Testus species",
            "query_species":"Testus species",
            "photo_id":90000000+i,
            "observation_id":5000000+i,
            "photo_url_large":f"https://example.org/photo/{i}.jpg",
            "photo_license":"cc-by",
            "latitude":3.2,"longitude":131.1,"observer_id":f"testperson{i}",
        })
    pics=source/"remaining_U100_authorized_metadata.csv.gz"
    pd.DataFrame(records).to_csv(pics,index=False,lineterminator="\n",
       compression={"method":"gzip","mtime":0})
    status=[]
    for i in range(4730):
        code=("original_149900_photo_measured" if i<1499 else
              "wave4_49800_photo_id_authorized" if i<1997 else
              "new_metadata_100_of_100_photo_IDs" if i==1997 else
              "new_H9_metadata_under_100_or_request_error")
        status.append({"inat_taxon_id":12345678 if i==1997 else i,
                       "u100_status":"synthetic",
                       "final_status":code})
    ledger=source/"full_U100_4730_status_after_H9.csv"
    pd.DataFrame(status).to_csv(ledger,index=False,lineterminator="\n")
    receipt={
        "schema":"fcp_u100_all4730_original_h9_only_v1",
        "stage":"all_4730_H9_metadata_terminal",
        "U100_identity_species":4730,
        "old_measured_species":1499,
        "wave4_original_photo_metadata_species":498,
        "new_attempted_species":2733,
        "new_complete_100_photo_species":1,
        "new_unmet_100_photo_species":2732,
        "new_original_photo_IDs":100,
        "verified_against_all_4_complete_shards":True,
        "metadata_gate_pass":True,
        "old_149900_method_and_results_changed":False,
        "photo_pixels_opened":False,
        "outputs_sha256":{
            pics.name:m.sha(pics),
            ledger.name:m.sha(ledger),
        },
    }
    (source/"all4730_metadata_result.json").write_text(json.dumps(receipt))
    return source


def test_all_species_complete_source_gate_ignores_climate_and_original_morph(tmp_path):
    source=make_source(tmp_path)
    meta,photo=m.authenticated_full_u100(source)
    assert len(photo)==100 and meta["new_complete_100_photo_species"]==1
    assert photo.inat_taxon_id.nunique()==1


def test_original_three_field_blind_worker_permits_no_species_coordinates_or_photo_urls(tmp_path):
    source=make_source(tmp_path)
    r=m.firewall(source,tmp_path/"blind")
    assert r["new_photo_rows"]==100
    assert r["partitions"]==256
    assert not r["photo_pixels_opened"]
    assert r["blind_worker_columns"]==["measurement_id","image_filename","photo_license"]
    for batch in (0,1):
        manifest=tmp_path/f"blind/batch_{batch}/worker_packet/measurement_manifest.csv"
        assert list(pd.read_csv(manifest).columns)==["measurement_id","image_filename","photo_license"]
        for forbidden in ("photo_url_large","latitude","longitude","inat_taxon_id","observer_id"):
            assert forbidden not in manifest.read_text()


def test_missing_exact_original_256_roi_partitions_fails_closed(tmp_path):
    source=make_source(tmp_path)
    m.firewall(source,tmp_path/"blind")
    (tmp_path/"parts").mkdir()
    with pytest.raises(RuntimeError,match="256 original classifier partitions"):
        m.reassemble(tmp_path/"blind",tmp_path/"parts",tmp_path/"result")


def test_source_shas_are_always_verified_before_pixels(tmp_path):
    source=make_source(tmp_path)
    with (source/"all4730_metadata_result.json").open() as f:
        r=json.load(f)
    r["outputs_sha256"]["remaining_U100_authorized_metadata.csv.gz"]="bad"
    (source/"all4730_metadata_result.json").write_text(json.dumps(r))
    with pytest.raises(RuntimeError,match="SHA drift"):
        m.authenticated_full_u100(source)


def test_new_original_h9_metadata_must_pass_completeness_gate(tmp_path):
    source=make_source(tmp_path)
    r=json.loads((source/"all4730_metadata_result.json").read_text())
    r["metadata_gate_pass"]=False
    (source/"all4730_metadata_result.json").write_text(json.dumps(r))
    with pytest.raises(RuntimeError,match="unverified"):
        m.firewall(source,tmp_path/"blind")


def test_original_inference_cannot_run_without_real_roi_terminal(tmp_path):
    source=make_source(tmp_path)
    pictures=source/"remaining_U100_authorized_metadata.csv.gz"
    record=tmp_path/"terminal.json"
    record.write_text(json.dumps({"schema":m.TERMINAL_SCHEMA,"measured_photo_file_sha256":"wrong"}))
    with pytest.raises(RuntimeError,match="not source-authenticated"):
        m.original_inference(pictures,record,tmp_path/"infer")


def test_orig_methods_and_partition_constants_not_changed():
    assert (m.N_BATCHES,m.N_SHARDS,m.N_PARTS,m.TOTAL_PARTS,m.TARGET)==(2,32,4,256,100)
    assert m.SOURCE_RUN_ID==38996743273
    assert m.SOURCE_HEAD=="61abc4f71e74db3a6633a790e37d1f30e8e3bfcc"
