"""Source-bound, blind, geographically matched original photo test ONLY.

Tests use synthetic metadata and a paired fake evaluator; no original photos
or original colour outcomes are opened in CI.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from fcp_pipeline.photo_polymorphism_screen import sha256
from fcp_pipeline.fcp_fresh_study import BIO
from scripts.acquisition.build_fcp_uncapped_blind_batches import new_id
from scripts.acquisition.plan_fcp_photo_first_screen import SELECTED_FIELDS
from scripts.acquisition.run_fcp_rgfca_first_geography_pair_roi import (
    PHOTO_SCHEMA,WORKER,ACQUISITION,geographical_cell,
    geo_balanced_first_photos,select_unit,measure_unit
)
from scripts.acquisition.stratify_fcp_rgfca_geographic_span import run as geo_run
from test_stratify_fcp_rgfca_geographic_span import fixture as geo_fixture


def setup(tmp_path):
    manifest,monthly_receipt,roster=geo_fixture(tmp_path)
    geo=tmp_path/"geo"
    geo_run(manifest,monthly_receipt,roster,geo,synthetic=True)
    src=tmp_path/"source"
    dated=src/"rgfca_1000_exact_source_photo_IDs"
    dated.mkdir(parents=True)
    shutil.copy(monthly_receipt,dated/"rgfca_1000_targeted_monthly_photo_ID_receipt.json")
    shutil.copy(manifest,dated/"rgfca_1000_complete_original_source_photo_IDs.csv")
    selection=src/"rgfca_1000_photo_colour_blind_selection"
    selection.mkdir()
    firewall=src/"rgfca_1000_original_photo_blind_worker_freeze"
    firewall.mkdir()

    with manifest.open("r",newline="",encoding="utf-8") as f:
        source=list(csv.DictReader(f))
    rows=[]
    seen_obs=set()
    for i,r in enumerate(source):
        if r["inat_taxon_id"] not in ("101","104"):
            continue
        if r["observation_id"] in seen_obs:
            continue
        seen_obs.add(r["observation_id"])
        n=sum(q["inat_taxon_id"]==r["inat_taxon_id"] for q in rows)
        row={
            "inat_taxon_id":r["inat_taxon_id"],"species":r["species"],
            "photo_id":r["photo_id"],"observation_id":r["observation_id"],
            "observer_id":r["observer_id"],"latitude":r["latitude"],
            "longitude":r["longitude"],"photo_url_large":r["photo_url_large"],
            "photo_license":r["photo_license"],"source_snapshot_date":r["source_snapshot_date"],
            "screen_fold":str(n%2),"selection_hash":f"{i+1:016x}",
        }
        rows.append(row)
    selected=selection/"photo_screen_original_metadata_selected.csv"
    with selected.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=SELECTED_FIELDS)
        w.writeheader();w.writerows(rows)
    selection_receipt=selection/"photo_screen_source_selection_receipt.json"
    selection_receipt.write_text(json.dumps({
        "source_snapshot_date":"2026-09-27",
        "photo_selection_csv_sha256":sha256(selected),
        "selected_photo_IDs":len(rows),
        "selected_sample_photo_pixels_measured":False,
    }))
    freeze=firewall/"photo_screen_freeze.json"
    freeze.write_text(json.dumps({
        "pre_colour_selection_csv_sha256":sha256(selected),
        "source_photo_screen_selection_receipt_sha256":sha256(selection_receipt),
        "image_worker_sees_only_three_fields":True,
        "no_old_photo_colour_rows_loaded":True,
        "original_image_pixels_opened":False,
    }))
    join=firewall/"metadata_join_key.csv"
    with join.open("w",newline="",encoding="utf-8") as f:
        cols=["measurement_id","inat_taxon_id","species","photo_id",
              "observation_id","observer_id","latitude","longitude","screen_fold"]
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow({"measurement_id":new_id(row["photo_id"]),**{
                k:row[k] for k in cols if k!="measurement_id"}})
    return src,geo,selected,selection_receipt,freeze


def test_first_matched_pair_is_source_blind_and_does_not_rewrite_1000_taxa(tmp_path):
    src,geo,selected,selection_receipt,freeze=setup(tmp_path)
    output=tmp_path/"prepared"
    receipt=select_unit(src,geo,output,_synthetic_source_gate=True)
    assert receipt["schema"]==PHOTO_SCHEMA
    assert receipt["status"]=="ORIGINAL_MONTHLY_SOURCE_GEOGRAPHIC_PAIR_PREPARED_NO_PIXELS"
    assert receipt["matched_pair_id"]=="1"
    assert receipt["broad_taxon_id"]=="101" and receipt["narrow_taxon_id"]=="104"
    assert receipt["n_original_photo_pixels_opened"]==0
    assert receipt["source_image_IDs_from_verified_monthly_20260927"] is False
    assert receipt["number_original_photo_IDs_in_first_unit"]==24
    assert receipt["per_taxon_per_observer_fold_photo_depth"]==6
    assert receipt["original_photo_first_unit_geographic_cell_round_robin"] is True
    assert receipt["broad_original_first_pixel_geographic_support"]["n_supported_selected_cells"]>=2
    assert receipt["broad_original_first_pixel_geographic_support"][
        "max_supported_selected_cell_centroid_distance_km"]>=1000
    assert receipt["source_location_accuracy_reverified_in_first_ROI_photo_subset"] is True
    assert receipt["n_source_geographic_cells_selected_by_group"]["broad"]>=2
    worker=list(csv.DictReader((output/"first_matched_pair_original_blind_worker.csv").open()))
    acq=list(csv.DictReader((output/"first_matched_pair_original_photo_source_acquisition.csv").open()))
    sidecar=list(csv.DictReader((output/"first_matched_pair_original_source_join_key.csv").open()))
    assert len(worker)==len(acq)==len(sidecar)==24
    assert list(worker[0])==list(WORKER)
    assert list(acq[0])==list(ACQUISITION)
    assert all("species" not in x and "photo_url_large" not in x for x in worker)
    assert {r["geographic_comparison_group"] for r in sidecar}=={"broad","narrow"}
    assert {r["inat_taxon_id"] for r in sidecar}=={"101","104"}
    assert all(acq[i]["measurement_id"]==worker[i]["measurement_id"] for i in range(24))
    assert receipt["original_blind_worker_sha256"]==sha256(
        output/"first_matched_pair_original_blind_worker.csv")
    with pytest.raises(ValueError,match="overwrite"):
        select_unit(src,geo,output,_synthetic_source_gate=True)


def test_original_roi_terminal_preserves_observer_blindness_and_source_IDs(tmp_path):
    source,geo,_,_,_=setup(tmp_path)
    prepared=tmp_path/"prepared"
    receipt=select_unit(source,geo,prepared,_synthetic_source_gate=True)
    exposure=[]
    def fake_loader(url):
        exposure.append(url)
        return b"fake_provider_original_pixel_bytes"
    def fake_evaluator(worker,images,model,detector,sam,dst):
        assert list(pd.read_csv(worker).columns)==list(WORKER)
        w=pd.read_csv(worker,dtype=str)
        result=[]
        for i,rec in w.iterrows():
            cols={f"flower_fraction_{key}":0.0 for key in BIO}
            if i%2==0:
                cols["flower_fraction_white"]=1.0
                result.append({"measurement_id":rec.measurement_id,
                               "measurement_status":"classified_four_state_morph",
                               "morph":"white",**cols})
            else:
                result.append({"measurement_id":rec.measurement_id,
                               "measurement_status":"not_evaluable_roi_or_flip_gate",
                               "morph":"mixed_uncertain",**cols})
        return pd.DataFrame(result)
    done=measure_unit(prepared,tmp_path/"measure",method=tmp_path/"source.py",
                      detector=tmp_path/"det.pt",sam=tmp_path/"sam",
                      loader=fake_loader,evaluator=fake_evaluator)
    assert done["status"]=="SYNTHETIC_ORIGINAL_GEOGRAPHIC_PAIR_TEST_NOT_BIOLOGY"
    assert done["new_original_photo_pixels_processed"]==24
    assert done["n_classifiable"]==12
    assert done["n_original_image_acquisition_failed"]==0
    assert done["by_geographic_observed_source_group"]["broad"][
        "n_same_original_roi_white_and_yellow_both_ge10pct"]==0
    assert done["by_geographic_observed_source_group"]["narrow"][
        "n_same_original_roi_white_and_yellow_both_ge10pct"]==0
    assert done["by_geographic_observed_source_group"]["broad"][
        "n_distinct_geographic_5degree_cells_in_measured_photos"]>=2
    assert done["old_measured_colours_imported"] is False
    assert done["full_RGFCAROI1000_complete"] is False
    assert len(exposure)==24
    rows=list(csv.DictReader((tmp_path/"measure"/"first_original_geographically_matched_photo_terminal.csv").open()))
    assert len(rows)==24
    assert {r["geographic_comparison_group"] for r in rows}=={"broad","narrow"}
    assert done["original_geographic_source_join_terminal_sha256"]==sha256(
        tmp_path/"measure"/"first_original_geographically_matched_photo_terminal.csv")
    with pytest.raises(ValueError,match="already sealed"):
        measure_unit(prepared,tmp_path/"measure",method=tmp_path/"source.py",
                     detector=tmp_path/"det.pt",sam=tmp_path/"sam",
                     loader=fake_loader,evaluator=fake_evaluator)


def test_geographic_pairs_reject_tampered_metadata_before_any_pixel(tmp_path):
    src,geo,_,_,_=setup(tmp_path)
    pair=geo/"rgfca_1000_original_source_broad_narrow_depthmatched_pairs.csv"
    pair.write_text(pair.read_text().replace("101,","999,",1))
    with pytest.raises(ValueError,match="not signed|ungrounded|unmatched"):
        select_unit(src,geo,tmp_path/"out",_synthetic_source_gate=True)


def test_real_production_cannot_use_fake_source_gate_or_synthetic_method(tmp_path):
    src,geo,_,_,_=setup(tmp_path)
    with pytest.raises(Exception):
        select_unit(src,geo,tmp_path/"production",_synthetic_source_gate=False)
    prepared=tmp_path/"prepared"
    select_unit(src,geo,prepared,_synthetic_source_gate=True)
    with pytest.raises(ValueError,match="Production original image pixels cannot use a synthetic provenance gate"):
        measure_unit(prepared,tmp_path/"false_production",method=tmp_path/"model.py",
                     detector=tmp_path/"det.pt",sam=tmp_path/"sam")



def test_first_roi_localities_round_robin_instead_of_one_city_photo_hashes():
    photos=[]
    for cell,longitude in enumerate([2.,32.,75.,121.]):
        for j in range(5):
            # All photos in first cell have lower original colour-blind hashes:
            # naive hash-order sample would consume just the first region.
            photos.append({
                "photo_id":str(10000+cell*100+j),
                "latitude":"10.5","longitude":str(longitude),
                "selection_hash":f"{cell*100+j:016x}",
            })
    got=geo_balanced_first_photos(photos,8)
    from collections import Counter
    counts=Counter(geographical_cell(r) for r in got)
    assert sorted(counts.values())==[2,2,2,2]
    assert len({r["photo_id"] for r in got})==8
    assert geographical_cell({"latitude":"0","longitude":"180"})==\
        geographical_cell({"latitude":"0","longitude":"-180"})
    with pytest.raises(ValueError,match="geographic coordinate invalid"):
        geographical_cell({"latitude":"invalid","longitude":"5"})


def test_original_photographic_yellow_white_same_roi_is_separate_structural_QC(tmp_path):
    src,geo,_,_,_=setup(tmp_path)
    unit=tmp_path/"prepared"
    select_unit(src,geo,unit,_synthetic_source_gate=True)

    def fake_evaluator(worker,images,model,detector,sam,dst):
        df=pd.read_csv(worker,dtype=str)
        output=[]
        for j,rec in df.iterrows():
            composition={f"flower_fraction_{key}":0.0 for key in BIO}
            composition["flower_fraction_white"]=0.7
            composition["flower_fraction_yellow"]=0.3
            output.append({"measurement_id":rec.measurement_id,
                           "measurement_status":"classified_four_state_morph",
                           "morph":"white",**composition})
        return pd.DataFrame(output)

    receipt=measure_unit(unit,tmp_path/"measured",method=tmp_path/"code.py",
                         detector=tmp_path/"model.pt",sam=tmp_path/"sam",
                         loader=lambda _:b"fake_byte_payload",
                         evaluator=fake_evaluator)
    assert receipt["n_classifiable"]==24
    for group in ("broad","narrow"):
        part=receipt["by_geographic_observed_source_group"][group]
        assert part["n_same_original_roi_white_and_yellow_both_ge10pct"]==12
        assert part["n_distinct_source_observations_in_measured_photos"]==12
    with (tmp_path/"measured"/"first_original_geographically_matched_photo_terminal.csv").open() as f:
        rows=list(csv.DictReader(f))
    assert all(r["same_roi_white_and_yellow_both_ge10pct"]=="True" for r in rows)
    assert receipt["measured_geographic_effect_estimates"]==0
    assert receipt["confirmed_genetic_FCP_species"]==0



def test_geographically_distant_single_original_photo_cannot_claim_broad_site_support():
    from scripts.acquisition.run_fcp_rgfca_first_geography_pair_roi import \
        supported_broad_photo_spread
    pictures=[]
    for i in range(12):
        lat,lon=(36.,139.) if i<11 else (-33.,151.)
        pictures.append({
            "latitude":str(lat),"longitude":str(lon),
            "photo_id":str(i+1),"observation_id":str(1000+i),
            "observer_id":str(2000+i//2),
            "positional_accuracy_m":"100",
        })
    result=supported_broad_photo_spread(pictures)
    assert result["n_supported_selected_cells"]==1
    assert result["max_supported_selected_cell_centroid_distance_km"]==0
    assert result["geographic_accuracy_of_selected_subset_independently_verified"] is True
    for new in (1,2,3):
        pictures.append({
            "latitude":str(-33.+new*.01),"longitude":"151",
            "photo_id":str(100+new),"observation_id":str(9000+new),
            "observer_id":str(7000+new),
            "positional_accuracy_m":"150000", # 150 km uncertainty; reject
        })
    result=supported_broad_photo_spread(pictures)
    assert result["n_supported_selected_cells"]==1
    assert result["max_supported_selected_cell_centroid_distance_km"]==0


def test_tampered_original_dated_photo_location_never_enters_real_pixel_roi(tmp_path):
    original,geo,_,_,_=setup(tmp_path)
    full=original/"rgfca_1000_exact_source_photo_IDs"
    path=full/"rgfca_1000_complete_original_source_photo_IDs.csv"
    with path.open(newline="") as file:
        rows=list(csv.DictReader(file))
    rows[0]["latitude"]="87.123" # Wrong original provider monthly location
    with path.open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    # Simulate a forged/corrupted monthly receipt that follows the altered
    # file checksum: the independent original photo-source selection still
    # disagrees with actual coordinates and must refuse it.
    receipt_file=full/"rgfca_1000_targeted_monthly_photo_ID_receipt.json"
    doc=json.loads(receipt_file.read_text())
    doc["selected_all_photo_IDs_csv_sha256"]=sha256(path)
    receipt_file.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="not signed|site differs|monthly"):
        select_unit(original,geo,tmp_path/"bad_geo",_synthetic_source_gate=True)
