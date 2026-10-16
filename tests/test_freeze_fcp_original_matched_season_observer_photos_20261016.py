"""Synthetic source-only matching tests; no original image pixels or flower labels."""
from __future__ import annotations
import csv,json
from pathlib import Path
import pytest

from scripts.acquisition.freeze_fcp_original_matched_season_observer_photos_20261016 import (
 SOURCE_FIELDS,BLIND,sha,freeze,local_year,measure_id
)

def csv_out(path,fields,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        w.writeheader();w.writerows(rows)

def inputs(tmp):
    src=tmp/"source.csv";receipt=tmp/"source.json"
    site=tmp/"sites.csv";calendar=tmp/"calendar.json"
    photos=[];monthly=[];pid=100000
    for season,month in (("spring",4),("summer",7)):
        for i in range(26):
            pid+=1
            owner=(999 if i==25 else (1000+i if month==4 else 2000+i))
            photos.append(dict(
                inat_taxon_id="101",species="Flower example",photo_id=str(pid),
                observation_id=f"original-{season}-{i}",observer_id=str(owner),
                latitude="35.125",longitude="135.125",
                observed_on=f"2024-{month:02d}-{i%27+1:02d}",
                positional_accuracy_m="125",photo_license="cc-by",photo_creator_id=str(owner),
                attribution="Example source attribution",
                photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                source_snapshot_date="2026-09-27"))
        monthly.append(dict(inat_taxon_id="101",species="Flower example",
            year_month=f"2024-{month:02d}",lat_0p25_cell="500",lon_0p25_cell="1260",
            all_source_photo_ids_with_5km_accuracy="26",
            distinct_source_observations="26",distinct_source_observers="26"))
    csv_out(src,SOURCE_FIELDS,photos)
    csv_out(site,list(monthly[0]),monthly)
    receipt.write_text(json.dumps({
        "status":"COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS",
        "source_snapshot_date":"2026-09-27",
        "selected_all_photo_IDs_csv_sha256":sha(src),
        "n_selected_original_photo_IDs":len(photos),
        "all_enrolled_species_match_independent_full_source_photo_ID_counts":True,
        "original_image_pixels_opened":False,
        "historical_flower_colour_CSV_imported":False}))
    calendar.write_text(json.dumps({
        "status":"COMPLETE_ORIGINAL_SOURCE_CALENDAR_CENSUS_NO_FLOWER_PIXELS",
        "original_source_site_month_file_sha256":sha(site),
        "source_photo_opportunity_is_not_flowering_detection":True,
        "new_flower_colour_or_white_frequency_estimated":False}))
    return src,receipt,site,calendar

def read(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def test_source_observer_disjoint_first_and_reserved_photo_work_waves(tmp_path):
    a=inputs(tmp_path)
    result=freeze(*a,tmp_path/"out",per_season_total=4,synthetic=True)
    assert result["status"]=="SYNTHETIC_NO_BIOLOGICAL_MEASUREMENT"
    assert result["selected_photo_ids"]==8 and result["photos_per_wave"]==4
    assert result["selected_source_pair"]["observers_shared"]==1
    assert result["selected_source_pair"]["disjoint_observers_a"]==25
    assert result["selected_source_pair"]["disjoint_observers_b"]==25
    rows=read(tmp_path/"out/original_source_season_frozen_128_photos.csv")
    assert len({r["photo_id"] for r in rows})==8
    assert len({r["observation_id"] for r in rows})==8
    assert len({r["observer_id"] for r in rows})==8
    assert {r["local_season"] for r in rows}=={"spring","summer"}
    assert {r["wave"] for r in rows}=={"first","reserve"}
    assert all(r["measurement_id"]==measure_id(r["photo_id"]) for r in rows)
    worker=read(tmp_path/"out/season_first_original_blind_worker.csv")
    assert len(worker)==4 and list(worker[0])==list(BLIND)
    assert all("species" not in r and "observed_on" not in r for r in worker)
    again=freeze(*a,tmp_path/"again",per_season_total=4,synthetic=True)
    assert result["file_sha256"]==again["file_sha256"]
    with pytest.raises(ValueError,match="overwrite"):
        freeze(*a,tmp_path/"out",per_season_total=4,synthetic=True)

def test_winter_rollover_and_tropical_rejection():
    assert local_year(12,2024,500)==("winter",2025)
    assert local_year(1,2025,500)==("winter",2025)
    assert local_year(7,2025,500)==("summer",2025)
    assert local_year(7,2025,200)==("winter",2025)

def test_all_originals_required_and_prior_colours_forbidden(tmp_path):
    src,receipt,site,calendar=inputs(tmp_path)
    src.write_text(src.read_text().replace("Flower example","New label",1))
    with pytest.raises(ValueError,match="authenticity"):
        freeze(src,receipt,site,calendar,tmp_path/"bad",per_season_total=4,synthetic=True)

def test_site_metadata_sha_guard_and_observer_support(tmp_path):
    a=inputs(tmp_path)
    a[2].write_text(a[2].read_text().replace("2024-07","2025-07"))
    with pytest.raises(ValueError,match="independently sealed"):
        freeze(*a,tmp_path/"bad",per_season_total=4,synthetic=True)
    b=inputs(tmp_path)
    with pytest.raises(ValueError,match="Insufficient observer-disjoint"):
        freeze(*b,tmp_path/"over",per_season_total=28,synthetic=True)

def test_no_partial_small_synthetic_source_can_pass_real(tmp_path):
    a=inputs(tmp_path)
    with pytest.raises(ValueError,match="authenticity"):
        freeze(*a,tmp_path/"bad",per_season_total=4,synthetic=False)
