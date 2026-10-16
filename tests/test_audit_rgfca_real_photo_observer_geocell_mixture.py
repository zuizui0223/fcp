"""No-network source/observer/geocell validation for six-taxon photo QC waves."""
from __future__ import annotations
import csv
import json
import pytest

from scripts.analysis.audit_rgfca_real_photo_observer_geocell_mixture import (
    PALETTE, ROI_STATUSES, VALID_SOURCE_SCHEMAS, VALID_ROI_SCHEMAS,
    sha256, summarize,
)


def build_fixture(tmp_path, wave):
    root=tmp_path/"original_pilot"
    root.mkdir()
    roster=tmp_path/"rgfca_roster.csv"
    with roster.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","source_photo_opportunities",
                    "opportunity_stratum","enrolment_hash"])
        w.writerows([[1000+i,f"Species {i}",500,"medium_256_1023","hash"]
                     for i in range(6)])
    prefix="rgfca_six" if wave==1 else "rgfca_six_wave2"
    suffix="" if wave==1 else "_wave2"
    source=root/f"{prefix}_source_photo_candidates_before_colour.csv"
    measured=root/f"rgfca_original_new_photo_six_taxa{suffix}_live_API_technical.csv"
    sources=[]
    outcomes=[]
    for i in range(6):
        for j in range(5):
            pid=str(50000+i*10+j)
            original={
                "inat_taxon_id":str(1000+i),"species":f"Species {i}",
                "photo_id":pid,"observation_id":str(80000+i*10+j),
                "observer_id":str(90000+i*10+j),
                "screen_fold":str(j%2),
                "latitude":str(35+j),"longitude":"139.5",
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.JPG",
                "photo_license":"cc-by","provider_photo_url":"technical",
                "source_opportunity_tier":"medium_256_1023","sample_hash":"blind"
            }
            good=j<4
            measurement={**original,
                         "global_classifiable":"true" if good else "false",
                         "morph":("yellow_orange" if j%2 else "white")
                                  if good else "mixed_uncertain",
                         "measurement_status":"classified_four_state_morph"
                               if good else "not_evaluable_roi_or_flip_gate"}
            for k in PALETTE:
                measurement[f"flower_fraction_{k}"]=(
                    ("0.52" if k=="white" else "0.48" if k=="yellow" else "0")
                    if good else "0")
            sources.append(original)
            outcomes.append(measurement)
    for path,rows in ((source,sources),(measured,outcomes)):
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    source_doc=root/f"{prefix}_source_precolour_receipt.json"
    source_receipt={
        "schema":next(k for k,v in VALID_SOURCE_SCHEMAS.items() if v==wave),
        "status":"SELECTED_ORIGINAL_PROVIDER_API_PHOTO_IDS_NO_IMAGE_PIXELS_YET",
        "NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS":True,
        "EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED":False,
        "photo_selection_used_literature_phenotype":False,
        "photo_selection_used_climate_soil_or_prior_photo_colours":False,
        "image_bytes_opened":False,
        "RGFCA_200_realization_photo_polymorphism_Ppoly_computed":False,
        "source_ROSTER_sha256":sha256(roster),
        "selected_metadata_sha256":sha256(source),
        "pilot_taxa":6,
        "technical_only_disjoint_previous_six_taxa":wave==2,
    }
    source_doc.write_text(json.dumps(source_receipt))
    roi_doc=root/f"rgfca_original_new_photo_six_taxa{suffix}_live_API_technical_receipt.json"
    roi_receipt={
        "schema":next(k for k,v in VALID_ROI_SCHEMAS.items() if v==wave),
        "status":ROI_STATUSES[wave],
        "source_kind":"iNaturalist_API_created_before_20260927_live_snapshot_NOT_dated_monthly_archive",
        "exact_dated_20260927_source_membership_admitted":False,
        "full_1000_species_ROI_study_complete":False,
        "prior_FCP_measured_colours_loaded":False,
        "new_terminal_photo_csv_sha256":sha256(measured),
        "source_image_ID_selection_sha256":sha256(source),
        "source_API_selection_receipt_sha256":sha256(source_doc),
        "number_taxa_selected":6,
        "selected_real_original_photo_IDs":30,
        "original_pixels_ROI_v4_evaluated":30,
        "n_classifiable_original_flower_photos":24,
        "technical_only_disjoint_previous_six_taxa":wave==2,
    }
    roi_doc.write_text(json.dumps(roi_receipt))
    return root,roster,source,measured,roi_doc


@pytest.mark.parametrize("wave",[1,2])
def test_synthetic_six_taxa_source_and_roi_match_across_observers(wave,tmp_path):
    root,roster,_,_,_=build_fixture(tmp_path,wave)
    out=tmp_path/"result"
    report=summarize(root,roster,out,synthetic=True)
    assert report["n_taxa"]==6
    assert report["pilot_wave"]==wave
    assert report["n_original_ROI_photos"]==30
    assert report["n_classifiable_flower_ROIs"]==24
    assert report["original_source_date_bulk_archive_membership_verified"] is False
    with (out/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.csv").open() as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==6
    assert all(int(x["n_classifiable_original_observers"])==4 for x in rows)
    assert all(int(x["within_one_ROI_white_and_yellow_both_ge10"])==4 for x in rows)
    assert all(float(x["depth_48_realized_photo_ID_overlap_fraction"])==1.0 for x in rows)
    assert all(int(x["depth_48_distinct_source_photo_ID_draw_sets"])==1 for x in rows)
    assert report["actual_pairwise_source_photo_ID_overlap_by_draw_depth_included"] is True
    assert all(x["genetic_or_population_colour_polymorphism_confirmed"]=="False" for x in rows)


def test_changed_original_photo_id_rejected_even_with_forced_result_hash(tmp_path):
    root,roster,source,measured,roi_doc=build_fixture(tmp_path,1)
    measured.write_text(measured.read_text().replace("50000,","99999,",1))
    d=json.loads(roi_doc.read_text())
    d["new_terminal_photo_csv_sha256"]=sha256(measured)
    roi_doc.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="Repeated or replaced"):
        summarize(root,roster,tmp_path/"tampered",synthetic=True)


def test_false_dated_snapshot_admission_never_enables_photo_FCP(tmp_path):
    root,roster,_,_,roi_doc=build_fixture(tmp_path,2)
    d=json.loads(roi_doc.read_text())
    d["exact_dated_20260927_source_membership_admitted"]=True
    roi_doc.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="negative FCP"):
        summarize(root,roster,tmp_path/"tampered",synthetic=True)
