"""Synthetic no-network, source-negative 12-taxon two original-photo RGFCA wave test."""
from __future__ import annotations
import csv
import json
import pytest

from scripts.analysis.compare_rgfca_two_disjoint_original_photo_roi_waves import (
    DEPTHS, read_wave, compare, sha256,
)


def make_two(tmp_path):
    roster=tmp_path/"rgfca_roster.csv"
    with roster.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","source_photo_opportunities",
                    "opportunity_stratum","enrolment_hash"])
        for wave in (1,2):
            for k in range(6):
                tid=1000+100*wave+k
                tier=["shallow_128_255","medium_256_1023","deep_ge1024"][k//2]
                w.writerow([tid,f"Species {tid}",128 if k<2 else 500 if k<4 else 1050,
                            tier,"frozen_source_hash"])
    roots=[]
    for wave in (1,2):
        root=tmp_path/f"wave{wave}"
        root.mkdir()
        rows=[]
        for k in range(6):
            tid=1000+100*wave+k
            valid=50 if k%2==0 else 20
            n=96
            cells=3
            data={
                "pilot_wave":wave,"inat_taxon_id":tid,"species":f"Species {tid}",
                "source_photo_ids_selected":n,"classifiable_ROI":valid,
                "classifiable_fraction":valid/n,
                "n_classifiable_original_observers":min(valid,12),
                "n_classifiable_5degree_cells":min(valid,cells),
                "photo_label_counts":json.dumps({"white":valid-5,"yellow_orange":5}),
                "within_one_ROI_white_and_yellow_fraction":.6 if valid>=40 else .1,
                "genetic_or_population_colour_polymorphism_confirmed":False,
            }
            for depth in DEPTHS:
                eligible=200 if valid>=40 and depth>=40 else 0
                positive=170 if eligible else 0
                data[f"depth_{depth}_evaluable_repetitions"]=eligible
                data[f"depth_{depth}_positive_photo_label_repetitions"]=positive
                data[f"depth_{depth}_conditional_photo_label_detection"]=(
                    positive/eligible if eligible else "")
            rows.append(data)
        table=root/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.csv"
        with table.open("w",newline="") as f:
            w=csv.DictWriter(f,list(rows[0]),lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        receipt={
            "schema":"fcp_real_rgfca_original_photo_technical_observer_geo_organ_QC_v1",
            "status":"TECHNICAL_REAL_ORIGINAL_ROI_ONLY_NOT_POLYMORPHISM_ESTIMATE",
            "pilot_wave":wave,"n_taxa":6,
            "n_original_ROI_photos":576,
            "n_classifiable_flower_ROIs":210,
            "source_pilot_random_photo_draws_per_taxon":200,
            "pilot_resampling_fold_depths":list(DEPTHS),
            "CSV_sha256":sha256(table),
            "original_source_date_bulk_archive_membership_verified":False,
            "main_1000_RGFCA_P_poly_estimated":False,
            "no_genetic_FCP_species_asserted":True,
            "no_old_colour_measurements_reused":True,
        }
        (root/f"wave{wave}_original_ROI_observer_geocell_mixed_colour_QC.json").write_text(
            json.dumps(receipt))
        roots.append(root)
    return roster,*roots


def test_photo_bias_two_original_technical_waves_remain_no_genetic_claim(tmp_path):
    roster,w1,w2=make_two(tmp_path)
    report=compare(w1,w2,roster,tmp_path/"comparison",synthetic=True)
    assert report["n_original_source_taxa_pilot"]==12
    assert report["source_taxon_overlaps_between_waves"]==0
    assert report["total_original_new_image_pixels_ROI_measured"]==1152
    assert report["total_classifiable_original_flower_ROIs"]==420
    assert report["n_species_with_at_least_40_classifiable_ROIs"]==6
    assert report["n_species_flagged_for_flower_organ_structure_review"]==6
    assert report["actual_1000_taxon_200x250_formal_Ppoly_estimated"] is False
    assert report["actual_20260927_original_source_photo_ID_membership_verified"] is False
    with (tmp_path/"comparison/rgfca_12_real_technical_taxa_photo_roi_observer_bias.csv").open() as f:
        data=list(csv.DictReader(f))
    assert len(data)==12
    assert all(x["phenotypic_polymorphism_confirmed"]=="False" for x in data)
    assert all(v["n_taxa"]==4
               for v in report["available_photo_opportunity_strata"].values())


def test_altered_RGFCA_redraw_denominator_or_biological_claim_is_rejected(tmp_path):
    roster,w1,w2=make_two(tmp_path)
    subject=w2/"wave2_original_ROI_observer_geocell_mixed_colour_QC.csv"
    source=subject.read_text()
    subject.write_text(source.replace("170,0.85","201,1.005",1))
    receipt=w2/"wave2_original_ROI_observer_geocell_mixed_colour_QC.json"
    data=json.loads(receipt.read_text())
    data["CSV_sha256"]=sha256(subject)
    receipt.write_text(json.dumps(data))
    with pytest.raises(ValueError,match="denominator|count"):
        compare(w1,w2,roster,tmp_path/"invalid",synthetic=True)


def test_two_waves_must_have_disjoint_original_taxon_IDs(tmp_path):
    roster,w1,w2=make_two(tmp_path)
    source=w2/"wave2_original_ROI_observer_geocell_mixed_colour_QC.csv"
    text=source.read_text()
    source.write_text(text.replace("2100,","1100,",1).replace("Species 2100","Species 1100",1))
    receipt=w2/"wave2_original_ROI_observer_geocell_mixed_colour_QC.json"
    data=json.loads(receipt.read_text())
    data["CSV_sha256"]=sha256(source)
    receipt.write_text(json.dumps(data))
    with pytest.raises(ValueError,match="not independent"):
        compare(w1,w2,roster,tmp_path/"invalid",synthetic=True)


def test_no_overwrite_for_verified_real_12_taxon_calibration(tmp_path):
    roster,w1,w2=make_two(tmp_path)
    out=tmp_path/"done"
    compare(w1,w2,roster,out,synthetic=True)
    with pytest.raises(ValueError,match="Cannot overwrite"):
        compare(w1,w2,roster,out,synthetic=True)
