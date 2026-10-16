"""Synthetic-only reproducibility tests for RGFCA photo-colour detection robustness.

200 Monte Carlo maps/draws resample the same original photo pool. These
tests deliberately do not diagnose biology or old legacy photo colours.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from fcp_pipeline.photo_polymorphism_screen import PALETTE, observer_split, sha256
from fcp_pipeline.rgfca_photo_repeatability import (
    load_new_roi, one_draw, repeat, validate_schedule
)


def original_photos(taxon: int, *, variable=True, quality=1.0):
    photos = []
    n = 0
    for fold in (0, 1):
        owners = []
        k = 0
        while len(owners) < 4:
            person = f"owner{taxon}_{k}"
            if observer_split(taxon,person) == fold:
                owners.append(person)
            k += 1
        for j in range(64):
            isgood = (j < int(64*quality))
            morph = ("red_pink" if (variable and j % 5 == 0)
                     else "white") if isgood else None
            photos.append({
                "photo_id":str(taxon*1000+n), "observer_id":owners[j % 4],
                "observation_id":f"{taxon}-{n}", "fold":fold,
                "classifiable":isgood, "morph":morph
            })
            n += 1
    return photos


def setup(tmp_path):
    cohort=tmp_path/"cohort.csv"
    with cohort.open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["inat_taxon_id","species","source_photo_opportunities",
                         "opportunity_stratum","enrolment_hash"])
        writer.writerows([[11,"Species redwhite",1000,"deep_ge1024","x"],
                          [12,"Species alwayswhite",1000,"deep_ge1024","y"]])
    sched=tmp_path/"schedule.csv"
    with sched.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["realization","cycle","slot","inat_taxon_id","species","source_opportunity_tier"])
        for n in range(1,51):
            w.writerow([n,n,1,11,"Species redwhite","deep_ge1024"])
            w.writerow([n,n,2,12,"Species alwayswhite","deep_ge1024"])
    provenance=tmp_path/"schedule_receipt.json"
    provenance.write_text(json.dumps({
        "schema":"fcp_rgfca_photo_polymorphism_balanced_sample_v1",
        "status":"ORIGINAL_SOURCE_BALANCED_REPETITION_PLAN_NO_COLOUR_MEASURED",
        "selected_species_csv_sha256":sha256(cohort),
        "realization_schedule_csv_sha256":sha256(sched),
        "selected_species":2,"realizations":50,
        "species_per_realization":2,"species_inclusions_each":50,
        "old_149900_flower_colours_used":False,
        "new_flower_images_opened":False,
        "photographic_poly_probabilities_measured":False
    }))
    measured=tmp_path/"fresh_original_ROI_screen.csv"
    with measured.open("w",newline="") as f:
        fields=["inat_taxon_id","species","photo_id","observer_id","observation_id",
                "screen_fold","morph","measurement_status","global_classifiable",
                *(f"flower_fraction_{x}" for x in PALETTE)]
        w=csv.DictWriter(f,fields,lineterminator="\n")
        w.writeheader()
        for taxon,name,variable in ((11,"Species redwhite",True),(12,"Species alwayswhite",False)):
            for row in original_photos(taxon,variable=variable):
                state=row["morph"]
                item={
                    "inat_taxon_id":taxon,"species":name,"photo_id":row["photo_id"],
                    "observer_id":row["observer_id"],"observation_id":row["observation_id"],
                    "screen_fold":row["fold"],"morph":state,
                    "measurement_status":"classified_four_state_morph",
                    "global_classifiable":"true",
                }
                item.update({f"flower_fraction_{x}":int(
                    x==("white" if state=="white" else "pink")) for x in PALETTE})
                w.writerow(item)
    receipt=tmp_path/"sealed_original_roi_receipt.json"
    receipt.write_text(json.dumps({
        "schema":"fcp_sealed_photo_screen_terminal_v1",
        "status":"SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL",
        "protocol":"PHOTO_FCP_OBSERVER_DISJOINT_SCREEN_V1",
        "source_date":"2026-09-27",
        "new_photo_table_sha256":sha256(measured),
        "historical_measured_flower_colours_reused":False,
        "photo_ID_selection_pre_colour":True,"source_photo_IDs_sealed":True,
        "screen_not_entire_original_93million_census":True,
        "full_source_image_measurement_complete":False,
        "n_terminal_measurements":256,
        "synthetic_test_only":True
    }))
    return cohort,sched,provenance,measured,receipt


def test_true_two_colour_photo_states_repeat_across_observer_disjoint_draws():
    photos=original_photos(11)
    draws=[one_draw(11,"Species redwhite",1000,photos,r,per_fold=40)
           for r in range(1,51)]
    assert all(x["n_selected"]==80 for x in draws)
    assert all(x["evaluable"] for x in draws)
    assert sum(x["photo_poly_detected"] for x in draws) >= 45
    assert draws[0]["fold_0_classifiable"] == 40
    assert draws[0]["fold_1_classifiable"] == 40


def test_one_photographic_state_is_not_monomorphic_evidence():
    draws=[one_draw(12,"Species alwayswhite",1000,original_photos(12,variable=False),r,per_fold=40)
           for r in range(1,21)]
    assert all(x["evaluable"] for x in draws)
    assert all(x["photo_poly_detected"] is False for x in draws)
    assert all("not_monomorphic" in x["coarse_morph_status"] for x in draws)


def test_unclassifiable_source_photos_are_not_replaced_or_dropped():
    photos=original_photos(11,quality=.25)
    r=one_draw(11,"Species redwhite",1000,photos,1,per_fold=40)
    assert r["n_selected"]==80
    assert r["n_classifiable"] < 40
    assert r["evaluable"] is False
    assert r["photo_poly_detected"] is False


def test_full_repeatability_with_sha_frozen_synthetic_screen(tmp_path):
    cohort,sched,provenance,measured,receipt=setup(tmp_path)
    species,assignments=validate_schedule(cohort,sched,provenance,synthetic=True)
    pools=load_new_roi(measured,receipt,species,synthetic=True)
    assert len(pools)==2
    report=repeat(species,assignments,pools,tmp_path/"output",synthetic=True)
    assert report["selected_species"]==2
    assert report["balanced_realizations"]==50
    assert report["number_photo_draw_assignments"]==100
    assert report["n_repeated_phenotypic_polymorphism_species"]==1
    assert report["number_independent_new_biological_worlds"]==1
    assert report["fcp_global_species_prevalence_estimated"] is False
    with (tmp_path/"output/species_rgfca_repeated_photo_polymorphism.csv").open() as f:
        rows=list(csv.DictReader(f))
    yes=next(x for x in rows if x["inat_taxon_id"]=="11")
    no=next(x for x in rows if x["inat_taxon_id"]=="12")
    assert float(yes["conditional_photo_colour_P_poly"])>=.8
    assert yes["repeatability_0_80_status"]=="repeatably_photo_polymorphic_phenotype"
    assert no["repeatability_0_80_status"]=="not_detected_in_evaluable_photo_resamples_not_monomorphic"
    assert no["genetic_polymorphism_established"]=="False"


def test_cannot_substitute_historical_or_tampered_photo_colour_receipt(tmp_path):
    cohort,sched,provenance,measured,receipt=setup(tmp_path)
    species,assignments=validate_schedule(cohort,sched,provenance,synthetic=True)
    doc=json.loads(receipt.read_text())
    doc["historical_measured_flower_colours_reused"]=True
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="No complete original fresh ROI-v4"):
        load_new_roi(measured,receipt,species,synthetic=True)


def test_fake_balanced_schedule_rejected(tmp_path):
    cohort,sched,provenance,measured,receipt=setup(tmp_path)
    lines=sched.read_text()
    sched.write_text(lines.replace("Species alwayswhite","Fake red flower species"))
    with pytest.raises(ValueError,match="nonoriginal or phenotype-adaptive"):
        validate_schedule(cohort,sched,provenance,synthetic=True)


def test_empty_photo_pool_never_yields_false_biological_negative(tmp_path):
    species={11:{"species":"Species redwhite","source_photo_opportunities":"1000"}}
    assignments=[{"realization":x,"inat_taxon_id":11} for x in range(1,51)]
    out=tmp_path/"missing"
    result=repeat(species,assignments,{},out,synthetic=True)
    assert result["n_repeated_phenotypic_polymorphism_species"]==0
    assert result["n_species_unresolved_or_not_evaluable"]==1


def test_after_39_of_96_empirical_roi_success_120photo_resampling_is_primary():
    photos=original_photos(11)
    primary=one_draw(11,"Species redwhite",1000,photos,1)
    original=one_draw(11,"Species redwhite",1000,photos,1,per_fold=40)
    assert primary["n_selected"]==120
    assert primary["fold_0_classifiable"]==60
    assert primary["fold_1_classifiable"]==60
    assert original["n_selected"]==80
    assert primary["evaluable"] and original["evaluable"]


def test_rgfca_40_80_120_128_sensitivities_do_not_redefine_taxon_universe(tmp_path):
    cohort,sched,provenance,measured,receipt=setup(tmp_path)
    species,assignments=validate_schedule(cohort,sched,provenance,synthetic=True)
    photos=load_new_roi(measured,receipt,species,synthetic=True)
    report=repeat(species,assignments,photos,tmp_path/"depths",synthetic=True)
    assert report["primary_120_photo_draw_and_60_photo_per_observer_fold"] is True
    assert report["sensitivity_photo_draw_sizes"]==[40,80,120,128]
    with (tmp_path/"depths/species_rgfca_repeated_photo_polymorphism.csv").open() as source:
        measured=list(csv.DictReader(source))
    assert len(measured)==2
    for record in measured:
        assert record["P_poly_primary_photo_draw"]=="120"
        counts=json.loads(record["draw_depth_sensitivity_counts"])
        assert set(counts)=={"40","80","120","128"}
        assert all(counts[t]["evaluable"]<=50 for t in counts)


def test_robust_photo_label_pair_can_be_single_inflorescence_and_must_be_held(tmp_path):
    taxon = 11
    photos = original_photos(taxon, variable=True)
    # A daisy-like flower repeatedly shows white ray + yellow disc tissue in
    # THE SAME ROI even though the original coarse photo class can alternate.
    for item in photos:
        if item["morph"] == "red_pink":
            item["morph"] = "yellow_orange"
        item["within_ROI_white_yellow"] = (.52, .48)
    assigned = [{"realization":r, "inat_taxon_id":taxon} for r in range(1,51)]
    registry = {taxon:{"species":"Synthetic composite-flower source",
                       "source_photo_opportunities":"1024"}}
    report = repeat(registry, assigned, {taxon:photos},
                    tmp_path / "hold_multicolour_photo_labels", synthetic=True)
    assert report["n_repeated_phenotypic_polymorphism_species"] == 0
    assert report["n_coarse_photo_label_recurring_but_structural_QC_HOLD"] == 1
    with (tmp_path / "hold_multicolour_photo_labels" /
          "species_rgfca_repeated_photo_polymorphism.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    subject = rows[0]
    assert subject["structural_organ_colour_review_hold"] == "True"
    assert float(subject["white_yellow_within_single_ROI_fraction"]) == 1.0
    assert float(subject["conditional_photo_colour_P_poly"]) >= .8
    assert subject["repeatability_0_80_status"] == (
        "repeatable_colour_LABELS_WITHIN_ROI_STRUCTURAL_QC_HOLD")
    assert subject["genetic_polymorphism_established"] == "False"


def test_experimental_original_roi_loader_preserves_within_photo_palette(tmp_path):
    cohort,sched,provenance,measured,receipt=setup(tmp_path)
    taxa,assignments=validate_schedule(cohort,sched,provenance,synthetic=True)
    photos=load_new_roi(measured,receipt,taxa,synthetic=True)
    assert all("within_ROI_white_yellow" in item for samples in photos.values()
               for item in samples)
    for taxon, samples in photos.items():
        for item in samples:
            if item["classifiable"]:
                assert len(item["within_ROI_white_yellow"]) == 2
