"""Source-blind pure tests for an actual six-taxon ROI photo technical audit.

Tests use synthetic ROI records: biological source result is checked separately
by the pinned GitHub Actions six-taxon original-image artifact.
"""
from __future__ import annotations
import json
import pytest

from scripts.analysis.audit_rgfca_real_six_taxa_photo_repeatability_20261016 import (
    DEPTHS, _draw, repeat_descriptive,
)
from fcp_pipeline.photo_polymorphism_screen import observer_split


def make_photos(taxon=117, *, mixed=True, fill=True):
    photos=[]
    for fold in (0,1):
        people=[]
        idx=0
        while len(people)<5:
            identity=f"observer_{taxon}_{idx}"
            if observer_split(taxon,identity)==fold:
                people.append(identity)
            idx+=1
        for j in range(48):
            state=("yellow_orange" if j%5==0 else "white") if mixed else "white"
            frac={"white":0.52,"yellow":0.48} if state in ("white","yellow_orange") else {}
            photos.append({
                "photo_id":str(taxon*1000+fold*50+j),
                "fold":fold,
                "observation_id":f"obs_{taxon}_{fold}_{j}",
                "observer_id":people[j%5],
                "classifiable":fill,
                "morph":state if fill else None,
                "fractions":frac,
            })
    return photos


def test_bounded_rgfca_real_photo_label_resampling_not_genetic_certification():
    photos=make_photos()
    rep=repeat_descriptive(117,"Synthetic multi-colour head",photos,repeat_count=30)
    assert rep["pilot_photos_selected"]==96
    assert rep["pilot_original_flower_ROI_classifiable"]==96
    assert rep["same_photo_white_and_yellow_both_ge_0_10"]==96
    assert rep["photo_colour_structural_confounding_high"] is True
    assert rep["within_flower_architecture_QC_required"] is True
    assert rep["strict_image_FCP_species_confirmed"] is False
    assert rep["depth_40_evaluable_repetitions"]==30
    assert rep["depth_40_positive_photo_label_repetitions"]>=22
    assert rep["depth_48_conditional_photo_label_detection"]==1.0
    assert rep["depth_48_realized_photo_ID_overlap_fraction"]==1.0
    assert rep["depth_48_distinct_source_photo_ID_draw_sets"]==1
    assert rep["depth_40_realized_photo_ID_overlap_fraction"]<1.0
    assert rep["depth_40_distinct_source_photo_ID_draw_sets"]>1


def test_single_coarse_morph_not_monomorphic():
    rep=repeat_descriptive(117,"One photo state",make_photos(mixed=False),repeat_count=25)
    assert rep["strict_image_FCP_species_confirmed"] is False
    assert rep["photo_colour_structural_confounding_high"] is False
    assert rep["depth_40_evaluable_repetitions"]==25
    assert rep["depth_40_conditional_photo_label_detection"]==0.0


def test_unclassifiable_images_keep_evaluable_denominator_zero():
    rep=repeat_descriptive(117,"Only failed ROI",make_photos(fill=False),repeat_count=10)
    assert rep["pilot_original_flower_ROI_classifiable"]==0
    assert rep["depth_40_evaluable_repetitions"]==0
    assert rep["depth_40_conditional_photo_label_detection"]==""
    assert rep["strict_image_FCP_species_confirmed"] is False


def test_repeatable_and_declared_depths_only():
    p=make_photos()
    a=repeat_descriptive(117,"Same original photos",p,repeat_count=20)
    b=repeat_descriptive(117,"Same original photos",p,repeat_count=20)
    assert a==b
    assert DEPTHS==(20,32,40,44,48)
    with pytest.raises(ValueError,match="predeclared"):
        _draw(117,"Same original photos",p,1,60)


def test_photo_class_repeatability_does_not_turn_api_snapshot_into_dated_bulk():
    r=repeat_descriptive(117,"S",make_photos(),repeat_count=4)
    assert r["unreconciled_date_snapshot"] is True
    assert r["genetic_or_population_FCP_proven"] is False
