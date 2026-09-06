import numpy as np
import pandas as pd

from fcp_pipeline.global_g5_matching import (
    N_CELLS,
    build_matched_sets,
    cell_center,
    dilate_cells,
    equal_area_cell_id,
    is_allopatric_after_dilation,
    is_sympatric,
    taxonomic_distance_class,
)


def test_equal_area_roundtrip_centres():
    for cell in [0, 1, 319, 320, 16000, N_CELLS - 1]:
        lat, lon = cell_center(cell)
        assert equal_area_cell_id(lat, lon) == cell


def test_longitude_wrap_dilation():
    west = 10 * 320
    east = 10 * 320 + 319
    assert east in dilate_cells({west})
    assert west in dilate_cells({east})


def test_sympatry_and_dilated_allopatry_rules():
    a = {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009}
    b = {1000, 1001, 1002, 1003, 1004, 1010, 1011, 1012, 1013, 1014}
    yes, shared, jaccard = is_sympatric(a, b)
    assert yes
    assert shared == 5
    assert jaccard >= 0.10
    assert not is_allopatric_after_dilation(a, b)
    far = {40000, 40001, 40002}
    assert is_allopatric_after_dilation(a, far)


def test_taxonomic_classes():
    assert taxonomic_distance_class("A", "F", "A", "F") == "congeneric"
    assert taxonomic_distance_class("A", "F", "B", "F") == "same_family_noncongeneric"
    assert taxonomic_distance_class("A", "F", "C", "G") == "different_family"


def test_build_matched_sets_obeys_exact_and_distance_matches():
    species = ["Focal", "Sym", "Ctl1", "Ctl2", "WrongRealm"]
    meta = pd.DataFrame(
        {
            "species": species,
            "accepted_genus": ["A", "B", "C", "D", "E"],
            "accepted_family": ["Fam", "Other", "Third", "Fourth", "Fifth"],
            "dominant_realm": ["R0", "R1", "R1", "R1", "R2"],
            "z_bio01": [0.0, 0.2, 0.21, 0.19, 0.2],
            "z_bio04": [0.0, -0.1, -0.11, -0.09, -0.1],
            "z_bio12": [0.0, 0.3, 0.31, 0.29, 0.3],
            "z_bio15": [0.0, 0.0, 0.01, -0.01, 0.0],
            "z_log1p_occupied_cells": [0.0, 0.2, 0.21, 0.19, 0.2],
            "z_log1p_gbif_occurrences": [0.0, 0.1, 0.11, 0.09, 0.1],
            "matching_eligible": [True] * 5,
        }
    )
    focal_cells = set(range(1000, 1020))
    sym_cells = set(range(1000, 1010)) | set(range(2000, 2010))
    occupancy = {
        "Focal": focal_cells,
        "Sym": sym_cells,
        "Ctl1": set(range(30000, 30020)),
        "Ctl2": set(range(35000, 35020)),
        "WrongRealm": set(range(40000, 40020)),
    }
    out = build_matched_sets(meta, occupancy, maximum_distance=1.5, controls_target=5, controls_minimum=2)
    focal_sym = out.loc[(out["focal_species"] == "Focal") & (out["sympatric_partner"] == "Sym")]
    assert len(focal_sym) == 2
    assert set(focal_sym["control_species"]) == {"Ctl1", "Ctl2"}
    assert np.all(focal_sym["matching_distance"].to_numpy() < 1.5)
