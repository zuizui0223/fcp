"""Synthetic and independent-design regression tests for PHOTO-FIRST RGFCA reuse.

No real flower-colour phenotype, older 149,900 colour outcomes, or source-image
pixels are provided to these tests. The expected result is a species/sample
INCLUSION schedule, not a biological flower polymorphism result.
"""
from __future__ import annotations
import csv
import json
from collections import Counter

import pytest

from scripts.analysis.freeze_fcp_rgfca_photo_repetition import (
    SALT, STRATA, choose_and_schedule, order_hash, output,
)


def frame():
    rows = []
    tid = 1000
    for key, lo, hi, quota in STRATA:
        for j in range(quota + 21):
            rows.append({"inat_taxon_id": tid, "species": f"Genus{tid} epithet{j}",
                         "n_photo": lo + min(j, 10)})
            tid += 1
    return rows


def test_rgfca_each_realization_250_and_each_taxon_exactly_50():
    species, schedule, diagnostics = choose_and_schedule(frame())
    assert len(species) == 1000
    assert len(schedule) == 200 * 250
    assert len({s["inat_taxon_id"] for s in species}) == 1000
    counts = Counter(s["inat_taxon_id"] for s in schedule)
    assert set(counts.values()) == {50}
    for rep in range(1, 201):
        sample = [x for x in schedule if x["realization"] == rep]
        assert len(sample) == 250
        assert len({x["inat_taxon_id"] for x in sample}) == 250
        assert len({x["slot"] for x in sample}) == 250
    assert diagnostics["coverage_checkpoints"]["100"] == {
        "n_species_seen": 1000, "min_inclusions": 25, "max_inclusions": 25}
    assert diagnostics["coverage_checkpoints"]["200"] == {
        "n_species_seen": 1000, "min_inclusions": 50, "max_inclusions": 50}


def test_species_equal_schedule_is_independent_of_source_order():
    selected_a, assignments_a, _ = choose_and_schedule(frame())
    selected_b, assignments_b, _ = choose_and_schedule(list(reversed(frame())))
    assert selected_a == selected_b
    assert assignments_a == assignments_b


def test_inclusion_is_metadata_only_and_not_colour_adaptive():
    baseline = frame()
    inserted_fake_colours = [{**r, "morph": "red_pink" if i % 2 else "white",
                              "temperature": 31 if i % 2 else -8,
                              "literature_C_S_state": "unknown"}
                             for i, r in enumerate(baseline)]
    a, aa, _ = choose_and_schedule(baseline)
    b, bb, _ = choose_and_schedule(inserted_fake_colours)
    assert [(x["inat_taxon_id"], x["tier"]) for x in a] == [
        (x["inat_taxon_id"], x["tier"]) for x in b]
    assert aa == bb
    assert "morph" not in aa[0] and "temperature" not in aa[0]


def test_underfilled_photo_opportunity_stratum_is_rejected():
    sparse = frame()[:500]
    with pytest.raises(ValueError, match="stratum underfilled"):
        choose_and_schedule(sparse)


def test_output_is_source_only_and_immutable(tmp_path):
    outdir = tmp_path / "rgfca"
    result = output(frame(), outdir)
    assert result["selected_species"] == 1000
    assert result["realizations"] == 200
    assert result["species_inclusions_each"] == 50
    assert result["species_sampling_imbalance"] == 0
    assert result["maximum_unique_original_image_pixels_if_pool_filled"] == 128000
    assert result["new_flower_images_opened"] is False
    assert result["photographic_poly_probabilities_measured"] is False
    assert result["species_not_enrolled_are_NOT_monomorphic_controls"] is True
    with (outdir / "rgfca_200x250_balanced_species_schedule.csv").open() as stream:
        assert sum(1 for _ in csv.DictReader(stream)) == 50000
    with pytest.raises(ValueError, match="Do not overwrite"):
        output(frame(), outdir)


def test_25_50_100_150_200_mc_convergence_bookkeeping():
    _, _, d = choose_and_schedule(frame())
    c = d["coverage_checkpoints"]
    assert {x: (c[str(x)]["min_inclusions"], c[str(x)]["max_inclusions"])
            for x in (25, 50, 100, 150, 200)} == {
                25: (6, 7), 50: (12, 13), 100: (25, 25),
                150: (37, 38), 200: (50, 50)}
    assert order_hash("enrolment", 1234, "Genus1234") == order_hash(
        "enrolment", 1234, "Genus1234")
