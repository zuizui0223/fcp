import pytest

from scripts.data.build_inaturalist_locked_image_packet import development_passed_species


def test_development_passed_species_is_exact_and_sorted():
    gate = {
        "protocol": "fcp-inaturalist-automated-colour-state-v2",
        "status": "complete_location_free_automated_colour_development_gate",
        "species_passed": 2,
        "species_results": [
            {"canonical_name": "B", "development_gate_status": "pass"},
            {"canonical_name": "A", "development_gate_status": "pass"},
            {"canonical_name": "C", "development_gate_status": "not_evaluable"},
        ],
    }
    assert development_passed_species(gate) == ["A", "B"]


def test_no_passing_species_keeps_locked_packet_closed():
    gate = {
        "protocol": "fcp-inaturalist-automated-colour-state-v2",
        "status": "complete_location_free_automated_colour_development_gate",
        "species_passed": 0,
        "species_results": [],
    }
    with pytest.raises(RuntimeError, match="must remain unopened"):
        development_passed_species(gate)
