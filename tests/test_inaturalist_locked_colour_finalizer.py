from scripts.data.finalize_inaturalist_locked_colour_table import locked_species_gate


def rows(admitted: int, background: int):
    values = []
    for index in range(120):
        is_admitted = index < admitted
        values.append(
            {
                "encounter_status": (
                    "automated_colour_state_admitted"
                    if is_admitted
                    else "automated_colour_state_not_evaluable"
                ),
                "background_control_status": (
                    "background_control_available"
                    if index < background
                    else "background_control_not_evaluable"
                ),
            }
        )
    return values


def test_locked_gate_passes_exact_70_percent_boundary():
    result = locked_species_gate(rows(84, 59))
    assert result["locked_gate_status"] == "pass"
    assert result["coordinates_joined"] is True


def test_locked_gate_withholds_coordinates_below_admission_threshold():
    result = locked_species_gate(rows(83, 83))
    assert result["locked_gate_status"] == "not_evaluable"
    assert result["coordinates_joined"] is False


def test_locked_gate_withholds_coordinates_below_background_threshold():
    result = locked_species_gate(rows(100, 69))
    assert result["locked_gate_status"] == "not_evaluable"
    assert result["coordinates_joined"] is False
