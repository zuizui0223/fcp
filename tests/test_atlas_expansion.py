from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.atlas_expansion import (
    BRANCHES,
    draw_disjoint_species_cohorts,
    joint_max_adjusted_p_values,
    monte_carlo_p_value,
    ordered_branch_decisions,
    validate_expansion_contract,
)


CONTRACT_PATH = Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json")


def contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def test_committed_expansion_contract_is_prospective_and_fail_closed() -> None:
    validate_expansion_contract(contract())

    changed = deepcopy(contract())
    changed["random_cohort_scaleout"]["early_stopping"] = True
    with pytest.raises(ValueError, match="early stopping"):
        validate_expansion_contract(changed)

    changed = deepcopy(contract())
    changed["ordered_inference"]["branches"] = list(reversed(BRANCHES))
    with pytest.raises(ValueError, match="branches changed"):
        validate_expansion_contract(changed)

    changed = deepcopy(contract())
    changed["estimator_qualification"]["independent_roi_benchmark"][
        "background_labels"
    ] = [2]
    with pytest.raises(ValueError, match="background palette decoding"):
        validate_expansion_contract(changed)


def test_random_cohorts_are_deterministic_disjoint_and_genus_capped() -> None:
    rows = [
        {"taxon_id": index, "species": f"Species {index}", "genus": f"Genus{index}"}
        for index in range(230)
    ]
    first = draw_disjoint_species_cohorts(rows, contract())
    second = draw_disjoint_species_cohorts(list(reversed(rows)), contract())

    assert first == second
    assert len(first) == 200
    assert len({row["taxon_id"] for row in first}) == 200
    assert len({row["genus"] for row in first}) == 200
    assert {row["cohort_id"] for row in first} == {f"C{i:02d}" for i in range(1, 9)}
    assert all(sum(row["cohort_id"] == cohort for row in first) == 25 for cohort in {row["cohort_id"] for row in first})
    assert all(row["target_observations"] == 300 for row in first)


def test_random_cohort_builder_rejects_outcome_leakage_and_shortfall() -> None:
    rows = [
        {"taxon_id": index, "species": f"Species {index}", "genus": f"Genus{index}"}
        for index in range(200)
    ]
    rows[0]["colour_effect"] = 0.4
    with pytest.raises(ValueError, match="outcome fields"):
        draw_disjoint_species_cohorts(rows, contract())

    clean_rows = [
        {"taxon_id": index, "species": f"Species {index}", "genus": f"Genus{index}"}
        for index in range(199)
    ]
    with pytest.raises(ValueError, match="not_evaluable"):
        draw_disjoint_species_cohorts(clean_rows, contract())


def test_monte_carlo_p_value_is_nonzero_and_uses_declared_tail() -> None:
    null = np.array([0.1, 0.2, 0.3, 0.4])
    assert monte_carlo_p_value(0.5, null) == pytest.approx(0.2)
    assert monte_carlo_p_value(0.05, null, alternative="less") == pytest.approx(0.2)


def test_joint_max_null_adjusts_every_evaluable_branch_together() -> None:
    observed = {"a": 0.8, "b": 0.7}
    null = {
        "a": np.array([0.1, 0.9, 0.4, 0.2]),
        "b": np.array([0.6, 0.3, 0.5, 0.1]),
    }
    adjusted = joint_max_adjusted_p_values(observed, null)
    assert adjusted["a"] == pytest.approx(0.4)
    assert adjusted["b"] == pytest.approx(0.4)


def test_ordered_decisions_preserve_negative_and_not_evaluable_branches() -> None:
    decisions = ordered_branch_decisions(
        {
            "shared_geographic_concentration": 0.3,
            "environmental_concordance": 0.02,
        },
        {
            "shared_geographic_concentration": True,
            "environmental_concordance": True,
            "pollinator_biogeographic_concordance": False,
        },
    )
    assert decisions["outcomes"] == {
        "shared_geographic_concentration": "not_supported",
        "environmental_concordance": "supported",
        "pollinator_biogeographic_concordance": "not_evaluable",
    }
    assert decisions["promoted_conclusion"] == "environmental_concordance"
