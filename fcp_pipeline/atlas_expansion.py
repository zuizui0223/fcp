"""Prospective repeated-cohort and ordered-inference contract for the FCP atlas.

The functions in this module operate on metadata, frozen branch statistics, or
simulation output.  They deliberately have no image-loading interface.  This keeps
cohort construction and multiplicity rules upstream of atlas colour measurement.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import math
from typing import Any

import numpy as np


PROTOCOL = "jbi-image-first-global-flower-colour-atlas-expansion-v2"
BRANCHES = (
    "shared_geographic_concentration",
    "environmental_concordance",
    "pollinator_biogeographic_concordance",
)
TERMINAL_STATES = ("supported", "not_supported", "not_evaluable")


def _stable_hash(salt: str, *parts: object) -> str:
    payload = "\x1f".join([salt, *(str(part) for part in parts)]).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_expansion_contract(contract: Mapping[str, Any]) -> None:
    """Fail closed if scale-out or pivot rules are outcome-adaptive."""

    if contract.get("protocol") != PROTOCOL:
        raise ValueError(f"unexpected expansion protocol: {contract.get('protocol')!r}")
    if contract.get("status") != (
        "frozen_after_three_species_negative_validation_and_before_any_atlas_candidate_pixels"
    ):
        raise ValueError("expansion contract is not frozen before atlas pixels")

    firewall = contract.get("outcome_firewall", {})
    forbidden_true = (
        "atlas_candidate_pixels_opened",
        "atlas_colour_fields_opened",
        "environmental_layers_joined_to_colour",
        "pollinator_layers_joined_to_colour",
        "cohort_count_selected_from_colour",
        "early_stopping_permitted",
    )
    opened = [key for key in forbidden_true if firewall.get(key) is not False]
    if opened:
        raise ValueError(f"outcome firewall is open for: {opened}")

    qualification = contract.get("estimator_qualification", {})
    if qualification.get("must_pass_before_atlas_pixels") is not True:
        raise ValueError("estimator qualification must precede atlas pixels")
    estimator = qualification.get("estimator", {})
    if estimator.get("retuning_on_benchmark") is not False:
        raise ValueError("benchmark retuning must be prohibited")
    benchmark = qualification.get("independent_roi_benchmark", {})
    if benchmark.get("selection") != (
        "all image IDs with an official trimap; no image-dependent subsampling"
    ):
        raise ValueError("ROI benchmark must score the complete official trimap set")
    if int(benchmark.get("foreground_label", -1)) != 1:
        raise ValueError("Oxford ROI foreground palette decoding changed")
    if list(benchmark.get("background_labels", [])) != [2, 3, 4]:
        raise ValueError("Oxford ROI background palette decoding changed")
    if int(benchmark.get("unlabelled_label", -1)) != 0:
        raise ValueError("Oxford ROI unlabelled palette decoding changed")
    if int(benchmark.get("minimum_scored_images", 0)) < 750:
        raise ValueError("ROI benchmark is too small for the frozen Oxford set")
    if benchmark.get("failure_rule", "").startswith("STOP") is False:
        raise ValueError("ROI failure must stop atlas pixel opening")

    recovery = qualification.get("signal_recovery", {})
    if int(recovery.get("simulation_repetitions", 0)) != 100:
        raise ValueError("signal recovery must use 100 repetitions")
    if int(recovery.get("permutations_per_repetition", 0)) != 999:
        raise ValueError("signal recovery must use 999 permutations per repetition")
    if list(recovery.get("effect_sizes", [])) != [0.0, 0.5, 1.0, 2.0]:
        raise ValueError("signal-recovery effect sizes changed")

    scaleout = contract.get("random_cohort_scaleout", {})
    expected = {
        "cohorts": 8,
        "species_per_cohort": 25,
        "observations_per_species": 300,
        "total_species": 200,
        "total_observations": 60_000,
    }
    changed = {
        key: scaleout.get(key)
        for key, value in expected.items()
        if int(scaleout.get(key, -1)) != value
    }
    if changed:
        raise ValueError(f"random-cohort dimensions changed: {changed}")
    if scaleout.get("species_overlap") != "none":
        raise ValueError("random cohorts must be species-disjoint")
    if scaleout.get("all_cohorts_required") is not True:
        raise ValueError("all frozen cohorts must be run")
    if scaleout.get("early_stopping") is not False:
        raise ValueError("early stopping is prohibited")
    if scaleout.get("single_nested_null") is not True:
        raise ValueError("cohort repetitions must enter one nested null")

    spatial = contract.get("spatial_design", {})
    if int(spatial.get("primary_scale_km", 0)) != 100:
        raise ValueError("100 km must remain the geometry-selected primary scale")
    if list(spatial.get("mandatory_sensitivity_scales_km", [])) != [250, 500]:
        raise ValueError("250/500-km sensitivities are mandatory")

    inference = contract.get("ordered_inference", {})
    if tuple(inference.get("branches", [])) != BRANCHES:
        raise ValueError("ordered inference branches changed")
    if tuple(inference.get("branch_outcomes", [])) != TERMINAL_STATES:
        raise ValueError("terminal outcome vocabulary changed")
    if "joint maximum statistic" not in str(inference.get("multiplicity", "")):
        raise ValueError("one joint maximum-statistic null is required")
    if int(inference.get("minimum_final_null_replicates", 0)) < 9_999:
        raise ValueError("final nested null must use at least 9,999 replicates")

    pollinator = contract.get("pollinator_overlay", {})
    if pollinator.get("ready_made_global_layer_available") is not False:
        raise ValueError("contract must not assert a ready-made global Bombus layer")
    if pollinator.get("freeze_before_colour_join") is not True:
        raise ValueError("pollinator regionalization must be colour-blind")

    publication = contract.get("publication_stop", {})
    if tuple(publication.get("terminal_states", [])) != TERMINAL_STATES:
        raise ValueError("publication must retain every terminal state")
    if publication.get("no_significance_chasing") is not True:
        raise ValueError("significance chasing must be explicitly prohibited")


def draw_disjoint_species_cohorts(
    eligible_species: Sequence[Mapping[str, Any]],
    contract: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Create all eight metadata-only species panels from one stable permutation."""

    validate_expansion_contract(contract)
    scaleout = contract["random_cohort_scaleout"]
    n_cohorts = int(scaleout["cohorts"])
    per_cohort = int(scaleout["species_per_cohort"])
    needed = n_cohorts * per_cohort
    salt = str(scaleout["stable_hash_salt"])

    forbidden_tokens = ("colour", "color", "pixel", "roi", "transition", "effect", "p_value")
    normalized: list[dict[str, Any]] = []
    seen_taxa: set[str] = set()
    seen_genera: set[str] = set()
    for raw in eligible_species:
        lowered = {str(key).casefold() for key in raw}
        leaked = sorted(
            key for key in lowered if any(token in key for token in forbidden_tokens)
        )
        if leaked:
            raise ValueError(f"colour/image outcome fields reached cohort sampling: {leaked}")
        taxon_id = str(raw.get("taxon_id", "")).strip()
        species = str(raw.get("species", "")).strip()
        genus = str(raw.get("genus", "")).strip()
        if not taxon_id or not species or not genus:
            raise ValueError("eligible rows require taxon_id, species, and genus")
        if taxon_id in seen_taxa:
            raise ValueError(f"duplicate eligible taxon_id: {taxon_id}")
        seen_taxa.add(taxon_id)
        if genus.casefold() in seen_genera:
            continue
        seen_genera.add(genus.casefold())
        normalized.append(
            {
                "taxon_id": taxon_id,
                "species": species,
                "genus": genus,
                "selection_hash": _stable_hash(salt, taxon_id, species),
            }
        )

    normalized.sort(key=lambda row: (row["selection_hash"], row["taxon_id"]))
    if len(normalized) < needed:
        raise ValueError(
            f"not_evaluable: need {needed} genus-distinct eligible species, found {len(normalized)}"
        )

    panels: list[dict[str, Any]] = []
    for index, row in enumerate(normalized[:needed]):
        cohort_index = index // per_cohort + 1
        panels.append(
            {
                **row,
                "cohort_id": f"C{cohort_index:02d}",
                "cohort_species_index": index % per_cohort + 1,
                "target_observations": int(scaleout["observations_per_species"]),
            }
        )
    return panels


def monte_carlo_p_value(
    observed: float,
    null: Sequence[float] | np.ndarray,
    *,
    alternative: str = "greater",
) -> float:
    """Return the non-zero Phipson-Smyth Monte Carlo p-value."""

    observed = float(observed)
    values = np.asarray(null, dtype=float)
    if not math.isfinite(observed) or values.ndim != 1 or values.size < 1:
        raise ValueError("observed and null statistics must be finite and non-empty")
    if not np.isfinite(values).all():
        raise ValueError("null statistics must be finite")
    if alternative == "greater":
        exceed = int(np.count_nonzero(values >= observed))
    elif alternative == "less":
        exceed = int(np.count_nonzero(values <= observed))
    else:
        raise ValueError("alternative must be 'greater' or 'less'")
    return float((exceed + 1) / (values.size + 1))


def joint_max_adjusted_p_values(
    observed: Mapping[str, float],
    null_statistics: Mapping[str, Sequence[float] | np.ndarray],
) -> dict[str, float]:
    """Familywise p-values from one joint, complete-pipeline maximum null."""

    if not observed:
        raise ValueError("at least one evaluable branch is required")
    names = tuple(observed)
    if set(names) != set(null_statistics):
        raise ValueError("observed and null branches must match exactly")
    arrays = [np.asarray(null_statistics[name], dtype=float) for name in names]
    sizes = {array.size for array in arrays}
    if len(sizes) != 1 or next(iter(sizes), 0) < 1:
        raise ValueError("all branch null vectors must have equal positive length")
    if any(array.ndim != 1 or not np.isfinite(array).all() for array in arrays):
        raise ValueError("branch null vectors must be finite one-dimensional arrays")
    if any(not math.isfinite(float(observed[name])) for name in names):
        raise ValueError("observed branch statistics must be finite")
    max_null = np.max(np.vstack(arrays), axis=0)
    return {
        name: monte_carlo_p_value(float(observed[name]), max_null, alternative="greater")
        for name in names
    }


def ordered_branch_decisions(
    adjusted_p: Mapping[str, float],
    evaluable: Mapping[str, bool],
    *,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Retain all branch outcomes and identify the first supported conclusion."""

    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")
    unknown = (set(adjusted_p) | set(evaluable)) - set(BRANCHES)
    if unknown:
        raise ValueError(f"unknown inference branches: {sorted(unknown)}")

    outcomes: dict[str, str] = {}
    promoted: str | None = None
    for branch in BRANCHES:
        is_evaluable = bool(evaluable.get(branch, False))
        if not is_evaluable:
            outcomes[branch] = "not_evaluable"
            continue
        p_value = float(adjusted_p.get(branch, math.nan))
        if not 0.0 <= p_value <= 1.0:
            raise ValueError(f"missing or invalid adjusted p-value for {branch}")
        outcome = "supported" if p_value <= alpha else "not_supported"
        outcomes[branch] = outcome
        if promoted is None and outcome == "supported":
            promoted = branch

    return {
        "alpha": alpha,
        "branch_order": list(BRANCHES),
        "outcomes": outcomes,
        "promoted_conclusion": promoted or "complete_negative_or_not_evaluable_tree",
    }
