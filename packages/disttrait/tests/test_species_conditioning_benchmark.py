import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
BENCHMARK = Path(__file__).resolve().parents[1] / "benchmarks" / "species_conditioning_benchmark.py"
RECEIPT = ROOT / "results" / "disttrait_species_conditioning_benchmark_v0_2_20260918" / "result.json"


def _load_benchmark():
    spec = importlib.util.spec_from_file_location("disttrait_species_conditioning_benchmark", BENCHMARK)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_species_conditioning_benchmark_matches_frozen_receipt() -> None:
    module = _load_benchmark()
    observed = module.run_benchmark()
    frozen = json.loads(RECEIPT.read_text(encoding="utf-8"))

    assert observed["worlds_per_condition"] == frozen["design"]["worlds_per_condition"]

    for section in ("null", "signal"):
        for key, value in observed[section].items():
            assert value == pytest.approx(frozen[section][key], abs=1e-12, rel=0)

    assert observed["null"]["naive_pooled_false_positive_fraction"] == 1.0
    assert observed["null"]["species_conditioned_false_positive_fraction"] == 0.0
    assert observed["signal"]["species_conditioned_detection_fraction"] == 1.0
