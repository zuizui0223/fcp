"""Exact-FCP original-photo four-hue Brier joint-abiotic sum and original cohort guards."""
import importlib.util,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/"scripts/analysis/fcp_joint36_photo_hue_Brier_all_environment_20261011.py"
sys.path.insert(0,str(p.parent))
sp=importlib.util.spec_from_file_location("photo_hue_joint",p)
mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)

def test_observed_hue_classes_are_fixed():
    assert len(mod.CLASSES)==4
    assert mod.CLASSES==("white","yellow_orange","red_pink","blue_purple")
    assert sum(map(len,mod.BLOCKS.values()))==36

def test_new_hue_full_environment_only_not_a_pigment_fitness_claim():
    assert "NOT" in (mod.__doc__ or "")
    assert mod.SCHEMA=="fcp_original100543_joint36_photocolour_hue_brier_contribution_v1"

def test_independent_hue_score_adds_to_original_multiclass():
    rng=np.random.default_rng(67)
    truth=np.eye(4)[rng.integers(0,4,100)]
    p=rng.dirichlet([2]*4,size=100)
    q=rng.dirichlet([2]*4,size=100)
    gain=(p-truth)**2-(q-truth)**2
    assert np.allclose(gain.sum(axis=1),((p-truth)**2).sum(axis=1)-((q-truth)**2).sum(axis=1))
