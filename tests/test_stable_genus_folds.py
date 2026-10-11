from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
from fcp_joint_environment_taxongeo_stable_20261011 import StableGenusKFold


def assignments(groups,k=3):
    groups=np.asarray(groups);out={}
    for f,(tr,te) in enumerate(StableGenusKFold(k).split(np.zeros((len(groups),1)),groups=groups)):
        assert not set(groups[tr])&set(groups[te])
        for g in groups[te]:out[str(g)]=f
    return out

def test_lexical_ties_and_lowest_fold_ties_are_explicit():
    assert assignments(['D','D','C','C','B','B','A','A'],3)=={'A':0,'B':1,'C':2,'D':0}

def test_input_order_and_outcome_values_never_change_genus_fold():
    rng=np.random.default_rng(81)
    groups=np.repeat([f'G{i:03}' for i in range(103)],np.tile([1,2,3,4],26)[:103])
    assert assignments(groups,5)==assignments(rng.permutation(groups),5)

def test_too_few_genera_is_explicit_failure():
    with pytest.raises(ValueError):list(StableGenusKFold(5).split(np.zeros((4,1)),groups=['A']*4))
