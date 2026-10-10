"""Guard original 42,111 species and all seven physical source covariate blocks."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import run_fcp_global42111_allenv_spatial_phylo_20261010 as M

def mock_original():
    d=pd.DataFrame({
        "inat_taxon_id":np.arange(1,42112),
        "measurement_status":["classified_four_state_morph"]*18457+["unclassified"]*23654
    })
    for features in M.predictors.BLOCKS.values():
        for feature in features:
            d[feature]=1.
    return d

def test_all_predesignated_seven_blocks_have_36_source_features():
    assert len(M.predictors.BLOCKS)==7
    assert sum(len(v) for v in M.predictors.BLOCKS.values())==36
    assert all(k in M.predictors.BLOCKS for k in
               ("elevation","temperature","precipitation","solar_radiation",
                "wind","vapor_pressure","soil"))

def test_baseline_uses_true_phylo_spatial_kernel_and_restores_block_defs(monkeypatch):
    original=mock_original()
    previous=M.baseline.BLOCKS
    def fake_baseline(source,ledger,trees,nboot):
        assert len(source)==42111
        assert nboot==19
        assert M.baseline.BLOCKS==M.predictors.BLOCKS
        assert set(trees)=={250,500}
        return {"full_872_1761_direct_tip_phylogenetic_coverage_HOLD":True,
                "cohorts":{}}
    monkeypatch.setattr(M.baseline,"run",fake_baseline)
    result=M.analyze(original,pd.DataFrame(),Path("fake_250.tre"),
                     Path("fake_500.tre"),nboot=19)
    assert result["schema"]==M.SCHEMA
    assert result["n_named_predictors"]==36
    assert result["original_spatial_and_phylogenetic_covariance_used"]
    assert M.baseline.BLOCKS is previous

def test_source_denominator_or_missing_col_never_imputed(monkeypatch):
    source=mock_original()
    with pytest.raises(ValueError,match="Original photo species"):
        M.analyze(source.iloc[:50],pd.DataFrame(),Path("t250"),Path("t500"))
    altered=source.drop(columns=["wc_srad_annual_mean"])
    with pytest.raises(ValueError,match="Missing a source-environment"):
        M.analyze(altered,pd.DataFrame(),Path("t250"),Path("t500"))
