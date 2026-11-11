"""Synthetic guards for classification-opportunity geographic negative-control limitations."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

script=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_polymorphism_classifiability_geography_20261010.py"
spec=importlib.util.spec_from_file_location("audit_classifiability",script)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def photos(classes,positions,dates=None):
    n=len(classes)
    dates=dates or ["2023-05-03"]*n
    return pd.DataFrame({"inat_taxon_id":[1]*n,"species":["A b"]*n,
       "photo_id":np.arange(n),"latitude":[0.]*n,"longitude":positions,
       "observed_on":dates,"morph":["white" if c else "" for c in classes],
       "_classifiable":np.asarray(classes,int)})

def test_spatial_failure_clustering_is_measured_separately_from_colour():
    g=photos([1]*50+[0]*50,[0]*50+[2]*50)
    row,null=m.measure(g,"discovery")
    assert row["technical_depletion"]>0.45
    assert row["colour_depletion"]==0
    assert abs(null["unconditional"].mean())<0.08

def test_fully_classifiable_species_has_zero_control_signal():
    row,null=m.measure(photos([1]*100,[0]*50+[2]*50),"discovery")
    assert row["technical_depletion"]==0
    assert np.max(np.abs(null["month"]))==0

def test_month_stratification_preserves_exact_photographic_status():
    g=photos([1]*50+[0]*50,[0]*50+[2]*50,
             ["2023-05-01"]*50+["2023-08-01"]*50)
    row,null=m.measure(g,"discovery")
    assert row["technical_month_exchangeable_photos"]==0
    assert np.allclose(null["month"],row["technical_depletion"])
    assert null["unconditional"].mean()<row["technical_depletion"]

def test_unstructured_status_yields_no_large_local_effect():
    row,null=m.measure(photos([1,0]*50,[0]*50+[2]*50),"validation")
    assert abs(row["technical_depletion"])<0.03
    assert abs(null["unconditional"].mean())<0.08
