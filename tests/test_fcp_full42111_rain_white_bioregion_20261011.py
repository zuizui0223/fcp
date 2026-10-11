"""Synthetic guards for original FCP global species and genus-conditional rainfall hue models."""
from pathlib import Path
import importlib.util,sys
import numpy as np,pandas as pd,pytest
p=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_fcp_full42111_rain_white_bioregion_20261011.py"
sys.path.insert(0,str(p.parent))
spec=importlib.util.spec_from_file_location("fcp_realm",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def frame():
    rng=np.random.default_rng(43)
    rows=[]
    for genus in range(35):
        for j in range(12):
            x=rng.normal(size=8)
            white=.35+.1*x[0]+rng.normal(scale=.01)
            y=np.array([white,(1-white)*.5,(1-white)*.3,(1-white)*.2])
            rows.append({"inat_taxon_id":genus*100+j,"genus":f"G{genus}","region":"0_tropical" if genus%2 else "3_extratropical",
                "photo_id":genus*100+j,"n_photos":1,**dict(zip(m.COLS,x)),**dict(zip(m.HUES,y))})
    return pd.DataFrame(rows)

def test_original_four_label_definitions():
    assert m.HUES==("white","yellow_orange","red_pink","blue_purple")
    assert m.COLS[0]=="wc_bio12"

def test_original_source_hash_fails_closed(tmp_path):
    q=tmp_path/"fake.zip";q.write_bytes(b"wrong")
    with pytest.raises(RuntimeError,match="SHA256"):m.get_source(q)

def test_genus_centered_model_detects_synthetic_photo_white_rainfall():
    x=frame()
    a=m.fit(x)
    assert a["status"]=="OK"
    assert a["n_species_or_taxon_region_units"]==420
    assert a["n_genera"]==35
    assert a["rainfall_hue"][0]["coef"]>.05
    assert a["rainfall_hue"][0]["CI95"][0]>0

def test_genus_equal_weight_and_colour_simplex():
    x=frame()
    assert np.allclose(x[list(m.HUES)].sum(axis=1),1)
    assert m.fit(x,True)["rainfall_hue"][0]["coef"]>0
    assert abs(sum(r["coef"] for r in m.fit(x)["rainfall_hue"]))<1e-10

def test_source_species_cannot_replicate_from_two_geographic_photo_cells():
    x=frame()
    copies=x.iloc[:17].copy()
    copies["region"]="3_tropical"
    d=pd.concat([x,copies],ignore_index=True)
    scale=np.ones(8);mean=np.zeros(8)
    global_species=m.species_means(d,mean,scale,by_region=False)
    within_region=m.species_means(d,mean,scale,by_region=True)
    assert len(global_species)==420
    assert len(within_region)==437
    assert not global_species.inat_taxon_id.duplicated().any()
    assert global_species.n_photos.sum()==len(d)
