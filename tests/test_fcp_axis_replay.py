import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fcp_pipeline import fcp_axis_replay as axis


def test_axial_fixed_white_mean_and_exact_source_identity(tmp_path,monkeypatch):
    source=tmp_path/'frozen.py'
    source.write_text('''
import pandas as pd
BIO=['white','yellow','orange','red','pink','magenta','purple','blue','bronze']
MIN_CLASSIFIABLE=40
def load_classifiable(path):
    return pd.read_csv(path)
def species_delta_vectors(work,threshold):
    delta={f'delta_{c}': 0.0 for c in BIO}
    delta['delta_white']=1.0;delta['delta_red']=-1.0
    delta.update(species='Sp1',n_classifiable=len(work))
    return pd.DataFrame([delta]), None, {'nclass_ge40':1,'H1_coarse_second_pass':1,'continuous_minor_fraction_pass':1,'nonzero_delta_pass':1}
''')
    h=axis.git_blob(source)
    monkeypatch.setattr(axis,'SOURCE_BLOB',h)
    photos=tmp_path/'p.csv';pd.DataFrame({'x':range(75)}).to_csv(photos,index=False)
    taxa=tmp_path/'s.csv';pd.DataFrame({'inat_taxon_id':[7],'species':['Sp1']}).to_csv(taxa,index=False)
    receipt=axis.run(photos,taxa,source,tmp_path/'out')
    assert receipt['n_classifiable_rows_used']==75
    assert receipt['n_eligible_W_species']==1
    assert receipt['status'].startswith('W_OBSERVED_ONLY')
    exp=(np.dot(np.array([1,0,0,-1,0,0,0,0,0],float)/np.sqrt(2),axis.TARGET))**2
    assert receipt['mean_species_equal_W']==pytest.approx(exp)
    with pytest.raises(ValueError,match='Git blob mismatch'):
        source.write_text(source.read_text()+'\n# altered')
        axis.run(photos,taxa,source,tmp_path/'bad')


def test_threshold_cannot_be_retuned():
    assert axis.TARGET.sum()==pytest.approx(0,abs=1e-12)
    assert axis.SOURCE_BLOB and axis.NULL_BLOB
