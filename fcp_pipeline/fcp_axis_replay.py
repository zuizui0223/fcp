"""Run *original frozen* nine-palette H2 geometry on all photos per species.

This is a separate expansion. The historical study code must be provided from
its immutable Git blob; disttrait is not bitwise equivalent for near ties. Full
construction-preserving W-null (999 worlds) is optional but REQUIRED for any
new cross-cohort axis-support claim, not implied by a descriptive per-species W.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_BLOB = '6bec7ee121170195f8795238f797ee7167406b48'
NULL_BLOB = '8af9e4a1b63c2ddae3240649f95babebf378a19b'
PALETTE = ('white','yellow','orange','red','pink','magenta','purple','blue','bronze')
TARGET = np.array([1.0]+[-1/8]*8,float)
TARGET /= np.linalg.norm(TARGET)


def git_blob(path:Path)->str:
    raw=path.read_bytes()
    return hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()


def _streaming_sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for data in iter(lambda: handle.read(1<<20),b''):
            digest.update(data)
    return digest.hexdigest()


def source_load(path:Path,expected_git_blob:str):
    if git_blob(path)!=expected_git_blob:
        raise ValueError('Original study-specific H2 source Git blob mismatch')
    spec=importlib.util.spec_from_file_location('_fcp_original_h2',path)
    if spec is None or spec.loader is None:
        raise ValueError('Unavailable original study H2 Python source')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(measured:Path,opportunity:Path,frozen:Path,output:Path,
        *, threshold:float=0.10, null_source:Path|None=None,
        null_count:int=999,seed:int=20260915,cohort:str='expansion')->dict:
    if threshold not in (.10,.20):
        raise ValueError('Only frozen 0.10 primary and 0.20 strict gates permitted')
    h2=source_load(frozen,SOURCE_BLOB)
    if tuple(h2.BIO)!=PALETTE or int(h2.MIN_CLASSIFIABLE)!=40:
        raise ValueError('Original phenotype palette / qualification drift')
    frame=pd.read_csv(opportunity,usecols=['inat_taxon_id','species'])
    if frame.inat_taxon_id.duplicated().any() or frame.species.duplicated().any():
        raise ValueError('Taxonomic identity / species names are not one to one')
    work=h2.load_classifiable(measured)
    # Full raw classifiable photo rows go to the frozen estimator, no per-species cap.
    vectors,_unused,gates=h2.species_delta_vectors(work,threshold)
    if vectors.empty:
        raise ValueError('No measured species passed the frozen H2 qualification')
    x=vectors[[f'delta_{c}' for c in PALETTE]].to_numpy(float)
    norms=np.linalg.norm(x,axis=1)
    if not np.all(np.isfinite(norms)) or np.any(norms<=1e-12):
        raise ValueError('Invalid original Delta vectors')
    w=((x/norms[:,None]) @ TARGET)**2
    join=vectors[['species']].merge(frame,on='species',how='left',validate='one_to_one',indicator=True)
    if not join['_merge'].eq('both').all():
        raise ValueError('Unknown or taxonomically unresolved species in H2 output')
    output.mkdir(parents=True,exist_ok=True)
    out=pd.DataFrame({'inat_taxon_id':join.inat_taxon_id.to_numpy(),
                      'cohort':cohort,'w_alignment_sq':w,
                      'method':'frozen_study_specific_H2',
                      'n_classifiable':vectors.n_classifiable.to_numpy(),
                      'threshold':threshold})
    out.to_csv(output/'species_W_fixed_axis.csv',index=False,float_format='%.17g')
    vectors.to_csv(output/'species_delta_vectors.csv',index=False,float_format='%.17g')
    status='W_OBSERVED_ONLY_NO_NEW_STRUCTURED_NULL_INFERENCE'
    null_summary=None
    if null_source is not None:
        if null_count!=999:
            raise ValueError('New W axis inference requires exactly 999 frozen construction-preserving null worlds')
        other=source_load(null_source,NULL_BLOB)
        prepared=other.prepare_selected(work,vectors.species.astype(str).tolist())
        rng=np.random.default_rng(seed)
        sims=[]
        for _ in range(null_count):
            permuted=other.permute_within_coarse_morph(prepared,rng)
            u=other.delta_units_from_compositions(permuted,prepared['species_indices'])
            sims.append(float(np.mean((u@TARGET)**2)))
        sims=np.asarray(sims)
        sims_obs=float(np.mean(w))
        null_summary={'replicates':len(sims),'null_mean':float(np.mean(sims)),
                      'null_median':float(np.median(sims)),
                      'p_upper':float((1+np.sum(sims>=sims_obs))/(len(sims)+1)),
                      'construction_preserves_coarse_state_counts':True,
                      'fixed_axis_not_fitted_from_new_data':True}
        pd.DataFrame({'W':sims}).to_csv(output/'W_structured_null.csv',index=False)
        status='EXPANSION_W_STRUCTURED_NULL_COMPUTED_RETROSPECTIVE'
    summary={'schema':'fcp_uncapped_H2_white_axis_replay_v1',
             'status':status,'n_classifiable_rows_used':len(work),
             'n_eligible_W_species':len(out),
             'mean_species_equal_W':float(np.mean(w)),
             'frozen_gates':{k:int(v) for k,v in gates.items()},
             'threshold':threshold,'null':null_summary,
             'reused_identical_Hellinger_two_means':True,
             'new_species_selection_not_prospective_confirmation':True,
             'original_paper_H2_decision_unchanged':True,
             'source_git_blob':SOURCE_BLOB,
             'input_sha256':_streaming_sha256(measured),
             'frozen_original_analysis_source_git_blob':SOURCE_BLOB}
    (output/'W_axis_receipt.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--measured',type=Path,required=True)
    p.add_argument('--opportunity',type=Path,required=True)
    p.add_argument('--frozen-h2-source',type=Path,required=True)
    p.add_argument('--frozen-null-source',type=Path)
    p.add_argument('--threshold',type=float,choices=[.1,.2],default=.1)
    p.add_argument('--cohort',default='expansion')
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args(argv)
    print(json.dumps(run(a.measured,a.opportunity,a.frozen_h2_source,a.outdir,
                         threshold=a.threshold,null_source=a.frozen_null_source,
                         cohort=a.cohort)))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
