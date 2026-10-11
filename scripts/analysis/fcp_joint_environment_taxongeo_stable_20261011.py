#!/usr/bin/env python3
"""Canonical entry point: make equal-sized genus fold ties version-independent.

The first cross-runtime replay had different equal-size group ordering between
NumPy/sklearn runtimes. This adapter replaces ONLY the core module's injected
fold builder with an explicit deterministic balanced assignment. No selection
of folds by outcome; no change to source data, fit, endpoint or buffer rule.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import fcp_joint_environment_taxongeo_20261011 as core

ALGORITHM='count_descending_genus_lexical_ascending_then_min_load_lowest_fold_v1'

class StableGenusKFold:
    def __init__(self,n_splits=5):
        self.n_splits=int(n_splits)
        if self.n_splits<2:raise ValueError('at least two folds')
    def split(self,X,y=None,groups=None):
        if groups is None:raise ValueError('groups are required')
        labels=np.asarray(groups,dtype=str)
        if len(labels)!=len(X):raise ValueError('group/data length mismatch')
        counts=Counter(labels.tolist())
        if len(counts)<self.n_splits:raise ValueError('fewer genera than folds')
        load=[0]*self.n_splits;mapping={}
        for genus,n in sorted(counts.items(),key=lambda item:(-item[1],item[0])):
            fold=min(range(self.n_splits),key=lambda f:(load[f],f))
            mapping[genus]=fold;load[fold]+=n
        assignment=np.asarray([mapping[g] for g in labels],int)
        for fold in range(self.n_splits):
            yield np.flatnonzero(assignment!=fold),np.flatnonzero(assignment==fold)


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--outdir',required=True,type=Path)
    args,_=parser.parse_known_args()
    # Explicit dependency injection; never modify sklearn's global implementation.
    core.GroupKFold=StableGenusKFold
    core.main()
    result_path=args.outdir/'result.json'
    result=json.loads(result_path.read_text())
    result['genus_fold_algorithm']=ALGORITHM
    result['first_library_default_fold_results_superseded']=True
    for name,case in result['cases'].items():
        rows=pd.read_csv(args.outdir/(name+'_oof_losses.csv.gz'),usecols=['genus','fold','inat_taxon_id'])
        if rows.groupby('genus').fold.nunique().max()!=1:raise RuntimeError('one genus in multiple folds')
        mapping=rows.groupby('genus',sort=True).fold.first().astype(int).to_dict()
        raw=json.dumps(mapping,sort_keys=True,ensure_ascii=True,separators=(',',':')).encode()
        case['genus_fold_map_sha256']=hashlib.sha256(raw).hexdigest()
        rows[['genus','fold']].drop_duplicates().sort_values('genus').to_csv(args.outdir/(name+'_genus_fold_map.csv'),index=False)
    result_path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('STABLE_FOLD_SOURCE_REPLAY_COMPLETE',flush=True)

if __name__=='__main__':main()
