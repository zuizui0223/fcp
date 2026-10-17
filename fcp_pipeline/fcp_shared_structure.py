"""Cross-species structural commonality of *within-species* flower-colour variation.

Only per-species photographed FCP properties are compared between species.
The species' mean white/yellow/red/blue frequency is *not* an interspecific
colour-on-climate response, which belongs to the separate mois repository.

The original frozen H2 fixed-axis construction null and the spatial vertex
permutation null remain the inferential authorities. This small module only
combines *reproduced and provenance-labelled* species metrics and checks which
cross-species FCP hypotheses are testable; it never invents new null p-values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

MIN_SPATIAL_PAIRS = 30
ENVIRONMENT_PHOTO_LOCATION = "direct_source_photo_coordinate"
W_SOURCE = "frozen_study_specific_H2"


def as_int(x: pd.Series) -> pd.Series:
    z = pd.to_numeric(x, errors="raise")
    if z.isna().any() or (z <= 0).any() or (z % 1 != 0).any():
        raise ValueError("Non-canonical taxon identity")
    return z.astype('int64')


def baseline(path: Path) -> pd.DataFrame:
    ledger = pd.read_csv(path, low_memory=False)
    required = {"inat_taxon_id", "species", "n_classifiable", "n_classifiable_observations",
                "n_classifiable_observers", "D_raw", "D_observer_equal", "measured_status"}
    if not required <= set(ledger):
        raise ValueError(f"Missing within-species evidence columns: {sorted(required-set(ledger))}")
    ledger['inat_taxon_id'] = as_int(ledger.inat_taxon_id)
    if ledger.inat_taxon_id.duplicated().any() or ledger.species.duplicated().any():
        raise ValueError("Source ledger duplicated species")
    if (ledger.n_classifiable < 0).any() or (ledger.n_classifiable_observations > ledger.n_classifiable).any():
        raise ValueError("Inconsistent photographed FCP measurement support")
    n = ledger.n_classifiable
    data = ledger[n.ge(40)]
    if not data.D_raw.between(0,1).all() or not data.D_observer_equal.between(0,1).all():
        raise ValueError("Invalid species FCP diversity")
    return ledger


def join_metric(source: pd.DataFrame, path: Path, what: str) -> pd.DataFrame:
    t = pd.read_csv(path, low_memory=False)
    t['inat_taxon_id'] = as_int(t.inat_taxon_id)
    if t.inat_taxon_id.duplicated().any():
        raise ValueError(f"Duplicate {what} by species")
    if not t.inat_taxon_id.isin(source.inat_taxon_id).all():
        raise ValueError(f"{what} has species not in photographic FCP source frame")
    if what == 'W':
        needed={"w_alignment_sq", "method"}
        if not needed <= set(t) or not t.method.eq(W_SOURCE).all():
            raise ValueError("W must be exported from the identical frozen study H2")
        if not t.w_alignment_sq.between(0,1).all():
            raise ValueError("Invalid W squared alignment")
        t = t[['inat_taxon_id','w_alignment_sq']]
    elif what == 'spatial':
        needed={"depletion_50km", "local_pairs_50km", "D_pair"}
        if not needed <= set(t):
            raise ValueError("Spatial metrics missing source-pair definition")
        if ((t.local_pairs_50km < MIN_SPATIAL_PAIRS) & t.depletion_50km.notna()).any():
            raise ValueError("Spatial allocation below frozen local-pair gate")
        t = t[['inat_taxon_id','depletion_50km','local_pairs_50km','D_pair']]
    else:
        raise ValueError(what)
    result = source.merge(t, on='inat_taxon_id', how='left', validate='one_to_one')
    result.loc[~result.n_classifiable.ge(40), t.columns.difference(['inat_taxon_id'])] = np.nan
    return result


def with_environment(source: pd.DataFrame, path: Path) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Use only species/colour environment contrasts from each photo's *own* coordinate.

    Nearest-neighbour mois environment transfers may be analyzed separately, but
    cannot silently become this direct-photo environmental gate.
    """
    raw = pd.read_csv(path, low_memory=False)
    needed={"inat_taxon_id", "variable", "white_minus_nonwhite_sd", "n_white", "n_nonwhite", "location_source"}
    if not needed <= set(raw):
        raise ValueError("Missing environment evidence columns")
    raw.inat_taxon_id = as_int(raw.inat_taxon_id)
    if raw.duplicated(['inat_taxon_id','variable']).any():
        raise ValueError("Duplicate environment effect per species and variable")
    if raw.location_source.ne(ENVIRONMENT_PHOTO_LOCATION).any():
        raise ValueError("Proxy or nearest-photo soil cannot be the primary source-photo environmental result")
    if (raw[['n_white','n_nonwhite']] < 5).any(axis=None):
        raise ValueError("Insufficient white/nonwhite environmental support")
    if not np.isfinite(raw.white_minus_nonwhite_sd).all():
        raise ValueError("Missing effects must remain missing, never zero")
    if not raw.inat_taxon_id.isin(source.inat_taxon_id).all():
        raise ValueError("Environment species absent from opportunity frame")
    # Keep long-form rows so 36 correlated predictors are not mistaken for
    # independent confirmatory tests.
    return source,raw


def with_environment_heterogeneity(source: pd.DataFrame, path: Path) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Read photo-coordinate matched, site-equal *within-species* spread.

    This is the primary ecological predictor for cross-species W, unlike the
    colour-state-specific signed environmental contrast (secondary).
    """
    raw=pd.read_csv(path,low_memory=False)
    needed={'inat_taxon_id','variable','environment_spread_site_equal',
            'n_environment_photos','n_distinct_spatial_cells','location_source'}
    if needed-set(raw):
        raise ValueError('Environment spread requires source-photo site cells and provenance')
    raw['inat_taxon_id']=as_int(raw['inat_taxon_id'])
    if raw.duplicated(['inat_taxon_id','variable']).any():
        raise ValueError('Duplicate species-variable environmental spread')
    if raw.location_source.ne(ENVIRONMENT_PHOTO_LOCATION).any():
        raise ValueError('Proxy or nearest-photo soil invalid for within-species FCP')
    if not raw.inat_taxon_id.isin(source.inat_taxon_id).all():
        raise ValueError('Environmental spread species missing from FCP frame')
    values=pd.to_numeric(raw['environment_spread_site_equal'],errors='coerce')
    if ((values.notna()) & ((values<0) | ~np.isfinite(values))).any():
        raise ValueError('Nonfinite or negative site environmental spread')
    bad=values.notna() & raw.n_distinct_spatial_cells.lt(3)
    if bad.any():
        raise ValueError('Cannot infer environmental heterogeneity from fewer than three cells')
    return source,raw


def _stats(x: pd.Series) -> dict:
    z=pd.to_numeric(x,errors='coerce').dropna()
    if not len(z):
        return {'n_species':0,'status':'not_estimable'}
    return {'n_species':int(len(z)),'median':float(z.median()),
            'mean_species_equal':float(z.mean()),'positive_species_fraction':float((z > 0).mean())}


def _rho(x: pd.DataFrame, a: str, b: str) -> dict:
    q=x[[a,b]].replace([np.inf,-np.inf],np.nan).dropna()
    if len(q)<20 or q[a].nunique()<2 or q[b].nunique()<2:
        return {'n_species':len(q),'status':'not_estimable'}
    return {'n_species':len(q),'spearman_rho':float(spearmanr(q[a],q[b]).statistic),
            'status':'posthoc_descriptive_not_null_test'}


def inspect_commonality(source: pd.DataFrame, env: pd.DataFrame | None=None,
                        environmental_spread: pd.DataFrame | None=None) -> tuple[dict,pd.DataFrame]:
    qualified=source.loc[source.n_classifiable.ge(40)].copy()
    # Only the species metrics, never species' average flower colours, are the
    # variables in the cross-species comparison.
    roles = {
        'within_species_W_direction':_stats(qualified.w_alignment_sq) if 'w_alignment_sq' in qualified else {'status':'not_provided'},
        'overall_sampled_FCP_D':_stats(qualified.D_raw),
        'observer_equal_FCP_D':_stats(qualified.D_observer_equal),
        'between_locality_50km_FCP_partition':_stats(qualified.depletion_50km) if 'depletion_50km' in qualified else {'status':'not_provided'},
    }
    rho={}
    for a,b in [('D_raw','D_observer_equal'),('D_raw','depletion_50km'),
                ('D_raw','w_alignment_sq'),('w_alignment_sq','depletion_50km')]:
        rho[f'{a}_vs_{b}']=_rho(qualified,a,b) if {a,b} <= set(qualified) else {'status':'not_provided'}
    effects=[]
    if env is not None:
        for var, sub in env.groupby('variable', sort=True):
            cols=sub[['inat_taxon_id','white_minus_nonwhite_sd']].copy()
            hit=qualified.merge(cols,on='inat_taxon_id',how='inner',validate='one_to_one')
            for metric in ('D_raw','D_observer_equal','w_alignment_sq','depletion_50km'):
                if metric not in hit: continue
                absdelta=hit.assign(abs_environment_difference=hit.white_minus_nonwhite_sd.abs())
                effects.append({'predictor':str(var),'fcp_metric':metric,
                                **_rho(absdelta,metric,'abs_environment_difference')})
    W_env=[]
    if environmental_spread is not None:
        for variable,sub in environmental_spread.groupby('variable',sort=True):
            q=qualified.merge(sub[['inat_taxon_id','environment_spread_site_equal']],
                              on='inat_taxon_id',how='inner',validate='one_to_one')
            if 'w_alignment_sq' in q:
                W_env.append({'variable':str(variable),'outcome':'W_white_chromatic_axis',
                              **_rho(q,'w_alignment_sq','environment_spread_site_equal')})
    result={
        'schema' :'fcp_uncapped_crossspecies_shared_structure_v1',
        'status':'EXPLORATORY_SPECIES_EQUAL_STRUCTURAL_COMPARISONS',
        'n_source_species':len(source),'n_D_evaluable_species':len(qualified),
        'species_equal_descriptive_summaries':roles,
        'cross_species_structural_relations':rho,
        'ecological_environment_associations':effects,
        'W_vs_withinspecies_environmental_heterogeneity':W_env,
        'W_environment_frozen_construction_null_completed':False,
        'frozen_primary_H2_W_null_recomputed':False,
        'frozen_primary_50km_null_recomputed':False,
        'environment_36_testing_not_prospective':True,
        'environment_confirmatory_support':False,
        'no_colour_prevalence_or_local_genetic_FCP_claim':True,
        'primary_next_tests':[
          'construction-preserving 999-world fixed white-axis H2 null over outcome-blind species',
          'species-exact 50-km local allocation versus composition-preserving null',
          'photo-specific direct environment join and observer/background/geography falsifications',
          'phylogenetic and sample-opportunity sensitivity by family/region',
        ],
    }
    return result,pd.DataFrame(effects)


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ledger',type=Path,required=True)
    p.add_argument('--w',type=Path)
    p.add_argument('--spatial',type=Path)
    p.add_argument('--environment',type=Path)
    p.add_argument('--environment-profile',type=Path)
    p.add_argument('--outdir',type=Path,required=True)
    args=p.parse_args(argv)
    src=baseline(args.ledger)
    if args.w:src=join_metric(src,args.w,'W')
    if args.spatial:src=join_metric(src,args.spatial,'spatial')
    env=None
    if args.environment:src,env=with_environment(src,args.environment)
    spread=None
    if args.environment_profile:src,spread=with_environment_heterogeneity(src,args.environment_profile)
    receipt,effects=inspect_commonality(src,env,spread)
    args.outdir.mkdir(parents=True,exist_ok=True)
    src.to_csv(args.outdir/'fcp_species_commonality_metrics.csv',index=False)
    effects.to_csv(args.outdir/'FCP_vs_environment_exploratory.csv',index=False)
    receipt['sources_sha256']={n:hashlib.sha256(path.read_bytes()).hexdigest() for n,path in
                               [('photo_ledger',args.ledger),('W',args.w),('spatial',args.spatial),('environment',args.environment),('environment_profile',args.environment_profile)] if path}
    (args.outdir/'fcp_commonality_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'n_species':receipt['n_source_species'],'n_D_evaluable':receipt['n_D_evaluable_species'],
                      'environment_36_candidate_associations':len(receipt['ecological_environment_associations'])}))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
