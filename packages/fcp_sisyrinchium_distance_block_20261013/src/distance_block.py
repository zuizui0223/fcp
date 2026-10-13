#!/usr/bin/env python3
"""Distance-bin conditioned flower colour allocation using unchanged public field counts.

Retrospective/source-exposed; synthetic random-label null NOT causal assignment of real plants.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np

ZIP_MD5='008ec196aa5acb769dbfd36e835e32ba'
INNER_SHA='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
PARENT_SHA='16c91defb8842ad52c63e96db02f2539f0dbbcf1ee33dbed3b2c08cf1d131a5f'
METHODS=('median','pairs')
SEEDS={'median':2026101312,'pairs':2026101313}


def source_rows(freq_path,zip_path):
    freq=Path(freq_path).read_bytes()
    zipped=Path(zip_path).read_bytes()
    if hashlib.sha256(freq).hexdigest()!=INNER_SHA or hashlib.md5(zipped).hexdigest()!=ZIP_MD5 or len(zipped)!=127502:
        raise ValueError('SOURCE_HASH_CHANGED')
    with zipfile.ZipFile(io.BytesIO(zipped)) as z:
        if z.testzip() is not None or z.read('data/frequency.txt')!=freq:
            raise ValueError('ZIP_OR_MEMBER_MISMATCH')
    rows=[]
    keyset=set()
    for r in csv.DictReader(freq.decode('utf-8-sig').splitlines(),delimiter='\t'):
        z=dict(r)
        for key in ['distance','line','v','w','l','s']:
            z[key]=int(z[key])
        if any(z[k]<0 for k in ['v','w','l','s']):raise ValueError('INVALID_NEGATIVE_PLANT_COUNTS')
        key=(z['year'],z['population'],z['line'],z['quad'])
        if key in keyset: raise ValueError('DUPLICATE_QUADRAT')
        keyset.add(key)
        z['n']=z['v']+z['w']
        rows.append(z)
    if len(rows)!=512 or {k:sum(z[k] for z in rows) for k in ['v','w','l','s']}!={'v':5081,'w':8723,'l':4463,'s':29}:
        raise ValueError('SISYRINCHIUM_ORIGINAL_COUNT_INVENTORY_CHANGED')
    return rows


def distance_bins(distances,n,method='median'):
    if method not in METHODS:raise ValueError('UNKNOWN_DISTANCE_BINNING')
    if len(distances)!=len(n):raise ValueError('SOURCE_LENGTH_MISMATCH')
    all_dist=sorted(set(int(d) for d,k in zip(distances,n) if k>0 and d>=0))
    if all_dist:
        if method=='median':
            median=all_dist[(len(all_dist)-1)//2]
            return ['NO_TARGET_PLANTS' if k<=0 else ('unknown' if d<0 else ('near' if d<=median else 'far')) for d,k in zip(distances,n)]
        order={d:i for i,d in enumerate(all_dist)}
        return ['NO_TARGET_PLANTS' if k<=0 else ('unknown' if d<0 else f"pair{order[d]//2:02d}") for d,k in zip(distances,n)]
    return ['NO_TARGET_PLANTS' if k<=0 else 'unknown' for d,k in zip(distances,n)]


def phenotype_partition(n,v,bins):
    if len(n)!=len(v) or len(n)!=len(bins):raise ValueError('SOURCE_LENGTH_MISMATCH')
    arrays=[(int(nn),int(vv),str(b)) for nn,vv,b in zip(n,v,bins)]
    if any(nn<0 or vv<0 or vv>nn for nn,vv,b in arrays):raise ValueError('SOURCE_COUNTS_INVALID')
    arrays=[a for a in arrays if a[0]>0]
    N=sum(a[0] for a in arrays);V=sum(a[1] for a in arrays)
    if N<2 or V==0 or V==N:raise ValueError('MONOMORPHIC_OR_INSUFFICIENT')
    p=V/N;den=N*p*(1-p)
    B=sum((vv-nn*p)**2/nn for nn,vv,b in arrays)/den
    groups=defaultdict(list)
    for nn,vv,b in arrays:groups[b].append((nn,vv))
    between=0.;within=0.
    for label,items in sorted(groups.items()):
        Nb=sum(x[0] for x in items);Vb=sum(x[1] for x in items);pb=Vb/Nb;m=len(items)
        between+=Nb*(pb-p)**2/den
        if Nb>1 and m>1:within+=Nb*pb*(1-pb)*(m-1)/(Nb-1)/den
    E=between+within
    return float(B),float(between),float(E),float(B-E)


def eligibility(groups):
    eligible=[]
    for label,data in groups.items():
        n=np.asarray(data['n'],dtype=int);v=np.asarray(data['v'],dtype=int)
        if len(n)>=2 and int(v.sum())>=1 and int((n-v).sum())>=1 and len(np.flatnonzero(n>0))>=2:
            eligible.append((label,data))
    distinct_lines={str(g['line']) for name,g in eligible}
    total_n=sum(int(sum(g['n'])) for name,g in eligible)
    status='ELIGIBLE' if len(distinct_lines)>=2 and total_n>=30 else 'HOLD_DISTANCE_BIN_SUPPORT'
    return {'status':status,'n_informative_bins':len(eligible),'n_informative_lines':len(distinct_lines),'n_informative_plants':total_n}


def hypergeom_trial(sizes,V,rng,R):
    n=np.asarray(sizes,dtype=int)
    if len(n)==0 or np.any(n<0) or V<0 or V>sum(n):raise ValueError('BAD_DRAW_INPUT')
    if len(n)==1:
        return np.full((R,1),V,dtype=int)
    return rng.multivariate_hypergeometric(n,int(V),size=R)


def _bins_for_site(rows,method):
    by_line=defaultdict(list)
    for r in rows:
        by_line[str(r['line'])].append(r)
    pairs=[]
    for line,rr in sorted(by_line.items()):
        bs=distance_bins([x['distance'] for x in rr],[x['n'] for x in rr],method=method)
        for x,b in zip(rr,bs):
            pairs.append((x,b, line+':'+b if b!='NO_TARGET_PLANTS' else 'NO_TARGET_PLANTS'))
    return pairs


def _site_stat(rows,method):
    occupied=[(x,b,name) for x,b,name in _bins_for_site(rows,method) if x['n']>0]
    grouping=defaultdict(list)
    for x,b,label in occupied:
        grouping[label].append(x)
    support=eligibility({label:dict(n=[x['n'] for x in items],v=[x['v'] for x in items],line=label.split(':')[0]) for label,items in grouping.items()})
    N=sum(x['n'] for x,b,name in occupied);V=sum(x['v'] for x,b,name in occupied)
    src_n=[x['n'] for x,b,name in occupied]
    src_v=[x['v'] for x,b,name in occupied]
    labels=[name for x,b,name in occupied]
    metric=None
    if V and N>V:
        B,bet,E,delta=phenotype_partition(src_n,src_v,labels)
        metric={'B_observed':B,'B_between_fixed_distance_bins':bet,'E_B_given_bin_inventory':E,
                'delta':delta,'n_occupied':len(occupied),'N':N,'V':V,'W':N-V,
                'n_bins':len(grouping), **support}
    else:
        metric={'N':N,'V':V,'W':N-V,'n_occupied':len(occupied),'status':'HOLD_NO_BOTH_M_COLOURS',
                'n_informative_bins':0,'n_informative_lines':0,'n_informative_plants':0}
    return metric,grouping


def _draws_for_site(metric,groups,rng,R):
    N=metric['N'];V=metric['V'];p=V/N
    den=N*p*(1-p)
    scores=np.zeros(R,dtype=float)
    for label,items in sorted(groups.items()):
        n=np.array([x['n'] for x in items],dtype=int)
        vv=sum(x['v'] for x in items)
        draws=hypergeom_trial(n,vv,rng,R)
        # Each bin retains observed M-morph total, quadrat opportunity and membership.
        assert np.all(draws.sum(axis=1)==vv) and np.all(draws>=0) and np.all(draws<=n)
        scores+=((draws-n[None,:]*p)**2/n[None,:]).sum(axis=1)
    return scores/den-metric['E_B_given_bin_inventory']


def _ensemble(records,method,R):
    rng=np.random.default_rng(SEEDS[method])
    eligible=[z for z in records if z[method]['status']=='ELIGIBLE']
    if len(eligible)<3:
        return {'status':'HOLD_INSUFFICIENT_SITEYEARS','n_eligible':len(eligible)}
    sim=[];obs=[]
    for z in eligible:
        x=z[method];
        obs.append(x['delta'])
        sim.append(_draws_for_site(x,z['_groups_'+method],rng,R))
    names=[z['siteyear'] for z in eligible]
    obs=np.array(obs,dtype=float);sim=np.array(sim,dtype=float)
    avg_sims=sim.mean(axis=0)
    mean_delta=float(obs.mean())
    pval=float((1+np.sum(avg_sims>=mean_delta))/(R+1))
    out={'status':'ELIGIBLE','n_eligible':len(eligible),'names':names,'n_named_populations':len({z['population'] for z in eligible}),
         'n_positive_siteyears':int((obs>0).sum()), 'mean_observed_B':float(np.mean([z[method]['B_observed'] for z in eligible])),
         'mean_B_fixed_between_bins':float(np.mean([z[method]['B_between_fixed_distance_bins'] for z in eligible])),
         'mean_exact_null_B':float(np.mean([z[method]['E_B_given_bin_inventory'] for z in eligible])),
         'mean_delta':mean_delta,'mean_delta_pp':100*mean_delta,'null_mean_delta':float(avg_sims.mean()),
         'null_CI95':list(map(float,np.quantile(avg_sims,[.025,.975]))),'p_one_sided':pval,
         'source_exposed_retrospective':True,'random_label_not_causal':True,
         'simulated_siteyear_delta_permutation_count':R,
         'siteyear_estimates':{k:float(v) for k,v in zip(names,obs)}}
    latest={}
    for i,z in enumerate(eligible):
        if z['population'] not in latest or z['year']>eligible[latest[z['population']]]['year']:
            latest[z['population']]=i
    chosen=sorted(latest.values());
    if len(chosen)>=3:
        latest_ref=sim[chosen,:].mean(axis=0);val=float(obs[chosen].mean())
        out['latest_year_per_site']={'status':'ELIGIBLE','n_named_populations':len(chosen),
            'names':[names[i] for i in chosen], 'mean_delta':val,'p_one_sided':float((1+np.sum(latest_ref>=val))/(R+1))}
    else:
        out['latest_year_per_site']={'status':'HOLD_LT3_POPULATIONS'}
    out['posthoc_leave_one_siteyear_out']=[]
    if len(obs)>=4:
        for i,name in enumerate(names):
            ref=(sim.sum(axis=0)-sim[i])/(len(obs)-1)
            v=(obs.sum()-obs[i])/(len(obs)-1)
            out['posthoc_leave_one_siteyear_out'].append({'dropped':name,'delta':float(v),
                   'p_one_sided':float((1+np.sum(ref>=v))/(R+1))})
    return out


def analyze(source,source_zip,parent,R=9999):
    rows=source_rows(source,source_zip)
    parent_bytes=Path(parent).read_bytes()
    if hashlib.sha256(parent_bytes).hexdigest()!=PARENT_SHA:raise ValueError('PARENT_RESULT_HASH_CHANGED')
    p=json.loads(parent_bytes)
    original_names=p['full']['results']['clumping']['names']
    if len(original_names)!=11:raise ValueError('PARENT_FROZEN_COHORT_CHANGED')
    by_site=defaultdict(list)
    for row in rows:
        by_site[row['year']+'|'+row['population']].append(row)
    support=[]; plot_ledger=[]
    for name, rr in sorted(by_site.items()):
        yr,pop=name.split('|',1)
        z={'siteyear':name,'year':yr,'population':pop,'parent_line_clumping_eligible':name in original_names,
           'n_source_quadrats':len(rr),'n_M_present_quadrats':sum(x['n']>0 for x in rr),
           'n_negative_distance_quadrats':sum(x['distance']<0 for x in rr)}
        for method in METHODS:
            metrics,groups=_site_stat(rr,method)
            if name not in original_names:
                metrics['status']='HOLD_PARENT_SITEYEAR_NOT_ELIGIBLE'
            z[method]=metrics
            z['_groups_'+method]=groups
        # ledger all 512 raw quadrats and their determined source-blind bins
        mapped={method:{(x['line'],x['quad']):b for x,b,label in _bins_for_site(rr,method)} for method in METHODS}
        for x in rr:
            plot_ledger.append({k:x[k] for k in ['year','population','line','quad','distance','v','w','l','s','n']}|
                               {'median_bin':mapped['median'][(x['line'],x['quad'])],
                                'pair_bin':mapped['pairs'][(x['line'],x['quad'])],
                                'no_M_plants':x['n']==0})
        support.append(z)
    if len(support)!=16 or len(plot_ledger)!=512:raise ValueError('SOURCE_GROUPING_INCOMPLETE')
    results={m:_ensemble(support,m,R) for m in METHODS}
    clean=[]
    for s in support:
        clean.append({k:v for k,v in s.items() if not k.startswith('_')})
    return {'schema':'fcp_sisyrinchium_water_distance_bin_conditioned_clumping_v1',
            'protocol_commit':'329099093265e46bca0c82b51804cd061a0a29a9',
            'source':{'doi':'10.5061/dryad.j8q34','zenodo_record':5003124,'md5_zip':ZIP_MD5,
                      'frequency_sha256':INNER_SHA,'parent_result_sha256':PARENT_SHA,
                      'source_quadrats':len(rows),'siteyears':len(support),'named_populations':len({z['population'] for z in support}),
                      'total_M_purple':5081,'total_M_white':8723,'other_L':4463,'other_S':29,
                      'negative_distance_quadrats':sum(x['distance']<0 for x in rows),
                      'zero_M_quadrats':sum(x['n']==0 for x in rows)},
            'sampling_design':'source-exposed single study; quadrats nested in lines/siteyears',
            'hypotheses':{'median':'coarse near-versus-far distance conditioned','pairs':'adjacent unique-distance pair conditioned'},
            'seeds':SEEDS,'R':R,'median':results['median'],'pairs':results['pairs'],'support':clean,'plot_ledger':plot_ledger,
            'caveats':['Conditional p benchmarks random morph labels among plots, not actual dispersal, field randomization or habitat selection.',
                'A distance-binned null is not a fitted environmental model or proof that water distance does not matter.',
                'Original study reported congener effects and water-distance collinearity, not newly discovered here.',
                'Within species conditional source results do not establish independent evolution or genetic F_ST.']}


def write_results(data,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    # Keep the complete 512-row auditable ledger in one canonical CSV, without
    # duplicating all original row-level records in the primary numeric JSON.
    primary={key:value for key,value in data.items() if key!='plot_ledger'}
    primary['plot_ledger_file']='results/quadrat_distance_bin_ledger.csv'
    primary['plot_ledger_record_count']=len(data['plot_ledger'])
    (out/'result.json').write_text(json.dumps(primary,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    flat=[]
    for z in data['support']:
        for method in METHODS:
            x=z[method]
            flat.append({k:z[k] for k in ['siteyear','year','population','n_source_quadrats','n_M_present_quadrats','n_negative_distance_quadrats']}|
                {'method':method, **x})
    cols=['siteyear','year','population','method','status','N','V','W','n_source_quadrats','n_M_present_quadrats',
        'n_negative_distance_quadrats','n_bins','n_informative_bins','n_informative_lines','n_informative_plants',
        'B_observed','B_between_fixed_distance_bins','E_B_given_bin_inventory','delta']
    with (out/'siteyear_distance_bin_ledger.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(flat)
    keys=['year','population','line','quad','distance','v','w','l','s','n','median_bin','pair_bin','no_M_plants']
    with (out/'quadrat_distance_bin_ledger.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(data['plot_ledger'])


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',default='source/frequency.txt');p.add_argument('--source-zip',default='source/Sisyrinchium_data.zip')
    p.add_argument('--parent',default='source/parent_line_result.json');p.add_argument('--output',default='results')
    p.add_argument('--draws',type=int,default=9999)
    opt=p.parse_args()
    result=analyze(Path(opt.source),Path(opt.source_zip),Path(opt.parent),R=opt.draws)
    write_results(result,opt.output)
    print(json.dumps({m:{k:result[m].get(k) for k in ['status','n_eligible','mean_delta','p_one_sided']} for m in METHODS},ensure_ascii=False))
