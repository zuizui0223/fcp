#!/usr/bin/env python3
"""Independent exact-Fraction original-source audit of distance-bin colour allocations.

Does NOT import the primary distance_block.py code.
"""
import argparse,csv,io,zipfile,hashlib,json,itertools
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import numpy as np

ZIP_MD5='008ec196aa5acb769dbfd36e835e32ba'
FREQ_SHA='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
PARENT_SHA='16c91defb8842ad52c63e96db02f2539f0dbbcf1ee33dbed3b2c08cf1d131a5f'
SEEDS={'median':2026101312,'pairs':2026101313}


def source(path,archive):
    src=Path(path).read_bytes(); bs=Path(archive).read_bytes()
    assert hashlib.sha256(src).hexdigest()==FREQ_SHA
    assert len(bs)==127502 and hashlib.md5(bs).hexdigest()==ZIP_MD5
    with zipfile.ZipFile(io.BytesIO(bs)) as z:
        assert z.testzip() is None and z.read('data/frequency.txt')==src
    rows=list(csv.DictReader(src.decode().splitlines(),delimiter='\t'))
    assert len(rows)==512
    for r in rows:
        for c in ['line','distance','v','w','l','s']:r[c]=int(r[c])
        r['n']=r['v']+r['w']
    assert len({(r['year'],r['population'],r['line'],r['quad']) for r in rows})==512
    return rows


def independent_bins(plots,mode):
    """Bin by unique physical-distance order within a line; no dependence on colour values."""
    nonneg=sorted({r['distance'] for r in plots if r['distance']>=0 and r['n']>0})
    out=[]
    for r in plots:
        if r['n']==0:
            out.append('NO_TARGET_PLANTS')
        elif r['distance']<0:
            out.append('unknown')
        elif mode=='median':
            cut=nonneg[(len(nonneg)-1)//2]
            out.append('near' if r['distance']<=cut else 'far')
        else:
            out.append(f"pair{nonneg.index(r['distance'])//2:02d}")
    return out


def exactly(src,mode):
    by_site=defaultdict(list)
    for r in src:
        by_site[r['year']+'|'+r['population']].append(r)
    computed={}
    for site,rows in sorted(by_site.items()):
        groups=defaultdict(list)
        for line in sorted(set(r['line'] for r in rows)):
            selection=[r for r in rows if r['line']==line]
            tags=independent_bins(selection,mode)
            for r,t in zip(selection,tags):
                if r['n']>0: groups[f'{line}:{t}'].append(r)
        N=sum(r['n'] for r in rows);V=sum(r['v'] for r in rows);W=N-V
        rmeta=dict(N=N,V=V,W=W,source_quadrats=len(rows),n_occupied=sum(r['n']>0 for r in rows))
        if N>=2 and V>0 and W>0:
            p=F(V,N);den=N*p*(1-p)
            B_num=sum((F(r['v'])-r['n']*p)**2/r['n'] for r in rows if r['n']>0)
            between=F(0);W_null=F(0)
            informative=[]
            for label,rr in sorted(groups.items()):
                m=len(rr);Nh=sum(r['n'] for r in rr);Vh=sum(r['v'] for r in rr)
                ph=F(Vh,Nh)
                between+=Nh*(ph-p)**2
                if m>1 and Nh>1:
                    W_null+=Nh*ph*(1-ph)*F(m-1,Nh-1)
                if m>=2 and Vh>0 and Vh<Nh:
                    informative.append((label,rr))
            E=(between+W_null)/den
            inf_lines={tag.split(':')[0] for tag,_ in informative}
            inf_n=sum(sum(r['n'] for r in rr) for _,rr in informative)
            rmeta.update(B=B_num/den, B_between=between/den,E=E, delta=B_num/den-E,
                 n_bins=len(groups),n_informative_bins=len(informative),
                 n_informative_lines=len(inf_lines),n_informative_plants=inf_n,
                 eligible_new=len(inf_lines)>=2 and inf_n>=30)
        else:
            rmeta.update(eligible_new=False)
        rmeta['groups']=groups
        computed[site]=rmeta
    return computed


def microcase_exact():
    # Under two disjoint bins of sizes 2+3 and 3+2, verify exact conditional expectation.
    plot_sizes=(2,3,3,2);bins=(0,0,1,1);V_bin=(2,3);N=sum(plot_sizes);V=sum(V_bin);p=F(V,N)
    draws=[]
    for left in itertools.combinations(range(5),2):
        for right in itertools.combinations(range(5),3):
            a=(sum(i<2 for i in left),sum(i>=2 for i in left),sum(i<3 for i in right),sum(i>=3 for i in right))
            B=sum((F(v)-nn*p)**2/nn for v,nn in zip(a,plot_sizes))/(N*p*(1-p))
            draws.append(B)
    expectation=sum(draws,F(0))/len(draws)
    between=F(0);within=F(0)
    for ns,vs in [(5,2),(5,3)]:
        pb=F(vs,ns)
        between+=ns*(pb-p)**2
        within+=ns*pb*(1-pb)*F(1,ns-1)
    assert expectation==(between+within)/(N*p*(1-p))
    return 1


def independent_monte_carlo(calcs,eligible_names,mode,R=9999):
    rng=np.random.default_rng(SEEDS[mode]);records=[];trials=[]
    for name in sorted(eligible_names):
        r=calcs[name];N=r['N'];V=r['V'];p=V/N;den=N*p*(1-p)
        score=np.zeros(R,dtype=float)
        for key,plots in sorted(r['groups'].items()):
            ns=np.array([q['n'] for q in plots],dtype=int);v=sum(q['v'] for q in plots)
            if len(ns)>1:
                draw=rng.multivariate_hypergeometric(ns,v,size=R)
            else:
                draw=np.full((R,1),v,dtype=int)
            assert np.all(draw.sum(axis=1)==v) and np.all((draw>=0)&(draw<=ns))
            score+=((draw-ns[None,:]*p)**2/ns[None,:]).sum(axis=1)
        trials.append(score/den-float(r['E']));records.append(float(r['delta']))
    ref=np.mean(np.array(trials),axis=0);observed=float(np.mean(records))
    return dict(observed=observed,p=(1+int(np.sum(ref>=observed)))/(R+1),
        null_mean=float(ref.mean()),ref=ref,trial_matrix=np.array(trials))


def verify(root):
    root=Path(root)
    src=source(root/'source/frequency.txt',root/'source/Sisyrinchium_data.zip')
    original_bytes=(root/'source/parent_line_result.json').read_bytes()
    assert hashlib.sha256(original_bytes).hexdigest()==PARENT_SHA
    parent=json.loads(original_bytes)
    results=json.loads((root/'results/result.json').read_text())
    assert results['schema']=='fcp_sisyrinchium_water_distance_bin_conditioned_clumping_v1'
    assert results['R']==9999 and results['source']['source_quadrats']==512
    tests=microcase_exact()+2;max_diff=0.0
    parent_names=set(parent['full']['results']['clumping']['names']);assert len(parent_names)==11
    for mode in ('median','pairs'):
        data=exactly(src,mode)
        assertions=[]
        for snap in results['support']:
            rec=data[snap['siteyear']]; metrics=snap[mode]
            parentok=snap['siteyear'] in parent_names
            assert snap['parent_line_clumping_eligible']==parentok;tests+=1
            expected_status='ELIGIBLE' if rec['eligible_new'] and parentok else (
                 'HOLD_PARENT_SITEYEAR_NOT_ELIGIBLE' if not parentok else 'HOLD_DISTANCE_BIN_SUPPORT')
            assert metrics['status']==expected_status,(mode,snap['siteyear'],metrics['status'],expected_status)
            tests+=1
            for k in ['N','V','W','n_occupied','n_bins','n_informative_bins','n_informative_lines','n_informative_plants']:
                if k in rec and k in metrics:
                    assert metrics[k]==rec[k],(mode,snap['siteyear'],k)
                    tests+=1
            for src_key,dest_key in [('B','B_observed'),('B_between','B_between_fixed_distance_bins'),
                                    ('E','E_B_given_bin_inventory'),('delta','delta')]:
                if src_key in rec:
                    diff=abs(float(rec[src_key])-metrics[dest_key]); max_diff=max(max_diff,diff)
                    assert diff<1e-12,(mode,snap['siteyear'],src_key,diff)
                    tests+=1
            if expected_status=='ELIGIBLE':assertions.append(snap['siteyear'])
        main=results[mode]
        assert main['names']==sorted(assertions) and main['n_eligible']==len(assertions);tests+=2
        sim=independent_monte_carlo(data,assertions,mode,R=9999)
        for key,src_value in [('mean_delta','observed'),('p_one_sided','p'),('null_mean_delta','null_mean')]:
            diff=abs(main[key]-sim[src_value]);max_diff=max(max_diff,diff)
            assert diff<1e-12,(mode,key,main[key],sim[src_value]);tests+=1
        for name in assertions:
            diff=abs(float(data[name]['delta'])-main['siteyear_estimates'][name])
            assert diff<1e-12;max_diff=max(max_diff,diff);tests+=1
        # SOURCE-EXPOSED robustness: LOO and latest site set are derived from the exact same replicate matrix.
        ob=np.array([float(data[name]['delta']) for name in assertions]);mat=sim['trial_matrix']
        for i,name in enumerate(assertions):
            record=next((x for x in main['posthoc_leave_one_siteyear_out'] if x['dropped']==name),None)
            assert record is not None
            delta=(ob.sum()-ob[i])/(len(ob)-1)
            pv=(1+int(np.sum((mat.sum(axis=0)-mat[i])/(len(ob)-1)>=delta)))/(9999+1)
            assert abs(delta-record['delta'])<1e-12 and abs(pv-record['p_one_sided'])<1e-12
            tests+=2
        latest={}
        for i,name in enumerate(assertions):
            y,p=name.split('|',1)
            if p not in latest or y>assertions[latest[p]].split('|',1)[0]:latest[p]=i
        ix=sorted(latest.values()); target=main['latest_year_per_site']
        assert target['names']==[assertions[i] for i in ix];tests+=1
        diff=abs(target['mean_delta']-np.mean(ob[ix]));max_diff=max(max_diff,diff);assert diff<1e-12;tests+=1
        pv=(1+int(np.sum(mat[ix].mean(axis=0)>=ob[ix].mean())))/10000
        assert target['p_one_sided']==pv;tests+=1
    out={'status':'PASS','independent_fraction_checks_and_null_checks':tests,'max_abs_diff':max_diff,
         'source_original_zip_md5':ZIP_MD5,'source_original_frequency_sha256':FREQ_SHA,
         'all_512_quadrats_retained':True,'no_import_of_primary_code':True,
         'recomputed_source_conditioned_9999_draws_per_method':True,
         'not_biological_independence_or_causation':True}
    (root/'results/independent_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out));return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='.');a=p.parse_args()
    verify(a.root)
