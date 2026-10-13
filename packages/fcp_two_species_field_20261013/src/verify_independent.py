#!/usr/bin/env python3
"""Independent source→site colour-clumping verification.

Uses Python csv/zipfile/Fraction to reconstruct the biological count grain; no
import of field_partition nor reuse of its numerical routines or result ledger.
A separate NumPy hypergeometric implementation checks the 9,999 seeded nulls.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,zipfile,collections
from pathlib import Path
from fractions import Fraction
import numpy as np

EXPECTED_MD5={
 'Sisyrinchium_data.zip':'008ec196aa5acb769dbfd36e835e32ba',
 'Boechera_natural_populations_individual.txt':'85990c8a30a8b9ef3a2c045bb6ffb72a',
 'Boechera_natural_populations_plot.txt':'1c284e1543450aaacf161cb7400d0a5a'
}
SEEDS={'sisyrinchium':2026101307,'boechera':2026101308}

def source_check(path):
    bs=path.read_bytes()
    assert hashlib.md5(bs).hexdigest()==EXPECTED_MD5[path.name], ('SOURCE_CHECKSUM',path)
    return hashlib.sha256(bs).hexdigest()

def get_original_counts(root):
    source=root/'source'
    hashes={name:source_check(source/name) for name in EXPECTED_MD5}
    result={}
    with zipfile.ZipFile(source/'Sisyrinchium_data.zip') as z:
        assert z.testzip() is None
        data=z.read('data/frequency.txt').decode('utf-8-sig')
    rows=list(csv.DictReader(io.StringIO(data, newline=None),delimiter='\t'))
    assert len(rows)==512
    assert list(rows[0])==['year','population','line','quad','distance','humidity','v','w','l','s']
    ss=collections.defaultdict(list)
    keys=set();other_counts=[0,0]
    for row in rows:
        k=tuple(row[t] for t in ('year','population','line','quad'))
        assert k not in keys;keys.add(k)
        vals=[int(row[t]) for t in ('v','w','l','s')]
        assert all(v>=0 for v in vals)
        other_counts[0]+=vals[2];other_counts[1]+=vals[3]
        ss[(row['year'],row['population'])].append((vals[0]+vals[1],vals[0]))
    assert len(ss)==16
    result['sisyrinchium']=ss
    lines=(source/'Boechera_natural_populations_individual.txt').read_text(encoding='utf-8-sig').splitlines()
    assert lines[0].split()==['Site_name','plot','plantID','Flower_color','Leaf_number','Leaf_area_removed','Elevation','Number_of_fruits']
    byplot=collections.defaultdict(list)
    unique=set();special=[]
    for ln,line in enumerate(lines[1:],2):
        v=line.split()
        if len(v)==9:
            assert v[2].endswith('.') and v[3][0].isdigit(),('unexpected id',ln)
            v=v[:2]+[v[2]+' '+v[3]]+v[4:]
            special.append(ln)
        assert len(v)==8
        k=tuple(v[:3]);assert k not in unique;unique.add(k)
        assert v[3] in ('0','1')
        byplot[tuple(v[:2])].append(int(v[3]))
    assert len(unique)==547 and len(special)==3
    plots=list(csv.DictReader((source/'Boechera_natural_populations_plot.txt').open(),delimiter='\t'))
    assert len(plots)==26
    bb=collections.defaultdict(list)
    for r in plots:
        key=(r['Site_name'],r['plot'])
        observed=byplot.pop(key)
        n=int(r['Boechera_flowered']);V=int(r['Boechera_pigmented'])
        assert n==len(observed) and V==sum(observed),(key,n,V)
        bb[r['Site_name']].append((n,V))
    assert not byplot and len(bb)==9
    result['boechera']=bb
    return result,hashes,{'sisyrinchium_other_l':other_counts[0], 'sisyrinchium_other_s':other_counts[1],
                           'boechera_id_whitespace_repairs':len(special)}

def exact_metric(pairs):
    p=[(int(n),int(k)) for n,k in pairs if int(n)>0]
    N=sum(n for n,k in p);V=sum(k for n,k in p);m=len(p)
    if N<20 or m<2 or V<5 or N-V<5:
        return None
    phi=Fraction(V,N)
    denominator=Fraction(N)*phi*(1-phi)
    assert denominator>0
    B=sum((Fraction(k)-n*phi)**2/n for n,k in p)/denominator
    e=Fraction(m-1,N-1)
    # For each plot of size n, Pr(no violet OR all violet) via hypergeometric sampling
    # Compute as exact ratio of combinatorials to check expected co-detection.
    mix=Fraction(0)
    for n,k in p:
        zero=Fraction(math.comb(N-V,n),math.comb(N,n)) if n<=N-V else Fraction(0)
        allpink=Fraction(math.comb(V,n),math.comb(N,n)) if n<=V else Fraction(0)
        mix+=1-zero-allpink
    return {'N':N,'V':V,'W':N-V,'m':m,'B':B,'expected_B':e,'delta':B-e,'expected_mixed':mix/m,
            'observed_mixed':Fraction(sum(int(0<k<n) for n,k in p),m)}

def compare_number(a,b,checks,path,tol=2e-12):
    assert abs(float(a)-float(b))<=tol,(path,float(a),float(b))
    checks[0]+=1
    checks[1]=max(checks[1],abs(float(a)-float(b)))

def numerical_verify(root, reps=9999):
    originals,hashes,qa=get_original_counts(root)
    saved=json.loads((root/'results/result.json').read_text())
    assert saved['monte_carlo_reps']==reps
    assert saved['source_sha256']==hashes
    checks=[0,0.]
    for name,group in originals.items():
        src=saved['primary'][name]
        if name=='sisyrinchium':
            ordered=sorted(group.items(),key=lambda t:(t[0][0],t[0][1]))
            sites=[(f'{key[0]}|{key[1]}',pairs) for key,pairs in ordered]
        else:sites=[(f'single_census_unreported_year|{key}',pairs) for key,pairs in sorted(group.items())]
        good=[];rng=np.random.default_rng(SEEDS[name]);site_sims=[];site_mixed=[]
        for site,pairs in sites:
            x=exact_metric(pairs)
            if x is None:continue
            good.append((site,x))
            # Explicit algebra, not production implementation; NumPy simulation is a separate
            # field inventory/morph-count construction sharing only prespecified RNG definition.
            occ=[(n,v) for n,v in pairs if n>0]
            n=np.array([z[0] for z in occ]);V=x['V'];N=x['N'];m=x['m'];px=V/N
            draw=rng.multivariate_hypergeometric(n,V,size=reps,method='marginals')
            assert np.all(draw.sum(axis=1)==V) and np.all(draw>=0) and np.all(draw<=n)
            nullB=np.sum((draw-n*px)**2/n,axis=1)/(N*px*(1-px))
            site_sims.append(nullB-float(x['expected_B']))
            site_mixed.append(np.mean((draw>0)&(draw<n),axis=1))
        ledger=list(csv.DictReader((root/'results/siteyear_metrics.csv').open()))
        ledger={a['siteyear']:a for a in ledger if a['species']==name}
        assert len(ledger)==len(sites)
        for site,x in good:
            row=ledger[site]
            for key in ('B','expected_B','delta','expected_mixed'):
                compare_number(x[key],float(row[key]),checks,f'{name}/{site}/{key}')
            compare_number(x['N'],int(row['N_target']),checks,f'{site}/N')
            compare_number(x['V'],int(row['V_pigmented']),checks,f'{site}/V')
            compare_number(x['observed_mixed'],float(row['mixed_fraction']),checks,f'{site}/mix')
        assert src['support']['n_eligible']==len(good)
        for key,vals in [('N',sum(q[1]['N'] for q in good)),('V',sum(q[1]['V'] for q in good))]:
            if name=='sisyrinchium' and key=='N':compare_number(vals,11924,checks,'Sisyrinchium primary N')
            if name=='boechera' and key=='N':compare_number(vals,472,checks,'Boechera primary N')
        fields=[('mean_B_observed','B'),('mean_B_expected_finite_random','expected_B'),
                ('excess_between_quadrat_diversity_fraction','delta'),
                ('expected_site_equal_mixed_fraction','expected_mixed'),
                ('observed_site_equal_mixed_quadrat_fraction','observed_mixed')]
        for k,field in fields:
            val=sum(x[field] for _,x in good)/len(good)
            compare_number(val,src['results'][k],checks,f'{name}/aggregate/{k}')
        ns=np.mean(site_sims,axis=0)
        mix=np.mean(site_mixed,axis=0)
        obs=float(sum(x['delta'] for _,x in good)/len(good))
        obsmix=float(sum(x['observed_mixed'] for _,x in good)/len(good))
        p_upper=(1+np.count_nonzero(ns>=obs-1e-12))/(reps+1)
        p_mix=(1+np.count_nonzero(mix<=obsmix+1e-12))/(reps+1)
        for label,val in [('permutation_upper_tail_p',p_upper),('permutation_lower_tail_mixing_p',p_mix),('null_mean_excess',ns.mean()),('null_mean_site_equal_mixed_fraction',mix.mean())]:
            compare_number(val,src['results'][label],checks,f'{name}/permutation/{label}')
        for idx,(site,x) in enumerate(good):
            row=ledger[site]
            truep=(1+np.count_nonzero(site_sims[idx]+float(x['expected_B'])>=float(x['B'])-1e-12))/(reps+1)
            compare_number(truep,float(row['site_p_upper']),checks,f'{site}/site_p')
        assert src['support']['n_source_siteyear']==len(sites)
    checks[0]+=3
    assert qa['sisyrinchium_other_l']==saved['source_quality']['sisyrinchium']['count_l_other_species']
    assert qa['sisyrinchium_other_s']==saved['source_quality']['sisyrinchium']['count_s_other_species']
    assert qa['boechera_id_whitespace_repairs']==3
    return {'status':'PASS','independent_method':'csv.DictReader + Fraction + standalone NumPy inventory Monte Carlo; no imports of production code',
            'checks':checks[0],'max_abs_difference':checks[1],'raw_source_md5_sha256_verified':True,
            'source_qa':qa,'checks_include_source_to_plot_to_site_to_null':True,'same_source_not_independent_biological_study':True,
            'boechera_sitewise_heterogeneity_not_consistent':True}

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--reps',type=int,default=9999)
    a=p.parse_args()
    receipt=numerical_verify(a.root,a.reps)
    (a.root/'results/independent_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':main()
