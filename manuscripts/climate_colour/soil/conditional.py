"""Paired secondary check: within distance bands, adjust soil contrast for both climate contrasts.
Usage: python conditional.py OUTPUT_DIR
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

out=Path(sys.argv[1]); d=pd.read_csv(out/'complete_sites.csv.gz')
main=pd.read_csv(out/'species_effects.csv').query("panel=='main'").set_index('taxon')
props=['phh2o','clay','nitrogen']; cols=props+['bio1','bio5','bio12','bio14']
rows=[]
for tax,g in d.groupby('inat_taxon_id'):
    if tax not in main.index or g.cell.nunique()<10:continue
    v=g[cols].to_numpy(dtype=float);v=(v-v.mean(0))/v.std(0)
    W=np.column_stack([v[:,:3],(v[:,3]+v[:,4])/2,(v[:,5]+v[:,6])/2])
    W=(W-W.mean(0))/W.std(0)
    i,j=np.triu_indices(len(g),1)
    keep=(g.cell.to_numpy()[i]!=g.cell.to_numpy()[j])&(g.observer_id.to_numpy()[i]!=g.observer_id.to_numpy()[j])
    i=i[keep];j=j[keep]
    la=np.radians(g.latitude.to_numpy());lo=np.radians(g.longitude.to_numpy())
    h=np.sin((la[i]-la[j])/2)**2+np.cos(la[i])*np.cos(la[j])*np.sin((lo[i]-lo[j])/2)**2
    km=12742.0176*np.arcsin(np.sqrt(np.clip(h,0,1)))
    b=np.floor(np.log1p(km)/np.log1p(.05)).astype(int)
    y=(g.morph.to_numpy()[i]!=g.morph.to_numpy()[j]).astype(float)
    n=np.bincount(b);m=np.bincount(b,weights=y)
    k=(m[b]>0)&(m[b]<n[b])
    if k.sum()<20:continue
    b=b[k];n=np.bincount(b);a=np.column_stack([y[k],np.abs(W[i[k]]-W[j[k]])])
    for t in range(a.shape[1]):
        a[:,t]-=np.bincount(b,weights=a[:,t])[b]/n[b]
    col=a[:,0]; climate=a[:,4:6]
    resid=col-climate@np.linalg.lstsq(climate,col,rcond=None)[0]
    if resid@resid<.10*(col@col):continue
    unadj=(col@a[:,1:4])/(col@col);adj=(resid@a[:,1:4])/(resid@resid)
    assert np.allclose(unadj,main.loc[tax,props].to_numpy(dtype=float),atol=1e-5)
    rows.append([tax,main.loc[tax,'genus'],*unadj,*adj])
df=pd.DataFrame(rows,columns=['taxon','genus',*props,*['adj_'+x for x in props]])
df.to_csv(out/'conditional_species.csv',index=False);cols=df.columns[2:]
s=df.groupby('genus')[cols].sum().to_numpy();n=df.groupby('genus').size().to_numpy()
rng=np.random.default_rng(20261111);boot=[]
for _ in range(50):
    idx=rng.integers(len(n),size=(200,len(n)))
    boot.append(s[idx].sum(1)/n[idx].sum(1)[:,None])
lo,hi=np.quantile(np.vstack(boot),[.025/3,1-.025/3],axis=0)
res=pd.DataFrame(dict(variable=cols,mean=df[cols].mean().values,lower=lo,upper=hi,species=len(df),genera=len(n)))
res.to_csv(out/'conditional_summary.csv',index=False);print(res.round(4))

