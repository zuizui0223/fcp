"""Usage: python range_soil.py measured.csv.gz ledger.csv Table_S2.csv.gz soil_dir bio_dir output_dir
bio_dir contains archived wc2.1_10m_bio.zip. soil_dir contains the three TIFFs and manifest.json.
"""
import hashlib,json,sys,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.warp import reproject,Resampling

photos,ledger,s2,soil,bio,out=map(Path,sys.argv[1:])
out.mkdir(exist_ok=True,parents=True)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
assert sha(photos)=='8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6'
arc=bio/'wc2.1_10m_bio.zip'
assert sha(arc)=='00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf'
d=pd.read_csv(ledger)
ids=d.loc[d.filter(regex='^miss_coarse_photo_').ge(5).sum(axis=1)>=2,'inat_taxon_id']
assert len(ids)==2400
cols='inat_taxon_id photo_id latitude longitude morph observed_on observer_id global_classifiable positional_accuracy_m'.split()
d=pd.read_csv(photos,usecols=cols,low_memory=False); d=d[d.inat_taxon_id.isin(ids)]
assert len(d)==240000 and d.global_classifiable.sum()==111059
d=d[d.global_classifiable & d.latitude.between(-90,90) & d.longitude.between(-180,180) & d.positional_accuracy_m.between(0,5000) & d.observer_id.notna()].copy()
d['month']=pd.to_datetime(d.observed_on,errors='coerce',format='mixed',utc=True).dt.month
ngeo=len(d)
with zipfile.ZipFile(arc) as z:
    for k in [1,5,12,14]:
        name=f'wc2.1_10m_bio_{k}.tif'
        z.extract(name,bio)
        with rio.open(bio/name) as src:
            transform=src.transform;crs=src.crs; shape=src.shape
            row,col=rio.transform.rowcol(transform,d.longitude.to_numpy(),d.latitude.to_numpy())
            row=np.clip(row,0,shape[0]-1);col=np.clip(col,0,shape[1]-1)
            d[f'bio{k}']=src.read(1,masked=True).astype(float).filled(np.nan)[row,col]
d['cell']=row*shape[1]+col
manifest=json.loads((soil/'manifest.json').read_text())
factors={'phh2o':10,'clay':10,'nitrogen':100}
for spec in manifest:
    key=spec['property']
    if key not in factors:continue
    path=soil/spec['file'];assert sha(path)==spec['sha256']
    with rio.open(path) as src:
        target=np.full(shape,np.nan,np.float32)
        reproject(rio.band(src,1),target,src_transform=src.transform,src_crs=src.crs,
                  src_nodata=src.nodata,dst_transform=transform,dst_crs=crs,dst_nodata=np.nan,
                  resampling=Resampling.average,num_threads=4,warp_mem_limit=128)
        d[key]=target[row,col]/factors[key]
    print('sampled',key,flush=True)
features=['phh2o','clay','nitrogen','temperature','precipitation']
base=['phh2o','clay','nitrogen','bio1','bio5','bio12','bio14']
d=d.dropna(subset=base)
d.to_csv(out/'complete_sites.csv.gz',index=False,compression='gzip')
genera=pd.read_csv(s2,usecols=['inat_taxon_id','genus']).set_index('inat_taxon_id').genus
def z(v):return (v-v.mean(axis=0))/v.std(axis=0,ddof=0)
def estimate(km,diff,E,width):
    bins=np.floor(np.log1p(km)/np.log1p(width)).astype(int)
    N=bins.max()+1
    same,other=np.bincount(bins[~diff],minlength=N),np.bincount(bins[diff],minlength=N)
    good=(same>0)&(other>0)
    if not good.any():return None
    # Within-bin contrasts weighted as the coefficient of mismatch in a bin fixed-effect regression.
    w=same[good]*other[good]/(same[good]+other[good])
    ss=[]
    for v in E.T:
        a=np.bincount(bins[~diff],weights=v[~diff],minlength=N)[good]/same[good]
        b=np.bincount(bins[diff],weights=v[diff],minlength=N)[good]/other[good]
        ss.append(np.average(b-a,weights=w))
    return int((same[good]+other[good]).sum()),int(other[good].sum()),ss
rows=[];support=[]
for tax,g in d.groupby('inat_taxon_id',sort=False):
    for thin in [False,True]:
        sub=g.sort_values('photo_id').drop_duplicates('observer_id') if thin else g
        raw=sub[base].to_numpy(dtype=float)
        if len(raw)<3 or np.any(raw.std(axis=0)<1e-10):continue
        b=z(raw); u=np.column_stack((b[:,:3],(b[:,3]+b[:,4])/2,(b[:,5]+b[:,6])/2))
        if np.any(u.std(axis=0)<1e-10):continue
        u=z(u)
        i,j=np.triu_indices(len(sub),1)
        sel=(sub.observer_id.to_numpy()[i]!=sub.observer_id.to_numpy()[j])&(sub.cell.to_numpy()[i]!=sub.cell.to_numpy()[j])
        i=i[sel];j=j[sel]
        if not len(i):continue
        lat=np.radians(sub.latitude.to_numpy()); lon=np.radians(sub.longitude.to_numpy())
        h=np.sin((lat[i]-lat[j])/2)**2+np.cos(lat[i])*np.cos(lat[j])*np.sin((lon[i]-lon[j])/2)**2
        km=12742.0176*np.arcsin(np.sqrt(np.clip(h,0,1)))
        diff=sub.morph.to_numpy()[i]!=sub.morph.to_numpy()[j]
        E=np.abs(u[i]-u[j]);mon=sub.month.to_numpy()
        for panel,mask,width in ([('one_per_observer',np.ones(len(i),bool),.05)] if thin else [
                ('main',np.ones(len(i),bool),.05),('narrow_distance',np.ones(len(i),bool),.025),
                ('same_calendar_month',mon[i]==mon[j],.05)]):
            if not mask.any():continue
            ans=estimate(km[mask],diff[mask],E[mask],width)
            if ans is None:continue
            pairs,mism,eff=ans
            rows.append([panel,tax,genera[tax],pairs,mism,*eff])
df=pd.DataFrame(rows,columns=['panel','taxon','genus','pairs','discordant',*features])
df.to_csv(out/'species_effects.csv',index=False)
results=[]
for panel,g in df.groupby('panel'):
    sums=g.groupby('genus')[features].sum().to_numpy(); n=g.groupby('genus').size().to_numpy()
    rng=np.random.default_rng(20261111);boots=[]
    for _ in range(50):
        draw=rng.integers(len(n),size=(200,len(n)))
        boots.append(sums[draw].sum(axis=1)/n[draw].sum(axis=1)[:,None])
    boot=np.vstack(boots);lo,hi=np.quantile(boot,[.025/5,1-.025/5],axis=0)
    for k,name in enumerate(features):
        results.append([panel,name,len(g),len(n),g.pairs.sum(),g.discordant.sum(),g[name].mean(),lo[k],hi[k],(g[name]>0).mean()])
res=pd.DataFrame(results,columns=['panel','property','species','genera','pairs','discordant','effect_sd','lower','upper','fraction_positive'])
res.to_csv(out/'results.csv',index=False);print(res.round(4).to_string(index=False))
print('geo',ngeo,'complete',len(d),d.inat_taxon_id.nunique())
(out/'coverage.json').write_text(json.dumps(dict(eligible_taxa=2400,photos=240000,classifiable=111059,coordinate_observer_eligible=ngeo,common_sites=len(d),common_taxa=d.inat_taxon_id.nunique()),indent=2))
