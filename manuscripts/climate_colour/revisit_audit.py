"""Usage: python revisit_audit.py MEASURED.csv.gz BIO.zip SRAD.zip OUTPUT_DIR"""
import hashlib
from pathlib import Path
import sys
import zipfile
import numpy as np
import pandas as pd
from rasterio.io import MemoryFile

src, bio, solar, out = map(Path, sys.argv[1:])
out.mkdir(exist_ok=True, parents=True)
for path, digest in [(src,'8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6'),
                     (bio,'00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf'),
                     (solar,'c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026')]:
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
root=Path(__file__).resolve().parent
species=pd.read_csv(root/'Table_S2_species_coefficients.csv.gz').set_index('inat_taxon_id')
freq=pd.read_csv(root/'Table_S6_white_frequencies.csv.gz').set_index('inat_taxon_id')
d=pd.read_csv(src,usecols='inat_taxon_id photo_id observer_id observed_on latitude longitude global_classifiable morph positional_accuracy_m'.split())
d=d[d.inat_taxon_id.isin(species.index)].copy()
assert len(d)==447500 and d.global_classifiable.sum()==227107 and d.photo_id.is_unique
date=pd.to_datetime(d.observed_on,format='mixed',utc=True,errors='coerce')
d['month'],d['year']=date.dt.month,date.dt.year
d=d[d.latitude.between(-90,90)&d.longitude.between(-180,180)&d.month.notna()].copy()
def raster(arc,name):
    with zipfile.ZipFile(arc) as z, MemoryFile(z.read(name)) as mem, mem.open() as ds:
        assert ds.width==2160 and ds.height==1080 and abs(ds.transform.a-1/6)<1e-9
        row=np.clip(np.floor((ds.transform.f-d.latitude)/(-ds.transform.e)).astype(int),0,1079)
        col=np.clip(np.floor((d.longitude-ds.transform.c)/ds.transform.a).astype(int),0,2159)
        d['cell']=row*2160+col
        return ds.read(1,masked=True).astype(float).filled(np.nan)[row,col]
for k in (1,5,12,14,4,15):
    d[f'b{k}']=raster(bio,f'wc2.1_10m_bio_{k}.tif')
d['sun']=np.mean([raster(solar,f'wc2.1_10m_srad_{m:02}.tif') for m in range(1,13)],axis=0)
climate=['b1','b5','b12','b14','b4','b15','sun']
rows=[]
for tax,g in d.dropna(subset=climate).groupby('inat_taxon_id'):
    if g.cell.nunique()<5 or g.observer_id.nunique()<20 or len(g)<40: continue
    x=g[climate].to_numpy()
    if np.any(x.std(axis=0)==0):continue
    x=(x-x.mean(axis=0))/x.std(axis=0)
    t=x[:,:2].mean(axis=1)
    lat,lon=np.deg2rad(g[['latitude','longitude']].to_numpy()).T
    design=np.column_stack([np.eye(12)[g.month.to_numpy(int)-1],np.cos(lat)*np.cos(lon),
        np.cos(lat)*np.sin(lon),np.sin(lat),x[:,2:4].mean(axis=1),x[:,4:]])
    resid=t-design@np.linalg.lstsq(design,t,rcond=1e-10)[0]
    if resid.var()<.05*t.var() or resid.var()<1e-12:continue
    rows.append((tax,len(g),g.global_classifiable.sum(),100*resid@g.global_classifiable.to_numpy(float)/(resid@resid)))
b=pd.DataFrame(rows,columns=['inat_taxon_id','n','classified','success_pp']).set_index('inat_taxon_id').join(species[['genus','thermal_level_white']])
b.to_csv(out/'Table_S7_classification.csv.gz')
for name,q in [('all',b),('thermal',b.dropna(subset=['thermal_level_white']))]:
    by=q.groupby('genus').success_pp.agg(['sum','count']).to_numpy()
    rng=np.random.default_rng(273)
    vals=[]
    for _ in range(40):
        draw=by[rng.integers(len(by),size=(100,len(by)))].sum(axis=1)
        vals.extend(draw[:,0]/draw[:,1])
    print(name,len(q),q.n.sum(),q.classified.sum(),q.success_pp.mean(),np.percentile(vals,[2.5,97.5]),flush=True)

# Choose earliest eligible pair of years without using colour. Within a year retain
# one photo per observer. The second year excludes all retained first-year observers.
q=d[d.inat_taxon_id.isin(freq.index[freq.repeated])&d.global_classifiable&
    d.positional_accuracy_m.between(0,5000)&d.observer_id.notna()].copy()
keys=['inat_taxon_id','month','cell','year']
q=q.sort_values('photo_id').drop_duplicates(keys+['observer_id'])
q['white']=q.morph.eq('white')
blocks=[]
for (tax,month),pool in q.groupby(['inat_taxon_id','month']):
    if len(pool)<4:continue
    pool=pool.reset_index(drop=True)
    pairs=[]
    for cell,g in pool.groupby('cell'):
        years=[(yr,h) for yr,h in g.groupby('year') if len(h)>=2]
        chosen=False
        for i,(yr,h) in enumerate(years):
            for yr2,h2 in years[i+1:]:
                follow=h2[~h2.observer_id.isin(h.observer_id)]
                if len(follow)>=2:
                    pairs.append((cell,yr,yr2,h.index.to_numpy(),follow.index.to_numpy()))
                    chosen=True
                    break
            if chosen:break
    if pairs:blocks.append((tax,month,pool.white.to_numpy(),pairs))
rows=[]; rng=np.random.default_rng(273); sims=np.zeros(9999,int)
for tax,month,white,pairs in blocks:
    # All index-year labels in these windows are fixed, including unmixed windows.
    fixed=np.concatenate([h for _,_,_,h,_ in pairs])
    remaining=np.setdiff1d(np.arange(len(white)),fixed)
    perm=white[remaining][np.argsort(rng.random((len(sims),len(remaining))),axis=1)]
    for cell,yr,yr2,h,f in pairs:
        baseline=white[h].sum(); actual=white[f].sum()
        draws=perm[:,np.searchsorted(remaining,f)].sum(axis=1)
        sim=(draws>0)&(draws<len(f))
        if 0<baseline<len(h):sims+=sim
        rows.append((tax,cell,month,yr,yr2,len(h),baseline,len(f),actual,sim.mean()))
r=pd.DataFrame(rows,columns='taxon cell month first_year later_year n0 white0 n1 white1 null_mixed'.split())
r.to_csv(out/'Table_S8_revisits.csv',index=False)
s=r[(r.white0>0)&(r.white0<r.n0)]
hit=s[(s.white1>0)&(s.white1<s.n1)]
print('eligible',len(r),r.taxon.nunique(),'initial mixed',len(s),s.taxon.nunique(),'redetected',len(hit),hit.taxon.nunique(),
      'null count mean/interval',sims.mean(),np.percentile(sims,[2.5,97.5]),'upper p',(1+sum(sims>=len(hit)))/(len(sims)+1))
