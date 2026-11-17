"""Paired thermal slopes. Supply complete_sites.csv.gz, measured photos, Table S2, WorldClim BIO and solar ZIPs, output CSV."""
from pathlib import Path
import sys
import zipfile

import numpy as np
import pandas as pd
from rasterio.io import MemoryFile

sites, photos, species, bio, solar, output = map(Path, sys.argv[1:])
d = pd.read_csv(sites).merge(pd.read_csv(photos, usecols=['photo_id','flower_fraction_white']), on='photo_id', validate='one_to_one')
for archive, fields in [(bio,[('bio4','wc2.1_10m_bio_4.tif'),('bio15','wc2.1_10m_bio_15.tif')]),
                        (solar,[(f'sun{i}',f'wc2.1_10m_srad_{i:02}.tif') for i in range(1,13)])]:
    with zipfile.ZipFile(archive) as z:
        for col, filename in fields:
            with MemoryFile(z.read(filename)) as mem, mem.open() as im:
                assert im.shape == (1080,2160) and abs(im.transform.a-1/6)<1e-7 and im.transform.c == -180 and im.transform.f == 90
                d[col] = im.read(1,masked=True).astype(float).filled(np.nan).ravel()[d.cell.to_numpy(int)]
d['sun'] = d[[f'sun{i}' for i in range(1,13)]].mean(axis=1,skipna=False)
cols=['bio1','bio5','bio12','bio14','bio4','bio15','sun','phh2o','clay','nitrogen']
d=d.dropna(subset=cols+['month','flower_fraction_white'])
rows=[]
for tax,g in d.groupby('inat_taxon_id'):
    if len(g)<40 or g.observer_id.nunique()<20 or g.cell.nunique()<5:continue
    x=g[cols].to_numpy(float)
    scale=x.std(axis=0)
    if (scale[:2]<1e-9).any():continue
    x=(x-x.mean(axis=0))/np.where(scale>1e-9,scale,1)
    t=x[:,:2].mean(axis=1)
    lat,lon=np.deg2rad(g[['latitude','longitude']].to_numpy()).T
    base=np.column_stack([np.eye(12)[g.month.to_numpy(int)-1],np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat),x[:,2:4].mean(axis=1),x[:,4:7]])
    slopes=[]
    for design in [base,np.column_stack([base,x[:,7:]])]:
        resid=t-design@np.linalg.lstsq(design,t,rcond=1e-10)[0]
        if t.var()<1e-9 or resid.var()<.05*t.var():break
        slopes.append(resid@g.flower_fraction_white/(resid@resid))
    if len(slopes)==2:rows.append((tax,len(g),*slopes))
out=pd.DataFrame(rows,columns=['inat_taxon_id','photos','climate','climate_soil'])
meta=pd.read_csv(species,usecols=['inat_taxon_id','genus','source_cohort'])
out=out.merge(meta,on='inat_taxon_id',validate='one_to_one')
out['change']=out.climate_soil-out.climate
out.to_csv(output,index=False)
