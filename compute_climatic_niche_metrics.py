#!/usr/bin/env python3
"""Extract WorldClim values and compute comparable occupied climatic niche metrics."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,rasterio
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

SELECT=[1,4,5,6,7,12,14,15,17]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--occurrences',required=True);ap.add_argument('--worldclim-dir',required=True);ap.add_argument('--out',required=True);ap.add_argument('--qc-out',required=True);ap.add_argument('--min-cells',type=int,default=10);a=ap.parse_args()
    d=pd.read_csv(a.occurrences); coords=list(zip(d.decimalLongitude,d.decimalLatitude)); vals=[]
    for b in SELECT:
        candidates=list(Path(a.worldclim_dir).glob(f'*bio_{b}.tif'))+list(Path(a.worldclim_dir).glob(f'*bio{b}.tif'))
        if not candidates: raise FileNotFoundError(f'BIO{b} raster not found')
        with rasterio.open(candidates[0]) as src: vals.append([x[0] if x and x[0]!=src.nodata else np.nan for x in src.sample(coords)])
    for b,v in zip(SELECT,vals):d[f'bio{b}']=v
    climate=[f'bio{x}' for x in SELECT]; d=d.dropna(subset=climate).copy()
    d=d.drop_duplicates(['canonical_name']+climate)
    scaler=StandardScaler(); z=scaler.fit_transform(d[climate]); pca=PCA(n_components=3).fit(z); pcs=pca.transform(z)
    for j,c in enumerate(climate):d[f'z_{c}']=z[:,j]
    for i in range(3):d[f'pc{i+1}']=pcs[:,i]
    result=[]
    for keys,g in d.groupby(['canonical_name','family','role','focal_species','match_level'],dropna=False):
        n=len(g); row=dict(zip(['canonical_name','family','role','focal_species','match_level'],keys));row['n_climate_cells']=n
        if n<a.min_cells:
            row['metric_status']='insufficient_cells';result.append(row);continue
        q=g[climate].quantile([.05,.95])
        row.update({f'{c}_q95q05':float(q.loc[.95,c]-q.loc[.05,c]) for c in climate})
        row['temperature_breadth']=float(np.mean([row[f'bio{x}_q95q05'] for x in [1,5,6,7]]))
        row['moisture_breadth']=float(np.mean([row[f'bio{x}_q95q05'] for x in [12,14,15,17]]))
        row['climatic_heterogeneity']=float(np.mean([g[f'z_{c}'].std(ddof=0) for c in climate]))
        p=g[['pc1','pc2','pc3']].to_numpy(); ctr=p.mean(axis=0)
        row['pca_dispersion']=float(np.mean(np.linalg.norm(p-ctr,axis=1)))
        row['pc1_q95q05']=float(g.pc1.quantile(.95)-g.pc1.quantile(.05));row['pc2_q95q05']=float(g.pc2.quantile(.95)-g.pc2.quantile(.05))
        try:row['pca_hull_area']=float(ConvexHull(g[['pc1','pc2']].to_numpy()).volume) if n>=3 else np.nan
        except Exception:row['pca_hull_area']=np.nan
        row['metric_status']='complete';result.append(row)
    out=pd.DataFrame(result);Path(a.out).parent.mkdir(parents=True,exist_ok=True);out.to_csv(a.out,index=False)
    qc={'occurrence_rows_with_climate':int(len(d)),'species_total':int(out.canonical_name.nunique()),'species_complete':int((out.metric_status=='complete').sum()),'min_cells':a.min_cells,'worldclim_bio_variables':SELECT,'global_pca_variance_explained':pca.explained_variance_ratio_.tolist(),'semantic_guard':'metrics describe occupied climatic breadth and heterogeneity, not physiological tolerance or causal niche expansion'}
    Path(a.qc_out).write_text(json.dumps(qc,indent=2)+'\n');print(json.dumps(qc))
if __name__=='__main__':main()
