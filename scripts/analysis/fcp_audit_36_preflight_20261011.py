import io, zipfile, os, numpy as np,pandas as pd
from pathlib import Path
ROOT=Path(os.environ.get('FCP_SOURCE_ROOT','/mnt/data'))
A=ROOT/'fcp_original_photo_site_climate_soil_85337_20261010.zip'
B=ROOT/'fcp_100543_all36_expanded_photo_site_environment_20261011.zip'
CLASSES=('white','yellow_orange','red_pink','blue_purple')
BLOCKS={"elevation":('wc_elevation_m',),"temperature":tuple('wc_bio'+str(i) for i in range(1,12)),"precipitation":tuple('wc_bio'+str(i) for i in range(12,20)),"solar_radiation":('wc_srad_annual_mean','wc_srad_monthly_cv'),"wind":('wc_wind_annual_mean','wc_wind_monthly_cv'),"vapor_pressure":('wc_vapr_annual_mean','wc_vapr_monthly_cv'),"soil":('soil_pH','soil_SOC','soil_N','soil_clay','soil_available_water_proxy','soil_cec_0_30cm_source_raw','soil_sand_0_30cm_source_raw','soil_silt_0_30cm_source_raw','soil_bdod_0_30cm_source_raw','soil_cfvo_0_30cm_source_raw')}
FEATURES=tuple(k for v in BLOCKS.values() for k in v)
assert len(FEATURES)==36

def gz(path, inner):
 with zipfile.ZipFile(path) as z: return pd.read_csv(io.BytesIO(z.read(inner)), compression='gzip',low_memory=False)

def run():
 a=gz(A,'fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz')
 b=gz(A,'fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz')
 ext=gz(B,'fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz')
 assert len(a)==42111 and len(b)==85337 and len(ext)==100543
 source=pd.concat([a[['inat_taxon_id','observation_id','photo_id','species','morph','measurement_status','flower_effective_pixels','nuisance_pixel_fraction']],b[['inat_taxon_id','observation_id','photo_id','species','morph','measurement_status','flower_effective_pixels','nuisance_pixel_fraction']]],ignore_index=True)
 duplicated=source[source.photo_id.duplicated(keep=False)]
 assert duplicated.photo_id.nunique()==26905
 assert duplicated.groupby('photo_id').inat_taxon_id.nunique().max()==1
 assert duplicated[duplicated.morph.isin(CLASSES)].groupby('photo_id').morph.nunique().max()==1
 source=source.drop_duplicates('photo_id')
 assert len(source)==100543 and source.inat_taxon_id.nunique()==42111
 flags=source.measurement_status.eq('classified_four_state_morph')&source.morph.isin(CLASSES)
 assert flags.sum()==45953
 env=ext[['inat_taxon_id','observation_id','photo_id','latitude','longitude',*FEATURES]]
 d=source.merge(env,on=['inat_taxon_id','observation_id','photo_id'],how='left',validate='one_to_one',indicator=True)
 flags=d.measurement_status.eq('classified_four_state_morph') & d.morph.isin(CLASSES)
 assert int(flags.sum())==45953
 assert d['_merge'].eq('both').all()
 d['abs_lat']=d.latitude.abs()
 d['lon_sin']=np.sin(np.deg2rad(d.longitude))
 d['lon_cos']=np.cos(np.deg2rad(d.longitude))
 d['genus']=d.species.astype(str).str.split().str[0]
 d['east']=d.longitude.ge(90)&d.longitude.le(180)
 climate=[x for k,v in BLOCKS.items() if k!='soil' for x in v]
 base=flags & d[['latitude','longitude','abs_lat','lon_sin','lon_cos',*climate]].notna().all(axis=1)
 soil=base & d[list(BLOCKS['soil'])].notna().all(axis=1)
 info={'all_rows':len(d),'classified':int(flags.sum()),'climate_complete':int(base.sum()),'soil_complete':int(soil.sum()),'climate_species':int(d.loc[base,'inat_taxon_id'].nunique()),'soil_species':int(d.loc[soil,'inat_taxon_id'].nunique()),'climate_genera':int(d.loc[base,'genus'].nunique()),'soil_genera':int(d.loc[soil,'genus'].nunique())}
 print(info,flush=True)
 for mode,mask in [('climate',base),('soil',soil)]:
  q=d[mask].copy()
  for label,case in [('world',q),('east',q[q.east]),('elsewhere',q[~q.east])]:
   sp=case.groupby('inat_taxon_id',sort=False).agg(n=('photo_id','size'),genus=('genus','first'))
   ct=sp.genus.value_counts()
   print(mode,label,'photo',len(case),'species',len(sp),'genus>=3',int((ct>=3).sum()),'species in >=3 genera',int(ct[ct>=3].sum()),flush=True)
 return d, info, {'climate':base,'soil':soil}
if __name__=='__main__':run()