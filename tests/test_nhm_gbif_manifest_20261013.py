import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from freeze_nhm_gbif_manifest_20261013 import choose,SALT

def fake(n_genera=110):
 return pd.DataFrame([{'taxon':f'Genus{g} {color}{k}','label':label}
  for g in range(n_genera) for label,color in [('white','alba'),('blue_purple','caerula')]
  for k in range(5)])
def test_deterministic_source_independent_selection():
 a,z=choose(fake());b,z2=choose(fake().sample(frac=1,random_state=42));assert a.equals(b);assert z==z2
 assert len(a.genus.unique())==96 and a.groupby('genus').size().eq(8).all()
def test_min_three_per_colour_entire_genus_required():
 d=fake(2);d=d[~((d.taxon.str.startswith('Genus1'))&(d.label=='white')&(d.taxon.str.endswith('3')|d.taxon.str.endswith('4')))];a,z=choose(d);assert len(a.genus.unique())==2
 d=d[~((d.taxon.str.startswith('Genus1'))&(d.label=='white')&(d.taxon.str.endswith('2')))];a,z=choose(d);assert set(a.genus)=={'Genus0'}
def test_single_morph_not_accepted():
 d=fake(3);d=d[~((d.taxon.str.startswith('Genus2'))&(d.label=='blue_purple'))];a,_=choose(d);assert set(a.genus)=={'Genus0','Genus1'}
def test_source_gap_not_replaced_by_other_hues():
 d=fake(2);d.loc[d.genus if 'genus' in d else d.taxon.str.startswith('Genus1'),'label']='yellow'
 a,_=choose(d);assert len(a)>0 and a.botanical_colour.isin(('white','blue_purple')).all()
