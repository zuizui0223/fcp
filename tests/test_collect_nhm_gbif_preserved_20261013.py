import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'/'analysis'))
from collect_nhm_gbif_preserved_20261013 import urls,valid

def example():
 return dict(basisOfRecord='PRESERVED_SPECIMEN',datasetKey='abcd',decimalLatitude=41.0,decimalLongitude=131.5,
        key=12345,species='Teucrium depressum',occurrenceStatus='PRESENT',year=1990,issues=[])
def test_query_hashed_offsets_are_disjoint_and_bounded():
 for n in (0,1,119,120,121,240,241,300,12300,100000):
  pairs=urls(42,n)
  assert sum(v for _,v in pairs)<=min(n,240)
  for off,lim in pairs: assert 0<=off and off+lim<=n
  if len(pairs)==2: assert pairs[0][0]+pairs[0][1]<=pairs[1][0]
def test_only_preserved_specimens():
 r=example();x,b=valid(r,'Teucrium depressum',546);assert b is None and x['basisOfRecord']=='PRESERVED_SPECIMEN'
 r['basisOfRecord']='HUMAN_OBSERVATION';assert valid(r,'Teucrium depressum',546)[1]=='not_preserved_specimen'
def test_inaturalist_dataset_excluded_even_if_mislabeled():
 r=example();r['datasetKey']='50c9509d-22c7-4a22-a47d-8c48425ef4a7'
 assert valid(r,'Teucrium depressum',546)[1]=='inat_dataset'
def test_no_genus_or_synonym_fallback():
 r=example();r['species']='Teucrium other';assert valid(r,'Teucrium depressum',546)[1]=='interpreted_species_mismatch'
def test_coordinates_and_uncertainty_hold():
 r=example();r['decimalLatitude']=0;r['decimalLongitude']=0
 assert valid(r,'Teucrium depressum',546)[1]=='invalid_coordinate'
 r=example();r['coordinateUncertaintyInMeters']=40000
 assert valid(r,'Teucrium depressum',546)[1]=='uncertainty_gt_20km'
def test_no_colour_estimation_from_gbif_data():
 r=example();r['flower_color']='white';x,b=valid(r,'Teucrium depressum',546)
 assert b is None and 'flower_color' not in x and 'botanical_colour' not in x
