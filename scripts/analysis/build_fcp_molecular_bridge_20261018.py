#!/usr/bin/env python3
"""Build an offline, provenance-preserving species bridge from frozen evidence.

No raw read reanalysis and no mechanism imputation. Exact normalized names only;
synonyms, author stripping, genus-only matches and fuzzy matching are forbidden.
The cross-species join supports an evidence audit, not an individual mechanism test.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
EVID=ROOT/'data/evidence/fcp_molecular_bridge_20261018'
ROSTER=ROOT/'results/fcp_real_447500_white_axis_decomposition_20261018/species_photo_composition.csv'
ENV=ROOT/'results/fcp_real_447500_environment_description_20261018/primary_ge5_species_environment.csv'
OUT=ROOT/'results/fcp_molecular_bridge_20261018'

def key(x):
    return re.sub(r'\s+', ' ',str(x).strip()).casefold()
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(d,p):
    p.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def table(d,p):
    d.to_csv(p,index=False,float_format='%.12g',lineterminator='\n')
def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--evidence-dir',type=Path,default=EVID)
    a.add_argument('--roster',type=Path,default=ROSTER)
    a.add_argument('--environment',type=Path,default=ENV)
    a.add_argument('--outdir',type=Path,default=OUT)
    z=a.parse_args();z.outdir.mkdir(parents=True,exist_ok=True)
    r=pd.read_csv(z.roster);e=pd.read_csv(z.environment)
    assert len(r)==4475 and r.inat_taxon_id.is_unique
    assert not e.duplicated(['inat_taxon_id','variable']).any()
    assert e.inat_taxon_id.nunique()==1574 and len(e)==11018
    assert set(e.inat_taxon_id)<=set(r.inat_taxon_id)
    assert e.groupby('inat_taxon_id').variable.nunique().eq(7).all()
    r['match_key']=r.species.map(key);assert r.match_key.is_unique
    lookup=r.set_index('match_key').inat_taxon_id.to_dict()
    panel=set(e.inat_taxon_id)
    c=json.loads((z.evidence_dir/'curated_contrasts.json').read_text())['contrasts']
    assert len({x['contrast_id'] for x in c})==len(c)
    records=[]
    for x in c:
        for tissue in x['tissues']:
            y={k:v for k,v in x.items() if k not in ['tissues','tissue_results']}
            y.update(tissue=tissue,tissue_result=x.get('tissue_results',{}).get(tissue,'see contrast-level evidence; no organ-specific quantitative result extracted'))
            y['match_key']=key(y['species']);tid=lookup.get(y['match_key'])
            y.update(inat_taxon_id=tid,roster_match_status='exact_normalized_name' if tid else 'unmatched_no_synonym_inference',environment_panel_match=tid in panel)
            for k,v in y.items():
                if isinstance(v,list):y[k]=';'.join(v)
            records.append(y)
    ev=pd.DataFrame(records)
    assert not ev.duplicated(['contrast_id','tissue']).any()
    table(ev,z.outdir/'contrast_tissue_evidence.csv')
    core=ev.drop_duplicates('species')[['species','inat_taxon_id','roster_match_status','environment_panel_match']]
    table(core,z.outdir/'curated_taxon_matches.csv')
    molecular_species=set(ev.loc[ev.molecular_tier.str.startswith(('T2','T3','T4')),'species'])
    annotated=set(ev.species)
    annotated_ids=set(ev.inat_taxon_id.dropna().astype(int))
    molecular_ids=set(ev.loc[ev.species.isin(molecular_species),'inat_taxon_id'].dropna().astype(int))
    r['in_environment_panel']=r.inat_taxon_id.isin(panel)
    r['curated_any_evidence']=r.inat_taxon_id.isin(annotated_ids)
    r['curated_molecular_tier_ge2']=r.inat_taxon_id.isin(molecular_ids)
    r['strict_same_FCP_individual_mechanism_link']=False
    r['missingness_meaning']=r.curated_any_evidence.map({True:'species_only_annotation_not_verified_for_FCP_morph',False:'not_curated_in_this_bounded_bridge_not_evidence_of_no_mechanism'})
    table(r.drop(columns='match_key'),z.outdir/'all_4475_species_coverage.csv')
    # One row per species x covariate. The 34 tissue rows never multiply this table.
    species_traits=ev.groupby('inat_taxon_id',dropna=True).agg(curated_contrasts=('contrast_id','nunique'),tissues_recorded=('tissue','nunique'),molecular_tiers=('molecular_tier',lambda v:';'.join(sorted(set(v)))),pigment_classes=('pigment_class',lambda v:';'.join(sorted(set(v)))),phenotype_origins=('phenotype_origin',lambda v:';'.join(sorted(set(v)))),source_dois=('source_dois',lambda v:';'.join(sorted({t for x in v for t in x.split(';')})))).reset_index()
    linked=e.merge(species_traits,on='inat_taxon_id',how='inner',validate='many_to_one')
    linked['linkage_level']='species_only_not_FCP_genotype_or_pigment_verified'
    table(linked,z.outdir/'species_environment_evidence_bridge.csv')
    screen=pd.read_csv(z.evidence_dir/'ena_title_screen_runs.tsv',sep='\t').fillna('')
    assert screen.run_accession.is_unique
    decisions=json.loads((z.evidence_dir/'screen_decisions.json').read_text())
    projects=[]
    for (pr,sp),d in screen.groupby(['study_accession','scientific_name'],sort=True):
        tid=lookup.get(key(sp));decision=decisions.get(pr,{})
        projects.append(dict(study_accession=pr,submitted_scientific_name=sp,study_title=d.study_title.iloc[0],run_count=d.run_accession.nunique(),experiment_count=d.experiment_accession.nunique(),biosample_count=d.sample_accession.nunique(),inat_taxon_id=tid,roster_exact_match=tid is not None,environment_panel_match=tid in panel,review_status=decision.get('status','pending_primary_review'),decision_reason=decision.get('reason','metadata alone does not establish a contrast or mechanism'),possible_title_taxon=decision.get('possible_title_taxon',''),source_url='https://www.ebi.ac.uk/ena/browser/view/'+pr))
    ps=pd.DataFrame(projects);table(ps,z.outdir/'ena_project_screen.csv')
    screen['inat_taxon_id']=screen.scientific_name.map(lambda x:lookup.get(key(x)))
    table(screen,z.outdir/'ena_run_roster_matches.csv')
    # Combine direct accession audit and screen metadata; deduplicate shared runs.
    allmeta=[]
    for p in sorted(z.evidence_dir.glob('ena_*_runs.tsv')):
        m=pd.read_csv(p,sep='\t').fillna('');m['retrieval_file']=p.name;allmeta.append(m)
    m=pd.concat(allmeta,ignore_index=True).drop_duplicates('run_accession').sort_values('run_accession')
    m['is_curated_silene_RNA']=m.scientific_name.eq('Silene littorea') & m.library_strategy.eq('RNA-Seq')
    table(m,z.outdir/'public_run_metadata_audit.csv')
    lib=[]
    for (pr,sp,sa,ex,ln),d in m.groupby(['study_accession','scientific_name','sample_accession','experiment_accession','library_name'],dropna=False,sort=True):
        lib.append(dict(study_accession=pr,scientific_name=sp,sample_accession=sa,experiment_accession=ex,library_name=ln,sample_title=d.sample_title.iloc[0],run_count=d.run_accession.nunique(),run_accessions=';'.join(d.run_accession),is_biological_replicate='not_determined_by_accession_uniqueness',is_curated_silene_RNA=bool(d.is_curated_silene_RNA.any())))
    table(pd.DataFrame(lib),z.outdir/'public_library_structure.csv')
    sl=m.loc[m.study_accession.eq('PRJNA229760')]
    slgood=sl.loc[sl.is_curated_silene_RNA]
    assert len(slgood)==18 and slgood.sample_accession.nunique()==9 and slgood.experiment_accession.nunique()==9
    assert slgood.groupby('sample_accession').size().eq(2).all()
    molids=set(ev.loc[ev.species.isin(molecular_species),'inat_taxon_id'].dropna().astype(int))
    stats=dict(schema='fcp_molecular_bridge_summary_v1',roster_species=len(r),environment_species=len(panel),environment_variables=sorted(e.variable.unique()),curated_species=len(annotated),curated_contrasts=len(c),curated_contrast_tissue_rows=len(ev),curated_roster_matched_species=int(core.inat_taxon_id.notna().sum()),curated_roster_unmatched_species=int(core.inat_taxon_id.isna().sum()),curated_environment_matched_species=int(core.environment_panel_match.sum()),curated_molecular_species=len(molecular_species),curated_molecular_roster_species=len(molids),curated_molecular_environment_species=len(molids&panel),strict_same_FCP_individual_mechanism_linked_species=0,species_environment_bridge_rows=len(linked),ena_title_screen_runs=len(screen),ena_title_screen_projects=screen.study_accession.nunique(),ena_title_screen_taxon_labels=screen.scientific_name.nunique(),ena_title_screen_exact_roster_species=ps.loc[ps.roster_exact_match,'inat_taxon_id'].nunique(),ena_title_screen_exact_environment_species=ps.loc[ps.environment_panel_match,'inat_taxon_id'].nunique(),ena_title_screen_unmatched_taxon_labels=ps.loc[~ps.roster_exact_match,'submitted_scientific_name'].nunique(),ena_title_screen_roster_matched_projects=int(ps.roster_exact_match.sum()),silene_reported_accession_total_runs=len(sl),silene_selected_RNA_runs=len(slgood),silene_selected_biosamples=slgood.sample_accession.nunique(),silene_selected_experiments=slgood.experiment_accession.nunique(),silene_excluded_non_Silene_runs=len(sl)-len(slgood),cross_species_mechanism_environment_model_fit=False,reason_no_model='only_four_molecular_annotated_environment_species;heterogeneous_nonrandom_evidence;published_mechanisms_not_confirmed_for_observed_photo_morphs')
    stats={k:int(v) if isinstance(v,(__import__('numpy').integer,)) else v for k,v in stats.items()}
    dump(stats,z.outdir/'summary.json')
    inputs=[z.roster,z.environment,*sorted(z.evidence_dir.glob('*'))]
    inputs=[p for p in inputs if p.is_file()]
    receipt=dict(script=str(Path(__file__).relative_to(ROOT)),python=platform.python_version(),pandas=pd.__version__,matching_rule='NFC-independent ASCII taxon names: trim+collapse whitespace+casefold; exact only; no synonyms or genus inference',inputs={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):digest(p) for p in inputs},outputs={p.name:digest(p) for p in sorted(z.outdir.iterdir()) if p.is_file() and p.name!='execution_receipt.json'},assertions=['4475 unique roster taxon IDs','1574 unique env species x 7 variables','no duplicate contrast x tissue rows','no environment join row multiplication','Silene: 18 RNA runs = 9 samples x 2 technical runs'])
    dump(receipt,z.outdir/'execution_receipt.json')
    print(json.dumps(stats,indent=2))
if __name__=='__main__':main()
