#!/usr/bin/env python3
"""Bounded ENA title expansion: frozen metadata, exact taxon joins, and contrasts.

Offline by default. --fetch writes one live, complete, small-metadata query into
a NEW directory. Follow-up accession metadata never alters title denominators.
No read download, RNA differential expression, causal adaptation, or imputation.
"""
from __future__ import annotations
import argparse, csv as csv_module, datetime, gzip, hashlib, io, json, math, platform, re
from pathlib import Path
import urllib.parse, urllib.request
import xml.etree.ElementTree as ET
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EVID = ROOT/'data/evidence/fcp_molecular_expansion_20261019'
OLD = ROOT/'data/evidence/fcp_molecular_bridge_20261018'
ROSTER = ROOT/'results/fcp_real_447500_white_axis_decomposition_20261018/species_photo_composition.csv'
ENV = ROOT/'results/fcp_real_447500_environment_description_20261018/primary_ge5_species_environment.csv'
OUT = ROOT/'results/fcp_molecular_expansion_20261019'
QUERY = 'tax_tree(33090) AND library_source="TRANSCRIPTOMIC" AND ((study_title="*flower*" OR study_title="*petal*") AND (study_title="*pigment*" OR study_title="*anthocyanin*" OR study_title="*carotenoid*" OR study_title="*color*" OR study_title="*colour*"))'
FIELDS = 'experiment_accession,run_accession,study_accession,secondary_study_accession,scientific_name,tax_id,study_title,sample_accession,sample_title,library_name,library_layout,library_strategy,library_source,instrument_platform'

def key(value):
    return re.sub(r'\s+', ' ', str(value).strip()).casefold()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(obj, path):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def csv(frame, path):
    frame.to_csv(path, index=False, lineterminator='\n', float_format='%.12g')

def read_runs(path):
    frame = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    assert frame.run_accession.is_unique, path
    return frame

def fetch(directory):
    directory.mkdir(parents=True, exist_ok=False)
    url = 'https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(
        dict(result='read_run', query=QUERY, fields=FIELDS, format='tsv', limit=0))
    data = urllib.request.urlopen(url, timeout=90).read(10000001)
    if len(data)>10000000:
        raise RuntimeError('Fixed 10 MB metadata cap exceeded')
    path = directory/'ena_expanded_title_screen_runs.tsv'
    path.write_bytes(data)
    dump(dict(query=QUERY, fields=FIELDS, results=[dict(file=path.name, url=url,
         bytes=len(data), sha256=sha(path), retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())]),
         directory/'retrieval_manifest.json')

def geo_samples(path):
    """Read primary GEO attributes; never decode phenotype from title initials."""
    text = gzip.decompress(path.read_bytes()).decode()
    records = []
    for section in text.split('^SAMPLE = ')[1:]:
        accession = section.splitlines()[0]
        fields = {}
        for line in section.splitlines()[1:]:
            if ' = ' in line:
                k, v = line.split(' = ', 1)
                fields.setdefault(k, []).append(v)
        attrs = dict(v.split(': ', 1) for v in fields.get('!Sample_characteristics_ch1', []) if ': ' in v)
        relations = ';'.join(fields.get('!Sample_relation', []))
        sample = re.search(r'SAM[NED]\d+', relations)
        experiment = re.search(r'[SED]RX\d+', relations)
        records.append(dict(geo_accession=accession, sample_accession=sample.group() if sample else '',
            experiment_accession=experiment.group() if experiment else '',
            sample_title=fields.get('!Sample_title', [''])[0],
            scientific_name=fields.get('!Sample_organism_ch1', [''])[0],
            phenotype=attrs.get('phenotype',''), tissue=attrs.get('tissue',''),
            stage=attrs.get('developmental stage',''),
            biological_replication='GEO series explicitly reports two biological replicates per sample class',
            growth_protocol=';'.join(fields.get('!Sample_growth_protocol_ch1', [])),
            processing_notes=';'.join(fields.get('!Sample_data_processing', []))))
    return pd.DataFrame(records)

def fpkm_descriptives(evidence_dir, geo, outdir):
    """Descriptive Arabidopsis-reference rows; never repair source sample names.

    Ratios of numeric-column means are conditional on submitted color labels.
    Matthiola's two white columns repeat W1: neither is relabelled as W2 or
    declared to be an independent biological replicate. No p-values, count
    model, ortholog validation, or inference of causal alleles is made here.
    """
    gene_map=pd.read_csv(evidence_dir/'pathway_gene_map.csv',dtype=str)
    genes=gene_map.set_index('reference_locus').to_dict('index')
    labels=geo.set_index('sample_title').to_dict('index')
    values, columns, ratios, checks=[],[],[],[]
    for species in ['Matthiola incana','Orychophragmus violaceus','Lobularia maritima']:
        path=evidence_dir/('GSE93879_'+species.replace(' ','_')+'.genes.fpkm_table.xls.gz')
        source=list(csv_module.reader(io.StringIO(gzip.decompress(path.read_bytes()).decode()),delimiter='\t'))
        headers=source[0]
        widths=sorted(set(len(row) for row in source[1:]))
        repeated=sorted({h for h in headers if headers.count(h)>1})
        checks.append(dict(species=species,file=path.name,rows=len(source)-1,header_fields=len(headers),
            numeric_row_fields=';'.join(map(str,widths)),duplicate_headers=';'.join(repeated),
            status='blocked_missing_header_separator' if widths!=[len(headers)] else
                   'conditional_duplicate_white_identity' if repeated else 'exact_sample_header_join'))
        if widths!=[len(headers)]:
            continue  # Lobularia's W1/W2 separator is malformed; no guessing.
        local=[]
        for index,label in enumerate(headers[1:],start=1):
            meta=labels[label]
            ambiguous=label in repeated
            col=dict(species=species,source_file=path.name,source_column_index=index,
                literal_source_header=label,source_header_occurrence=headers[:index+1].count(label),
                phenotype=meta['phenotype'],stage=meta['stage'],
                sample_accession='' if ambiguous else meta['sample_accession'],
                geo_accession='' if ambiguous else meta['geo_accession'],
                identity_status='duplicate_W1_header_cannot_assign_W1_vs_W2' if ambiguous else 'exact_GEO_sample_title_join')
            columns.append(col);local.append(col)
        for row in source[1:]:
            if row[0] not in genes:
                continue
            gene=genes[row[0]]
            nums=[float(v) for v in row[1:]]
            assert all(math.isfinite(v) and v>=0 for v in nums)
            for col,value in zip(local,nums):
                values.append(dict(col,reference_locus=row[0],gene_label=gene['gene_label'],fpkm=value,
                    annotation_status='submitted_Arabidopsis_reference_assignment; focal_species_orthology_not_revalidated'))
            phenotypes=sorted({col['phenotype'] for col in local})
            white='white flower' if species=='Matthiola incana' else 'white-violet flower'
            w=[v for col,v in zip(local,nums) if col['phenotype']==white]
            for color in [v for v in phenotypes if v!=white]:
                c=[v for col,v in zip(local,nums) if col['phenotype']==color]
                wm,cm=sum(w)/len(w),sum(c)/len(c)
                ratio=math.log2(wm/cm) if wm>0 and cm>0 else None
                logs=[math.log2(v/cm) for v in w if v>0 and cm>0]
                ratios.append(dict(species=species,reference_locus=row[0],gene_label=gene['gene_label'],
                    white_label=white,colored_label=color,white_column_values=';'.join(map(str,w)),
                    colored_column_values=';'.join(map(str,c)),white_numeric_column_mean=wm,
                    colored_numeric_column_mean=cm,log2_white_over_colored_numeric_means=ratio,
                    individual_white_column_log2_ratios=';'.join(map(str,logs)),
                    minimum_individual_white_log2_ratio=min(logs) if logs else None,
                    maximum_individual_white_log2_ratio=max(logs) if logs else None,
                    white_numeric_columns=len(w),colored_numeric_columns=len(c),
                    biological_replication='GEO reports n2/class; source numeric white-column identity unresolved' if repeated else 'GEO source reports n2 biological/class; exact source-column join',
                    estimand='descriptive_ratio_of_submitted_FPKM_column_means; no statistical_DE_or_causality'))
    csv(pd.DataFrame(columns),outdir/'processed_expression_column_map.csv')
    csv(pd.DataFrame(values),outdir/'pathway_fpkm_per_source_column.csv')
    csv(pd.DataFrame(ratios),outdir/'pathway_fpkm_label_conditional_ratios.csv')
    csv(pd.DataFrame(checks),outdir/'processed_expression_validation.csv')
    assert len(values)==36 and len(ratios)==12
    return dict(genes=len(genes),species_with_descriptive_values=2,per_column_values=len(values),
        label_conditional_gene_contrasts=len(ratios),count_DE_tests=0,
        source_header_repaired=False,replication_and_annotation_limitations_retained=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-dir', type=Path, default=EVID)
    parser.add_argument('--outdir', type=Path, default=OUT)
    parser.add_argument('--fetch', type=Path, help='new unreviewed metadata snapshot directory')
    args = parser.parse_args()
    if args.fetch:
        fetch(args.fetch)
        return
    ev, out = args.evidence_dir, args.outdir
    out.mkdir(parents=True, exist_ok=True)
    roster, env = pd.read_csv(ROSTER), pd.read_csv(ENV)
    roster['match_key'] = roster.species.map(key)
    assert len(roster)==4475 and roster.match_key.is_unique
    assert env.inat_taxon_id.nunique()==1574
    assert not env.duplicated(['inat_taxon_id','variable']).any()
    lookup = roster.set_index('match_key').inat_taxon_id.to_dict()
    panel = set(env.inat_taxon_id)
    old = read_runs(OLD/'ena_title_screen_runs.tsv')
    new = read_runs(ev/'ena_expanded_title_screen_runs.tsv')
    oldruns, newruns = set(old.run_accession), set(new.run_accession)
    # New rows take precedence for changed metadata; differences are saved below.
    union = pd.concat([new,old],ignore_index=True).drop_duplicates('run_accession')
    changed=[]
    for run in sorted(oldruns & newruns):
        a=old.set_index('run_accession').loc[run]
        b=new.set_index('run_accession').loc[run]
        for field in old.columns.drop('run_accession'):
            if a[field]!=b[field]: changed.append(dict(run_accession=run, field=field, old_value=a[field], new_value=b[field]))
    csv(pd.DataFrame(changed,columns=['run_accession','field','old_value','new_value']),out/'overlapping_run_metadata_changes.csv')
    for frame in [old,new,union]:
        frame['inat_taxon_id']=frame.scientific_name.map(lambda x:lookup.get(key(x)))
        frame['environment_panel_match']=frame.inat_taxon_id.isin(panel)
    union['in_original_title_query']=union.run_accession.isin(oldruns)
    union['in_expanded_title_query']=union.run_accession.isin(newruns)
    csv(union.sort_values('run_accession'),out/'union_run_metadata.csv')
    stats=[]
    for name,frame in [('original',old),('expanded',new),('union',union),('expanded_only',new.loc[~new.run_accession.isin(oldruns)])]:
        stats.append(dict(screen=name,runs=len(frame),projects=frame.study_accession.nunique(),
            submitted_taxon_labels=frame.scientific_name.nunique(),experiments=frame.experiment_accession.nunique(),
            biosamples=frame.sample_accession.nunique(),exact_roster_species=frame.inat_taxon_id.nunique(),
            exact_environment_species=frame.loc[frame.environment_panel_match,'inat_taxon_id'].nunique()))
    csv(pd.DataFrame(stats),out/'screen_denominators.csv')
    decisions=json.loads((ev/'screen_decisions.json').read_text())
    assert set(decisions)==set(new.study_accession), 'Every expanded project needs an explicit disposition'
    olddec=json.loads((OLD/'screen_decisions.json').read_text())
    projects=[]
    for (project,species),frame in union.groupby(['study_accession','scientific_name'],sort=True):
        d=decisions.get(project,olddec.get(project))
        assert d is not None
        tid=lookup.get(key(species))
        projects.append(dict(study_accession=project,submitted_scientific_name=species,
            study_title=frame.study_title.iloc[0],run_count=len(frame),experiment_count=frame.experiment_accession.nunique(),
            biosample_count=frame.sample_accession.nunique(),inat_taxon_id=tid,
            exact_roster_match=tid is not None,environment_panel_match=tid in panel,
            original_query=project in set(old.study_accession),expanded_query=project in set(new.study_accession),
            disposition=d['status'],reason=d['reason'],
            source_url='https://www.ebi.ac.uk/ena/browser/view/'+project))
    csv(pd.DataFrame(projects),out/'all_project_dispositions.csv')
    samples=[]
    for path in sorted(ev.glob('ena_sample_*.xml')):
        tree=ET.parse(path)
        attrs={x.findtext('TAG'):x.findtext('VALUE') for x in tree.findall('.//SAMPLE_ATTRIBUTE')}
        sample=path.stem.removeprefix('ena_sample_')
        frame=new.loc[new.sample_accession.eq(sample)]
        assert len(frame)>0
        samples.append(dict(study_accession=frame.study_accession.iloc[0],scientific_name=frame.scientific_name.iloc[0],
            sample_accession=sample,experiment_accessions=';'.join(sorted(frame.experiment_accession.unique())),
            run_accessions=';'.join(sorted(frame.run_accession.unique())),tissue=attrs.get('tissue',''),
            stage=attrs.get('dev_stage',''),collection_date=attrs.get('collection_date',''),
            phenotype=attrs.get('flower color',attrs.get('phenotype','')),treatment=attrs.get('treatment',''),
            cultivar=attrs.get('cultivar',''),isolate=attrs.get('isolate',''),submitted_label=attrs.get('tmp',''),
            independent_plant_unit='not_confirmed_by_accession_uniqueness',source_file=path.name))
    sampleframe=pd.DataFrame(samples)
    csv(sampleframe,out/'sample_attribute_map.csv')
    drought=sampleframe.loc[sampleframe.study_accession.eq('PRJNA1475932')].copy()
    drought['cell']=drought.phenotype+' | '+drought.treatment
    assert len(drought)==18 and drought.groupby(['phenotype','treatment']).size().eq(3).all()
    assert drought[['phenotype','treatment']].drop_duplicates().shape[0]==6
    drought['replicate_within_cell']=drought.submitted_label.str.extract(r'(\d+)')[0].astype(int).sub(1).mod(3).add(1)
    drought['biological_unit_status']='three_distinct_submitted_samples_per_cell; plant/pool identity and randomization unverified'
    csv(drought.sort_values(['phenotype','treatment','replicate_within_cell']),out/'drought_factorial_sample_sheet.csv')
    geo=geo_samples(ev/'geo_GSE93879_family.soft.gz')
    csv(geo,out/'geo_GSE93879_sample_sheet.csv')
    assert len(geo)==44 and geo.sample_accession.is_unique
    followup=read_runs(ev/'ena_PRJNA883049_followup_runs.tsv')
    followup_rows=[]
    for path in sorted(ev.glob('ena_followup_sample_*.xml')):
        sample=path.stem.removeprefix('ena_followup_sample_')
        tree=ET.parse(path)
        attrs={x.findtext('TAG'):x.findtext('VALUE') for x in tree.findall('.//SAMPLE_ATTRIBUTE')}
        frame=followup.loc[followup.sample_accession.eq(sample)]
        followup_rows.append(dict(sample_accession=sample,run_accessions=';'.join(frame.run_accession),
            experiment_accessions=';'.join(frame.experiment_accession),study_accession='PRJNA883049',
            submitted_replicate_label=attrs.get('replicate',''),tissue=attrs.get('tissue',''),stage=attrs.get('dev_stage',''),
            source_file=path.name,query_scope='direct_accession_followup_outside_title_denominators',
            source_paper_mapping='unresolved: metadata has five color labels and bud stage; main paper three varieties freshly opened flowers'))
    csv(pd.DataFrame(followup_rows),out/'Matthiola_2025_followup_sample_sheet.csv')
    expression_summary=fpkm_descriptives(ev,geo,out)
    manifest=json.loads((ev/'retrieval_manifest.json').read_text())
    for record in manifest['results']:
        if 'sha256' in record: assert sha(ev/record['file'])==record['sha256']
    curated=json.loads((ev/'curated_contrasts.json').read_text())['contrasts']
    oldcurated=json.loads((OLD/'curated_contrasts.json').read_text())['contrasts']
    assert not {v['contrast_id'] for v in curated}&{v['contrast_id'] for v in oldcurated}
    rows=[]
    for origin,records in [('original',oldcurated),('expanded',curated)]:
        for record in records:
            for tissue in record['tissues']:
                row={k:v for k,v in record.items() if k not in ['tissues','tissue_results']}
                row.update(tissue=tissue,curation_snapshot=origin,inat_taxon_id=lookup.get(key(record['species'])))
                row['environment_panel_match']=row['inat_taxon_id'] in panel
                for k,v in row.items():
                    if isinstance(v,(list,dict)):row[k]=json.dumps(v,ensure_ascii=False)
                rows.append(row)
    evidence=pd.DataFrame(rows)
    assert not evidence.duplicated(['contrast_id','tissue']).any()
    csv(evidence,out/'union_contrast_tissue_evidence.csv')
    candidates=pd.read_csv(ev/'reanalysis_contrasts.csv').fillna('')
    candidates['inat_taxon_id']=candidates.species.map(lambda x:lookup.get(key(x)))
    candidates['environment_panel_match']=candidates.inat_taxon_id.isin(panel)
    csv(candidates,out/'reanalysis_contrasts.csv')
    species=evidence.groupby('inat_taxon_id').agg(curated_contrasts=('contrast_id','nunique'),
        molecular_tiers=('molecular_tier',lambda x:';'.join(sorted(set(x)))),
        source_dois=('source_dois',lambda x:';'.join(sorted(set(x))))).reset_index()
    bridge=env.merge(species,on='inat_taxon_id',how='inner',validate='many_to_one')
    bridge['linkage_level']='species_annotation; published material not equated to observed FCP morph'
    csv(bridge,out/'species_environment_evidence_bridge.csv')
    annotation_ids=set(evidence.inat_taxon_id.dropna())
    molecular_ids=set(evidence.loc[evidence.molecular_tier.str.startswith(('T2','T3','T4')),'inat_taxon_id'].dropna())
    roster['environment_panel_match']=roster.inat_taxon_id.isin(panel)
    roster['curated_any_evidence']=roster.inat_taxon_id.isin(annotation_ids)
    roster['curated_molecular_tier_ge2']=roster.inat_taxon_id.isin(molecular_ids)
    roster['evidence_missingness']='unannotated is unknown; fixed title terms and selective primary follow-up'
    csv(roster.drop(columns='match_key'),out/'all_4475_species_coverage.csv')
    summary=dict(schema='fcp_molecular_expansion_v1',query=QUERY,screen_denominators=stats,
        overlapping_runs=len(oldruns&newruns),overlapping_projects=len(set(old.study_accession)&set(new.study_accession)),
        expanded_projects_dispositioned=len(decisions),new_curated_contrasts=len(curated),
        new_curated_species=len({x['species'] for x in curated}),union_curated_contrasts=len(oldcurated)+len(curated),
        union_curated_species=evidence.species.nunique(),union_curated_tissue_rows=len(evidence),
        union_molecular_environment_species=len(molecular_ids&panel),
        new_species_with_unpatterned_white_chromatic_primary_evidence=['Matthiola incana','Centaurea cyanus'],
        new_species_with_partial_white_pattern_primary_evidence=['Antirrhinum majus'],
        drought_factorial_submitted_samples=18,drought_factorial_cells=6,
        drought_independent_plant_replication_verified=False,drought_primary_paper_located=False,
        raw_read_downloads=0,rna_differential_expression_fit=False,processed_expression_descriptive=expression_summary,cross_species_mechanism_environment_fit=False,
        reason_no_cross_species_fit='Seven heterogeneous selectively documented environment species, including partial pattern and commercial contrasts; unsuitable for representative mechanism-frequency or environment-mechanism regression. Species annotation is permitted without same-individual data.')
    dump(summary,out/'summary.json')
    inputs=[ROSTER,ENV,OLD/'curated_contrasts.json',OLD/'ena_title_screen_runs.tsv',OLD/'screen_decisions.json',*sorted(p for p in ev.iterdir() if p.is_file())]
    dump(dict(script=str(Path(__file__).relative_to(ROOT)),script_sha256=sha(Path(__file__)),python=platform.python_version(),
         pandas=pd.__version__,inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs},
         outputs={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='execution_receipt.json'},
         checks=['4475 exact roster names','1574 environment species','all 65 expanded projects dispositioned',
         'overlapping runs deduplicated','39 sample XML records preserved','18 drought samples in 6 cells x 3',
         '44 GEO samples linked by accessions','source snapshot hashes verified','no contrast-tissue duplication',
         'species environment many-to-one join']),out/'execution_receipt.json')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
