#!/usr/bin/env python3
"""Bounded public-file retrieval; no login, challenge solving, or invented rows."""
import argparse,csv,hashlib,io,json
from pathlib import Path
from urllib.request import Request,urlopen

FILES={'README.md':4398699,'Infection_between_sites.csv':4398689,
       'Within_stand_density.csv':4398690,'Fitness_data.csv':4398691}


def profile_bytes(name,body):
    text=body.decode('utf-8-sig')
    start=text.lstrip()[:500].lower()
    if '<html' in start or '<!doctype' in start or '<title>validating' in start:
        raise ValueError('HTML/challenge page is not a data file')
    if name.endswith('.md'):
        if not text.strip():raise ValueError('empty README')
        return {'kind':'README','text':text,'census_semantics_verified':False}
    rows=list(csv.reader(io.StringIO(text)))
    if len(rows)<2 or len(rows[0])<2 or len(set(rows[0]))!=len(rows[0]):
        raise ValueError('CSV needs unique headers and data rows')
    nonempty=[r for r in rows[1:] if r]
    if any(len(r)!=len(rows[0]) for r in nonempty):raise ValueError('ragged CSV')
    return {'kind':'CSV','headers':rows[0],'rows':len(nonempty),
            'first_three_rows':nonempty[:3], 'census_semantics_verified':False}


def recover(outdir):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    entries=[]
    for name,fileid in FILES.items():
        url='https://datadryad.org/downloads/file_stream/'+str(fileid)
        row={'name':name,'url':url,'status':'NOT_ATTEMPTED'}
        try:
            request=Request(url,headers={'User-Agent':'FCP-reproducible-research/1.0'})
            with urlopen(request,timeout=20) as response:
                body=response.read(2_000_001)
                row.update(http_status=response.status,final_url=response.url,
                           content_type=response.headers.get('Content-Type'))
            if len(body)>2_000_000:raise ValueError('public-file size cap exceeded')
            row.update(bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
            profile=profile_bytes(name,body)
            (outdir/name).write_bytes(body)
            row.update(status='BYTES_RECOVERED_SEMANTICS_UNVERIFIED',profile=profile)
        except Exception as exc:row.update(status='NOT_RECOVERED',error=str(exc),error_type=type(exc).__name__)
        entries.append(row); print(json.dumps(row),flush=True)
    result={'schema':'fcp_hesperis_public_recovery_v1','dataset_doi':'10.5061/dryad.5x69p8dj1',
        'dataset_role':'companion viral-infection study, not yet verified as target flower-colour census',
        'files':entries,'n_byte_recovered':sum(r['status'].startswith('BYTES') for r in entries),
        'census_semantics_verified':False,'new_biological_results_computed':False}
    (outdir/'recovery_status.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    recover(p.parse_args().outdir)
