#!/usr/bin/env python3
"""Optional read-only refresh of official ENA/Europe PMC metadata, never FASTQ.

Offline reproduction uses the frozen files. Refresh deliberately changes the
search snapshot and hashes; keep the resulting manifest with any updated screen.
"""
from __future__ import annotations
import argparse,concurrent.futures,datetime,hashlib,json
from pathlib import Path
import urllib.parse,urllib.request
ROOT=Path(__file__).resolve().parents[2]
QUERY='tax_tree(33090) AND library_source="TRANSCRIPTOMIC" AND ((study_title="*flower*" AND (study_title="*white*" OR study_title="*polymorph*")) OR (study_title="*petal*" AND study_title="*white*"))'
FIELDS='experiment_accession,run_accession,study_accession,secondary_study_accession,scientific_name,tax_id,study_title,sample_accession,sample_title,library_name,library_layout,library_strategy,library_source,instrument_platform'
def main():
 a=argparse.ArgumentParser(description=__doc__)
 a.add_argument('--outdir',type=Path,required=True,help='new snapshot directory; avoid replacing a reviewed snapshot')
 a.add_argument('--curation',type=Path,default=ROOT/'data/evidence/fcp_molecular_bridge_20261018/curated_contrasts.json')
 z=a.parse_args();z.outdir.mkdir(parents=True,exist_ok=True)
 jobs={'ena_title_screen_runs.tsv':'https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(dict(result='read_run',query=QUERY,fields=FIELDS,format='tsv',limit=0))}
 for acc in ['SRP033277','SRA028419','PRJNA504609','PRJNA1136702','SRR998575','SRR998853']:
  jobs[f'ena_{acc}_runs.tsv']='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode(dict(accession=acc,result='read_run',fields=FIELDS,format='tsv'))
 for doi in sorted({doi for x in json.loads(z.curation.read_text())['contrasts'] for doi in x['source_dois']}):
  jobs['epmc_'+doi.split('/')[-1]+'.json']='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode(dict(query='DOI:'+doi,format='json'))
 def get(item):
  name,url=item
  try:
   data=urllib.request.urlopen(url,timeout=60).read();(z.outdir/name).write_bytes(data)
   return dict(file=name,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
  except Exception as err:return dict(file=name,url=url,error=str(err))
 with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:results=list(ex.map(get,jobs.items()))
 (z.outdir/'retrieval_manifest.json').write_text(json.dumps(dict(query=QUERY,fields=FIELDS,results=results),indent=2)+'\n')
 failures=[r for r in results if 'error' in r]
 print(json.dumps(dict(files=len(results),failures=failures),indent=2))
 if failures:raise SystemExit(1)
if __name__=='__main__':main()
