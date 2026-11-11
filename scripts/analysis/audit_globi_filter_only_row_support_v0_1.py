#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import duckdb

DOIS=[
    "10.1126/science.257.5073.1107",
    "10.1086/674445",
    "10.1007/s10641-010-9606-0",
]
MANGAL_NS="globalbioticinteractions/mangal"

def sha_lines(values:list[str])->str:
    h=hashlib.sha256()
    for x in sorted(values):
        h.update(x.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--url",required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()

    con=duckdb.connect(database=":memory:")
    con.execute("INSTALL httpfs")
    con.execute("LOAD httpfs")
    con.execute("INSTALL parquet")
    con.execute("LOAD parquet")
    con.execute("SET threads=4")
    con.execute("SET memory_limit='6GB'")

    # Flatten the frozen DOI × field family deterministically.
    clauses=[]
    for d in DOIS:
        for field in ("referenceDoi","sourceDOI","referenceCitation","sourceCitation"):
            clauses.append(f"lower(coalesce(\"{field}\",'')) LIKE '%{d.lower()}%'")
    doi_pred=" OR ".join(clauses)
    ns_pred=f"lower(trim(coalesce(\"sourceNamespace\",''))) = '{MANGAL_NS}'"

    scan=f"""
      SELECT
        "sourceNamespace","referenceDoi","sourceDOI","referenceCitation","sourceCitation",
        "sourceTaxonSpeciesName","sourceTaxonSpeciesId",
        "targetTaxonSpeciesName","targetTaxonSpeciesId","interactionTypeName"
      FROM read_parquet(?)
    """
    escaped_url=a.url.replace("'","''")
    con.execute("CREATE OR REPLACE TEMP VIEW projected AS "+scan.replace("read_parquet(?)",f"read_parquet('{escaped_url}')"))
    con.execute(f"""
      CREATE OR REPLACE TEMP VIEW flagged AS
      SELECT *,
        ({ns_pred}) AS excluded_namespace,
        ({doi_pred}) AS excluded_doi,
        (({ns_pred}) OR ({doi_pred})) AS excluded_union,
        coalesce(nullif(trim("sourceTaxonSpeciesId"),''), nullif(trim("sourceTaxonSpeciesName"),'')) AS focal_id,
        coalesce(nullif(trim("targetTaxonSpeciesId"),''), nullif(trim("targetTaxonSpeciesName"),'')) AS partner_id,
        coalesce(nullif(trim("interactionTypeName"),''),'__MISSING__') AS interaction_type
      FROM projected
    """)

    total=con.execute("""
      SELECT
        count(*) AS total_rows,
        count(*) FILTER (WHERE excluded_namespace) AS excluded_namespace_rows,
        count(*) FILTER (WHERE excluded_doi) AS excluded_doi_rows,
        count(*) FILTER (WHERE excluded_namespace AND excluded_doi) AS excluded_both_rows,
        count(*) FILTER (WHERE excluded_union) AS excluded_union_rows,
        count(*) FILTER (WHERE NOT excluded_union) AS remaining_rows,
        count(DISTINCT focal_id) FILTER (WHERE NOT excluded_union AND focal_id IS NOT NULL) AS distinct_focal_taxa,
        count(DISTINCT partner_id) FILTER (WHERE NOT excluded_union AND partner_id IS NOT NULL) AS distinct_partner_taxa
      FROM flagged
    """).fetchone()

    groups=con.execute("""
      SELECT interaction_type, focal_id,
             count(*) AS n_records,
             count(DISTINCT partner_id) AS n_partners
      FROM flagged
      WHERE NOT excluded_union
        AND focal_id IS NOT NULL
        AND partner_id IS NOT NULL
      GROUP BY interaction_type, focal_id
      HAVING count(*) >= 10 AND count(DISTINCT partner_id) >= 2
    """).fetchall()

    by_type={}
    eligible_ids=set()
    for itype,focal,nrec,npartners in groups:
        by_type[str(itype)]=by_type.get(str(itype),0)+1
        eligible_ids.add(str(focal))

    status=(
      "FILTER_ONLY_SUPPORT_PASS_SPECIALIZATION_VALUE_STILL_SEALED"
      if len(eligible_ids)>=20 else
      "HOLD_INSUFFICIENT_INDEPENDENT_GLOBI_PARTNER_SUPPORT"
    )

    result={
      "version":"v0.1",
      "status":status,
      "source":{
        "zenodo_record":22691479,
        "version":"0.11",
        "file":"interactions.parquet",
        "bytes":2375721521,
        "md5":"8212776b7ba0c4fa8de26a00b5ce8595",
        "url":a.url
      },
      "exact_projection":[
        "sourceNamespace","referenceDoi","sourceDOI","referenceCitation","sourceCitation",
        "sourceTaxonSpeciesName","sourceTaxonSpeciesId",
        "targetTaxonSpeciesName","targetTaxonSpeciesId","interactionTypeName"
      ],
      "exclusion_rule":{
        "namespace_exact_case_insensitive":MANGAL_NS,
        "response_source_dois":DOIS
      },
      "counts":{
        "total_rows":int(total[0]),
        "excluded_namespace_rows":int(total[1]),
        "excluded_doi_rows":int(total[2]),
        "excluded_both_rows":int(total[3]),
        "excluded_union_rows":int(total[4]),
        "remaining_rows":int(total[5]),
        "distinct_focal_taxa_remaining":int(total[6]),
        "distinct_partner_taxa_remaining":int(total[7]),
        "eligible_focal_taxa_total":len(eligible_ids),
        "eligible_focal_taxa_by_interaction_type":dict(sorted(by_type.items()))
      },
      "eligible_focal_identity_set_sha256":sha_lines(list(eligible_ids)),
      "eligibility_rule":{
        "minimum_records":10,
        "minimum_distinct_partners":2,
        "minimum_eligible_focal_taxa_total":20
      },
      "persisted_focal_taxon_names":False,
      "persisted_partner_taxon_names":False,
      "latitude_opened":False,
      "longitude_opened":False,
      "event_date_opened":False,
      "specialization_values_computed":False,
      "partner_entropy_computed":False,
      "turnover_computed":False,
      "next_gate":(
        "Freeze focal-panel crosswalk/support against the future response panel before any partner entropy is computed."
        if status.startswith("FILTER_ONLY_SUPPORT_PASS") else
        "STOP_SPECIALIZATION_PREDICTOR_SUPPORT"
      )
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
