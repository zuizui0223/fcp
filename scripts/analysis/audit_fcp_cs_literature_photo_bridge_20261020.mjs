#!/usr/bin/env node
/**
 * Strict exact-name bridge between literature-documented C/S flower-colour
 * organisation and frozen FCP photographed-colour results.
 * No evidence is inferred from missing C/S labels. Photo-only colour states do
 * not establish inherited local polymorphism; this script counts direct coverage.
 *
 * Usage:
 * node scripts/analysis/audit_fcp_cs_literature_photo_bridge_20261020.mjs \
 *   --root . --outdir results/fcp_cs_literature_photo_bridge_20261020
 */
import {readFileSync, writeFileSync, mkdirSync} from "node:fs";
import {join} from "node:path";
import {createHash} from "node:crypto";

const argv=process.argv.slice(2);
const flags={};
for(let i=0;i<argv.length;i+=2){if(!argv[i]?.startsWith("--"))throw Error("Flag expected");flags[argv[i].slice(2)]=argv[i+1];}
const root=flags.root||".";
const outdir=flags.outdir||"results/fcp_cs_literature_photo_bridge_20261020";
const files={
  C_S_LITERATURE:"data/frozen/jbi_cs_evidence_freeze_v22.csv",
  PHOTOS:"results/fcp_real_447500_white_axis_decomposition_20261018/species_photo_composition.csv",
  SPATIAL:"results/fcp_real_447500_white_axis_decomposition_20261018/species_spatial_decomposition.csv",
  ENVIRONMENT:"results/fcp_real_447500_environment_model_20261018/species_slopes.csv"
};

function parseCsv(source) {
  const rows=[]; let field="",row=[],inQuotes=false;
  for(let i=0;i<source.length;i++) {
    const c=source[i];
    if(inQuotes){
      if(c==='"'&&source[i+1]==='"'){field+='"';i++;}
      else if(c==='"')inQuotes=false;
      else field+=c;
    }else if(c==='"'){inQuotes=true;}
    else if(c===","){row.push(field);field="";}
    else if(c==="\n"){row.push(field.replace(/\r$/,""));rows.push(row);row=[];field="";}
    else field+=c;
  }
  if(inQuotes)throw Error("Unclosed quote in CSV");
  if(field||row.length){row.push(field.replace(/\r$/,""));rows.push(row);}
  const head=rows.shift();
  return rows.filter(x=>x.length!==1||x[0]!=="").map((values,i)=>{
    if(values.length!==head.length)throw Error("CSV columns disagree on data row "+(i+2));
    return Object.fromEntries(head.map((k,j)=>[k,values[j]]));
  });
}
const raw={},tbl={};
for(const [key,p] of Object.entries(files)){
  raw[key]=readFileSync(join(root,p),"utf8");
  tbl[key]=parseCsv(raw[key]);
}
const expected={C_S_LITERATURE:34,PHOTOS:4475,SPATIAL:1761};
for(const [key,n] of Object.entries(expected))
  if(tbl[key].length!==n)throw Error("Incorrect frozen row count "+key+": "+tbl[key].length+" != "+n);
const norm=x=>String(x||"").normalize("NFKC").replace(/_/g," ").trim().replace(/\s+/g," ").toLowerCase();
const index=(rows,k)=>{
  const out=new Map();
  for(const r of rows){const name=norm(r[k]);if(out.has(name))throw Error("Duplicate species "+name);out.set(name,r);}
  return out;
};
const ph=index(tbl.PHOTOS,"species"),sp=index(tbl.SPATIAL,"species");
const env=tbl.ENVIRONMENT.filter(x=>x.model==="primary"&&x.variable==="BIO5"&&x.status==="ok");
const en=index(env,"species");
const entries=tbl.C_S_LITERATURE.map(x=>{
  const name=norm(x.canonical_name),p=ph.get(name),s=sp.get(name),e=en.get(name);
  return {
    scientific_name:x.canonical_name,
    literature_state:x.organization_state,
    literature_C_positive:Number(x.C_local_coexistence_documented),
    literature_S_positive:Number(x.S_spatial_segregation_documented),
    literature_sources:Number(x.n_resolved_sources),
    photo_roster_match:Number(!!p),
    n_classifiable:p?Number(p.n_classifiable):"",
    n_white:p?Number(p.n_white):"",
    n_chromatic:p?Number(p.n_nonwhite):"",
    n_photo_states:p?Number(p.n_states):"",
    spatial_50km_match:Number(!!s),
    spatial_depletion:s?Number(s.depletion_total):"",
    n_local_50km_pairs:s?Number(s.local_pairs_50km):"",
    BIO5_model_match:Number(!!e),
    BIO5_slope_per_species_SD:e?Number(e.slope):"",
    photo_color_is_wild_morph_validated:0,
    publication_C_S_is_not_negative_when_unreported:1
  };
});
const classes=["local_coexistence_only","spatial_segregation_only","coexistence_and_segregation"];
const summary=classes.map(k=>{
  const rows=entries.filter(x=>x.literature_state===k);
  return {
    literature_state:k,n_literature_species:rows.length,
    n_photo_roster:rows.filter(x=>x.photo_roster_match===1).length,
    n_spatial_50km:rows.filter(x=>x.spatial_50km_match===1).length,
    n_BIO5:rows.filter(x=>x.BIO5_model_match===1).length,
    n_BIO5_and_spatial:rows.filter(x=>x.BIO5_model_match===1&&x.spatial_50km_match===1).length
  };
});
const sums=keys=>summary.reduce((total,row)=>total+row[keys],0);
if(sums("n_literature_species")!==34||sums("n_photo_roster")!==16||
   sums("n_spatial_50km")!==5||sums("n_BIO5")!==11)
  throw Error("Historical exact-name overlap counts no longer match reviewed freeze");
const sourceHash=Object.fromEntries(Object.entries(raw).map(([k,v])=>[files[k],createHash("sha256").update(v).digest("hex")]));
const output={
  schema:"fcp_literature_CS_photo_bridge_v1",date:"2026-10-20",
  matching:"exact NFKC normalized accepted name only; no guessed synonyms",
  study_units:{literature_n:34,photo_roster_n:4475,photo_spatial_50km_n:1761,
    fcp_environment_BIO5_n:env.length},
  summary,totals:{
    literature_species:34,photo_roster:16,spatial_50km:5,BIO5_models:11,
    BIO5_and_spatial:entries.filter(r=>r.spatial_50km_match&&r.BIO5_model_match).length,
    spatial_segregation_only_spatial_matches:summary.find(r=>r.literature_state==="spatial_segregation_only").n_spatial_50km
  },source_sha256:sourceHash,
  inference:"no literature-C/S versus photographic-spatial response model is estimable: no S-only species in the strict five-species spatial intersection",
  warnings:[
    "C/S literature evidence is positive-only; zero means not documented, never absence.",
    "50-km photographed mismatch is not a within-population C/S indicator.",
    "Photo labels have no linked individual-genotype or independently audited morph identity.",
    "The 74-species exploratory JBI core is a different selection from this 34-species strict evidence freeze.",
    "The FCP sample has 100 originally measured photo rows per species and 227107 classifiable overall; cannot represent global genetically validated FCP prevalence.",
    "Exact-name join can miss legitimate synonyms; do not silently substitute or inflate overlap."
  ]
};
function quote(v){let s=String(v==null?"":v);return /[",\r\n]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s;}
function toCsv(records){
  const h=Object.keys(records[0]);
  return h.join(",")+"\n"+records.map(r=>h.map(k=>quote(r[k])).join(",")).join("\n")+"\n";
}
mkdirSync(join(root,outdir),{recursive:true});
writeFileSync(join(root,outdir,"summary.json"),JSON.stringify(output,null,2)+"\n","utf8");
writeFileSync(join(root,outdir,"strict34_species_crosswalk.csv"),toCsv(entries),"utf8");
writeFileSync(join(root,outdir,"class_coverage.csv"),toCsv(summary),"utf8");
console.log(JSON.stringify(output.totals));
