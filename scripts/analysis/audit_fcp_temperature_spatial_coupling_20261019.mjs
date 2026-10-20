#!/usr/bin/env node
/**
 * Exploratory audit: does species-level temperature-white association magnitude
 * track photographic geographic white/chromatic mismatch depletion?
 * No new Firth fitting, no causal interpretation, no joint photo-label null.
 *
 * Run: node scripts/analysis/audit_fcp_temperature_spatial_coupling_20261019.mjs [repo-root]
 */
import {readFileSync,writeFileSync,mkdirSync} from "node:fs";
import {join} from "node:path";
import {createHash} from "node:crypto";
const root=process.argv[2]??".";
const pSlope=join(root,"results/fcp_real_447500_environment_model_20261018/species_slopes.csv");
const pSpatial=join(root,"results/fcp_real_447500_white_axis_decomposition_20261018/species_spatial_decomposition.csv");
const outdir=join(root,"results/fcp_temperature_spatial_coupling_20261019");
function csv(buf){
  const [head,...lines]=buf.toString("utf8").trim().split(/\r?\n/),keys=head.split(",");
  return lines.map((l,j)=>{
    const v=l.split(",");
    if(v.length!==keys.length)throw Error("Malformed CSV line "+(j+2));
    return Object.fromEntries(keys.map((k,i)=>[k,v[i]]));
  });
}
function sha(buf){return createHash("sha256").update(buf).digest("hex");}
const slopeRaw=readFileSync(pSlope),spatialRaw=readFileSync(pSpatial);
const slopes=csv(slopeRaw).filter(r=>r.model==="primary"&&r.variable==="BIO5"&&r.status==="ok");
const spatial=csv(spatialRaw);
if(slopes.length!==1574||spatial.length!==1761)throw Error("Frozen row counts changed");
const byId=new Map(spatial.map(r=>[r.inat_taxon_id,r]));
const rows=slopes.filter(r=>byId.has(r.inat_taxon_id)).map(r=>{
  const s=byId.get(r.inat_taxon_id);
  const x={
    taxon_id:r.inat_taxon_id,species:r.species,genus:r.species.split(" ")[0],cohort:r.source_cohort,
    beta:Number(r.slope), se:Number(r.se_approx), n:Number(r.n_retained),
    nWhite:Number(r.n_white),nNon:Number(r.n_nonwhite),
    tempSD:Number(r.primary_within_species_sd_native),
    localPairs:Number(s.local_pairs_50km),globalWhite:Number(s.global_white_nonwhite)
  };
  x.whiteFrac=x.nWhite/x.n;
  x.absBeta=Math.abs(x.beta);
  x.normSpatial=(Number(s.depletion_white_nonwhite)-Number(s.null_mean_depletion_white_nonwhite))/x.globalWhite;
  if(!Number.isFinite(x.normSpatial)||!Number.isFinite(x.absBeta))throw Error("Non-finite joined data");
  return x;
});
if(rows.length!==787||new Set(rows.map(r=>r.taxon_id)).size!==rows.length)throw Error("Frozen joint eligibility changed");
function ranks(a){
  const p=a.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
  const r=new Array(a.length);
  for(let i=0;i<p.length;){
    let j=i+1;
    while(j<p.length&&p[j].v===p[i].v)j++;
    const v=(i+j+1)/2;
    for(let k=i;k<j;k++)r[p[k].i]=v;
    i=j;
  }
  return r;
}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;}
function Pearson(a,b){
  const ma=a.reduce((s,x)=>s+x,0)/a.length,mb=b.reduce((s,x)=>s+x,0)/b.length;
  const x=a.map(v=>v-ma),y=b.map(v=>v-mb);
  return dot(x,y)/Math.sqrt(dot(x,x)*dot(y,y));
}
function Spearman(rr,a,b){return Pearson(ranks(rr.map(r=>r[a])),ranks(rr.map(r=>r[b])));}
const COHORTS=[...new Set(rows.map(r=>r.cohort))].sort();
function partialSpearman(rr,xName,yName){
  const n=rr.length;
  let design=[Array(n).fill(1)];
  for(const k of ["se","whiteFrac","n","tempSD","localPairs","globalWhite"])design.push(ranks(rr.map(r=>r[k])));
  for(const cohort of COHORTS.slice(1))design.push(rr.map(r=>r.cohort===cohort?1:0));
  const Q=[];
  for(const col of design){
    const v=[...col];
    for(const q of Q){
      const scalar=dot(v,q);
      for(let i=0;i<n;i++)v[i]-=scalar*q[i];
    }
    const norm=Math.sqrt(dot(v,v));
    if(norm>1e-9)Q.push(v.map(x=>x/norm));
  }
  const residual=k=>{
    const v=ranks(rr.map(r=>r[k]));
    for(const q of Q){
      const scalar=dot(v,q);
      for(let i=0;i<n;i++)v[i]-=scalar*q[i];
    }
    return v;
  };
  const x=residual(xName),y=residual(yName);
  return dot(x,y)/Math.sqrt(dot(x,x)*dot(y,y));
}
const baseline={
  n_species:rows.length,n_genera:new Set(rows.map(r=>r.genus)).size,
  raw_abs_temperature_slope_vs_spatial:Spearman(rows,"absBeta","normSpatial"),
  raw_signed_temperature_slope_vs_spatial:Spearman(rows,"beta","normSpatial"),
  covariate_adjusted_abs_vs_spatial:partialSpearman(rows,"absBeta","normSpatial"),
  covariate_adjusted_signed_vs_spatial:partialSpearman(rows,"beta","normSpatial"),
};
const higher=rows.filter(r=>r.n>=60&&r.nWhite>=10&&r.nNon>=10&&r.localPairs>=100);
const highSupport={n_species:higher.length,n_genera:new Set(higher.map(r=>r.genus)).size,
  adjusted_abs_vs_spatial:partialSpearman(higher,"absBeta","normSpatial"),
  adjusted_signed_vs_spatial:partialSpearman(higher,"beta","normSpatial")
};
let state=2026101921;
function rand(){state^=state<<13;state^=state>>>17;state^=state<<5;return(state>>>0)/4294967296;}
const grouped=new Map();
for(const r of rows){if(!grouped.has(r.genus))grouped.set(r.genus,[]);grouped.get(r.genus).push(r);}
const genera=[...grouped.values()],draws=[];
for(let b=0;b<500;b++){
  const sample=[];
  for(let i=0;i<genera.length;i++)sample.push(...genera[Math.floor(rand()*genera.length)]);
  draws.push([partialSpearman(sample,"absBeta","normSpatial"),partialSpearman(sample,"beta","normSpatial")]);
}
function percentile(x,pr){const z=[...x].sort((a,b)=>a-b);return z[Math.floor(z.length*pr)];}
baseline.genus_bootstrap_95_abs=[percentile(draws.map(x=>x[0]),.025),percentile(draws.map(x=>x[0]),.975)];
baseline.genus_bootstrap_95_signed=[percentile(draws.map(x=>x[1]),.025),percentile(draws.map(x=>x[1]),.975)];
const byCohort=[...new Set(rows.map(r=>r.cohort))].map(cohort=>{
  const rr=rows.filter(r=>r.cohort===cohort);
  return{cohort,n_species:rr.length,spearman_abs:Spearman(rr,"absBeta","normSpatial"),
    spearman_signed:Spearman(rr,"beta","normSpatial")};
});
const result={
  schema:"fcp_temperature_vs_geographic_colour_sorting_exploratory_v1",
  source_sha256:{species_slopes:sha(slopeRaw),spatial_decomposition:sha(spatialRaw)},
  baseline,high_support:highSupport,by_cohort:byCohort,
  controls:["first-stage slope SE","photographed white fraction","photo count",
            "sampled BIO5 climate SD","50-km photo-pair count","global white mismatch","source cohort"],
  bootstrap:{unit:"genus",n_draws:500,seed:2026101921,conditional_on_fitted_slopes:true},
  limitations:[
    "Both slope and geographic sorting are computed from the same photographic colour labels; mathematically coupled under geographic climate autocorrelation.",
    "No matched null jointly re-fitting slopes and spatial effects after within-species label permutation was executed.",
    "Species inclusion is the overlap of two previously selected analyses, not a random sample.",
    "Partial rank correlation is descriptive, not a mechanism adjustment or causal selection coefficient.",
    "The sign of temperature slopes is not systematically associated with spatial sorting strength.",
    "Distinct biological colour morph frequencies, genetic isolation and pigment costs are unmeasured."
  ]
};
mkdirSync(outdir,{recursive:true});
writeFileSync(join(outdir,"result.json"),JSON.stringify(result,null,2)+"\n");
const fields=Object.keys(byCohort[0]);
writeFileSync(join(outdir,"cohort_correlations.csv"),
  [fields.join(","),...byCohort.map(r=>fields.map(k=>r[k]).join(","))].join("\n")+"\n");
console.log(JSON.stringify({baseline,highSupport}));
