#!/usr/bin/env node
/**
 * Species-equal descriptive seven-axis environment screen, using actual
 * white-versus-chromatic within-species measurements (no photo reclassification).
 *
 * Run from repo root:
 * node scripts/analysis/audit_fcp_seven_environment_axes_20261020.mjs --root . \
 *   --outdir results/fcp_seven_environment_axes_20261020
 *
 * Outputs genus bootstrap simultaneous max-|t| intervals across seven climate
 * axes. These are post-hoc descriptive summaries, not causal selection tests.
 */
import {readFileSync,writeFileSync,mkdirSync} from "node:fs";
import {join} from "node:path";
import {createHash} from "node:crypto";

const args=process.argv.slice(2),flags={};
for(let i=0;i<args.length;i+=2){if(!args[i]?.startsWith("--"))throw Error("bad CLI flags");flags[args[i].slice(2)]=args[i+1];}
const root=flags.root||".",outdir=flags.outdir||"results/fcp_seven_environment_axes_20261020";
const infile="results/fcp_real_447500_environment_description_20261018/primary_ge5_species_environment.csv";
const source=readFileSync(join(root,infile),"utf8");
const vars=["BIO1","BIO4","BIO5","BIO12","BIO14","BIO15","SRADmean"];
const lines=source.trim().split(/\r?\n/);
const head=lines.shift().split(",");
const records=lines.map((line,i)=>{
  const a=line.split(","); if(a.length!==head.length)throw Error("CSV width issue line "+i);
  return Object.fromEntries(head.map((k,j)=>[k,a[j]]));
});
const perSpecies=new Map();
for(const r of records){
  if(!vars.includes(r.variable))throw Error("unknown climate variable "+r.variable);
  if(r.sd_status!=="estimable")throw Error("unestimable row");
  const id=r.inat_taxon_id, genus=r.species.trim().split(" ")[0];
  if(!perSpecies.has(id))perSpecies.set(id,{
    id,species:r.species,genus,cohort:r.source_cohort,n:Number(r.n_complete_classifiable),
    white:Number(r.n_white),chromatic:Number(r.n_nonwhite),value:{}});
  const s=perSpecies.get(id);
  if(s.species!==r.species||s.n!==Number(r.n_complete_classifiable)||s.cohort!==r.source_cohort)
    throw Error("species-level join inconsistency "+id);
  if(s.value[r.variable]!==undefined)throw Error("duplicate species-variable "+id+" "+r.variable);
  s.value[r.variable]=Number(r.standardized_difference);
  if(!Number.isFinite(s.value[r.variable]))throw Error("nonfinite contrast");
}
const species=[...perSpecies.values()];
if(species.length!==1574||records.length!==1574*7)throw Error("frozen input count drift");
for(const s of species)if(Object.keys(s.value).length!==7||s.white<5||s.chromatic<5)throw Error("incorrect frozen species coverage");
const photos=species.reduce((z,s)=>z+s.n,0);
if(photos!==86114)throw Error("frozen 86114 photograph count drift");
const mean=a=>a.reduce((x,y)=>x+y,0)/a.length;
const point=vars.map(v=>mean(species.map(s=>s.value[v])));
const genera=[...Map.groupBy(species,s=>s.genus).entries()].sort((a,b)=>a[0].localeCompare(b[0]))
 .map(([genus,rows])=>({genus,n:rows.length,sum:vars.map(v=>rows.reduce((z,s)=>z+s.value[v],0))}));
const seed=2026102019;let state=seed;
function rand(){state^=state<<13;state^=state>>>17;state^=state<<5;return(state>>>0)/4294967296;}
const B=5000, draws=vars.map(()=>[]);
const bio5minusBio1=[];
for(let b=0;b<B;b++){
  let sums=Array(7).fill(0),n=0;
  for(let j=0;j<genera.length;j++){
    const g=genera[Math.floor(rand()*genera.length)];n+=g.n;
    for(let k=0;k<7;k++)sums[k]+=g.sum[k];
  }
  for(let k=0;k<7;k++)draws[k].push(sums[k]/n);
  bio5minusBio1.push(sums[2]/n-sums[0]/n);
}
function pct(a,p){const v=[...a].sort((x,y)=>x-y);return v[Math.floor((v.length-1)*p)];}
const se=draws.map(a=>Math.sqrt(a.reduce((s,x)=>s+(x-mean(a))**2,0)/(a.length-1)));
const maxima=[];
for(let b=0;b<B;b++){
  let v=0;
  for(let k=0;k<vars.length;k++)v=Math.max(v,Math.abs((draws[k][b]-point[k])/se[k]));
  maxima.push(v);
}
const familyQuantile=pct(maxima,.95),rows=[];
for(let k=0;k<vars.length;k++){
  const v=vars[k],values=species.map(s=>s.value[v]);
  rows.push({variable:v,n_species:species.length,n_genera:genera.length,white_minus_chromatic_mean_SD:point[k],
    positive_species_fraction:values.filter(x=>x>0).length/values.length,
    genus_boot_ci_low:pct(draws[k],.025),genus_boot_ci_high:pct(draws[k],.975),
    simultaneous_family7_low:point[k]-familyQuantile*se[k],
    simultaneous_family7_high:point[k]+familyQuantile*se[k],
    family7_ci_excludes_zero: (point[k]-familyQuantile*se[k]>0)||(point[k]+familyQuantile*se[k]<0)});
}
const cohortRows=[];
for(const [cohort,subset] of Map.groupBy(species,s=>s.cohort)){
  for(const v of vars)cohortRows.push({cohort,variable:v,n_species:subset.length,mean_difference_SD:mean(subset.map(s=>s.value[v]))});
}
const contrast={BIO5_minus_BIO1_mean_SD:point[2]-point[0],
  genus_boot_CI95:[pct(bio5minusBio1,.025),pct(bio5minusBio1,.975)],
  note:"difference between species-specific within-species SD units; not the physical 1-degree slopes or a multivariate conditional comparison"};
const result={
 schema:"fcp_seven_climate_axis_mean_difference_v1",date:"2026-10-20",
 source:infile,source_sha256:createHash("sha256").update(source).digest("hex"),
 n_source_rows:records.length,n_species:species.length,n_genera:genera.length,
 n_original_photos_on_species:photos,bootstrap_draws:B,bootstrap_seed:seed,
 simultaneous_method:"genus-cluster fixed-effect seven-axis centered max-|t| bootstrap intervals, 95%",
 family95_maxT_quantile:familyQuantile,axis_results:rows,cohort_results:cohortRows,
 descriptive_contrast:contrast,limitations:[
   "Exploratory and conditional on already selected classifiable white and chromatic photographs.",
   "Mean standardized climate differences are univariate descriptive offsets, not Firth joint-model slopes.",
   "Confidence bands resample genus groupings, not photos, and condition on fitted species-level offsets.",
   "Current WorldClim climate normals are not flowering-date weather or historical range reconstructions.",
   "BIO4 climatological temperature seasonality is not year-specific temporal fluctuations.",
   "Species selection, image exposure, observation conditions, geography and habitat may confound recorded white.",
   "Positive association does not demonstrate pigment production costs, fitness tradeoffs, adaptation or phylogeographic history."
 ]};
mkdirSync(join(root,outdir),{recursive:true});
function csv(arr){
  const h=Object.keys(arr[0]);return h.join(",")+"\n"+arr.map(r=>h.map(k=>String(r[k])).join(",")).join("\n")+"\n";
}
writeFileSync(join(root,outdir,"summary.json"),JSON.stringify(result,null,2)+"\n");
writeFileSync(join(root,outdir,"axis_summary.csv"),csv(rows));
writeFileSync(join(root,outdir,"cohort_summary.csv"),csv(cohortRows));
console.log(JSON.stringify({photos,n_species:species.length,maxT:familyQuantile,axis_results:rows,contrast}));
