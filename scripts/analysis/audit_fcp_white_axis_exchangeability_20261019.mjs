#!/usr/bin/env node
/**
 * Reproducible post-hoc audit of four-state white-label frequency and
 * white/chromatic vs chromatic/chromatic spatial depletion.
 * Usage:
 *   node scripts/analysis/audit_fcp_white_axis_exchangeability_20261019.mjs [repoRoot]
 *
 * No reclassification or model fitting. Four colour labels are NOT assumed
 * ecologically equivalent; the label-exchangeability model is a diagnostic,
 * not a biological evolutionary null.
 */
import {readFileSync, writeFileSync, mkdirSync} from "node:fs";
import {join} from "node:path";
import {createHash} from "node:crypto";

const root=process.argv[2] ?? ".";
const photoFile=join(root,"results/fcp_real_447500_white_axis_decomposition_20261018/species_photo_composition.csv");
const spatialFile=join(root,"results/fcp_real_447500_white_axis_decomposition_20261018/species_spatial_decomposition.csv");
const out=join(root,"results/fcp_white_axis_exchangeability_20261019");
const nBootstrap=5000;
const rngSeed=20261019;

function csv(source) {
  const lines=source.trim().split(/\r?\n/);
  const head=lines.shift().split(",");
  return lines.map((line,j)=>{
    const cells=line.split(",");
    if(cells.length!==head.length)throw Error("CSV width mismatch, line "+(j+2));
    const r=Object.fromEntries(head.map((k,i)=>[k,cells[i]]));
    for(const k of head)
      if(!["inat_taxon_id","species","source_cohort"].includes(k))
        r[k]=r[k]===""?null:Number(r[k]);
    r.genus=r.species.split(" ")[0];
    return r;
  });
}
const sourcePhoto=readFileSync(photoFile);
const sourceSpatial=readFileSync(spatialFile);
const photos=csv(sourcePhoto.toString("utf8"));
const spatial=csv(sourceSpatial.toString("utf8"));
if(process.env.FCP_ALLOW_FIXTURE!=="1"&&(photos.length!==4475||spatial.length!==1761))
  throw Error("frozen-data species count contract changed");
const sha=buf=>createHash("sha256").update(buf).digest("hex");
function generator(seed) {
  let v=seed|0;
  return ()=>{v^=v<<13;v^=v>>>17;v^=v<<5;return (v>>>0)/4294967296;};
}
const random=generator(rngSeed);
function percentile(values,p) {
  const v=[...values].sort((a,b)=>a-b);
  return v[Math.floor((v.length-1)*p)];
}
function bootCI(rows, effect) {
  const grouped=Object.values(Object.groupBy(rows,r=>r.genus))
    .map(arr=>({n:arr.length,s:arr.reduce((a,r)=>a+effect(r),0)}));
  const draws=[];
  for(let b=0;b<nBootstrap;b++){
    let n=0,sum=0;
    for(let j=0;j<grouped.length;j++){
      const g=grouped[Math.floor(random()*grouped.length)];
      n+=g.n;sum+=g.s;
    }
    draws.push(sum/n);
  }
  return [percentile(draws,.025),percentile(draws,.975)];
}
const mean=(rows,f)=>rows.reduce((sum,r)=>sum+f(r),0)/rows.length;
const eligible=photos.filter(r=>r.n_classifiable>=40&&r.n_states>=2);
if(process.env.FCP_ALLOW_FIXTURE!=="1"&&eligible.length!==3029)
  throw Error("colour-eligible species count drift");
const semanticWhite=r=>(r.n_white>0?1:0);
const presenceExcess=r=>semanticWhite(r)-r.n_states/4;
const filters={
  all:eligible,
  classifiable_ge60:eligible.filter(r=>r.n_classifiable>=60),
  classifiable_ge80:eligible.filter(r=>r.n_classifiable>=80),
  two_or_three_states:eligible.filter(r=>r.n_states<=3),
};
const presence=Object.entries(filters).map(([name,rows])=>{
  const ci=bootCI(rows,presenceExcess);
  return {
    group:name,n_species:rows.length,n_genera:new Set(rows.map(r=>r.genus)).size,
    observed_white_species:rows.filter(semanticWhite).length,
    observed_fraction:mean(rows,semanticWhite),
    exchangeable_label_expected_fraction:mean(rows,r=>r.n_states/4),
    excess_over_label_exchangeability:mean(rows,presenceExcess),
    genus_bootstrap_low:ci[0],genus_bootstrap_high:ci[1]
  };
});
const byStates=[2,3,4].map(k=>{
  const arr=eligible.filter(r=>r.n_states===k);
  return {n_states:k,n_species:arr.length,n_white:arr.filter(semanticWhite).length,
    observed_fraction:mean(arr,semanticWhite),label_expected_fraction:k/4};
});
const cohort=Object.entries(Object.groupBy(eligible,r=>r.source_cohort))
  .map(([source_cohort,arr])=>({source_cohort,n_species:arr.length,
    white_species:arr.filter(semanticWhite).length,
    observed_fraction:mean(arr,semanticWhite),
    expected_exchangeable:mean(arr,r=>r.n_states/4)}));
function spatialDiff(r){
  return (r.depletion_white_nonwhite-r.null_mean_depletion_white_nonwhite)/r.global_white_nonwhite
       -(r.depletion_nonwhite_nonwhite-r.null_mean_depletion_nonwhite_nonwhite)/r.global_nonwhite_nonwhite;
}
function whiteNorm(r){return (r.depletion_white_nonwhite-r.null_mean_depletion_white_nonwhite)/r.global_white_nonwhite;}
function hueNorm(r){return (r.depletion_nonwhite_nonwhite-r.null_mean_depletion_nonwhite_nonwhite)/r.global_nonwhite_nonwhite;}
const thresholds=[0,.005,.01,.02,.05];
const normalized=thresholds.map(min=>{
  const arr=spatial.filter(r=>r.global_white_nonwhite>min&&r.global_nonwhite_nonwhite>min);
  const ci=bootCI(arr,spatialDiff);
  return {min_global_mismatch: min,n_species:arr.length,
    n_genera:new Set(arr.map(r=>r.genus)).size,
    normalized_white_mean:mean(arr,whiteNorm),
    normalized_nonwhite_mean:mean(arr,hueNorm),
    difference_white_minus_nonwhite:mean(arr,spatialDiff),
    fraction_white_relatively_stronger:mean(arr,r=>Number(spatialDiff(r)>0)),
    genus_bootstrap_95_low:ci[0],genus_bootstrap_95_high:ci[1]};
});
const output={
  schema:"fcp_white_axis_exchangeability_and_normalized_spatial_v1",
  date:"2026-10-19",model:"post-hoc descriptive, not an evolutionary or causal null",
  sources:{photo_sha256:sha(sourcePhoto),spatial_sha256:sha(sourceSpatial),
    photo_rows:photos.length,spatial_rows:spatial.length},
  n_bootstrap:nBootstrap,rng_seed:rngSeed,
  presence,by_states:byStates,by_cohort:cohort,spatial_normalized:normalized,
  limitations:[
    "White is 1 of 4 photo-label classes; label exchangeability is hypothetical and not physically justified.",
    "Four-state species include white mechanically, so 2-3 state subsets are reported.",
    "No independent local-population genotype or fitness information.",
    "Normalized spatial comparison conditions on both baseline components exceeding a chosen threshold.",
    "199-permutation component null means treated as fixed; resampling genus does not propagate uncertainty in those means.",
    "Bootstrap resamples named genera and is not a phylogenetic covariance model.",
    "No held-out confirmation: all estimates are post-hoc exploratory."
  ]
};
mkdirSync(out,{recursive:true});
writeFileSync(join(out,"result.json"),JSON.stringify(output,null,2)+"\n");
function outputCSV(path,rs){
  if(!rs.length)return;
  const head=Object.keys(rs[0]);
  const text=[head.join(","),...rs.map(r=>head.map(k=>String(r[k])).join(","))].join("\n")+"\n";
  writeFileSync(join(out,path),text);
}
outputCSV("presence_vs_label_exchangeability.csv",presence);
outputCSV("by_states.csv",byStates);
outputCSV("cohort_presence.csv",cohort);
outputCSV("normalized_spatial_components.csv",normalized);
console.log(JSON.stringify({presence:presence[0],spatial:normalized.find(x=>x.min_global_mismatch===.01)}));
