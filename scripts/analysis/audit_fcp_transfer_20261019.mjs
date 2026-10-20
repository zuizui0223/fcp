// Recompute exploratory FCP transfer diagnostics from source species-level CSV files.
// Node 20+, standard library only. Run: node scripts/analysis/audit_fcp_transfer_20261019.mjs
import {readFileSync, mkdirSync, writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
const PATHS={
  genus:'results/fcp_environment_transport_20261019/genus_species_metrics.csv',
  spatial:'results/fcp_environment_transport_20261019/spatial_species_metrics.csv',
  output:'results/fcp_transfer_independent_audit_20261019/summary.json'
};
function csv(path){
 const raw=readFileSync(path,'utf8').trim();
 const ls=raw.split(/\r?\n/),h=ls.shift().split(',');
 return {hash:createHash('sha256').update(raw+'\n').digest('hex'),rows:ls.map(l=>{
   const vs=l.split(',');if(vs.length!==h.length)throw Error('Unexpected CSV fields: '+path);
   return Object.fromEntries(h.map((name,i)=>[name,vs[i]]));
 })};
}
const mean=(xs,field)=>xs.reduce((s,x)=>s+Number(x[field]),0)/xs.length;
const group=(xs,key)=>{const m=new Map();for(const x of xs){const v=key(x);if(!m.has(v))m.set(v,[]);m.get(v).push(x)}return m};
const frac=(xs,fn)=>xs.filter(fn).length/xs.length;
const quant=(xs,p)=>{const s=[...xs].sort((a,b)=>a-b);return s[Math.floor((s.length-1)*p)]};
function rng(seed){let st=seed|0;return()=>{st^=st<<13;st^=st>>>17;st^=st<<5;return(st>>>0)/4294967296}};
function seedFor(str){let s=20261019;for(const c of str)s=Math.imul(s^c.charCodeAt(0),16777619)>>>0;return s||1}
function clusterCI(rows,field,seed,n=5000){
 const gs=[...group(rows,r=>r.genus).values()].map(a=>({
   n:a.length,total:a.reduce((s,r)=>s+Number(r[field]),0)
 }));
 const next=rng(seed),values=[];
 for(let b=0;b<n;b++){let n=0,sum=0;for(let i=0;i<gs.length;i++){const g=gs[Math.floor(next()*gs.length)];sum+=g.total;n+=g.n}values.push(sum/n)}
 return [quant(values,.025),quant(values,.975)];
}
function metricSummary(xs){return {
 species:xs.length,genera:group(xs,r=>r.genus).size,
 mean_gain:mean(xs,'pair_logloss_gain'),
 mean_auc:mean(xs,'pair_auc'),
 fraction_gain_positive:frac(xs,r=>Number(r.pair_logloss_gain)>0)
}};
const G=csv(PATHS.genus),S=csv(PATHS.spatial);
const bio=G.rows.filter(r=>r.model==='BIO5_component'),joint=G.rows.filter(r=>r.model==='joint_three');
if(bio.length!==1574||joint.length!==1574)throw Error('Unexpected genus panel');
const lookup=new Map(joint.map(r=>[r.inat_taxon_id,r]));
const paired=bio.map(r=>{
 const other=lookup.get(r.inat_taxon_id);if(!other||other.genus!==r.genus)throw Error('Missing joint pair '+r.inat_taxon_id);
 return {genus:r.genus,delta:Number(other.pair_logloss_gain)-Number(r.pair_logloss_gain)};
});
const results={
 provenance:{generated_with:'scripts/analysis/audit_fcp_transfer_20261019.mjs',
  source_csv_sha256:{genus:G.hash,spatial:S.hash},
  input_rows:{genus:G.rows.length,spatial:S.rows.length},
  inference:'Exploratory post hoc audit of saved species-level cross-validation results; not retrained models.'},
 genus:{
  bio5:{...metricSummary(bio),genus_bootstrap_95:clusterCI(bio,'pair_logloss_gain',seedFor('BIO5_component'))},
  joint:{...metricSummary(joint),genus_bootstrap_95:clusterCI(joint,'pair_logloss_gain',seedFor('joint_three'))},
  paired_joint_minus_bio5_mean:mean(paired,'delta'),
  paired_joint_minus_bio5_genus_bootstrap_95:clusterCI(paired,'delta',seedFor('joint_minus_bio5')),
  by_cohort:[...group(bio,r=>r.source_cohort)].map(([cohort,a])=>({cohort,...metricSummary(a),genus_bootstrap_95:clusterCI(a,'pair_logloss_gain',seedFor('cohort_'+cohort))})),
  by_white_frequency:[
    {label:'0-20%',rows:bio.filter(r=>Number(r.n_white)/Number(r.n_photos)<=.2)},
    {label:'20-50%',rows:bio.filter(r=>Number(r.n_white)/Number(r.n_photos)>.2&&Number(r.n_white)/Number(r.n_photos)<=.5)},
    {label:'50-80%',rows:bio.filter(r=>Number(r.n_white)/Number(r.n_photos)>.5&&Number(r.n_white)/Number(r.n_photos)<=.8)},
    {label:'80-100%',rows:bio.filter(r=>Number(r.n_white)/Number(r.n_photos)>.8)}
  ].map(v=>({bin:v.label,...metricSummary(v.rows)}))
 },
 spatial:[...group(S.rows,r=>r.buffer_km+'|'+r.comparison)].map(([key,a])=>{
  const sorted=[...a].sort((x,y)=>Number(x.logloss_gain)-Number(y.logloss_gain));
  const byCohort=[...group(a,r=>r.source_cohort)].map(([cohort,vs])=>({
    cohort,n:vs.length,mean_gain:mean(vs,'logloss_gain'),fraction_positive:frac(vs,x=>Number(x.logloss_gain)>0)
  }));
  const noExtreme=a.filter(x=>Number(x.augmented_abs_train_z_gt_3_photos)===0);
  return {key,species:a.length,genera:group(a,r=>r.genus).size,
   photos:a.reduce((s,x)=>s+Number(x.n_photos),0),
   mean_gain:mean(a,'logloss_gain'),mean_brier_gain:mean(a,'brier_gain'),
   median_gain:quant(a.map(r=>Number(r.logloss_gain)),.5),
   fraction_positive:frac(a,x=>Number(x.logloss_gain)>0),
   genus_bootstrap_95:clusterCI(a,'logloss_gain',seedFor(key)),
   drop_most_negative_10percent_mean:mean(sorted.slice(Math.floor(a.length*.1)),'logloss_gain'),
   no_absolute_train_z_over_3_n:noExtreme.length,
   no_absolute_train_z_over_3_mean:noExtreme.length?mean(noExtreme,'logloss_gain'):null,
   by_cohort:byCohort};
 })
};
if(S.rows.length!==5607)throw Error('Unexpected spatial panel');
mkdirSync(PATHS.output.split('/').slice(0,-1).join('/'),{recursive:true});
writeFileSync(PATHS.output,JSON.stringify(results,null,2)+'\n');
console.log(JSON.stringify({output:PATHS.output,counts:results.provenance.input_rows,genus:results.genus.bio5.mean_gain,spatial:results.spatial.map(x=>[x.key,x.mean_gain])},null,2));
