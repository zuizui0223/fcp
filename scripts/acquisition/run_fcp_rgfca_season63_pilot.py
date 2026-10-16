#!/usr/bin/env python3
"""Independent source-only 63-original-photo seasonal FCP pilot; NOT main P_poly."""
import argparse,csv,hashlib,json,math,re
from collections import defaultdict
from datetime import date
from pathlib import Path

SOURCE_SHA="58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01"
ID_SHA="870b08609f1192c41e36982a31c063fd223e3cfbeba0c5c8dec4b55ba2427f64"
PHOTO_RE=re.compile(r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/(\d+)/large\.(jpg|jpeg|png)",re.I)
# Fixed *before* any pixels were opened. Each group is one taxon/site/year/two months.
FROZEN={
"740971":("Ranunculus formosa-montanus",(456,1205),{
"2026-06":"673316202 682069264 685787050 685070217 681484568 683321331 675444667",
"2026-07":"691275753 710820680 690489146 693813185 694347131 693243155"}),
"870121":("Guzmania donnellsmithii",(401,381),{
"2026-01":"665007719 606517505 614094540 607052318 609074861",
"2026-02":"614890888 623094941 618854509 629177388 617539721"}),
"1023613":("Pedicularis refracta",(456,1205),{
"2026-07":"690490176 700018618 693232222 702719688 705208864 694357554 697652075",
"2026-06":"681446694 690517687 672583223 680029268 685789336"}),
"586210":("Geissorhiza confusa",(225,796),{
"2024-12":"461160175 455706328 457285131 455736419 458454796 458229466",
"2024-11":"451160782 449584967 451823046 501282527"}),
"586370":("Gladiolus angustus",(223,793),{
"2025-10":"584924576 584778828 583284361 585979770",
"2025-12":"597976075 605642267 599723384 600494818"}),
"861076":("Caladenia aurantiaca",(204,1305),{
"2024-10":"457220363 442824285 440829251 441202742 439535020 441151438",
"2024-09":"432910453 440152680 435326735 435660112"})}
SIDE=("measurement_id","inat_taxon_id","species","photo_id","observation_id",
"observer_id","latitude","longitude","screen_fold","positional_accuracy_m",
"observed_on","year_month","lat_cell","lon_cell")

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda:f.read(1024*1024),b""):h.update(buf)
    return h.hexdigest()
def photo_ids():
    ids={}
    for tax,(_,_,months) in FROZEN.items():
        if len(months)!=2 or len({x[:4] for x in months})!=1:
            raise ValueError("Not two months of the same year")
        for ym,txt in months.items():
            for pid in txt.split():
                if pid in ids:raise ValueError("Source photo repeated")
                ids[pid]=(tax,ym)
    digest=hashlib.sha256(("\n".join(sorted(ids,key=int))+"\n").encode()).hexdigest()
    if len(ids)!=63 or len(FROZEN)!=6 or digest!=ID_SHA:
        raise ValueError("Original pre-colour 63 photo IDs changed")
    return ids
def cell(r):
    lat=float(r["latitude"]);lon=float(r["longitude"])
    if not all(map(math.isfinite,(lat,lon))) or abs(lat)>90 or abs(lon)>180:
        raise ValueError("Invalid original coordinates")
    return min(719,math.floor((lat+90)*4)),math.floor(((lon+180)%360)*4)
def verify_row(r,tax,month):
    species,expected,_=FROZEN[tax]
    try:day=date.fromisoformat(r["observed_on"]);acc=float(r["positional_accuracy_m"])
    except (TypeError,ValueError) as exc:raise ValueError("Bad true observation date or location accuracy") from exc
    if (r["inat_taxon_id"]!=tax or r["species"]!=species
        or r["source_snapshot_date"]!="2026-09-27"
        or len(r["observed_on"])!=10 or r["observed_on"][:7]!=month
        or day>date(2026,9,27) or cell(r)!=expected
        or not math.isfinite(acc) or not 0<acc<=5000
        or not r["observation_id"] or not r["observer_id"]):
        raise ValueError("Original same-site same-year independent observation mismatch")
def csvout(path,columns,rows):
    with Path(path).open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=columns,lineterminator="\n");w.writeheader();w.writerows(rows)
def prepare(root,out):
    from scripts.acquisition.prepare_first_rgfca_original_monthly_roi_unit import prepare as gate
    from scripts.acquisition.build_fcp_uncapped_blind_batches import new_id
    from scripts.acquisition.freeze_fcp_photo_first_blind_workers import WORKER,ACQUISITION
    if out.exists() and any(out.iterdir()):raise ValueError("Do not overwrite pre-colour evidence")
    gate(root,out/"full_original_source_authentication")
    manifest=root/"rgfca_1000_exact_source_photo_IDs/rgfca_1000_complete_original_source_photo_IDs.csv"
    if sha(manifest)!=SOURCE_SHA:raise ValueError("Original complete 3.43m photo-ID SHA changed")
    wanted=photo_ids(); selected={}; original={};join={}
    path=root/"rgfca_1000_photo_colour_blind_selection/photo_screen_original_metadata_selected.csv"
    with path.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            pid=r["photo_id"]
            if pid not in wanted:continue
            if pid in selected:raise ValueError("Duplicate blind original selection")
            tax,_=wanted[pid];m=PHOTO_RE.fullmatch(r["photo_url_large"])
            if (r["inat_taxon_id"]!=tax or r["screen_fold"] not in ("0","1")
                or m is None or m.group(1)!=pid or not r["photo_license"]):
                raise ValueError("Source-only selected original identity/license invalid")
            selected[pid]=r
    if set(selected)!=set(wanted):raise ValueError("63 photo IDs not in existing 106376-photo pool")
    with manifest.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            pid=r["photo_id"]
            if pid not in wanted:continue
            if pid in original:raise ValueError("Duplicate original photo ID")
            tax,month=wanted[pid];verify_row(r,tax,month)
            s=selected[pid]
            if any(r[k]!=s[k] for k in ("inat_taxon_id","species","photo_id","observation_id",
                "observer_id","photo_license","photo_url_large")):
                raise ValueError("Selected photo differs from original provider archive")
            if any(abs(float(r[k])-float(s[k]))>1e-7 for k in ("latitude","longitude")):
                raise ValueError("Frozen site coordinate changed")
            original[pid]=r
    if set(original)!=set(wanted):raise ValueError("Full provider source missing frozen photo")
    mids={new_id(pid):pid for pid in wanted}
    if len(mids)!=63:raise ValueError("Measurement ID collision")
    path=root/"rgfca_1000_original_photo_blind_worker_freeze/metadata_join_key.csv"
    with path.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mid=r["measurement_id"]
            if mid not in mids:continue
            if mid in join:raise ValueError("Duplicate frozen sidecar")
            pid=mids[mid]
            if any(r[k]!=(selected[pid][k] if k=="screen_fold" else original[pid][k])
                   for k in ("inat_taxon_id","species","photo_id","observation_id","observer_id","screen_fold")):
                raise ValueError("128-worker source metadata differs from immutable original")
            join[mid]=r
    if set(join)!=set(mids):raise ValueError("Missing original frozen worker join")
    blind=[];acq=[];side=[];seen_obs=set();month_users=defaultdict(set);month_folds=defaultdict(set)
    for tax,(species,site,months) in FROZEN.items():
        for month,words in months.items():
            for pid in words.split():
                r=original[pid];s=selected[pid];mid=new_id(pid)
                if r["observation_id"] in seen_obs:raise ValueError("Same observation photographed twice")
                seen_obs.add(r["observation_id"])
                month_users[(tax,month)].add(r["observer_id"]);month_folds[(tax,month)].add(s["screen_fold"])
                file=mid+".jpg"
                blind.append(dict(measurement_id=mid,image_filename=file,photo_license=r["photo_license"]))
                acq.append(dict(measurement_id=mid,image_filename=file,
                     photo_url_large=r["photo_url_large"],photo_license=r["photo_license"]))
                side.append(dict(measurement_id=mid,inat_taxon_id=tax,species=species,
                     photo_id=pid,observation_id=r["observation_id"],observer_id=r["observer_id"],
                     latitude=r["latitude"],longitude=r["longitude"],screen_fold=s["screen_fold"],
                     positional_accuracy_m=r["positional_accuracy_m"],observed_on=r["observed_on"],
                     year_month=month,lat_cell=site[0],lon_cell=site[1]))
    if len(month_users)!=12 or any(len(u)<2 for u in month_users.values()) or any(
            folds!={"0","1"} for folds in month_folds.values()):
        raise ValueError("Each original month must have two distinct observers/folds")
    out.mkdir(parents=True,exist_ok=True)
    files=[("season63_blind.csv",WORKER,blind),("season63_acquisition.csv",ACQUISITION,acq),
           ("season63_join.csv",SIDE,side)]
    for n,c,rows in files:csvout(out/n,c,rows)
    report={"status":"FROZEN_ORIGINAL_PRECOLOUR_SEASON63_NO_PIXELS",
      "n_photo_ids":63,"n_taxa":6,"n_taxon_months":12,"n_independent_original_observations":63,
      "original_source_sha256":SOURCE_SHA,"photo_id_set_sha256":ID_SHA,
      "all_128_original_worker_shas_verified":True,"source_selection_not_morph_based":True,
      "no_photo_pixels_opened":True,"no_genetic_FCP_or_climate_effect_estimated":True,
      "file_sha256":{n:sha(out/n) for n,_,_ in files}}
    (out/"precolour_receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report
class SourceIndex:
    def __init__(self,root):
        import pandas as pd
        from scripts.acquisition.freeze_fcp_photo_first_blind_workers import WORKER,ACQUISITION
        self.worker=pd.read_csv(root/"season63_blind.csv",dtype=str).fillna("")
        a=pd.read_csv(root/"season63_acquisition.csv",dtype=str).fillna("")
        if (tuple(self.worker.columns)!=WORKER or tuple(a.columns)!=ACQUISITION
            or len(self.worker)!=63 or len(a)!=63
            or self.worker.measurement_id.duplicated().any() or a.measurement_id.duplicated().any()):
            raise ValueError("63-photo ROI evaluator input not blind and unique")
        self.acq=self.worker.merge(a,on=list(WORKER),validate="one_to_one",how="left")
        if self.acq.photo_url_large.isna().any():raise ValueError("Blind worker lacks source URL")
    def get_batch(self,part,batch,size):
        if (part,batch,size)!=(0,0,64):raise ValueError("Only bounded 63-photo batch admitted")
        return self.worker,self.acq,1
def month_qc(records):
    by=defaultdict(list)
    for r in records:by[(r["inat_taxon_id"],r["year_month"])].append(r)
    months=[]
    for (tax,month),rows in sorted(by.items()):
        obs=defaultdict(list); good=mixed=failed=0
        for r in rows:
            if r["measurement_status"]=="image_acquisition_failed":failed+=1
            if r["measurement_status"]!="classified_four_state_morph":continue
            if r["morph"] not in {"white","yellow_orange","red_pink","blue_purple"}:
                raise ValueError("Invalid original photo colour category")
            good+=1
            w=float(r["flower_fraction_white"]);y=float(r["flower_fraction_yellow"])
            if not all(map(math.isfinite,(w,y))) or not(0<=w<=1 and 0<=y<=1):
                raise ValueError("Invalid original photo palette fractions")
            if w>=.1 and y>=.1:mixed+=1;continue
            obs[r["observation_id"]].append((r["observer_id"],r["morph"]))
        valid=[];contradictions=0
        for vals in obs.values():
            if len({o for o,m in vals})!=1 or len({m for o,m in vals})!=1:contradictions+=1
            else:valid.append(vals[0])
        n=len(valid);white=sum(m=="white" for o,m in valid)
        months.append(dict(taxon=tax,species=FROZEN[tax][0],month=month,
            n_source_photo_ids=len(rows),n_classifiable_rois=good,n_provider_failures=failed,
            n_same_roi_white_yellow_mixed_excluded=mixed,n_discordant_observations=contradictions,
            n_usable_independent_observations=n,n_usable_observers=len({o for o,m in valid}),
            n_white_obs=white,n_nonwhite_obs=n-white,descriptive_white_rate=f"{white/n:.8f}" if n else ""))
    contrasts=[]
    for tax,(_,_,expected) in FROZEN.items():
        groups=sorted([r for r in months if r["taxon"]==tax],key=lambda r:r["month"])
        if len(groups)!=2 or set(x["month"] for x in groups)!=set(expected):
            raise ValueError("Missing frozen original species/month after photo measurement")
        eligible=all(r["n_usable_independent_observations"]>=2 and
                     r["n_usable_observers"]>=2 for r in groups)
        contrasts.append(dict(taxon=tax,month1=groups[0]["month"],month2=groups[1]["month"],
            n_month1=groups[0]["n_usable_independent_observations"],
            n_month2=groups[1]["n_usable_independent_observations"],
            both_months_two_obs_and_observers=eligible,
            descriptive_white_rate_change=(
                f"{float(groups[1]['descriptive_white_rate'])-float(groups[0]['descriptive_white_rate']):.8f}"
                if eligible else ""),causal_or_confirmatory_inference=False))
    return months,contrasts
def measure(root,out,method,detector,sam):
    import pandas as pd
    from scripts.acquisition.run_fcp_fresh_image_partition import run_batch
    a=json.loads((root/"precolour_receipt.json").read_text())
    if a.get("status")!="FROZEN_ORIGINAL_PRECOLOUR_SEASON63_NO_PIXELS" or a.get("n_photo_ids")!=63:
        raise ValueError("No sealed original 63-image sampling certificate")
    if any(sha(root/n)!=h for n,h in a["file_sha256"].items()):
        raise ValueError("Original source blind worker tampered")
    if out.exists() and any(out.iterdir()):raise ValueError("Photo output already exists")
    result=run_batch(root,0,0,64,out,source=method,detector=detector,sam=sam,indexed_source=SourceIndex(root))
    original=pd.read_csv(root/"season63_join.csv",dtype=str).fillna("")
    terminal=pd.read_csv(out/"terminal_p000_b0000000.csv",dtype=str).fillna("")
    j=original.merge(terminal,on="measurement_id",validate="one_to_one")
    if len(j)!=63 or result["n_source_photo_ids"]!=63:raise ValueError("Loss of source photo in blind ROI join")
    jf=out/"season63_measured_postblind.csv";j.to_csv(jf,index=False)
    months,contrasts=month_qc(j.to_dict("records"))
    if len(months)!=12 or len(contrasts)!=6:raise ValueError("Month/source outcome lost")
    mf=out/"season63_month_white_QC.csv";cf=out/"season63_siteyear_month_contrasts.csv"
    csvout(mf,list(months[0]),months);csvout(cf,list(contrasts[0]),contrasts)
    rec={"status":"REAL_ORIGINAL_EXPLORATORY_SEASON63_ROI_MEASURED",
        "n_original_photo_pixels_attempted":63,"n_four_state_classifiable":result["n_classifiable"],
        "n_taxa":6,"n_month_groups":12,
        "n_species_with_two_usable_months":sum(q["both_months_two_obs_and_observers"] for q in contrasts),
        "n_same_roi_organ_mixtures_excluded":sum(m["n_same_roi_white_yellow_mixed_excluded"] for m in months),
        "real_photo_source_date_added_only_after_blind_ROI":True,
        "not_Ppoly_climate_soil_or_genetic_FCP":True,
        "files_sha256":{f.name:sha(f) for f in (jf,mf,cf)}}
    (out/"season63_real_roi_receipt.json").write_text(json.dumps(rec,indent=2,sort_keys=True)+"\n")
    return rec
def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest="action",required=True)
    a=s.add_parser("prepare");a.add_argument("--source-root",required=True,type=Path)
    a.add_argument("--outdir",required=True,type=Path)
    b=s.add_parser("measure");b.add_argument("--prepared",required=True,type=Path)
    b.add_argument("--outdir",required=True,type=Path)
    b.add_argument("--method",required=True,type=Path)
    b.add_argument("--detector",required=True,type=Path);b.add_argument("--sam",required=True,type=Path)
    args=p.parse_args()
    result=prepare(args.source_root,args.outdir) if args.action=="prepare" else measure(
        args.prepared,args.outdir,args.method,args.detector,args.sam)
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
