#!/usr/bin/env python3
"""Outcome-blind 2026-09-27 original-photo season-matched measurement queue.

Source: SHA-sealed all 3,430,502 original RGFCA photo IDs plus independently
sealed original observed-on calendar/site source counts. This code neither opens
photographs nor consumes image-colour, literature or phenology annotations.
The 2x2 32-photo work units do not cap the full source population.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from collections import defaultdict
from datetime import date
from pathlib import Path

SNAPSHOT="2026-09-27"
SALT="FCP_SEASONAL_SOURCE_OBSERVER_DISJOINT_20261016_V1"
PHOTO_SALT="FCP_ALL_OBSERVABLE_PHOTO_BLIND_20261014_V1"
VALID_SOURCE="COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"
VALID_CALENDAR="COMPLETE_ORIGINAL_SOURCE_CALENDAR_CENSUS_NO_FLOWER_PIXELS"
VALID_URL=re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/([0-9]+)/large\.(?:jpg|jpeg|png)$",re.I)
SOURCE_FIELDS=("inat_taxon_id","species","photo_id","observation_id","observer_id",
 "latitude","longitude","observed_on","positional_accuracy_m","photo_license",
 "photo_creator_id","attribution","photo_url_large","source_snapshot_date")
FIELDS=("measurement_id","inat_taxon_id","species","photo_id","observation_id",
 "observer_id","observed_on","observed_year_month","latitude","longitude",
 "positional_accuracy_m","photo_url_large","photo_license","attribution",
 "source_snapshot_date","local_season","season_year","site_cell_lat_0p25",
 "site_cell_lon_0p25","wave","fold","selection_hash")
BLIND=("measurement_id","image_filename","photo_license")
ACQ=("measurement_id","image_filename","photo_url_large","photo_license")

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()

def ranking(*v):
    return hashlib.sha256(("|".join([SALT,*map(str,v)])).encode()).hexdigest()

def measure_id(pid):
    return "FCPR-"+hashlib.sha256(f"{PHOTO_SALT}\x1fphoto\x1f{pid}".encode()).hexdigest().upper()[:24]

def local_year(month,year,latcell):
    latitude=-90+(latcell+.5)/4
    if abs(latitude)<23.5:return None
    s="winter" if month in (12,1,2) else "spring" if month in (3,4,5) else "summer" if month in (6,7,8) else "autumn"
    if latitude<0:s={"winter":"summer","spring":"autumn","summer":"winter","autumn":"spring"}[s]
    return s,year+int(month==12)

def precolour_sites(csvpath,calendar,*,synthetic=False):
    d=json.loads(Path(calendar).read_text())
    if (d.get("status")!=VALID_CALENDAR
        or d.get("original_source_site_month_file_sha256")!=sha(csvpath)
        or d.get("source_photo_opportunity_is_not_flowering_detection") is not True
        or d.get("new_flower_colour_or_white_frequency_estimated") is not False
        or (not synthetic and (d.get("n_original_photo_ids")!=3430502 or d.get("n_taxa")!=1000))):
        raise ValueError("Not exact independently sealed full original site-month opportunity census")
    accum=defaultdict(lambda:[0,0])
    names={}
    with Path(csvpath).open(newline="",encoding="utf-8") as h:
        reader=csv.DictReader(h)
        need={"inat_taxon_id","species","year_month","lat_0p25_cell","lon_0p25_cell",
              "all_source_photo_ids_with_5km_accuracy","distinct_source_observations","distinct_source_observers"}
        if not need.issubset(reader.fieldnames or []):raise ValueError("Source site-month schema changed")
        for r in reader:
            tid=r["inat_taxon_id"]; name=r["species"]
            if not tid.isdecimal() or not name or (tid in names and names[tid]!=name):raise ValueError("Original taxon identity drift")
            names[tid]=name
            ym=r["year_month"]
            try:
                if len(ym)!=7 or ym[4]!="-":raise ValueError()
                yr=int(ym[:4]);mo=int(ym[5:7]);date(yr,mo,1)
                lat=int(r["lat_0p25_cell"]);lon=int(r["lon_0p25_cell"])
                photos=int(r["all_source_photo_ids_with_5km_accuracy"])
                obs=int(r["distinct_source_observations"]);reporters=int(r["distinct_source_observers"])
            except (TypeError,ValueError) as e:raise ValueError("Invalid source month/site cell or count") from e
            if not 0<=lat<=719 or not 0<=lon<=1439 or not 0<=reporters<=obs<=photos:
                raise ValueError("Invalid source position/count bound")
            if obs<2 or reporters<2 or yr<1958 or yr>2025:continue
            sy=local_year(mo,yr,lat)
            if sy is None or sy[1]>2025:continue
            s,year=sy
            agg=accum[(tid,lat,lon,year,s)]
            agg[0]+=photos;agg[1]+=obs
    by_site=defaultdict(dict)
    for (tax,lat,lon,year,s),(photos,observations) in accum.items():
        if photos>=24 and observations>=24:by_site[(tax,lat,lon,year)][s]=(photos,observations)
    pairs=[]
    for key,months in by_site.items():
        if len(months)<2:continue
        options=[]
        for i,(a,A) in enumerate(sorted(months.items())):
            for b,B in list(sorted(months.items()))[i+1:]:
                options.append(((min(A[1],B[1]),min(A[0],B[0]),A[1]+B[1]),a,b))
        _,a,b=max(options,key=lambda x:x[0])
        pairs.append((key,names[key[0]],a,b))
    pairs.sort()
    if not synthetic and (len(pairs)!=268 or len({p[0][0] for p in pairs})!=39):
        raise ValueError("Frozen 268-site/39-taxon outcome-blind source opportunity frame drift")
    wanted={(tax,lat,lon,year,s) for (tax,lat,lon,year),_,a,b in pairs for s in (a,b)}
    return pairs,wanted

def date_season(day,lat):
    try:
        if len(day)!=10 or day[4]!="-" or day[7]!="-":return None
        d=date.fromisoformat(day)
        if d>date.fromisoformat(SNAPSHOT) or d.year<1958 or d.year>2025:return None
    except ValueError:return None
    if abs(lat)<23.5:return None
    s="winter" if d.month in (12,1,2) else "spring" if d.month in (3,4,5) else "summer" if d.month in (6,7,8) else "autumn"
    if lat<0:s={"winter":"summer","spring":"autumn","summer":"winter","autumn":"spring"}[s]
    return s,d.year+int(d.month==12)

def one_per_observer(original):
    rep={}
    for r in original.values():
        owner=r["observer_id"];key=(ranking("observer-representative",r["observation_id"],r["photo_id"]),r["photo_id"])
        if owner not in rep or key<rep[owner][0]:rep[owner]=(key,r)
    return {k:v[1] for k,v in rep.items()}

def month_round_robin(reps,n):
    bins=defaultdict(list)
    for r in reps.values():bins[r["observed_on"][:7]].append(r)
    for m in bins:bins[m].sort(key=lambda x:(ranking("selected-observer",x["observer_id"],x["photo_id"]),x["photo_id"]))
    sel=[]
    while len(sel)<n:
        made=False
        for month in sorted(bins):
            if bins[month]:
                sel.append(bins[month].pop(0));made=True
                if len(sel)==n:break
        if not made:raise ValueError("Observer reservation not available")
    return sel

def freeze(manifest,source_receipt,site,calendar,out,*,per_season_total=64,synthetic=False):
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError("Do not overwrite original-season source photo freeze")
    if per_season_total<4 or per_season_total>256 or per_season_total%4:raise ValueError("Invalid bounded seasonal work depth")
    doc=json.loads(Path(source_receipt).read_text())
    if (doc.get("status")!=VALID_SOURCE or doc.get("source_snapshot_date")!=SNAPSHOT
        or doc.get("selected_all_photo_IDs_csv_sha256")!=sha(manifest)
        or doc.get("all_enrolled_species_match_independent_full_source_photo_ID_counts") is not True
        or doc.get("original_image_pixels_opened") is not False
        or doc.get("historical_flower_colour_CSV_imported") is not False
        or (not synthetic and doc.get("n_selected_original_photo_IDs")!=3430502)):
        raise ValueError("Original complete source photo-ID authenticity gate failed")
    pairs,wanted=precolour_sites(site,calendar,synthetic=synthetic)
    photo_by_observation=defaultdict(dict)
    n=0;n_eligible_rows=0
    with Path(manifest).open(newline="",encoding="utf-8") as h:
        reader=csv.DictReader(h)
        if tuple(reader.fieldnames or ())!=SOURCE_FIELDS:raise ValueError("Original source schema was changed")
        for row in reader:
            n+=1
            if row["source_snapshot_date"]!=SNAPSHOT or not row["photo_id"].isascii() or not row["photo_id"].isdecimal():raise ValueError("Original source ID/date invalid")
            try:
                lat=float(row["latitude"]);lon=float(row["longitude"]);accuracy=float(row["positional_accuracy_m"])
                if not (math.isfinite(lat) and math.isfinite(lon) and math.isfinite(accuracy) and -90<=lat<=90 and -180<=lon<=180 and 0<accuracy<=5000):continue
            except (ValueError,TypeError):continue
            period=date_season(row["observed_on"],lat)
            if period is None:continue
            season,year=period
            ilat=min(719,math.floor((lat+90)*4));ilon=math.floor(((lon+180)%360)*4)
            key=(row["inat_taxon_id"],ilat,ilon,year,season)
            if key not in wanted:continue
            match=VALID_URL.fullmatch(row["photo_url_large"])
            if not match or match.group(1)!=row["photo_id"] or not row["observer_id"] or not row["observation_id"] or not row["photo_license"]:raise ValueError("Original selected photo URL or observation identity invalid")
            n_eligible_rows+=1
            ident=row["observation_id"];old=photo_by_observation[key].get(ident)
            if old:
                if any(old[k]!=row[k] for k in ("species","observer_id","observed_on","latitude","longitude","positional_accuracy_m")):raise ValueError("Same original observation has conflicting metadata")
                if ranking("photo-within-observation",row["photo_id"]) < ranking("photo-within-observation",old["photo_id"]):
                    photo_by_observation[key][ident]=row
            else:photo_by_observation[key][ident]=row
    if n!=doc["n_selected_original_photo_IDs"]:raise ValueError("Original full source truncation")
    candidates=[]
    for key,name,a,b in pairs:
        A=photo_by_observation.get((*key,a),{});B=photo_by_observation.get((*key,b),{})
        if any(r["species"]!=name for r in [*A.values(),*B.values()]):raise ValueError("Source taxon name drift")
        aa=one_per_observer(A);bb=one_per_observer(B)
        disjoint_a={u:r for u,r in aa.items() if u not in bb}
        disjoint_b={u:r for u,r in bb.items() if u not in aa}
        evidence=dict(inat_taxon_id=key[0],species=name,site_cell_lat_0p25=key[1],site_cell_lon_0p25=key[2],season_year=key[3],
             local_season_a=a,local_season_b=b,source_observations_a=len(A),source_observations_b=len(B),
             observers_a=len(aa),observers_b=len(bb),observers_shared=len(set(aa)&set(bb)),
             disjoint_observers_a=len(disjoint_a),disjoint_observers_b=len(disjoint_b))
        candidates.append((evidence,disjoint_a,disjoint_b))
    candidates.sort(key=lambda x:(-min(x[0]["disjoint_observers_a"],x[0]["disjoint_observers_b"]),
        -min(x[0]["source_observations_a"],x[0]["source_observations_b"]),
        int(x[0]["inat_taxon_id"]),x[0]["site_cell_lat_0p25"],x[0]["site_cell_lon_0p25"],x[0]["season_year"]))
    adequate=[x for x in candidates if min(x[0]["disjoint_observers_a"],x[0]["disjoint_observers_b"])>=per_season_total]
    if not adequate:raise ValueError("Insufficient observer-disjoint source support; never replace on colour")
    pair,A,B=adequate[0]
    selected=[]
    for season,reps in ((pair["local_season_a"],A),(pair["local_season_b"],B)):
        for index,r in enumerate(month_round_robin(reps,per_season_total)):
            wave="first" if index<per_season_total//2 else "reserve"
            mid=measure_id(r["photo_id"])
            selected.append(dict(
                measurement_id=mid,**{k:r[k] for k in (
                    "inat_taxon_id","species","photo_id","observation_id","observer_id",
                    "observed_on","latitude","longitude","positional_accuracy_m",
                    "photo_url_large","photo_license","attribution","source_snapshot_date")},
                observed_year_month=r["observed_on"][:7],local_season=season,
                season_year=pair["season_year"],site_cell_lat_0p25=pair["site_cell_lat_0p25"],
                site_cell_lon_0p25=pair["site_cell_lon_0p25"],wave=wave,fold=index%2,
                selection_hash=ranking("photo-within-observation",r["photo_id"])[:20]))
    if any(len({r[k] for r in selected})!=len(selected) for k in ("photo_id","observation_id","observer_id","measurement_id")):
        raise ValueError("Repeated independent observation/observer/photo in frozen sample")
    out.mkdir(parents=True,exist_ok=True)
    hashes={}
    def save(name,columns,rows):
        path=out/name
        with path.open("x",newline="",encoding="utf-8") as h:
            w=csv.DictWriter(h,fieldnames=columns,lineterminator="\n")
            w.writeheader();w.writerows(rows)
        hashes[name]=sha(path)
    save("original_source_season_frozen_128_photos.csv",FIELDS,selected)
    save("all_268_precolour_site_observer_support.csv",tuple(candidates[0][0]),[x[0] for x in candidates])
    for wave in ("first","reserve"):
        rows=[r for r in selected if r["wave"]==wave]
        save("season_"+wave+"_original_blind_worker.csv",BLIND,[dict(
            measurement_id=r["measurement_id"],image_filename=r["measurement_id"]+".jpg",
            photo_license=r["photo_license"]) for r in rows])
        save("season_"+wave+"_original_acquisition.csv",ACQ,[dict(
            measurement_id=r["measurement_id"],image_filename=r["measurement_id"]+".jpg",
            photo_url_large=r["photo_url_large"],photo_license=r["photo_license"]) for r in rows])
    report=dict(
        schema="fcp_observer_disjoint_original_seasonal_photo_freeze_v1",
        status="SYNTHETIC_NO_BIOLOGICAL_MEASUREMENT" if synthetic else "REAL_ORIGINAL_SEASONAL_SOURCE_PHOTOS_FROZEN_NO_PIXELS",
        original_source_snapshot=SNAPSHOT,original_source_photo_ids_count=n,
        precolour_site_year_pairs=len(pairs),precolour_pairs_with_enough_independent_observers=len(adequate),
        source_photos_matching_prior_site_months=n_eligible_rows,
        selected_source_pair=pair,selected_photo_ids=len(selected),
        per_season_total_not_global_inclusion_cap=per_season_total,
        photos_per_wave=per_season_total,observer_ids_unique_across_seasons=True,
        prior_random_sample_106376_stays_unchanged=True,
        original_pixels_opened=0,flowering_or_fruiting_from_photo_not_annotated=True,
        no_climate_or_colour_outcome_used_to_select_photos=True,
        original_manifest_sha256=sha(manifest),
        original_source_receipt_sha256=sha(source_receipt),
        original_calendar_site_month_sha256=sha(site),
        original_calendar_receipt_sha256=sha(calendar),
        file_sha256=hashes)
    (out/"original_season_photo_freeze_receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("complete-original-source","original-source-receipt","source-site-month","original-calendar-receipt","outdir"):
        p.add_argument("--"+flag,type=Path,required=True)
    p.add_argument("--per-season-total",type=int,default=64)
    x=p.parse_args()
    print(json.dumps(freeze(x.complete_original_source,x.original_source_receipt,
        x.source_site_month,x.original_calendar_receipt,x.outdir,
        per_season_total=x.per_season_total),sort_keys=True))
if __name__=="__main__":main()
