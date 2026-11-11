"""Audited U100 follow-up bridge; never re-estimates frozen H2 white-axis tests.

Input data are explicit: the 4,730-species opportunity ledger, disjoint frozen or
new photo tables, and optionally photo-level environment plus study-specific H2
and spatial species metrics. Missing photographs do not imply monomorphism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
ENVIRONMENTS = ("bio5", "bio14", "srad_mean")
MIN_D = 40
MIN_ENV_PER_COLOUR = 5
MIN_PAIRED_OBSERVERS = 3
COHORT_SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "prospective": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def require(df: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(df.columns))
    if missing:
        raise ValueError(f"{label}: missing columns: {missing}")


def bool_values(col: pd.Series) -> pd.Series:
    return col.fillna("").astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def identifiers(df: pd.DataFrame, label: str) -> pd.DataFrame:
    df = df.copy()
    df["inat_taxon_id"] = pd.to_numeric(df["inat_taxon_id"], errors="raise").astype("int64")
    if (df["inat_taxon_id"] <= 0).any():
        raise ValueError(f"{label}: invalid taxon IDs")
    return df


def opportunity(path: Path, expected: int) -> pd.DataFrame:
    frame = identifiers(pd.read_csv(path, low_memory=False), "U100")
    require(frame, {"inat_taxon_id", "species"}, "U100")
    if frame["inat_taxon_id"].duplicated().any() or frame["species"].isna().any():
        raise ValueError("U100: taxon IDs must be unique and species names present")
    if "after_observer_cap" in frame:
        cap = pd.to_numeric(frame["after_observer_cap"], errors="coerce")
        if cap.isna().any() or (cap < 100).any():
            raise ValueError("U100: candidate has missing or <100 observer-capped photos")
    if len(frame) != expected:
        raise ValueError(f"U100: expected {expected} species, found {len(frame)}")
    frame = frame.sort_values("inat_taxon_id", kind="stable").reset_index(drop=True)
    return frame


def parse_cohorts(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise ValueError("--measured must be cohort=path")
        label, path = item.split("=", 1)
        if not label or not path or label in out or "=" in label:
            raise ValueError(f"duplicate or invalid cohort: {label}")
        out[label] = Path(path)
    return out


def photo_tables(paths: dict[str, Path], frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    sheets, hashes = [], {}
    reference = dict(zip(frame["inat_taxon_id"], frame["species"].astype(str)))
    seen_species: set[int] = set()
    seen_photos: set[str] = set()
    for cohort, path in paths.items():
        actual_hash = sha256(path)
        if cohort in COHORT_SHA256 and actual_hash != COHORT_SHA256[cohort]:
            raise ValueError(f"{cohort}: immutable measured table SHA256 mismatch")
        d = identifiers(pd.read_csv(path, low_memory=False), cohort)
        require(d, {"inat_taxon_id", "species", "photo_id", "morph", "global_classifiable"}, cohort)
        if d["photo_id"].isna().any():
            raise ValueError(f"{cohort}: missing photo_id")
        d["photo_id"] = d["photo_id"].astype(str)
        if d["photo_id"].duplicated().any():
            raise ValueError(f"{cohort}: duplicate photo_id")
        taxa = set(map(int, d["inat_taxon_id"].unique()))
        if taxa - reference.keys():
            raise ValueError(f"{cohort}: measured taxa missing from U100")
        if taxa & seen_species:
            raise ValueError("A species occurs in multiple disjoint measurement cohorts")
        seen_species.update(taxa)
        photos = set(d["photo_id"])
        if photos & seen_photos:
            raise ValueError("Photo ID reused across supposedly disjoint cohorts")
        seen_photos.update(photos)
        for ident, name in d[["inat_taxon_id", "species"]].drop_duplicates().itertuples(index=False, name=None):
            if str(name).strip() != str(reference[int(ident)]).strip():
                raise ValueError(f"{cohort}: taxon-name mismatch for {ident}")
        d["cohort"] = cohort
        d["usable_colour"] = bool_values(d["global_classifiable"]) & d["morph"].isin(MORPHS)
        sheets.append(d)
        hashes[cohort] = {"sha256": actual_hash, "rows": int(len(d)), "species": len(taxa)}
    if not sheets:
        return pd.DataFrame(), hashes
    return pd.concat(sheets, ignore_index=True), hashes


def registry(frame: pd.DataFrame, photos: pd.DataFrame, nshards: int) -> pd.DataFrame:
    if nshards < 1:
        raise ValueError("shards must be >=1")
    out = frame[[x for x in ("inat_taxon_id", "species", "family", "after_observer_cap") if x in frame]].copy()
    hvals = [hashlib.sha256(f"FCP_U100_FOLLOWUP_V1|{i}|{s}".encode()).hexdigest()
             for i, s in out[["inat_taxon_id", "species"]].itertuples(index=False, name=None)]
    out["order_hash"] = hvals
    out["shard"] = [int(h[:16], 16) % nshards for h in hvals]
    out["cohort"] = pd.NA
    out["n_measured"] = 0
    out["n_classifiable"] = 0
    out["n_states_observed"] = pd.NA
    out["white_fraction"] = np.nan
    out["D"] = np.nan
    out["D_pair"] = np.nan
    out["measurement_status"] = "not_measured"
    if photos.empty:
        return out
    rows = []
    for taxon, g in photos.groupby("inat_taxon_id", sort=True):
        sub = g.loc[g["usable_colour"]]
        n = len(sub)
        states = sub["morph"].value_counts()
        proportions = (states / n).to_numpy(float) if n else np.array([], dtype=float)
        D = float(1.0 - (proportions ** 2).sum()) if n >= MIN_D else np.nan
        rows.append({
            "inat_taxon_id": taxon,
            "cohort": str(g["cohort"].iloc[0]),
            "n_measured": int(len(g)),
            "n_classifiable": int(n),
            "n_states_observed": int(len(states)) if n >= MIN_D else pd.NA,
            "white_fraction": float(states.get("white", 0) / n) if n >= MIN_D else np.nan,
            "D": D,
            "D_pair": D * n / (n - 1) if n >= MIN_D else np.nan,
            "measurement_status": "D_evaluable" if n >= MIN_D else "insufficient_classifiable",
        })
    from_photos = pd.DataFrame(rows)
    out = out.drop(columns=[c for c in from_photos if c != "inat_taxon_id"])
    out = out.merge(from_photos, how="left", on="inat_taxon_id", validate="one_to_one")
    out["measurement_status"] = out["measurement_status"].fillna("not_measured")
    out[["n_measured", "n_classifiable"]] = out[["n_measured", "n_classifiable"]].fillna(0).astype(int)
    return out


def attach_species_metrics(ledger: pd.DataFrame, path: Path, *, kind: str) -> pd.DataFrame:
    extra = pd.read_csv(path, low_memory=False)
    extra = identifiers(extra, kind)
    if kind == "spatial":
        require(extra, {"inat_taxon_id", "cohort", "D_pair", "local_pairs_50km", "depletion_50km"}, kind)
        selected = ["inat_taxon_id", "cohort", "D_pair", "local_pairs_50km", "depletion_50km"]
        extra = extra[selected].rename(columns={"D_pair": "spatial_D_pair"})
    elif kind == "white_axis":
        require(extra, {"inat_taxon_id", "cohort", "w_alignment_sq", "method"}, kind)
        if not extra["method"].eq("frozen_study_specific_H2").all():
            raise ValueError("W must be exported from the frozen study-specific H2 implementation")
        extra = extra[["inat_taxon_id", "cohort", "w_alignment_sq"]]
        values = pd.to_numeric(extra["w_alignment_sq"], errors="coerce")
        if ((values < 0) | (values > 1) | ~np.isfinite(values)).any():
            raise ValueError("out-of-range W alignment values")
    else:
        raise ValueError(kind)
    if extra.duplicated(["inat_taxon_id", "cohort"]).any():
        raise ValueError(f"{kind}: duplicate species-cohort keys")
    measured_keys = set(ledger.dropna(subset=["cohort"])[["inat_taxon_id", "cohort"]].itertuples(index=False, name=None))
    data_keys = set(extra[["inat_taxon_id", "cohort"]].itertuples(index=False, name=None))
    if data_keys - measured_keys:
        raise ValueError(f"{kind}: data for unmeasured species or different cohort")
    ledger = ledger.merge(extra, on=["inat_taxon_id", "cohort"], how="left", validate="one_to_one")
    if kind == "spatial":
        both = ledger["spatial_D_pair"].notna() & ledger["D_pair"].notna()
        if not np.allclose(ledger.loc[both, "spatial_D_pair"], ledger.loc[both, "D_pair"], atol=1e-9):
            raise ValueError("spatial metric D_pair differs from the measured input")
        bad = (ledger["local_pairs_50km"] < 30) & ledger["depletion_50km"].notna()
        if bad.any():
            raise ValueError("spatial metric provided despite <30 local pairs")
    return ledger


def environment_contrasts(photos: pd.DataFrame, env_path: Path) -> pd.DataFrame:
    """Return species effects only with externally frozen exposure controls.

    This layer never overwrites the September 2026 confirmatory environment gates.
    """
    if photos.empty:
        raise ValueError("environment join needs measured photos")
    d = pd.read_csv(env_path, low_memory=False)
    require(d, {"inat_taxon_id", "photo_id", "near_clip_fraction", "high_clip_excluded", *ENVIRONMENTS}, "environment")
    d = identifiers(d, "environment")
    if d["photo_id"].isna().any() or d["high_clip_excluded"].isna().any():
        raise ValueError("environment missing photo ID or high-clip exclusion flag")
    d["photo_id"] = d["photo_id"].astype(str)
    if d.duplicated(["inat_taxon_id", "photo_id"]).any():
        raise ValueError("environment duplicate (taxon,photo) keys")
    cols = ["inat_taxon_id", "photo_id", "near_clip_fraction", "high_clip_excluded", *ENVIRONMENTS]
    valid_n = photos.loc[photos["usable_colour"]].groupby("inat_taxon_id").size()
    eligible_ids = set(valid_n[valid_n >= MIN_D].index)
    rows = photos.loc[photos["usable_colour"] & photos["inat_taxon_id"].isin(eligible_ids)].merge(d[cols], how="left", on=["inat_taxon_id", "photo_id"], validate="one_to_one", indicator=True)
    # An unmatched geography/clip record is not a zero environmental effect.
    rows = rows.loc[rows["_merge"].eq("both")].copy()
    rows["near_clip_fraction"] = pd.to_numeric(rows["near_clip_fraction"], errors="coerce")
    quality = rows["near_clip_fraction"].between(0, 1) & ~bool_values(rows["high_clip_excluded"])
    rows = rows.loc[quality].copy()
    result = []
    for (cohort, taxon), g in rows.groupby(["cohort", "inat_taxon_id"], sort=True):
        for env in ENVIRONMENTS:
            gx = g.copy()
            gx[env] = pd.to_numeric(gx[env], errors="coerce")
            gx = gx.loc[np.isfinite(gx[env])].copy()
            nw = int(gx["morph"].eq("white").sum())
            nc = int(len(gx) - nw)
            if min(nw, nc) < MIN_ENV_PER_COLOUR:
                continue
            scale = float(gx[env].std(ddof=0))
            if not np.isfinite(scale) or scale <= 0:
                continue
            gx["z"] = (gx[env] - gx[env].mean()) / scale
            effect = float(gx.loc[gx["morph"].eq("white"), "z"].mean() -
                           gx.loc[gx["morph"].ne("white"), "z"].mean())
            row = {"cohort": cohort, "inat_taxon_id": taxon, "variable": env,
                   "n_white": nw, "n_nonwhite": nc, "white_minus_nonwhite_sd": effect,
                   "paired_observers": 0, "paired_observer_delta": np.nan}
            if "observer_id" in gx:
                paired = []
                for _, obs in gx.dropna(subset=["observer_id"]).groupby("observer_id"):
                    w = obs.loc[obs["morph"].eq("white"), "z"]
                    c = obs.loc[obs["morph"].ne("white"), "z"]
                    if len(w) and len(c):
                        paired.append(float(w.mean() - c.mean()))
                row["paired_observers"] = len(paired)
                if len(paired) >= MIN_PAIRED_OBSERVERS:
                    row["paired_observer_delta"] = float(np.mean(paired))
            result.append(row)
    return pd.DataFrame(result, columns=["cohort", "inat_taxon_id", "variable", "n_white", "n_nonwhite",
                                         "white_minus_nonwhite_sd", "paired_observers", "paired_observer_delta"])


def stats(x: pd.Series) -> dict:
    a = pd.to_numeric(x, errors="coerce").to_numpy(float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"n": 0, "status": "not_estimable"}
    return {"n": len(a), "median": float(np.median(a)), "mean": float(np.mean(a)),
            "fraction_positive": float(np.mean(a > 0))}


def correlation(frame: pd.DataFrame, x: str, y: str) -> dict:
    if x not in frame or y not in frame:
        return {"status": "not_estimable", "reason": "metric_missing"}
    tmp = frame[[x, y, "family"]].replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y]) if "family" in frame else frame[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(tmp) < 20 or tmp[x].nunique() < 2 or tmp[y].nunique() < 2:
        return {"status": "not_estimable", "reason": "n<20_or_constant", "n": len(tmp)}
    rho = float(spearmanr(tmp[x], tmp[y]).statistic)
    out = {"status": "exploratory", "n": len(tmp), "rho": rho}
    if "family" in tmp and tmp["family"].nunique() >= 5:
        leave = []
        for fam in sorted(tmp["family"].dropna().unique()):
            sub = tmp.loc[tmp["family"].ne(fam)]
            if len(sub) >= 20 and sub[x].nunique() >= 2 and sub[y].nunique() >= 2:
                leave.append(float(spearmanr(sub[x], sub[y]).statistic))
        out["family_loo_range"] = [min(leave), max(leave)] if leave else None
    else:
        out["family_loo_range"] = None
    return out


def summary(ledger: pd.DataFrame, env: pd.DataFrame) -> dict:
    eligible = ledger.loc[ledger["measurement_status"].eq("D_evaluable")].copy()
    answer = {"scope": "exploratory_followup_not_frozen_H2_retest",
              "n_u100": len(ledger), "n_measured": int((ledger["n_measured"] > 0).sum()),
              "n_D_evaluable": len(eligible),
              "measurement_status": {str(k): int(v) for k, v in ledger["measurement_status"].value_counts().items()},
              "cohorts": {}, "correlations": {}}
    for cohort, g in eligible.groupby("cohort"):
        answer["cohorts"][str(cohort)] = {
            "n_D_evaluable": len(g), "D": stats(g["D"]),
            "white_fraction": stats(g["white_fraction"]),
            "W_alignment_squared": stats(g["w_alignment_sq"]) if "w_alignment_sq" in g else {"n": 0, "status": "not_provided"},
            "local_depletion_50km": stats(g["depletion_50km"]) if "depletion_50km" in g else {"n": 0, "status": "not_provided"},
        }
    for x, y in (("D", "depletion_50km"), ("D", "w_alignment_sq"), ("w_alignment_sq", "depletion_50km")):
        answer["correlations"][f"{x}_vs_{y}"] = correlation(eligible, x, y)
    answer["white_environment"] = {}
    for var in ENVIRONMENTS:
        sub = env.loc[env["variable"].eq(var)]
        merged = sub.merge(eligible, how="inner", on=["inat_taxon_id", "cohort"], validate="one_to_one")
        answer["correlations"][f"D_vs_white_delta_{var}"] = correlation(merged, "D", "white_minus_nonwhite_sd")
        answer["correlations"][f"W_vs_white_delta_{var}"] = correlation(merged, "w_alignment_sq", "white_minus_nonwhite_sd")
        answer["correlations"][f"spatial_vs_white_delta_{var}"] = correlation(merged, "depletion_50km", "white_minus_nonwhite_sd")
        answer["white_environment"][var] = {
            "species_equal_weight": stats(sub["white_minus_nonwhite_sd"]),
            "paired_observer_species": stats(sub["paired_observer_delta"]),
            "cohorts": {str(c): stats(g["white_minus_nonwhite_sd"]) for c, g in sub.groupby("cohort")},
        }
    answer["limits"] = [
        "U100 is observer-capped photographic opportunity, not population polymorphism prevalence",
        "photo colour is not pigment biochemistry or genetic morph identity",
        "W is sign-invariant and valid only when imported from frozen study-specific H2",
        "environment values are long-term site covariates, not weather on photo dates",
        "photo-level white associations require independent exposure/observer/background checks",
        "no pooling of original retrospective and prospective null p-values",
        "no original manuscript/confirmatory decisions changed",
    ]
    answer["confirmatory_decisions_changed"] = False
    return answer


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--u100", type=Path, required=True)
    p.add_argument("--measured", action="append", default=[], help="cohort=path (repeat for species-disjoint tables)")
    p.add_argument("--environment", type=Path, help="exact photo/taxon joined WorldClim + frozen exposure controls")
    p.add_argument("--w-metrics", type=Path, help="study-specific per-species W squared projection, no generic refit")
    p.add_argument("--spatial-metrics", type=Path, help="study-specific 50-km per-species D_pair/depletion")
    p.add_argument("--expected-u100", type=int, default=4730)
    p.add_argument("--shards", type=int, default=24)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args(argv)
    paths = parse_cohorts(args.measured)
    frame = opportunity(args.u100, args.expected_u100)
    photos, hashes = photo_tables(paths, frame)
    led = registry(frame, photos, args.shards)
    if args.w_metrics:
        led = attach_species_metrics(led, args.w_metrics, kind="white_axis")
    if args.spatial_metrics:
        led = attach_species_metrics(led, args.spatial_metrics, kind="spatial")
    env = environment_contrasts(photos, args.environment) if args.environment else pd.DataFrame(columns=["cohort", "inat_taxon_id", "variable", "white_minus_nonwhite_sd", "paired_observer_delta"])
    report = summary(led, env)
    report["input_sha256"] = {"u100": sha256(args.u100),
        "measured": hashes,
        **({"environment": sha256(args.environment)} if args.environment else {}),
        **({"w_metrics": sha256(args.w_metrics)} if args.w_metrics else {}),
        **({"spatial_metrics": sha256(args.spatial_metrics)} if args.spatial_metrics else {})}
    args.outdir.mkdir(parents=True, exist_ok=True)
    led.to_csv(args.outdir / "u100_species_registry.csv", index=False)
    env.to_csv(args.outdir / "white_environment_species_contrasts.csv", index=False)
    (args.outdir / "comparison_summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"n_u100": report["n_u100"], "n_measured": report["n_measured"],
                      "n_D_evaluable": report["n_D_evaluable"], "output": str(args.outdir)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
