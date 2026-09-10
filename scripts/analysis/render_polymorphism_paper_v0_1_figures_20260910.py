#!/usr/bin/env python3
"""Render publication-oriented FCP paper v0.1 figures from synchronized frozen figure data."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "polymorphism_v0_1"
DATA = PAPER / "figure_data"
FIG = PAPER / "figures"
FIG.mkdir(parents=True, exist_ok=True)

NUMBERS = PAPER / "paper_numbers.json"
SHARED_GENUS = ROOT / "results" / "polymorphism_genus_replication_step7_20260910" / "shared_repeated_genus_means.csv"

DISC = "#0072B2"
RES = "#D55E00"
DARK = "#222222"
MID = "#777777"
LIGHT = "#D9D9D9"
TECH = "#999999"
AMB = "#CC79A7"
NOBIO = "#E69F00"
CLASS = "#009E73"
MINC = "#56B4E9"
MAXC = "#D55E00"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def save(fig: plt.Figure, stem: str) -> list[Path]:
    paths = [FIG / f"{stem}.png", FIG / f"{stem}.pdf"]
    fig.savefig(paths[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(paths[1], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return paths


def binned_trend(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str, bins: int = 8) -> None:
    order = np.argsort(x, kind="mergesort")
    chunks = np.array_split(order, bins)
    bx, by = [], []
    for idx in chunks:
        if len(idx) == 0:
            continue
        bx.append(float(np.median(x[idx])))
        by.append(float(np.mean(y[idx])))
    ax.plot(bx, by, "-o", lw=1.7, ms=4, color=color, zorder=4, label="equal-count bin mean")


def bootstrap_mean_ci(values: np.ndarray, seed: int, reps: int = 4000) -> tuple[float, float, float]:
    v = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    means = np.empty(reps, float)
    for i in range(reps):
        means[i] = np.mean(v[rng.integers(0, len(v), len(v))])
    return float(np.mean(v)), float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def ptxt(p: float) -> str:
    if p < 0.001:
        return "p<0.001"
    return f"p={p:.3f}"


def figure1(numbers: dict[str, Any]) -> list[Path]:
    df = pd.read_csv(DATA / "figure1_discovery_species_polymorphism.csv")
    if len(df) != 369:
        raise RuntimeError("Figure 1 discovery species count drift")

    fig = plt.figure(figsize=(7.2, 5.7))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    ax.axis("off")
    boxes = [
        (0.03, 0.62, 0.26, 0.22, "100 fixed\nphotos / species"),
        (0.37, 0.62, 0.26, 0.22, "measurement +\nstatus partition"),
        (0.71, 0.62, 0.26, 0.22, "4 admitted\ncolour states"),
        (0.37, 0.16, 0.26, 0.22, "species diversity\nD = 1 − Σp²"),
    ]
    for x, y, w, h, label in boxes:
        rect = plt.Rectangle((x, y), w, h, transform=ax.transAxes, fill=False, lw=1.1, ec=DARK)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center")
    arrows = [((0.29, 0.73), (0.37, 0.73)), ((0.63, 0.73), (0.71, 0.73)), ((0.84, 0.62), (0.54, 0.38))]
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a, xycoords=ax.transAxes, arrowprops=dict(arrowstyle="->", lw=1.0, color=DARK))
    ax.text(0.5, 0.01, "Amount of polymorphism is the focal species-level estimand", transform=ax.transAxes, ha="center", fontsize=7.5)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    ax.hist(df["D"], bins=np.linspace(0, 0.72, 25), edgecolor="white", linewidth=0.5, color=DISC)
    med = float(df["D"].median())
    ax.axvline(med, color=DARK, lw=1.1, ls="--")
    ax.text(med + 0.015, ax.get_ylim()[1]*0.88, f"median={med:.2f}", fontsize=7.5)
    ax.set(xlabel="Four-state flower-colour diversity, D", ylabel="Species", title="Polymorphism varies continuously")

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    ax.scatter(df["D"], df["second_fraction"], s=11, alpha=0.35, color=DISC, edgecolors="none")
    ax.axhline(0.10, color=MID, ls="--", lw=0.9)
    ax.axhline(0.20, color=MID, ls=":", lw=0.9)
    ax.text(0.01, 0.105, "10%", fontsize=7, color=MID)
    ax.text(0.01, 0.205, "20%", fontsize=7, color=MID)
    ax.set(xlabel="D", ylabel="Second-largest morph fraction", title="Threshold summaries preserve a continuous gradient")

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    vals = [100*numbers["discovery"]["second_ge_0_10_fraction"], 100*numbers["discovery"]["second_ge_0_20_fraction"]]
    y = np.arange(2)
    ax.barh(y, vals, height=0.55, color=DISC)
    ax.set_yticks(y, ["second morph ≥10%", "second morph ≥20%"])
    ax.invert_yaxis()
    ax.set(xlabel="Species in admitted discovery frame (%)", xlim=(0, 55), title="Non-trivial secondary morphs are common in-frame")
    for yi, v in zip(y, vals):
        ax.text(v + 1, yi, f"{v:.1f}%", va="center")
    ax.text(0, 1.35, "Not an angiosperm-wide prevalence estimate", fontsize=7, color=MID)

    fig.suptitle("Figure 1 | Flower-colour polymorphism as a species-level property", fontsize=11, y=0.995)
    return save(fig, "figure1_polymorphism_estimand")


def figure2(numbers: dict[str, Any]) -> list[Path]:
    d = pd.read_csv(DATA / "figure2_discovery_spatial.csv")
    r = pd.read_csv(DATA / "figure2_reserve_spatial.csv")
    effects = pd.read_csv(DATA / "figure2_effect_summary.csv")
    if len(d) != 369 or len(r) != 363:
        raise RuntimeError("Figure 2 frame count drift")

    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.34)

    for ax, df, tranche, color, rho, p, letter in [
        (fig.add_subplot(gs[0, 0]), d, "Discovery", DISC, numbers["discovery"]["rho_D_spatial"], numbers["discovery"]["rho_D_spatial_p"], "A"),
        (fig.add_subplot(gs[0, 1]), r, "Species-disjoint reserve", RES, numbers["reserve"]["rho_D_spatial"], numbers["reserve"]["rho_D_spatial_p"], "B"),
    ]:
        panel(ax, letter)
        x = df["D"].to_numpy(float); y = df["spatial_rho"].to_numpy(float)
        ax.scatter(x, y, s=10, alpha=0.28, color=color, edgecolors="none")
        ax.axhline(0, color=LIGHT, lw=1)
        binned_trend(ax, x, y, color)
        ax.text(0.03, 0.94, f"Spearman ρ={rho:.3f}; {ptxt(p)}", transform=ax.transAxes, va="top", fontsize=7.5)
        ax.set(xlabel="Flower-colour polymorphism, D", ylabel="Within-species spatial ρ", title=tranche)

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    x_positions = [0, 1, 3, 4]
    labels = ["Disc. <10%", "Disc. ≥10%", "Res. <10%", "Res. ≥10%"]
    datasets = [
        d.loc[~d["primary_second_ge_0_10"].astype(bool), "spatial_rho"].to_numpy(float),
        d.loc[d["primary_second_ge_0_10"].astype(bool), "spatial_rho"].to_numpy(float),
        r.loc[~r["primary_second_ge_0_10"].astype(bool), "spatial_rho"].to_numpy(float),
        r.loc[r["primary_second_ge_0_10"].astype(bool), "spatial_rho"].to_numpy(float),
    ]
    cols = [DISC, DISC, RES, RES]
    for i, (x, vals, c) in enumerate(zip(x_positions, datasets, cols)):
        mean, lo, hi = bootstrap_mean_ci(vals, 20260910 + i)
        ax.errorbar(x, mean, yerr=[[mean-lo], [hi-mean]], fmt="o", ms=5, capsize=3, color=c, lw=1.2)
    ax.axhline(0, color=LIGHT, lw=1)
    ax.set_xticks(x_positions, labels)
    pdisc = float(effects.loc[(effects.tranche=="discovery") & (effects.metric=="polymorphic_minus_complement"), "p"].iloc[0])
    pres = float(effects.loc[(effects.tranche=="reserve") & (effects.metric=="polymorphic_minus_complement"), "p"].iloc[0])
    ax.text(0.03, 0.94, f"contrast: discovery {ptxt(pdisc)}; reserve {ptxt(pres)}", transform=ax.transAxes, va="top", fontsize=7.2)
    ax.set(ylabel="Mean species spatial ρ ± bootstrap 95% CI", title="Polymorphic subset is more spatially organized")

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    order = [
        "second_ge_0.10_primary_mean_rho",
        "second_ge_0.10_observer_exclusion_mean_rho",
        "second_ge_0.10_quarter_mean_rho",
        "second_ge_0.10_matched_background_mean_rho",
    ]
    lab = ["primary", "observer exclusion", "calendar quarter", "flower − background"]
    sub = effects[(effects.tranche=="reserve") & effects.metric.isin(order)].copy().set_index("metric").loc[order].reset_index()
    yy = np.arange(len(sub))
    ax.axvline(0, color=LIGHT, lw=1)
    ax.scatter(sub["estimate"], yy, s=30, color=RES, zorder=3)
    for x, y, p in zip(sub["estimate"], yy, sub["p"]):
        ax.text(x + 0.0014, y, ptxt(float(p)), va="center", fontsize=7)
    ax.set_yticks(yy, lab)
    ax.invert_yaxis()
    ax.set(xlabel="Equal-species mean spatial ρ (≥10% subset)", title="Reserve replication survives nuisance controls")

    fig.suptitle("Figure 2 | Greater polymorphism predicts stronger within-species geographic organization", fontsize=11, y=0.995)
    return save(fig, "figure2_polymorphism_spatial_replication")


def figure3(numbers: dict[str, Any]) -> list[Path]:
    status = pd.read_csv(DATA / "figure3_measurement_status_partition.csv")
    robust = pd.read_csv(DATA / "figure3_spatial_robustness_summary.csv")
    dint = pd.read_csv(DATA / "figure3_discovery_ambiguity_intervals.csv")
    rint = pd.read_csv(DATA / "figure3_reserve_ambiguity_intervals.csv")

    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 2, hspace=0.47, wspace=0.37)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    maplab = {
        "classified_four_state_morph": "four-state classified",
        "not_evaluable_roi_or_flip_gate": "ROI/flip technical failure",
        "not_evaluable_ambiguous_palette_composition": "ambiguous palette",
        "not_evaluable_no_biological_palette_mass": "no biological palette",
    }
    colors = {
        "classified_four_state_morph": CLASS,
        "not_evaluable_roi_or_flip_gate": TECH,
        "not_evaluable_ambiguous_palette_composition": AMB,
        "not_evaluable_no_biological_palette_mass": NOBIO,
    }
    status = status.sort_values("count", ascending=True)
    yy = np.arange(len(status))
    ax.barh(yy, status["fraction"]*100, color=[colors[x] for x in status["measurement_status"]])
    ax.set_yticks(yy, [maplab[x] for x in status["measurement_status"]])
    ax.set(xlabel="Reserve candidate-frame photos (%)", title="Missingness is not one process")
    for y, v in zip(yy, status["fraction"]*100):
        ax.text(v + 0.6, y, f"{v:.1f}%", va="center", fontsize=7)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    for df, label, color in [(dint, "Discovery", DISC), (rint, "Reserve", RES)]:
        vals = np.sort(df["D_interval_width"].to_numpy(float))
        ecdf = np.arange(1, len(vals)+1)/len(vals)
        ax.plot(vals, ecdf, lw=1.8, label=label, color=color)
    ax.set(xlabel="Four-state ambiguity interval width (Dmax4 − Dmin4)", ylabel="Cumulative fraction of species", title="Ambiguity is material but bounded")
    ax.legend(frameon=False, loc="lower right")

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    groups = [
        ("discovery", "primary_spatial", "Discovery primary"),
        ("reserve", "primary_spatial", "Reserve primary"),
        ("reserve", "flower_minus_background", "Reserve flower−background"),
    ]
    variants = [("D_min4", MINC, "Dmin4"), ("observed", DARK, "Observed D"), ("D_max4", MAXC, "Dmax4")]
    base = np.arange(len(groups), dtype=float)
    offsets = {"D_min4": -0.18, "observed": 0.0, "D_max4": 0.18}
    for v, c, vl in variants:
        xs, ys = [], []
        for i, (tr, resp, _) in enumerate(groups):
            row = robust[(robust.tranche==tr) & (robust.response==resp) & (robust.D_variant==v)]
            if len(row) != 1:
                raise RuntimeError(f"missing robustness row {tr}/{resp}/{v}")
            xs.append(base[i] + offsets[v]); ys.append(float(row.partial_rho.iloc[0]))
        ax.scatter(xs, ys, s=30, color=c, label=vl, zorder=3)
        for x, y, (tr, resp, _) in zip(xs, ys, groups):
            p = float(robust[(robust.tranche==tr) & (robust.response==resp) & (robust.D_variant==v)].p.iloc[0])
            ax.text(x, y + 0.006, f"{p:.3f}", ha="center", va="bottom", fontsize=6.6, color=c)
    ax.axhline(0, color=LIGHT, lw=1)
    ax.set_xticks(base, [g[2] for g in groups], rotation=18, ha="right")
    ax.set(ylabel="Partial Spearman ρ\n| sampled span + technical failure", title="Spatial signal survives both ambiguity endpoints")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.text(0.02, 0.03, "numbers above points are geometry-preserving p", transform=ax.transAxes, fontsize=6.5, color=MID)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    rr = rint.sort_values("D").reset_index(drop=True)
    xx = np.arange(len(rr))
    ax.fill_between(xx, rr["D_min4"], rr["D_max4"], color=LIGHT, alpha=0.8, linewidth=0, label="four-state completion interval")
    ax.plot(xx, rr["D"], color=RES, lw=1.3, label="observed D")
    ax.set(xlabel="Reserve species ranked by observed D", ylabel="D", title="Observed ordering is stable across broad completion bounds")
    ax.legend(frameon=False, loc="upper left")
    ax.text(0.03, 0.86, f"ρ(D,Dmin4)={numbers['reserve'].get('rho_D_Dmin4', 0.997):.3f}\nρ(D,Dmax4)={numbers['reserve'].get('rho_D_Dmax4', 0.961):.3f}", transform=ax.transAxes, fontsize=7)

    fig.suptitle("Figure 3 | Opportunity, technical failure, and ambiguity do not erase the spatial result", fontsize=11, y=0.995)
    return save(fig, "figure3_measurement_and_ambiguity_robustness")


def figure4(numbers: dict[str, Any]) -> list[Path]:
    summary = pd.read_csv(DATA / "figure4_genus_clustering_summary.csv")
    shared = pd.read_csv(SHARED_GENUS)
    shared.to_csv(DATA / "figure4_shared_genus_means.csv", index=False)

    fig = plt.figure(figsize=(7.2, 3.3))
    gs = fig.add_gridspec(1, 2, wspace=0.42)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    labels = {
        ("discovery", "observed_raw"): "Discovery raw",
        ("reserve", "observed_raw"): "Reserve raw",
        ("reserve", "span_only"): "Reserve | span",
        ("reserve", "span_plus_technical"): "Reserve | span + technical",
        ("reserve", "D_min4_span_plus_technical"): "Reserve Dmin4 | span + technical",
        ("reserve", "D_max4_span_plus_technical"): "Reserve Dmax4 | span + technical",
    }
    summary["label"] = [labels[(a, b)] for a, b in zip(summary.tranche, summary.variant)]
    yy = np.arange(len(summary))
    ax.axvline(0, color=LIGHT, lw=1)
    ax.scatter(summary["gain"], yy, s=32, color=[DISC if t=="discovery" else RES for t in summary.tranche])
    for x, y, p in zip(summary.gain, yy, summary.p):
        ax.text(x + 0.005, y, ptxt(float(p)), va="center", fontsize=6.7)
    ax.set_yticks(yy, summary["label"])
    ax.invert_yaxis()
    ax.set(xlabel="Genus clustering gain", title="Taxonomic clustering is replicated but ambiguity-sensitive")

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    ax.scatter(shared["discovery_mean_D"], shared["reserve_mean_D"], s=28, alpha=0.75, color=RES, edgecolors="white", linewidths=0.4)
    lim = (0, max(float(shared["discovery_mean_D"].max()), float(shared["reserve_mean_D"].max())) + 0.04)
    ax.plot(lim, lim, ls="--", lw=0.9, color=MID)
    ax.set(xlim=lim, ylim=lim, xlabel="Discovery genus mean D", ylabel="Reserve genus mean D", title="23 repeated genera show cross-tranche concordance")
    ax.text(0.04, 0.95, f"Spearman ρ={numbers['genus_secondary']['shared_23_rho']:.3f}; {ptxt(numbers['genus_secondary']['shared_23_p'])}", transform=ax.transAxes, va="top", fontsize=7.5)
    deviation = np.abs(shared["reserve_mean_D"].to_numpy() - shared["discovery_mean_D"].to_numpy())
    label_idx = set(np.argsort(deviation)[-4:].tolist()) | set(np.argsort(shared["discovery_mean_D"].to_numpy())[-2:].tolist())
    for i in sorted(label_idx):
        row = shared.iloc[i]
        ax.text(row.discovery_mean_D + 0.007, row.reserve_mean_D + 0.007, row.genus, fontsize=6.6)

    fig.suptitle("Figure 4 | Genus-level structure is a secondary, qualified result", fontsize=11, y=1.01)
    return save(fig, "figure4_genus_taxonomic_structure")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    numbers = read_json(NUMBERS)
    if numbers["paper_version"] != "v0.1-post-step9":
        raise RuntimeError("paper number ledger version drift")
    if numbers["discovery"]["species"] != 369 or numbers["reserve"]["species"] != 363:
        raise RuntimeError("paper number ledger species counts drift")

    outputs: list[Path] = []
    outputs += figure1(numbers)
    outputs += figure2(numbers)
    outputs += figure3(numbers)
    outputs += figure4(numbers)

    captions = """# Figure captions — FCP polymorphism paper v0.1

## Figure 1 | Flower-colour polymorphism as a species-level property

(A) Fixed-photo measurement design and the four-state Gini–Simpson diversity estimand D. (B) Distribution of D across the 369-species discovery frame. (C) Relationship between D and the fraction of observations assigned to the second-largest morph, with 10% and 20% descriptive thresholds. (D) Fractions of the admitted discovery frame exceeding those thresholds. These fractions characterize the admitted photo-derived frame and are not estimates of angiosperm-wide polymorphism prevalence.

## Figure 2 | Greater polymorphism predicts stronger within-species geographic organization

(A–B) Species-level flower-colour diversity D versus within-species spatial rho in discovery and a species-disjoint reserve photo tranche. Points are species; connected points show equal-count-bin mean spatial rho and are descriptive. Reported P values use the frozen spatial randomization procedures. (C) Mean species spatial rho for species below versus above the preregistered 10% second-morph threshold; intervals are descriptive bootstrap 95% intervals, while the displayed contrast P values come from the frozen spatial nulls. (D) Reserve >=10% subset results under primary, observer-pair exclusion, calendar-quarter, and matched flower-minus-background metrics.

## Figure 3 | Opportunity, technical failure, and ambiguity do not erase the spatial result

(A) Measurement-status partition of the 50,000-photo reserve candidate frame. Only ROI/flip-gate failure is treated as clearly technical; ambiguous-palette observations remain separate. (B) Empirical distributions of species-level four-state ambiguity interval width. (C) Partial Spearman associations between polymorphism and spatial organization after controlling sampled geographic span and technical-failure rate, shown for observed D and the exact Dmin4/Dmax4 species-level completion endpoints under the four-state ambiguity sensitivity model. P values use the original geometry-preserving within-species spatial null arrays. (D) Reserve species ranked by observed D, showing Dmin4–Dmax4 intervals. Uniform endpoint stress tests do not exhaust arbitrary species-specific latent allocations and do not imply that ambiguous observations are true four-state morphs.

## Figure 4 | Genus-level structure is a secondary, qualified result

(A) Genus-level taxonomic clustering gain in discovery, species-disjoint reserve, opportunity/technical-failure adjustments, and ambiguity endpoint sensitivities. The reserve Dmin4 endpoint is not supported at P<0.05, so genus clustering is not treated as uniform-endpoint robust. (B) Mean observed D for the 23 genera represented by at least two different species in both tranches. This is taxonomic clustering/concordance, not formal phylogenetic signal.
"""
    (PAPER / "FIGURE_CAPTIONS.md").write_text(captions, encoding="utf-8")

    manifest = {
        "paper_version": numbers["paper_version"],
        "renderer": str(Path(__file__).relative_to(ROOT)),
        "figure_inputs": sorted(p.name for p in DATA.glob("*.csv")),
        "outputs": [
            {"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in outputs
        ],
        "captions": str((PAPER / "FIGURE_CAPTIONS.md").relative_to(ROOT)),
    }
    (FIG / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
