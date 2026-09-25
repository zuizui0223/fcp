#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

NODE_PATTERNS = {
    "MYB": re.compile(r"(?:\b(?:MYB90|MYB75|MYB113|MYB114|PAP1|PAP2)\b|R2R3[- _]?MYB)", re.I),
    "CHS": re.compile(r"(?:\bCHS\b|chalcone synthase)", re.I),
    "F3H": re.compile(r"(?:\bF3H\b|flavanone[- _]?3[- _]?hydroxylase)", re.I),
    "ANS": re.compile(r"(?:\bANS\b|anthocyanidin synthase)", re.I),
    "UFGT_UGT": re.compile(r"(?:\bUFGT\b|\bU78D2\b|\bU75C1\b)", re.I),
    "FLS": re.compile(r"(?:\bFLS\b|flavonol synthase)", re.I),
    "bHLH": re.compile(r"(?:\bbHLH\b|basic helix[- _]?loop[- _]?helix)", re.I),
    "WD40_TTG1": re.compile(r"(?:\bWD40\b|\bTTG1\b|transparent testa glabra\s*1)", re.I),
}

CONTEXT_PATTERNS = {
    "PAL": re.compile(r"(?:\bPAL\b|phenylalanine ammonia[- _]?lyase)", re.I),
    "4CL": re.compile(r"(?:\b4CL\b|4[- _]?coumarate.*CoA ligase)", re.I),
    "CHI": re.compile(r"(?:\bCHI\b|chalcone isomerase)", re.I),
    "DFR": re.compile(r"(?:\bDFR\b|dihydroflavonol.*reductase)", re.I),
}

LOGFC_RX = re.compile(r"(?:log\s*2?\s*(?:fold.?change|fc)|fold.?change|logfc)", re.I)
FDR_RX = re.compile(r"(?:fdr|adj(?:usted)?[ ._-]*p|padj|q[ ._-]*value)", re.I)
P_RX = re.compile(r"(?:^p$|p[ ._-]*value|pval)", re.I)


def safe_float(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip().replace(",", "")
    if not s:
        return np.nan
    s = re.sub(r"^[<>~= ]+", "", s)
    try:
        return float(s)
    except Exception:
        return np.nan


def stringify_row(row):
    return ["" if pd.isna(x) else str(x) for x in row.tolist()]


def nearest_header(df: pd.DataFrame, row_idx: int):
    for h in range(row_idx - 1, -1, -1):
        vals = stringify_row(df.iloc[h])
        joined = " | ".join(vals)
        score = int(bool(LOGFC_RX.search(joined))) + int(bool(FDR_RX.search(joined))) + int(bool(P_RX.search(joined)))
        if score > 0:
            return h, vals
    return None


def locate_col(headers, rx):
    hits = [i for i, h in enumerate(headers) if rx.search(str(h))]
    return hits[0] if hits else None


def load_workbook(path: Path):
    engine = "xlrd" if path.suffix.lower() == ".xls" else "openpyxl"
    return pd.read_excel(path, sheet_name=None, header=None, engine=engine)


def scan_workbook(path: Path):
    source_label = path.stem
    sheets = load_workbook(path)
    manifest = []
    matches = []
    all_patterns = {**NODE_PATTERNS, **CONTEXT_PATTERNS}
    for sheet, df in sheets.items():
        manifest.append({
            "source": source_label,
            "file_name": path.name,
            "sheet": str(sheet),
            "rows": int(df.shape[0]),
            "cols": int(df.shape[1]),
        })
        for ridx in range(len(df)):
            vals = stringify_row(df.iloc[ridx])
            joined = " | ".join(vals)
            hit_nodes = [name for name, rx in all_patterns.items() if rx.search(joined)]
            if not hit_nodes:
                continue
            hdr = nearest_header(df, ridx)
            header_row = None
            headers = [""] * len(vals)
            if hdr is not None:
                header_row, headers = hdr
            log_col = locate_col(headers, LOGFC_RX)
            fdr_col = locate_col(headers, FDR_RX)
            p_col = locate_col(headers, P_RX)
            logfc = safe_float(vals[log_col]) if log_col is not None and log_col < len(vals) else np.nan
            fdr = safe_float(vals[fdr_col]) if fdr_col is not None and fdr_col < len(vals) else np.nan
            pval = safe_float(vals[p_col]) if p_col is not None and p_col < len(vals) else np.nan
            for node in hit_nodes:
                matches.append({
                    "source": source_label,
                    "file_name": path.name,
                    "sheet": str(sheet),
                    "excel_row": int(ridx + 1),
                    "header_excel_row": int(header_row + 1) if header_row is not None else np.nan,
                    "node": node,
                    "panel": "FCP" if node in NODE_PATTERNS else "context",
                    "logFC_summer_vs_spring": logfc,
                    "FDR_or_adjusted_p": fdr,
                    "p_value": pval,
                    "row_text": joined[:8000],
                    "header_text": " | ".join(headers)[:5000],
                })
    return pd.DataFrame(manifest), pd.DataFrame(matches)


def summarize(matches: pd.DataFrame):
    out = []
    for node in NODE_PATTERNS:
        q = matches[(matches.panel == "FCP") & (matches.node == node)].copy() if len(matches) else pd.DataFrame()
        vals = pd.to_numeric(q.logFC_summer_vs_spring, errors="coerce").dropna().to_numpy(float) if len(q) else np.array([])
        if len(vals):
            med = float(np.median(vals))
            direction = "summer_down" if med < 0 else "summer_up" if med > 0 else "zero"
        else:
            med = np.nan
            direction = "direction_not_machine_readable" if len(q) else "not_found"
        sig_down = False
        if len(q):
            lv = pd.to_numeric(q.logFC_summer_vs_spring, errors="coerce")
            fv = pd.to_numeric(q.FDR_or_adjusted_p, errors="coerce")
            sig_down = bool(((lv < 0) & (fv < 0.05)).fillna(False).any())
        out.append({
            "node": node,
            "matched_rows": int(len(q)),
            "machine_readable_logFC_rows": int(len(vals)),
            "median_logFC_summer_vs_spring": med,
            "node_direction": direction,
            "any_significant_summer_down": sig_down,
        })
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workbooks", nargs="+", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    manifests = []
    match_frames = []
    for raw in a.workbooks:
        path = Path(raw)
        m, x = scan_workbook(path)
        manifests.append(m)
        if len(x):
            match_frames.append(x)

    manifest = pd.concat(manifests, ignore_index=True)
    matches = pd.concat(match_frames, ignore_index=True) if match_frames else pd.DataFrame(
        columns=["source","file_name","sheet","excel_row","header_excel_row","node","panel","logFC_summer_vs_spring","FDR_or_adjusted_p","p_value","row_text","header_text"]
    )

    manifest.to_csv(out / "workbook_manifest.csv", index=False)
    matches.to_csv(out / "all_pathway_matches.csv", index=False)

    summary = summarize(matches)
    summary.to_csv(out / "fcp_node_summary.csv", index=False)

    found = int((summary.matched_rows > 0).sum())
    readable = int((summary.machine_readable_logFC_rows > 0).sum())
    down = int(summary.node_direction.eq("summer_down").sum())
    sig_down = int(summary.any_significant_summer_down.sum())

    result = {
        "schema": "fcp_moricandia_public_pathway_reanalysis_v1",
        "status": "complete",
        "role": "post_publication_reproducibility_and_mechanistic_triangulation",
        "diagnostic_note": "The initial two-workbook pass found no FCP nodes because Supplementary Dataset 2 was a GO-enrichment table. The all-workbook rerun located the DEG table, and this corrected pass fixes two parser issues without changing the frozen 8-node panel: persistent header detection and restriction of broad MYB/UGT aliases to pigment-relevant classes.",
        "fcp_node_panel": list(NODE_PATTERNS),
        "fcp_nodes_total": len(NODE_PATTERNS),
        "fcp_nodes_found_in_public_workbooks": found,
        "fcp_nodes_with_machine_readable_logFC": readable,
        "fcp_nodes_with_median_summer_down": down,
        "fcp_nodes_with_any_significant_summer_down_transcript": sig_down,
        "comparison_semantics": "summer relative to spring; negative logFC means lower in summer-white flowers when the workbook exposes a fold-change column",
        "workbooks": [Path(x).name for x in a.workbooks],
        "bioproject": "PRJNA604514",
        "hard_nonclaims": [
            "not an outcome-blind confirmation because the paper already reports anthocyanin-pathway repression",
            "temperature and photoperiod changed together in the paired RNA-seq comparison",
            "does not rescue the failed cross-cohort BIO5 replication",
            "does not establish identical orthology or causal nucleotide changes across FCP systems",
        ],
    }
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
