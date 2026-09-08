#!/usr/bin/env python3
"""Spatial recurrence + flower/background decomposition for recurrent RGFCA modes.

Requires the upstream recurrence runner/results frozen on 2026-09-09.
No ecological covariates are used here.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

NREP_DEFAULT = 200
OBS_SEED_DEFAULT = 20260909
NULL_SEED_DEFAULT = 20260910
REPORT_CELL_M = 500_000.0
MIN_REPORT_SPECIES = 10
MAIN_SUPPORT_FRACTION = 0.80


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--recurrence-runner', type=Path, default=Path('scripts/analysis/run_rgfca_global_signal_recurrence_20260909.py'))
    p.add_argument('--input', type=Path, default=Path('data/derived/rgfca_reserve_replication_measured_photos_v1.csv'))
    p.add_argument('--reference-loadings', type=Path, default=Path('results/rgfca_global_signal_recurrence_20260909/reference_loadings.csv'))
    p.add_argument('--upstream-realizations', type=Path, default=Path('results/rgfca_global_signal_recurrence_20260909/recurrence_realizations.csv'))
    p.add_argument('--outdir', type=Path, default=Path('results/rgfca_global_signal_spatial_component_20260909'))
    p.add_argument('--reps', type=int, default=NREP_DEFAULT)
    p.add_argument('--seed', type=int, default=OBS_SEED_DEFAULT)
    p.add_argument('--null-seed', type=int, default=NULL_SEED_DEFAULT)
    return p.parse_args()


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location('rgfca_recurrence_runner', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def weighted_corr(a, b, w):
    sw = w.sum()
    ma = (w * a).sum() / sw
    mb = (w * b).sum() / sw
    da = a - ma
    db = b - mb
    den = np.sqrt((w * da * da).sum() * (w * db * db).sum())
    return float((w * da * db).sum() / den) if den > 0 else np.nan


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    mod = load_module(args.recurrence_runner)
    df, sp, obs, d_values, x, y = mod.load_measurements(args.input)
    f_values = np.column_stack([
        df[f'palette_count_{c}'].to_numpy(float) / df['flower_total12'].to_numpy(float)
        for c in mod.PALETTE
    ])
    b_values = np.column_stack([
        df[f'background_palette_count_{c}'].to_numpy(float) / df['background_total12'].to_numpy(float)
        for c in mod.PALETTE
    ])
    if np.max(np.abs(d_values - (f_values - b_values))) > 1e-12:
        raise RuntimeError('photo-level D=F-B identity failed')

    groups_for_shift, _, _ = mod.make_engine(sp, obs, d_values, x, y)
    reference = pd.read_csv(args.reference_loadings, index_col=0).to_numpy(float)[:, :5]
    permutations = list(itertools.permutations(range(5)))
    obs_rng = np.random.default_rng(args.seed)
    null_rng = np.random.default_rng(args.null_seed)

    def state_vectors_triplet(inv, eligible):
        order = np.argsort(inv, kind='stable')
        sorted_inv = inv[order]
        starts = np.r_[0, np.flatnonzero(np.diff(sorted_inv)) + 1]
        ends = np.r_[starts[1:], len(order)]
        state_ids, ds, fs, bs = [], [], [], []
        for start, end in zip(starts, ends):
            state_id = int(sorted_inv[start])
            if not eligible[state_id]:
                continue
            idx = order[start:end]
            state_obs = obs[idx]
            unique_obs = np.unique(state_obs)
            drawn_obs = obs_rng.choice(unique_obs, size=len(unique_obs), replace=True)
            u = obs_rng.random((len(unique_obs), len(idx)))
            u[drawn_obs[:, None] != state_obs[None, :]] = np.inf
            chosen = idx[np.argmin(u, axis=1)]
            state_ids.append(state_id)
            ds.append(d_values[chosen].mean(axis=0))
            fs.append(f_values[chosen].mean(axis=0))
            bs.append(b_values[chosen].mean(axis=0))
        return np.asarray(state_ids, np.int32), np.vstack(ds), np.vstack(fs), np.vstack(bs)

    def decompose(state_ids, dm, fm, bm, xc, yc, state_species):
        species_id = state_species[state_ids]
        order = np.argsort(species_id, kind='stable')
        species_id = species_id[order]
        state_ids = state_ids[order]
        dm, fm, bm = dm[order], fm[order], bm[order]
        starts = np.r_[0, np.flatnonzero(np.diff(species_id)) + 1]
        ends = np.r_[starts[1:], len(species_id)]
        covariance = np.zeros((12, 12))
        admitted = []
        n_species = 0
        n_states = 0
        for start, end in zip(starts, ends):
            if end - start < 2:
                continue
            ids = state_ids[start:end]
            xx, yy = xc[ids], yc[ids]
            dx = xx[:, None] - xx[None, :]
            dy = yy[:, None] - yy[None, :]
            if np.sqrt(dx * dx + dy * dy).max() < mod.SEPARATION_M:
                continue
            zd = dm[start:end] - dm[start:end].mean(axis=0, keepdims=True)
            zf = fm[start:end] - fm[start:end].mean(axis=0, keepdims=True)
            zb = bm[start:end] - bm[start:end].mean(axis=0, keepdims=True)
            if np.max(np.abs(zd - (zf - zb))) > 1e-12:
                raise RuntimeError('state-level D=F-B identity failed')
            k = end - start
            covariance += (zd.T @ zd) / k
            admitted.append((int(species_id[start]), ids.copy(), zd, zf, zb))
            n_species += 1
            n_states += k
        covariance /= n_species
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        explained = eigenvalues / eigenvalues.sum()
        return n_species, n_states, eigenvectors, explained, admitted

    def aggregate(frame, cols):
        a = frame.groupby(['rx', 'ry', 'spid'], as_index=False)[cols].mean()
        n = a.groupby(['rx', 'ry'])['spid'].nunique().rename('n_species').reset_index()
        m = a.groupby(['rx', 'ry'], as_index=False)[cols].mean()
        out = m.merge(n, on=['rx', 'ry'])
        return out[out['n_species'] >= MIN_REPORT_SPECIES].copy()

    observed_cell_rows = []
    null_cell_rows = []
    decomposition_rows = []
    upstream_check_rows = []

    for rep in range(1, args.reps + 1):
        sx = float(obs_rng.uniform(0, mod.CELL_M))
        sy = float(obs_rng.uniform(0, mod.CELL_M))
        inv, eligible, xc, yc, state_species = groups_for_shift(sx, sy)
        state_ids, dm, fm, bm = state_vectors_triplet(inv, eligible)
        n_species, n_states, evec, explained, admitted = decompose(
            state_ids, dm, fm, bm, xc, yc, state_species
        )
        first5 = evec[:, :5]
        similarity = np.abs(reference.T @ first5)
        perm = max(permutations, key=lambda p: sum(similarity[i, p[i]] for i in range(5)))
        matched = np.column_stack([first5[:, perm[i]] for i in range(5)])
        cosines = []
        for i in range(5):
            dot = float(reference[:, i] @ matched[:, i])
            if dot < 0:
                matched[:, i] *= -1
                dot = -dot
            cosines.append(dot)
        upstream_check_rows.append({
            'rep': rep,
            'n_species': n_species,
            'n_states': n_states,
            **{f'M{i+1}_cosine': cosines[i] for i in range(5)},
            **{f'M{i+1}_matched_run_rank': perm[i] + 1 for i in range(5)},
            **{f'M{i+1}_explained': explained[perm[i]] for i in range(5)},
        })

        all_d, all_f, all_b, all_w = [], [], [], []
        spatial_obs = []
        spatial_null = []
        for spid, ids, zd, zf, zb in admitted:
            sd = zd @ matched
            sf = zf @ matched
            sb = zb @ matched
            k = len(ids)
            all_d.append(sd)
            all_f.append(sf)
            all_b.append(sb)
            all_w.append(np.full(k, 1.0 / k))
            null_scores = sd[null_rng.permutation(k)]
            for j, state_id in enumerate(ids):
                rx = int(np.floor(xc[state_id] / REPORT_CELL_M))
                ry = int(np.floor(yc[state_id] / REPORT_CELL_M))
                ro = {'spid': spid, 'rx': rx, 'ry': ry}
                rn = {'spid': spid, 'rx': rx, 'ry': ry}
                for m in range(5):
                    ro[f'D{m+1}'] = sd[j, m]
                    rn[f'D{m+1}'] = null_scores[j, m]
                    if m < 3:
                        ro[f'F{m+1}'] = sf[j, m]
                        ro[f'B{m+1}'] = sb[j, m]
                spatial_obs.append(ro)
                spatial_null.append(rn)

        sd = np.vstack(all_d)
        sf = np.vstack(all_f)
        sb = np.vstack(all_b)
        w = np.concatenate(all_w)
        for m in range(3):
            d, f, b = sd[:, m], sf[:, m], sb[:, m]
            decomposition_rows.append({
                'rep': rep,
                'mode': f'M{m+1}',
                'corr_D_F': weighted_corr(d, f, w),
                'corr_D_negB': weighted_corr(d, -b, w),
                'mean_abs_D': np.average(np.abs(d), weights=w),
                'mean_abs_F': np.average(np.abs(f), weights=w),
                'mean_abs_B': np.average(np.abs(b), weights=w),
                'frac_absF_gt_absB': np.average(np.abs(f) > np.abs(b), weights=w),
                'frac_absB_gt_absF': np.average(np.abs(b) > np.abs(f), weights=w),
                'max_identity_error': float(np.max(np.abs(d - (f - b)))),
            })

        obs_frame = pd.DataFrame(spatial_obs)
        null_frame = pd.DataFrame(spatial_null)
        obs_cols = [f'D{i}' for i in range(1, 6)] + [f'F{i}' for i in range(1, 4)] + [f'B{i}' for i in range(1, 4)]
        null_cols = [f'D{i}' for i in range(1, 6)]
        aobs = aggregate(obs_frame, obs_cols)
        anull = aggregate(null_frame, null_cols)
        for row in aobs.itertuples(index=False):
            for m in range(1, 6):
                observed_cell_rows.append({
                    'rep': rep, 'rx': int(row.rx), 'ry': int(row.ry),
                    'n_species': int(row.n_species), 'mode': f'M{m}',
                    'D': getattr(row, f'D{m}'),
                    'F': getattr(row, f'F{m}') if m <= 3 else np.nan,
                    'B': getattr(row, f'B{m}') if m <= 3 else np.nan,
                })
        for row in anull.itertuples(index=False):
            for m in range(1, 6):
                null_cell_rows.append({
                    'rep': rep, 'rx': int(row.rx), 'ry': int(row.ry),
                    'n_species': int(row.n_species), 'mode': f'M{m}',
                    'D': getattr(row, f'D{m}'),
                })

    check = pd.DataFrame(upstream_check_rows)
    upstream = pd.read_csv(args.upstream_realizations)
    check_cols = ['rep', 'n_species', 'n_states'] + [
        x for i in range(1, 6)
        for x in (f'M{i}_cosine', f'M{i}_matched_run_rank', f'M{i}_explained')
    ]
    merged = check[check_cols].merge(upstream[check_cols], on='rep', suffixes=('_new', '_up'))
    max_diff = 0.0
    for col in check_cols[1:]:
        max_diff = max(max_diff, float(np.max(np.abs(
            merged[f'{col}_new'].to_numpy(float) - merged[f'{col}_up'].to_numpy(float)
        ))))
    if max_diff > 1e-12:
        raise RuntimeError(f'upstream recurrence mismatch {max_diff}')

    obs_long = pd.DataFrame(observed_cell_rows)
    null_long = pd.DataFrame(null_cell_rows)

    def spatial_summary(frame, value_col, prefix):
        rows = []
        for (rx, ry, mode), g in frame.groupby(['rx', 'ry', 'mode']):
            ppos = float((g[value_col] > 0).mean())
            rows.append({
                'rx': int(rx), 'ry': int(ry), 'mode': mode,
                f'{prefix}_support_reps': int(g['rep'].nunique()),
                f'{prefix}_support_fraction': float(g['rep'].nunique() / args.reps),
                f'{prefix}_species_median': float(g['n_species'].median()),
                f'{prefix}_score_median': float(g[value_col].median()),
                f'{prefix}_score_q025': float(g[value_col].quantile(.025)),
                f'{prefix}_score_q975': float(g[value_col].quantile(.975)),
                f'{prefix}_P_positive': ppos,
                f'{prefix}_R_sign': max(ppos, 1.0 - ppos),
                f'{prefix}_direction': 'positive' if ppos > .5 else ('negative' if ppos < .5 else 'tie'),
            })
        return pd.DataFrame(rows)

    obs_sum = spatial_summary(obs_long, 'D', 'obs')
    null_sum = spatial_summary(null_long, 'D', 'null')
    spatial = obs_sum.merge(null_sum, on=['rx', 'ry', 'mode'])
    spatial['R_sign_excess_over_null'] = spatial['obs_R_sign'] - spatial['null_R_sign']
    spatial['main_map'] = spatial['obs_support_fraction'] >= MAIN_SUPPORT_FRACTION
    spatial.to_csv(args.outdir / 'spatial_recurrence_summary.csv', index=False)

    global_spatial = []
    for mode in [f'M{i}' for i in range(1, 6)]:
        g = spatial[(spatial['mode'] == mode) & spatial['main_map']]
        global_spatial.append({
            'mode': mode,
            'n_main_map_cells': len(g),
            'n_cells_Rsign_ge_0_8': int((g['obs_R_sign'] >= .8).sum()),
            'n_cells_Rsign_ge_0_9': int((g['obs_R_sign'] >= .9).sum()),
            'n_cells_excess_Rsign_ge_0_1': int((g['R_sign_excess_over_null'] >= .1).sum()),
            'max_obs_Rsign': float(g['obs_R_sign'].max()) if len(g) else np.nan,
            'max_Rsign_excess': float(g['R_sign_excess_over_null'].max()) if len(g) else np.nan,
        })
    pd.DataFrame(global_spatial).to_csv(args.outdir / 'spatial_recurrence_global_summary.csv', index=False)

    decomp = pd.DataFrame(decomposition_rows)
    global_decomp = []
    for mode, g in decomp.groupby('mode'):
        row = {'mode': mode}
        for col in ['corr_D_F','corr_D_negB','mean_abs_D','mean_abs_F','mean_abs_B','frac_absF_gt_absB','frac_absB_gt_absF','max_identity_error']:
            row[f'{col}_median'] = g[col].median()
            row[f'{col}_q025'] = g[col].quantile(.025)
            row[f'{col}_q975'] = g[col].quantile(.975)
        global_decomp.append(row)
    pd.DataFrame(global_decomp).to_csv(args.outdir / 'global_flower_background_decomposition_summary.csv', index=False)

    component_rows = []
    for (rx, ry, mode), g in obs_long[obs_long['mode'].isin(['M1','M2','M3'])].groupby(['rx','ry','mode']):
        rec = {'rx': int(rx), 'ry': int(ry), 'mode': mode, 'support_reps': int(g['rep'].nunique())}
        for label in ['D','F','B']:
            p = float((g[label] > 0).mean())
            rec[f'{label}_P_positive'] = p
            rec[f'{label}_R_sign'] = max(p, 1-p)
            rec[f'{label}_median'] = float(g[label].median())
        component_rows.append(rec)
    component = pd.DataFrame(component_rows)
    main_cells = spatial[spatial['main_map'] & spatial['mode'].isin(['M1','M2','M3'])][
        ['rx','ry','mode','obs_R_sign','obs_P_positive','obs_direction','R_sign_excess_over_null']
    ]
    component = main_cells.merge(component, on=['rx','ry','mode'], how='left')
    max_spatial_diff = float(np.max(np.abs(component['obs_R_sign'] - component['D_R_sign'])))
    if max_spatial_diff > 1e-12:
        raise RuntimeError(f'D spatial recurrence mismatch {max_spatial_diff}')
    component.to_csv(args.outdir / 'main_map_flower_background_decomposition.csv', index=False)
    check.to_csv(args.outdir / 'upstream_mode_reproduction_check.csv', index=False)

    print(f'upstream max diff={max_diff:.3g}; D spatial max diff={max_spatial_diff:.3g}')


if __name__ == '__main__':
    main()
