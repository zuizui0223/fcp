"""Render frozen summaries only; no photo access, model fitting or new tests."""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
SOURCES = {
    'h1_direct': 'results/polymorphism_h1_observer_disjoint_d_20260913/result.json',
    'h1_strict': 'results/polymorphism_h1_observer_disjoint_D_reliability_20260913/result.json',
    'h1_repeated': 'results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json',
    'h2': 'results/polymorphism_white_axis_targeted_test_20260912/result.json',
    'h3a': 'docs/supporting/h123/h3a_signal_by_scenario.csv',
    'h3b': 'docs/supporting/h123/h3b_span_summary.csv',
}


def figure_data(root=ROOT):
    sources, hashes = {}, {}
    for key, path in SOURCES.items():
        raw = (root / path).read_bytes().replace(b'\r\n', b'\n')
        hashes[path] = hashlib.sha256(raw).hexdigest()
        sources[key] = (json.loads(raw) if path.endswith('.json')
                        else list(csv.DictReader(raw.decode().splitlines())))
    direct = sources['h1_direct']['analyses']['reserve']['D']
    strict = sources['h1_strict']['reserve']
    repeated = sources['h1_repeated']['reserve']['primary20']
    h1 = [
        dict(label='Direct split', estimate=direct['spearman_rho'],
             interval=direct['bootstrap_95_percentile_ci'], interval_type='bootstrap 95% CI', n=direct['n']),
        dict(label='Later strict split', estimate=strict['spearman_D_A_D_B'],
             interval=strict['spearman_bootstrap_95_percentile_ci'], interval_type='bootstrap 95% CI', n=strict['reliability_eligible_n']),
        dict(label='200 partitions', estimate=repeated['rho_median'],
             interval=[repeated['rho_q05'], repeated['rho_q95']], interval_type='partition 5th–95th percentiles', n=repeated['paired_n_median']),
    ]
    h2 = []
    for threshold, short in [('primary_0_10', 'Primary'), ('strict_0_20', 'Strict')]:
        for cohort in ('discovery', 'reserve'):
            r = sources['h2']['thresholds'][threshold][cohort]
            h2.append(dict(label=f'{short}: {cohort}', n=r['species'],
                           observed=r['observed_mean_squared_white_axis_alignment'],
                           median=r['structured_null_summary']['q50'],
                           interval=[r['structured_null_summary']['q025'], r['structured_null_summary']['q975']],
                           p=r['structured_null_upper_p']))
    return dict(h1=h1, h2=h2, h3a=sources['h3a'], h3b=sources['h3b'],
                source_sha256_lf=hashes, scientific_results_recomputed=False,
                p500_opened_by_renderer=False)


def render(output_dir, root=ROOT):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    data = figure_data(root)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42})
    blue, gold, grey = '#235789', '#986F19', '#777777'
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.4))
    fig.subplots_adjust(left=.17, right=.97, bottom=.19, top=.91, wspace=.68, hspace=.78)
    a, b, c, d = axes.ravel()
    a.set_title('A  H1: reserve sampling repeatability', loc='left', pad=15)
    for y, r in enumerate(data['h1']):
        color = gold if y == 2 else blue
        a.plot(r['interval'], [y, y], color=color, lw=2, linestyle='--' if y == 2 else '-')
        a.plot(r['estimate'], y, 'D' if y == 2 else 'o', color=color, ms=6)
        a.text(.91, y, f"{r['estimate']:.3f}", ha='right', va='center')
    a.axvline(.8, color=grey, ls=':', lw=1)
    a.set_yticks(range(3), ['Direct; N=363\n95% bootstrap CI',
                           'Later strict; N=363\n95% bootstrap CI',
                           'Repeated; median N=329\n5–95% partition range'])
    a.set_ylim(2.6, -.7)
    a.set_xlim(.65, .92)
    a.set_xlabel('Spearman rho (focused scale)')
    a.text(0, -0.34, 'Dotted line: 0.80 criterion for direct / later strict only.\nLater strict fails; repeated partitions have a different rule.',
           transform=a.transAxes, fontsize=9)

    b.set_title('B  H2: targeted white-axis alignment', loc='left', pad=15)
    for y, r in enumerate(data['h2']):
        b.plot(r['interval'], [y, y], color=grey, lw=3)
        b.plot(r['median'], y, 's', color=grey, ms=5)
        b.plot(r['observed'], y, 'o', color=blue, ms=6)
        b.text(.60, y, f"p={r['p']:.3f}", ha='right', va='center', fontsize=9)
    b.set_yticks(range(4), [f"{r['label']}\nN={r['n']}" for r in data['h2']])
    b.set_xlim(.38, .61)
    b.set_ylim(3.7, -.7)
    b.set_xlabel('Mean squared alignment W (focused scale)')
    b.text(0, -0.34, 'Circle: observed. Grey: null median and 2.5–97.5% range.\nNull range is NOT an observed-effect confidence interval.',
           transform=b.transAxes, fontsize=9)

    c.set_title('C  H3a: raw phylogenetic signal', loc='left', pad=15)
    for r in data['h3a']:
        y = int(r['scenario'][1:]) - 1 + (-.13 if r['cohort'] == 'discovery' else .13)
        reserve = r['cohort'] == 'reserve'
        c.plot(float(r['K_raw']), y, 's' if reserve else 'o', color=blue if reserve else gold, ms=6)
        c.text(.153, y, f"p={float(r['p_K_raw']):.4f}", ha='right', va='center', fontsize=9)
    c.set_yticks(range(3), ['S1', 'S2', 'S3'])
    c.set_ylim(2.6, -.6)
    c.set_xlim(0, .16)
    c.set_xlabel('Blomberg K (raw D)')
    c.text(0, -.34, 'Circle: discovery, 368 tips. Square: reserve, 341 tips.\nScenarios are not confidence bounds; no CI is supplied.',
           transform=c.transAxes, fontsize=9)

    d.set_title('D  H3b: sampled photographic span', loc='left', pad=15)
    for y, r in enumerate(data['h3b']):
        d.plot(float(r['rho_D_span']), y, 'o' if y == 0 else 's', color=gold if y == 0 else blue, ms=6)
        d.text(.27, y, f"p={float(r['p_D_span']):.4f}", ha='right', va='center', fontsize=9)
    d.axvline(0, color=grey, lw=1, ls=':')
    d.set_yticks(range(2), [f"{r['cohort'].capitalize()}\nN={r['n_species']}" for r in data['h3b']])
    d.set_ylim(1.7, -.7)
    d.set_xlim(-.08, .28)
    d.set_xlabel('Spearman rho: D versus sampled span')
    d.text(0, -.34, 'No CI is supplied by the retained summary.\nSampled span is not verified biological range size.', transform=d.transAxes, fontsize=9)
    for ax in axes.ravel():
        ax.grid(axis='x', alpha=.16)
        ax.set_axisbelow(True)
    fig.suptitle('Frozen H1–H3 evidence: repeatability, targeted geometry and explanatory limits', fontsize=14)
    fig.text(.17, .025, 'Image-derived descriptors, not validated morph frequencies. H2 is retrospective. Measurement qualification remains unresolved.', fontsize=9)
    fig.savefig(output_dir / 'h123_evidence.png', dpi=180)
    fig.savefig(output_dir / 'h123_evidence.pdf', metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    (output_dir / 'h123_evidence_data.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    return data


if __name__ == '__main__':
    render(ROOT / 'docs/figures/h123')
