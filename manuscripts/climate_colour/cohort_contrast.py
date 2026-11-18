"""Build Table S3 from paired species coefficients and the archived primary estimates."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
HUES = 'white yellow orange red pink magenta purple blue bronze'.split()
DOMAINS = ['thermal_level', 'precip_amount']
Q = np.array([8] + [-1] * 8) / np.sqrt(72)
COHORTS = {'discovery': 'Discovery', 'validation': 'Validation',
           'prospective': 'Prospective', 'wave4': 'Fourth wave',
           'new_U100': 'Largest cohort'}

def alignment(vectors):
    lengths = (vectors ** 2).sum(axis=-1)
    if np.any(lengths == 0):
        raise ValueError('Alignment undefined for a zero vector')
    return 100 * (vectors @ Q) ** 2 / lengths

def omission_result(frame, cohort):
    kept = frame.loc[frame.source_cohort != cohort]
    vectors = np.stack([kept[[f'{domain}_{h}' for h in HUES]].to_numpy()
                        for domain in DOMAINS], axis=1)
    genus, names = pd.factorize(kept.genus, sort=True)
    genus_count = len(names)
    counts = np.bincount(genus)
    sums = np.zeros((genus_count, 2, 9))
    np.add.at(sums, genus, vectors)
    rng = np.random.default_rng(20261111)
    differences = []
    for _ in range(20):
        draw = rng.integers(genus_count, size=(200, genus_count))
        means = sums[draw].sum(axis=1) / counts[draw].sum(axis=1)[:, None, None]
        temperature, precipitation = alignment(means).T
        differences.extend(temperature - precipitation)
    lo, hi = np.quantile(differences, [.025, .975])
    temp, precip = alignment(vectors.mean(axis=0))
    return [f'Omit {COHORTS[cohort]}', len(kept), genus_count,
            temp, precip, temp - precip, lo, hi, 4000]

if __name__ == '__main__':
    species = pd.read_csv(ROOT/'Table_S2_species_coefficients.csv.gz')
    cols = [[f'{d}_{h}' for h in HUES] for d in DOMAINS]
    paired = species.dropna(subset=cols[0] + cols[1])
    individual = alignment(np.stack([paired[c].to_numpy() for c in cols], axis=1))
    assert len(paired) == 2689 and np.sum(individual[:, 0] > individual[:, 1]) == 1345
    assert set(paired.source_cohort) == set(COHORTS) and paired.genus.notna().all()
    source = ROOT.parents[1]/'results/fcp4475_shared_thermal_moisture_photo_hues_verified_20261024.json'
    archive = json.loads(source.read_text(encoding='utf8'))
    records = []
    for key, label in [('original', 'All matched species'), ('observer_thinned', 'One photo per observer')]:
        r = archive['results'][key]['matched_climatic_photo_colour_geometry']
        records.append([label, r['n_species'], r['n_genera'],
                        100*r['thermal_white_axis_alignment'], 100*r['precip_white_axis_alignment'],
                        100*r['difference_white_axis_alignment'],
                        *(100*np.array(r['paired_genus_bootstrap95'])), archive['n_whole_genus_resamples']])
    records.extend(omission_result(paired, cohort) for cohort in COHORTS)
    fields = ['panel','species','genera','temperature_pct','precipitation_pct',
              'difference_pp','lower95_pp','upper95_pp','resamples']
    pd.DataFrame(records, columns=fields).to_csv(ROOT/'Table_S3_cohort_contrasts.csv',
                                                  index=False, float_format='%.8f')
