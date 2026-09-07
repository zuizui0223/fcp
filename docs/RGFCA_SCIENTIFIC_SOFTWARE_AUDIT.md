# Scientific software citations and runtime boundaries

Checked 7 September 2026 against the projects' own citation guidance and the
implemented FCP entry points. This is a bounded software-credit audit, not a
complete dependency bill of materials, licence clearance or ecological validation.
No reserve outcome file was opened. No frozen analysis was rerun or upgraded.

## Verified citation records

- **NumPy:** Harris, C. R., et al. (2020), *Array programming with NumPy*,
  *Nature* 585:357–362. [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2).
  Read the project's [citation and bibliographic record](https://numpy.org/citing-numpy/).
  NumPy implements array operations and randomization support in the inspected
  FCP inference code; this paper is a software credit, not evidence that our
  ecological null or sample is valid.
- **SciPy:** Virtanen, P., et al. (2020), *SciPy 1.0: fundamental algorithms for
  scientific computing in Python*, *Nature Methods* 17:261–272.
  [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).
  Read the project's [citation record](https://scipy.org/citing-scipy/).
  The inspected omnibus imports `rankdata`; reserve inference additionally
  imports `spearmanr` for equivalence checks. The citation does not imply use
  of SciPy's asymptotic pair-count p-values as the RGFCA test.
- **pandas:** McKinney, W. (2010), *Data Structures for Statistical Computing
  in Python*, *Proceedings of the 9th Python in Science Conference*, pp. 56–61.
  [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a).
  Read [official citation guidance](https://pandas.pydata.org/about/citing.html)
  and the linked [six-page proceedings PDF](https://pub.curvenote.com/01908378-3686-7168-a380-d82bbf21c799/public/mckinney-57fc0d4e8a08cd7f26a4b8bf468a71f4.pdf)
  for title, author, year and page range. The DOI resolver returned HTTP 403 in
  this audit, while the official linked PDF was accessible. pandas is used for
  tables and census checks. The official guide also requests a version-specific
  software citation; its generic `latest` example is not a release identifier.
- **Matplotlib:** Hunter, J. D. (2007), *Matplotlib: A 2D graphics environment*,
  *Computing in Science & Engineering* 9(3):90–95.
  [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55).
  Read [official citation guidance](https://matplotlib.org/stable/project/citing.html).
  This credits static figure generation, not measurement, model fitting or
  inference. Current documentation is newer than the pinned rendering runtime;
  the code and output manifests retain the actual version used.

These records establish bibliographic metadata and the projects' requested
credits. They are not claims to have audited every algorithm in the four papers.
Version-specific software archive identifiers remain to be reconciled with the
final package; no all-version DOI is labelled as an exact release DOI.

## Implementation and version evidence

The inspected entry points are
`scripts/analysis/run_global_rgfca_within_species_spatial_omnibus.py`,
`fcp_pipeline/rgfca_reserve_inference.py`,
`scripts/analysis/make_rgfca_publication_figures.py`, and
`scripts/analysis/make_rgfca_photo_bar.py`. The frozen omnibus workflow and
prospective reserve-inference workflow pin Python 3.11, NumPy 2.3.5,
SciPy 1.17.0 and pandas 2.2.3. This does not imply reserve inference has run.

The first-two-figure rendering workflow uses NumPy 2.2.6, pandas 2.3.3 and
Matplotlib 3.10.8. The discovery photo-bar acquisition receipt records its
separate actual runtime, including NumPy 2.3.5, pandas 2.2.3 and Pillow 12.1.1.
Do not copy display runtime versions into the scientific inference provenance.
Re-rendering existing lossless crops is not a new model measurement.

The model/data-provider audit is separate from these four software citations.
Training data and model rights, observation-metadata redistribution, all other
dependencies, environmental-layer references and final release citations remain
submission-package checks. Publication figures and a passing test suite do not
close those gates or the independent ecological validation gate.
