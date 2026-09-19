# disttrait release checklist

Current validated package version: **0.12.0**

## Already executable

- [x] package source isolated under `packages/disttrait/`;
- [x] PEP 517 / setuptools `pyproject.toml`;
- [x] deterministic unit and regression tests;
- [x] frozen FCP/RGFCA algorithm-equivalence fixtures;
- [x] non-flower empirical transports;
- [x] synthetic calibration/power/estimand benchmarks;
- [x] sdist + wheel build gate;
- [x] Twine metadata validation gate;
- [x] fresh-virtualenv wheel install smoke test;
- [x] versioned changelog.

## Requires an explicit ownership/release decision

These items are intentionally **not guessed by automation**:

- [ ] choose the package software licence;
- [ ] confirm package authors / maintainers and ordering;
- [ ] add final `CITATION.cff`;
- [ ] choose the standalone GitHub repository name/visibility;
- [ ] create the standalone repository (the current GitHub connector does not
      expose repository creation);
- [ ] decide whether to publish to PyPI/TestPyPI;
- [ ] tag the first external release;
- [ ] archive the tagged release and record a DOI (for example through Zenodo).

## Recommended first external release sequence

1. resolve licence and author metadata;
2. copy `packages/disttrait/` plus its tests/examples/benchmarks into the
   standalone repository;
3. preserve the FCP equivalence manifest and source Git SHAs;
4. run the release-package check in the standalone repository;
5. create tag `v0.12.0`;
6. archive the tag and record the DOI;
7. only then enable an authenticated PyPI publishing workflow.

No package-index credentials should be added before the licence/authorship and
repository decisions are frozen.
