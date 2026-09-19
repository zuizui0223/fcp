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

- [x] package software licence fixed as **MIT**;
- [ ] confirm package authors / maintainers and ordering;
- [ ] choose the standalone GitHub repository name/visibility;
- [ ] create the standalone repository (the current GitHub connector does not
      expose repository creation);
- [ ] decide whether to publish to PyPI/TestPyPI;
- [ ] tag the first external release;
- [ ] archive the tagged release and record a DOI (for example through Zenodo).

Metadata rendering itself is now automated. Copy
`release_metadata.example.toml` to `release_metadata.toml`, fill only explicit
ownership/release decisions, place the selected licence text in `LICENSE`, set
`release.release_ready = true`, and run:

```bash
python scripts/release/finalize_disttrait_release_metadata.py \
  --metadata packages/disttrait/release_metadata.toml \
  --package-dir packages/disttrait

python scripts/release/check_disttrait_public_release.py \
  --metadata packages/disttrait/release_metadata.toml \
  --package-dir packages/disttrait
```

The finalizer generates the final `CITATION.cff`, author/maintainer/license
fields and standalone URLs in `pyproject.toml`, plus
`RELEASE_METADATA.json`. It refuses placeholders; it does not decide those
values.

## Recommended first external release sequence

1. finalize author / maintainer metadata;
2. copy `packages/disttrait/` plus its tests/examples/benchmarks into the
   standalone repository;
3. preserve the FCP equivalence manifest and source Git SHAs;
4. run the release-package check in the standalone repository;
5. create tag `v0.12.0`;
6. archive the tag and record the DOI;
7. only then enable an authenticated PyPI publishing workflow.

No package-index credentials should be added before authorship and repository
decisions are frozen.
