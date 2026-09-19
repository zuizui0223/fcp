from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_package_uses_mit_license_metadata():
    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    release_template = (ROOT / "release_metadata.example.toml").read_text(
        encoding="utf-8"
    )

    assert license_text.startswith("MIT License")
    assert "Copyright (c) 2026 disttrait contributors" in license_text
    assert 'license = { file = "LICENSE" }' in pyproject
    assert '"License :: OSI Approved :: MIT License"' in pyproject
    assert 'license_spdx = "MIT"' in release_template
