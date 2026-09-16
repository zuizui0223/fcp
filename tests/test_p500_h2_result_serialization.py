import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py"


def _load_repo_display_path():
    source = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(source)
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "repo_display_path"
    )
    isolated = ast.Module(body=[function], type_ignores=[])
    namespace = {"Path": Path, "ROOT": ROOT}
    exec(compile(isolated, str(SCRIPT), "exec"), namespace)
    return namespace["repo_display_path"], source


def test_relative_measurement_result_path_does_not_require_relative_to_root():
    display, source = _load_repo_display_path()
    relative = Path("results/polymorphism_h2_p500_prospective_execution_20260915/measurement_result.json")

    assert display(relative) == str(relative)
    assert "args.measurement_result.relative_to(ROOT)" not in source
    assert "repo_display_path(args.measurement_result)" in source


def test_absolute_repo_path_is_serialized_repo_relative():
    display, _ = _load_repo_display_path()
    relative = Path("results/example/measurement_result.json")

    assert display(ROOT / relative) == str(relative)


def test_absolute_external_path_is_preserved_for_auditability(tmp_path):
    display, _ = _load_repo_display_path()
    external = tmp_path / "measurement_result.json"

    assert display(external) == str(external)
