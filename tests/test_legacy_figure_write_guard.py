"""Replay the real write step without permitting repository or network writes."""

import os
from pathlib import Path
import shutil
import subprocess
import textwrap

import pytest


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/34species-paper.yml"
WRITE_STEP = "      - name: Commit generated figures and manuscript integration to PR branch"


def replay_write_step(actor, *, changed=True, triggering_actor="zuizui0223"):
    workflow = WORKFLOW.read_text(encoding="utf-8")
    step = workflow.split(WRITE_STEP, 1)[1].split("      - ", 1)[0]
    script = textwrap.dedent(step.split("        run: |\n", 1)[1])
    script = script.replace("${{ github.head_ref }}", "fixture/figure-guard")
    # Execute the production shell body. Every git call is intercepted here,
    # including commit/push; no checkout, object, index, or remote is changed.
    stub = """
    git() {
      printf 'GIT_CALL %s\\n' "$*"
      if [ "$1" = diff ]; then return "$FCP_FIXTURE_DIFF_EXIT"; fi
      return 0
    }
    """
    env = os.environ.copy()
    env.pop("GITHUB_ACTOR", None)
    if actor is not None:
        env["GITHUB_ACTOR"] = actor
    env["GITHUB_TRIGGERING_ACTOR"] = triggering_actor
    env["FCP_FIXTURE_DIFF_EXIT"] = "1" if changed else "0"
    bash = (Path(os.environ.get("ProgramFiles", "C:/Program Files"))
            / "Git/bin/bash.exe") if os.name == "nt" else shutil.which("bash")
    assert bash and Path(bash).is_file(), "Bash is required to replay the CI write step"
    return subprocess.run(
        [str(bash), "--noprofile", "--norc", "-c", textwrap.dedent(stub) + script],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=15,
    )


@pytest.mark.parametrize("triggering_actor", ["github-actions[bot]", "zuizui0223"])
def test_bot_generated_update_cannot_write_again_even_after_maintainer_approval(triggering_actor):
    result = replay_write_step("github-actions[bot]", triggering_actor=triggering_actor)
    assert result.returncode == 0, result.stderr
    assert "GIT_CALL" not in result.stdout, result.stdout


def test_human_change_still_commits_and_pushes_generated_figures():
    result = replay_write_step("zuizui0223")
    assert result.returncode == 0, result.stderr
    assert "GIT_CALL commit -m Integrate expanded JBI figures into manuscript" in result.stdout
    assert "GIT_CALL push origin HEAD:fixture/figure-guard" in result.stdout


def test_unchanged_human_run_does_not_commit_or_push():
    result = replay_write_step("zuizui0223", changed=False)
    assert result.returncode == 0, result.stderr
    assert "GIT_CALL commit" not in result.stdout
    assert "GIT_CALL push" not in result.stdout


def test_missing_original_actor_fails_closed_before_git_writes():
    result = replay_write_step(None)
    assert result.returncode != 0
    assert "GIT_CALL" not in result.stdout


def test_bot_guard_does_not_skip_scientific_or_figure_validation():
    workflow = WORKFLOW.read_text(encoding="utf-8")
    validation, write = workflow.split(WRITE_STEP, 1)
    assert "github.actor" not in validation
    assert "GITHUB_ACTOR" not in validation
    assert "github.event_name == 'pull_request'" in validation
    for name in (
        "Run five symmetric niche models",
        "Run design-based power and precision",
        "Validate frozen scope and numerical regression",
        "Verify bit-reproducible PNG and PDF figures",
        "Validate figure outputs and manuscript calls",
    ):
        assert f"- name: {name}" in validation
    assert "GITHUB_ACTOR" in write
