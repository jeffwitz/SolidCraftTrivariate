#!/usr/bin/env python3
"""Rebuild and audit standalone geometry recipes and recipe-bearing corpora.

Corpus admission fields (expected, accepted, or admitted) define intended
refusals. Cases without one are treated as positive recipes. Refused builds
must preserve a pre-existing output file. This runner does not export or
publish generated geometry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
SEARCH_ROOTS = (ROOT / "tests/references", ROOT / "tests/regressions")
STATUS_KEYS = ("expected", "accepted", "admitted")
NON_CLI_PATHS = {
    "tests/references/g2_plate_jets.input.json": (
        "historical S1 curve-analysis input; no production geometry builder is registered",
        "g2_plate_v1",
    ),
    "tests/regressions/pocket_wall_internal_face.json": (
        "audit-selection regression fixture, not a build recipe",
        "pocket_depth_v1",
    ),
    "tests/regressions/rational_sector_endpoint_y_control.json": (
        "coefficient-mutation regression fixture, not a build recipe",
        "rational_sector_prism_v1",
    ),
}


def corpus_records(path: Path, value):
    """Yield each embedded recipe record once, with its explicit status."""
    if isinstance(value, list):
        for index, item in enumerate(value):
            if isinstance(item, dict) and isinstance(item.get("recipe"), dict):
                recipe = item["recipe"]
                if isinstance(recipe.get("family"), str):
                    status = next((item[k] for k in STATUS_KEYS if isinstance(item.get(k), bool)), True)
                    yield f"{path.relative_to(ROOT)}[{index}]", recipe, status
            else:
                yield from corpus_records(path, item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from corpus_records(path, item)


def inventory():
    rows, excluded = [], []
    for root in SEARCH_ROOTS:
        for path in sorted(root.rglob("*.json")):
            try:
                value = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError) as error:
                raise RuntimeError(f"cannot read recipe inventory JSON {path.relative_to(ROOT)}") from error
            rel = str(path.relative_to(ROOT))
            if isinstance(value, dict) and isinstance(value.get("family"), str):
                if rel in NON_CLI_PATHS:
                    reason, expected_family = NON_CLI_PATHS[rel]
                    if value["family"] != expected_family:
                        raise RuntimeError(f"non-build fixture family changed at {rel}: {value['family']}")
                    excluded.append({"source": rel, "family": value["family"], "reason": reason})
                else:
                    rows.append((rel, value, True))
            else:
                rows.extend(corpus_records(path, value))
    # The same physical recipe can appear as a corpus case and as a standalone
    # file; keep both records because they represent separately named fixtures.
    return rows, excluded


def run(cli: Path, work: Path):
    cases, excluded = inventory()
    inventory_payload = json.dumps(
        {"cases": cases, "excluded": excluded}, sort_keys=True, separators=(",", ":")
    ).encode()
    inventory_sha256 = hashlib.sha256(inventory_payload).hexdigest()
    work.mkdir(parents=True, exist_ok=True)
    rows = []
    canonical_parent = None
    for index, (source, recipe, expected_success) in enumerate(cases):
        recipe_path = work / f"recipe-{index:04}.json"
        model_path = work / f"model-{index:04}.json"
        recipe_path.write_text(json.dumps(recipe, separators=(",", ":"), allow_nan=False))
        sentinel = b"preserve-on-refusal\n"
        model_path.write_bytes(sentinel)
        build = subprocess.run([str(cli), "build", str(recipe_path), str(model_path)],
                               capture_output=True, text=True)
        accepted = build.returncode == 0
        preserved = model_path.read_bytes() == sentinel if not accepted else None
        audit_ok = None
        audit_status = None
        audit_stderr = ""
        parent_match = None
        if accepted:
            if isinstance(recipe.get("parent_model"), dict):
                if canonical_parent is None:
                    canonical_recipe = ROOT / "tests/references/g2_cross_sketches.recipe.json"
                    canonical_input = work / "canonical-g2-cross.recipe.json"
                    canonical_output = work / "canonical-g2-cross.model.json"
                    canonical_input.write_bytes(canonical_recipe.read_bytes())
                    canonical_build = subprocess.run(
                        [str(cli), "build", str(canonical_input), str(canonical_output)],
                        capture_output=True, text=True,
                    )
                    if canonical_build.returncode:
                        raise RuntimeError("canonical G2 parent recipe failed to build")
                    canonical_parent = json.loads(canonical_output.read_text())
                parent_match = recipe["parent_model"] == canonical_parent
            audit = subprocess.run([str(cli), "audit", str(model_path)],
                                   capture_output=True, text=True)
            try:
                audit_json = json.loads(audit.stdout)
                audit_status = audit_json.get("status")
                audit_ok = audit.returncode == 0 and audit_status == "valid"
            except json.JSONDecodeError:
                audit_stderr = (audit.stderr or audit.stdout)[-500:]
                audit_ok = False
        passed = (accepted == expected_success and preserved is not False
                  and audit_ok is not False and parent_match is not False)
        rows.append({
            "source": source,
            "family": recipe["family"],
            "expected_success": expected_success,
            "build_success": accepted,
            "audit_status": audit_status,
            "refusal_preserved_output": preserved,
            "embedded_parent_matches_current_canonical": parent_match,
            "passed": passed,
            "diagnostic": "" if passed else (build.stderr or build.stdout or audit_stderr)[-700:],
        })
        if not passed:
            print(f"FAIL {source}: expected_success={expected_success} build={accepted} audit={audit_status}", flush=True)
    return rows, excluded, inventory_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    before = hashlib.sha256(cli.read_bytes()).hexdigest()
    rows, excluded, inventory_sha256 = run(cli, args.work.resolve())
    after = hashlib.sha256(cli.read_bytes()).hexdigest()
    if before != after:
        raise RuntimeError("CLI binary changed during replay")
    report = {
        "schema": "solidcraft_geometry_recipe_replay_v1",
        "cli": os.path.relpath(cli, ROOT),
        "cli_sha256": before,
        "checkout_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "inventory_sha256": inventory_sha256,
        "case_count": len(rows),
        "positive_count": sum(r["expected_success"] for r in rows),
        "intentional_refusal_count": sum(not r["expected_success"] for r in rows),
        "passed_count": sum(r["passed"] for r in rows),
        "failed_count": sum(not r["passed"] for r in rows),
        "excluded_non_cli_fixture_count": len(excluded),
        "excluded_non_cli_fixtures": excluded,
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"{report['passed_count']}/{report['case_count']} recipes passed; report: {args.output}")
    if report["failed_count"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
