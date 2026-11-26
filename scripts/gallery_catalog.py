#!/usr/bin/env python3
"""Regenerate all gallery geometry bundles from tracked recipes and the Rust CLI."""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from fractions import Fraction
from pathlib import Path
from check_gallery_geometry import step_content


ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/site"
MANIFEST = SITE / "experiments.json"
POLYGON_RESULTS = ROOT / "docs/results/convex-host-rings.json"
BRACKET_RESULTS = ROOT / "docs/results/propagated-bores.json"
EXPECTED_IDS = 68


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def format_exact_volume(value) -> str:
    if isinstance(value, (list, tuple)) and len(value) == 2:
        numerator, denominator = map(int, value)
        amount = str(numerator) if denominator == 1 else f"{numerator}/{denominator}"
        return f"{amount} mm³"
    return f"{value} mm³" if value is not None else "see audit"


def sync_volume_metric(entry, audit):
    native = audit["model_audit"]
    exact = native.get("volume_exact_mm3")
    if exact is not None:
        value = format_exact_volume(exact)
    elif "exact_map_volume_mm3" in audit:
        value = format_exact_volume(audit["exact_map_volume_mm3"])
    else:
        source = native.get("volume_exact") or native.get("source_certificate", {}).get("volume_exact")
        if not source:
            entry["metrics"] = [row for row in entry["metrics"] if row["label"] != "Exact volume"]
            return
        scale = Fraction(abs(native.get("determinant_numerator", 1)), native.get("determinant_denominator", 1))
        constant = scale * Fraction(source["twice_rational_part_mm3"], 2)
        coefficient = scale * Fraction(*source["pi_coefficient_mm3"])
        sign = "−" if coefficient < 0 else "+"
        value = f"({constant} {sign} {abs(coefficient)}π) mm³"
    for row in entry["metrics"]:
        if row["label"] == "Exact volume":
            row["value"] = value


def format_bound(value) -> str:
    if value is None:
        return "see audit"
    return f"{float(value):.8g} mm³"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def site_path(url: str) -> Path:
    return SITE / url.split("?", 1)[0]


def run(cli: Path, *arguments: str) -> str:
    result = subprocess.run(
        [str(cli), *map(str, arguments)], check=True, text=True,
        capture_output=True, timeout=120,
    )
    return result.stdout.strip()


def sample_condition(cli: Path, model: Path, temp: Path, block_count: int) -> dict:
    vtk = temp / "geometry.vtu"
    run(cli, "export-geometry-vtu", model, "4", vtk)
    arrays = ET.parse(vtk).findall('.//PointData/DataArray[@Name="condition_number"]')
    per_block = [list(map(float, (array.text or "").split())) for array in arrays]
    if len(per_block) != block_count or any(len(values) != 125 for values in per_block):
        raise ValueError(f"incomplete closed 5x5x5 condition sample: {model}")
    values = [value for block in per_block for value in block]
    if not values or any(not (1.0 <= value < float("inf")) for value in values):
        raise ValueError(f"invalid sampled condition number: {model}")
    return {
        "method": "maximum spectral Jacobian condition number on the closed 5x5x5 parameter grid per block; diagnostic only",
        "points_per_block": 125,
        "sample_count": len(values),
        "maximum_condition_2": max(values),
    }


def polygon_latest() -> tuple[dict[str, dict], dict]:
    report = read_json(POLYGON_RESULTS)
    latest = report["latest_run_details"]
    by_id = {case["id"]: case for case in latest["cases"]}
    if len(by_id) != 45:
        raise ValueError(f"Expected 45 latest polygon records; found {len(by_id)}")
    return by_id, latest


def bracket_latest() -> tuple[dict[str, dict], dict]:
    report = read_json(BRACKET_RESULTS)
    cases = report["candidate"]["cases"]
    if len(cases) != 10:
        raise ValueError(f"Expected 10 latest bracket records; found {len(cases)}")
    return cases, report


def preserve_independent_record(case_id: str, old: dict, bracket_cases: dict, bracket_report: dict) -> dict:
    name = case_id.removeprefix("bench-bracket-").replace("-", "_")
    if name in bracket_cases:
        return {
            "scope": "Frozen bracket comparison at the recorded candidate and baseline revisions.",
            "candidate_revision": bracket_report.get("candidate_revision"),
            "baseline_revision": bracket_report.get("baseline_revision"),
            "candidate": bracket_cases[name],
            "reference": bracket_report.get("baseline", {}).get("cases", {}).get(name),
        }
    return old.get("independent_reference_validation", {})


def prepare_audit(case_id: str, old: dict, kernel_audit: dict, build_report: dict,
                  metadata: dict, kernel_sha: str, revision: str,
                  recipe: Path, model_bytes: bytes, step_bytes: bytes,
                  condition: dict, bracket_cases: dict, bracket_report: dict) -> dict:
    if case_id.startswith("bench_poly_"):
        audit = {
            "schema": "polygon_gallery_case_audit_v1",
            "case": metadata if metadata else {"id": case_id, "diagnostic": "projection-only control"},
        }
    elif case_id.startswith("bench-bracket-"):
        audit = {
            "schema": "bracket_propagation_gallery_audit_v1",
            "case": case_id.removeprefix("bench-bracket-").replace("-", "_"),
            "independent_reference_validation": preserve_independent_record(case_id, old, bracket_cases, bracket_report),
        }
    elif case_id.startswith("fieldcraft_"):
        audit = {
            "schema": "gallery_kernel_rebuild_audit_v1",
            "case": case_id,
            "independent_reference_validation": {
                "report": "https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/results/fieldcraft-compact-current-oracle.json",
                "scope": "Independent host B-rep comparison for the five FieldCraft reference geometries.",
            },
        }
    else:
        audit = dict(old)
        prior_readback = audit.pop("step_readback", None)
        if prior_readback is not None:
            audit["historical_independent_validation"] = {
                "scope": "Historical independent readback for the former gallery STEP; it does not validate the regenerated STEP.",
                "record": prior_readback,
            }
        # Do not nest an earlier regeneration record on repeated runs.
        audit.pop("kernel_rebuild", None)

    audit["model_audit"] = kernel_audit
    audit["kernel_rebuild"] = {
        "command": "build recipe.json model.json; audit model.json; export-step model.json reference.step",
        "kernel_revision": revision,
        "kernel_binary_sha256": kernel_sha,
        "recipe_sha256": sha256(recipe.read_bytes()),
        "model_sha256": sha256(model_bytes),
        "step_sha256": sha256(step_bytes),
        "build_report": build_report,
        "sampled_condition": condition,
    }
    return audit


def refresh_entry_metadata(entry: dict, build_report: dict, kernel_audit: dict,
                           model: dict, polygon_record: dict | None,
                           bracket_record: dict | None, condition: dict) -> None:
    method = model.get("partition_plan", {}).get("method", build_report.get("method", kernel_audit.get("method", model.get("schema", "native Rust builder"))))
    blocks = int(kernel_audit.get("blocks", len(model.get("blocks", []))))
    interfaces = int(kernel_audit.get("shared_faces", len(model.get("topology", {}).get("interfaces", []))))
    sampled_kappa = condition["maximum_condition_2"]
    if polygon_record:
        old_kappa = max(float(row["max_condition_2"]) for row in polygon_record["condition_grid_5cube_per_block"])
        if int(polygon_record["block_count"]) != blocks or method != polygon_record["construction_method"]:
            raise ValueError(f"latest polygon ledger does not match current kernel for {entry['id']}")
        latest_build = json.loads(polygon_record["build_stdout"])
        if int(latest_build["shared_faces"]) != interfaces:
            raise ValueError(f"latest polygon interface count does not match current kernel for {entry['id']}")
        if abs(sampled_kappa - old_kappa) > 1e-10 * max(1.0, old_kappa):
            raise ValueError(f"latest polygon condition measurement differs for {entry['id']}: {sampled_kappa} != {old_kappa}")
        current_det = float(kernel_audit["determinant_lower_bound_mm3"])
        latest_det = float(polygon_record["closed_domain_certificate"]["determinant_lower_bound_mm3"])
        if abs(current_det - latest_det) > 1e-10 * max(1.0, abs(latest_det)):
            raise ValueError(f"latest polygon determinant bound differs for {entry['id']}")
        entry["metrics"] = [
            {"label": "Partition method", "value": method},
            {"label": "Rust model blocks", "value": blocks},
            {"label": "Rust shared interfaces", "value": interfaces},
            {"label": "Sampled κ₂ max", "value": f"{sampled_kappa:.6g}"},
            {"label": "Certified determinant lower bound", "value": format_bound(kernel_audit.get("determinant_lower_bound_mm3", polygon_record.get("closed_domain_certificate", {}).get("determinant_lower_bound_mm3")))},
            {"label": "Exact volume", "value": format_exact_volume(kernel_audit.get("volume_exact_mm3", build_report.get("volume_exact_mm3")))},
        ]
        entry["description"] = (
            f"Latest tracked Rust release: {method}, {blocks} blocks and {interfaces} shared interfaces. "
            "Exact layer coverage and the closed-domain determinant certificate passed. "
            "The 5×5×5 Jacobian condition sample is a diagnostic, not a continuous bound."
        )
    if bracket_record:
        if int(bracket_record["blocks"]) != blocks or int(bracket_record["interfaces"]) != interfaces:
            raise ValueError(f"latest bracket ledger does not match current kernel for {entry['id']}")
        old_kappa = float(bracket_record["max_sampled_kappa"])
        if abs(sampled_kappa - old_kappa) > 1e-10 * max(1.0, old_kappa):
            raise ValueError(f"latest bracket condition measurement differs for {entry['id']}: {sampled_kappa} != {old_kappa}")
        entry["metrics"] = [
            {"label": "Partition method", "value": bracket_record.get("method", method)},
            {"label": "Rust model blocks", "value": blocks},
            {"label": "Rust shared interfaces", "value": interfaces},
            {"label": "Sampled κ₂ max", "value": f"{sampled_kappa:.6g}"},
            {"label": "Certified determinant lower bound", "value": format_bound(kernel_audit.get("determinant_lower_bound_mm3", bracket_record.get("det_lower_mm3")))},
            {"label": "Benchmark material probes", "value": f'{bracket_record.get("probes", "not reported")} probes; {bracket_record.get("probe_mismatches", "not reported")} disagreements'},
            {"label": "Exact volume", "value": format_exact_volume(kernel_audit.get("volume_exact_mm3", build_report.get("volume_exact_mm3")))},
        ]
    elif entry.get("id", "").startswith("bench-bracket-"):
        entry["metrics"] = [
            {"label": "Partition method", "value": method},
            {"label": "Rust model blocks", "value": blocks},
            {"label": "Rust shared interfaces", "value": interfaces},
            {"label": "Sampled κ₂ max", "value": f"{sampled_kappa:.6g}"},
            {"label": "Certified determinant lower bound", "value": format_bound(kernel_audit.get("determinant_lower_bound_mm3"))},
            {"label": "Exact volume", "value": format_exact_volume(kernel_audit.get("volume_exact_mm3", build_report.get("volume_exact_mm3")))},
        ]
    if not (polygon_record or bracket_record or entry.get("id", "").startswith("bench-bracket-")):
        # Geometry facts are regenerated from the native report. Keep only the
        # projection overlap annotations that are independent of the partition.
        retained = []
        if entry.get("id") == "bench_poly_projection_diamonds":
            retained = [row for row in entry.get("metrics", []) if row.get("label") in {
                "X projection overlap", "Y projection overlap", "Through-holes",
            }]
        entry["metrics"] = retained + [
            {"label": "Rust kernel method", "value": method},
            {"label": "Rust model blocks", "value": blocks},
            {"label": "Rust shared interfaces", "value": interfaces},
            {"label": "Sampled κ₂ max", "value": f"{sampled_kappa:.6g}"},
            {"label": "Certified determinant lower bound", "value": format_bound(kernel_audit.get("determinant_lower_bound_mm3"))},
            {"label": "Exact volume", "value": format_exact_volume(kernel_audit.get("volume_exact_mm3", build_report.get("volume_exact_mm3")))},
        ]
    if entry.get("id", "").startswith("bench-bracket-"):
        entry["stepLabel"] = "OCCT-checked boundary export (STEP)"
        entry["description"] = re.sub(
            r" The viewer renders the trivariate model;.*$",
            " The linked STEP is exported from this same native model; independent benchmark checks are scoped separately in the audit.",
            entry.get("description", ""),
        )
    else:
        entry["stepLabel"] = "OCCT-checked boundary export (STEP)"


def rebuild(entries: list[dict], cli: Path) -> None:
    if not cli.is_file():
        raise FileNotFoundError(f"Rust CLI not found: {cli}; run cargo build --release --locked --offline")
    kernel_sha = sha256(cli.read_bytes())
    revision = subprocess.run(["git", "log", "-1", "--format=%H", "--", "src", "Cargo.toml", "Cargo.lock"], cwd=ROOT, check=True,
                              text=True, capture_output=True).stdout.strip()
    polygon_records, polygon_run = polygon_latest()
    bracket_records, bracket_report = bracket_latest()

    for entry in entries:
        case_id = entry["id"]
        recipe, model, audit, step = (site_path(entry[key]) for key in ("source", "model", "audit", "step"))
        for path in (recipe, model, audit, step):
            if not path.is_file():
                raise FileNotFoundError(f"gallery asset is missing: {path}")
        old_audit = read_json(audit)
        old_step = step.read_bytes()
        polygon_record = polygon_records.get(case_id.removeprefix("bench_poly_")) if case_id.startswith("bench_poly_") else None
        bracket_name = case_id.removeprefix("bench-bracket-").replace("-", "_")
        bracket_record = bracket_records.get(bracket_name)

        with tempfile.TemporaryDirectory(prefix="gallery-regenerate-", dir=model.parent) as temp_name:
            temp = Path(temp_name)
            staged_model, staged_step, staged_audit = temp / "model.json", temp / "reference.step", temp / "audit.json"
            build_report = json.loads(run(cli, "build", recipe, staged_model))
            kernel_audit = json.loads(run(cli, "audit", staged_model))
            run(cli, "export-step", staged_model, staged_step)
            model_bytes = staged_model.read_bytes()
            fresh_step = staged_step.read_bytes()
            same_step_except_timestamp = step_content(old_step) == step_content(fresh_step)
            if same_step_except_timestamp:
                published_step = old_step
            else:
                published_step = fresh_step

            model_doc = read_json(staged_model)
            condition = sample_condition(cli, staged_model, temp, len(model_doc.get("blocks", [])))
            audit_doc = prepare_audit(
                case_id, old_audit, kernel_audit, build_report,
                polygon_record, kernel_sha, revision, recipe,
                model_bytes, published_step, condition, bracket_records, bracket_report,
            )
            if polygon_record:
                audit_doc["measurement_release"] = {
                    "tag": polygon_run["tag"], "cli_sha256": polygon_run["cli_sha256"],
                }
            refresh_entry_metadata(entry, build_report, kernel_audit, model_doc, polygon_record, bracket_record, condition)
            sync_volume_metric(entry, audit_doc)
            write_json(staged_audit, audit_doc)

            # Do not replace a fresh STEP whose only difference is its generated timestamp.
            os.replace(staged_model, model)
            if not same_step_except_timestamp:
                os.replace(staged_step, step)
            os.replace(staged_audit, audit)

        for key in ("model", "source", "step", "audit"):
            path = entry[key].split("?", 1)[0]
            entry[key] = f"{path}?v={sha256(site_path(path).read_bytes())[:12]}"


def remove_unreferenced_data(entries: list[dict]) -> int:
    referenced = set()
    for entry in entries:
        for value in entry.values():
            if isinstance(value, str) and value.startswith("data/"):
                referenced.add(site_path(value).resolve())
    root = (SITE / "data").resolve()
    removed = 0
    for path in root.rglob("*"):
        if path.is_file() and path.resolve() not in referenced:
            path.unlink()
            removed += 1
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_dir():
            try:
                path.rmdir()
            except OSError:
                pass
    return removed


def version_modules(filename: str, dependencies: list[str]) -> None:
    path = SITE / filename
    source = path.read_text(encoding="utf-8")
    for dependency in dependencies:
        url = "./" + dependency
        source, hits = re.subn(
            re.escape(url) + r"\?v=[a-f0-9]+",
            url + "?v=" + sha256((SITE / dependency).read_bytes())[:12],
            source,
        )
        if hits != 1:
            raise ValueError(f"Expected one versioned dependency {dependency} in {filename}; got {hits}")
    path.write_text(source, encoding="utf-8")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate",
                        help="Current release binary produced by cargo build --release --locked --offline")
    args = parser.parse_args()
    entries = read_json(MANIFEST)
    if len(entries) != EXPECTED_IDS or len({entry["id"] for entry in entries}) != EXPECTED_IDS:
        raise ValueError(f"Expected {EXPECTED_IDS} unique published gallery entries")
    rebuild(entries, args.cli.resolve())
    removed = remove_unreferenced_data(entries)
    write_json(MANIFEST, entries)
    version_modules("app.js", ["renderer.js"])
    version_modules("index.html", ["app.js", "experiments.json", "style.css"])
    print(f"Regenerated {len(entries)} gallery bundles from tracked recipes with {args.cli}; removed {removed} unreferenced assets")


if __name__ == "__main__":
    main()
