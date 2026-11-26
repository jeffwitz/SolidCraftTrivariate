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
from evaluate_polygon_extrusion_quality import summarize_model


ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/site"
MANIFEST = SITE / "experiments.json"
POLYGON_RESULTS = ROOT / "docs/results/convex-host-rings.json"
POLYGON_COMPARISON = ROOT / "docs/results/sketch-nurbs-comparison.json"
POLYGON_EVALUATION = ROOT / "build/sketch-nurbs-evaluation.final.json"
POLYGON_MODELS = ROOT / "build/sketch-nurbs-evaluation"
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


def polygon_selections() -> tuple[dict[str, dict], dict[str, dict]]:
    """Select the smallest exact NURBS assembly from the frozen two-mode study."""
    comparison = read_json(POLYGON_COMPARISON)
    evaluation = read_json(POLYGON_EVALUATION) if POLYGON_EVALUATION.is_file() else {}
    evaluated = {case["id"]: case for case in evaluation.get("cases", [])}
    selections = {}
    for row in comparison["cases"]:
        modes = row["modes"]
        selected = min(modes, key=lambda mode: (
            int(modes[mode]["patches"]), int(modes[mode]["cells"]), mode != "assemble"
        ))
        detail = evaluated.get(row["case"])
        model_path = POLYGON_MODELS / f"{row['case']}.{selected}.model.json"
        model_sha = sha256(model_path.read_bytes()) if model_path.is_file() else None
        if detail is not None and selected in detail and model_sha:
            if model_sha != detail[selected]["model"]["model_sha256"]:
                raise ValueError(f"frozen NURBS model digest mismatch: {row['case']} / {selected}")
        selections[row["case"]] = {
            "selected_mode": selected,
            "selected": modes[selected],
            "alternatives": modes,
            "volume_mm3": row["volume_mm3"],
            "selected_evaluation_model_sha256": model_sha,
            "evaluation_cli_sha256": comparison["kernel_binary_sha256"],
            "comparison_report_sha256": sha256(POLYGON_COMPARISON.read_bytes()),
        }
    if len(selections) != 46:
        raise ValueError(f"Expected 46 NURBS partition comparisons; found {len(selections)}")
    return selections, evaluated


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
            "scope": "Latest frozen bracket benchmark record. Its STEP checks concern the benchmark candidate STEP, not the regenerated gallery STEP.",
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
                           polygon_selection: dict | None,
                           bracket_record: dict | None, condition: dict, height_factor: int = 1) -> None:
    method = model.get("partition_plan", {}).get("method", build_report.get("method", kernel_audit.get("method", model.get("schema", "native Rust builder"))))
    blocks = int(kernel_audit.get("blocks", len(model.get("blocks", []))))
    interfaces = int(kernel_audit.get("shared_faces", len(model.get("topology", {}).get("interfaces", []))))
    sampled_kappa = condition["maximum_condition_2"]
    if polygon_selection:
        selected = polygon_selection["selected"]
        if (int(selected["patches"]) != int(kernel_audit["nurbs_patch_count"])
                or int(selected["cells"]) != blocks
                or int(selected["interfaces_verified"]) != interfaces):
            raise ValueError(f"selected NURBS partition differs from the native model: {entry['id']}")
        if abs(sampled_kappa - float(selected["kappa_3d"])) > 1e-10 * max(1.0, sampled_kappa):
            raise ValueError(f"selected NURBS condition measurement differs for {entry['id']}")
        current_volume = Fraction(*kernel_audit["volume_exact_mm3"])
        if current_volume != Fraction(polygon_selection["volume_mm3"]):
            raise ValueError(f"scaled polygon volume differs from the frozen NURBS comparison: {entry['id']}")
        retained = [row for row in entry.get("metrics", []) if row.get("label") in {
            "X projection overlap", "Y projection overlap", "Through-holes",
        }] if entry.get("id") == "bench_poly_projection_diamonds" else []
        entry["metrics"] = retained + [
            {"label": "NURBS partition choice", "value": polygon_selection["selected_mode"]},
            {"label": "Partition method", "value": method},
            {"label": "Rust model blocks", "value": blocks},
            {"label": "Rust shared interfaces", "value": interfaces},
            {"label": "Sampled κ₂ max", "value": f"{sampled_kappa:.6g}"},
            {"label": "Certified determinant lower bound", "value": format_bound(kernel_audit.get("determinant_lower_bound_mm3"))},
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
    if not (polygon_selection or polygon_record or bracket_record or entry.get("id", "").startswith("bench-bracket-")):
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


def rebuild(entries: list[dict], cli: Path, kernel_revision: str | None = None) -> None:
    if not cli.is_file():
        raise FileNotFoundError(f"Rust CLI not found: {cli}; run cargo build --release --locked --offline")
    kernel_sha = sha256(cli.read_bytes())
    revision = kernel_revision or subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", "src", "Cargo.toml", "Cargo.lock"],
        cwd=ROOT, check=True, text=True, capture_output=True,
    ).stdout.strip()
    polygon_records, polygon_run = polygon_latest()
    polygon_choice, _ = polygon_selections()
    frozen = read_json(ROOT / "docs/results/sketch-planner-polygon-evaluation/cases.json")
    original_recipes = {case["id"]: case["recipe"] for case in frozen["cases"]}
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
        polygon_selection = polygon_choice.get(case_id) if case_id.startswith("bench_poly_") else None
        bracket_name = case_id.removeprefix("bench-bracket-").replace("-", "_")
        bracket_record = bracket_records.get(bracket_name)

        with tempfile.TemporaryDirectory(prefix="gallery-regenerate-", dir=model.parent) as temp_name:
            temp = Path(temp_name)
            staged_model, staged_step, staged_audit = temp / "model.json", temp / "reference.step", temp / "audit.json"
            staged_recipe = temp / "recipe.json"
            build_recipe = recipe
            if polygon_selection:
                expected_mode = polygon_selection["selected_mode"]
                current = read_json(recipe)
                if polygon_record:
                    original = original_recipes[polygon_record["id"]]
                    height_factor, remainder = divmod(current["height_mm"], original["height_mm"])
                    expected = dict(original, nurbs_partition=expected_mode,
                                    height_mm=original["height_mm"] * height_factor,
                                    cuts=[dict(cut, z_interval_mm=[z * height_factor for z in cut["z_interval_mm"]])
                                          for cut in original["cuts"]])
                    legacy = dict(expected, nurbs_partition=current.get("nurbs_partition"))
                    if remainder or height_factor < 1 or current not in (expected, legacy):
                        raise ValueError(f"gallery recipe differs from the frozen planar benchmark: {case_id}")
                else:
                    height_factor = 1
                    expected = dict(current, nurbs_partition=expected_mode)
                    legacy = dict(expected, nurbs_partition=current.get("nurbs_partition"))
                    if current not in (expected, legacy):
                        raise ValueError(f"projection-only recipe differs from its tracked comparison: {case_id}")
                write_json(staged_recipe, expected)
                build_recipe = staged_recipe
            build_report = json.loads(run(cli, "build", build_recipe, staged_model))
            kernel_audit = json.loads(run(cli, "audit", staged_model))
            run(cli, "export-step", staged_model, staged_step)
            model_bytes = staged_model.read_bytes()
            if polygon_selection:
                selected_path = POLYGON_MODELS / f"{case_id}.{polygon_selection['selected_mode']}.model.json"
                if selected_path.is_file():
                    if model_bytes != selected_path.read_bytes():
                        raise ValueError(f"native model differs byte-for-byte from frozen {polygon_selection['selected_mode']} evaluation: {case_id}")
                    expected_sha = polygon_selection["selected_evaluation_model_sha256"]
                    if expected_sha and sha256(model_bytes) != expected_sha:
                        raise ValueError(f"native model digest differs from frozen evaluation: {case_id}")
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
                polygon_record, kernel_sha, revision, build_recipe,
                model_bytes, published_step, condition, bracket_records, bracket_report,
            )
            height_factor = 1
            if polygon_record:
                original = original_recipes[polygon_record["id"]]
                current = read_json(build_recipe)
                height_factor, remainder = divmod(current["height_mm"], original["height_mm"])
                old_volume = Fraction(*polygon_record["volume_exact_cli_mm3"])
                if Fraction(*kernel_audit["volume_exact_mm3"]) != height_factor * old_volume:
                    raise ValueError(f"scaled polygon volume mismatch: {case_id}")
                audit_doc["case"] = {"id": polygon_record["id"], "family": polygon_record["family"]}
                audit_doc["frozen_reference"] = {
                    "height_mm": original["height_mm"], "gallery_height_factor": height_factor,
                    "tag": polygon_run["tag"], "report": "docs/results/convex-host-rings.json",
                }
            elif polygon_selection:
                if Fraction(*kernel_audit["volume_exact_mm3"]) != Fraction(polygon_selection["volume_mm3"]):
                    raise ValueError(f"projection-control volume differs from the frozen comparison: {case_id}")
            if polygon_selection:
                audit_doc["partition_selection"] = dict(polygon_selection)
                audit_doc["partition_selection"].update({
                    "criterion": "lexicographic (NURBS patches, Bézier cells), with assemble preferred on ties",
                    "native_model_sha256": sha256(model_bytes),
                    "evaluation_report": "docs/results/sketch-nurbs-comparison.json",
                    "optional_evaluation_model_sha256": polygon_selection["selected_evaluation_model_sha256"],
                })
            refresh_entry_metadata(entry, build_report, kernel_audit, model_doc, polygon_record,
                                   polygon_selection, bracket_record, condition, height_factor)
            if case_id.startswith("bench_poly_"):
                height = read_json(recipe)["height_mm"]
                entry["metrics"].insert(0, {"label": "Stock thickness", "value": f"{height} mm"})
                patches = kernel_audit["nurbs_patch_count"]
                spans = kernel_audit["bezier_span_count"]
                entry["metrics"].insert(1, {"label": "NURBS volumes", "value": patches})
                for metric in entry["metrics"]:
                    if metric["label"] == "Rust model blocks":
                        metric["label"] = "Bézier cells"
                entry["description"] = (
                    f"{height} mm extrusion. {spans} Bézier cells form {patches} NURBS volumes."
                    if patches < spans else
                    f"{height} mm extrusion. Each of the {spans} Bézier cells forms a single-span NURBS volume."
                )
                planar = summarize_model(model_doc)
                audit_doc["extrusion_quality"] = planar
                for metric in entry["metrics"]:
                    if metric["label"] == "Sampled κ₂ max":
                        metric["label"] = "Sampled 3D κ₂ max"
                entry["metrics"].extend([
                    {"label": "Sampled in-plane κ₂ max", "value": f"{planar['planar_spectral_kappa']:.6g}"},
                    {"label": "Minimum corner sine", "value": f"{planar['minimum_corner_sine']:.6g}"},
                    {"label": "Maximum in-plane span / normal thickness", "value": f"{planar['thickness_direction_aspect']:.6g}"},
                ])
            sync_volume_metric(entry, audit_doc)
            write_json(staged_audit, audit_doc)

            # Do not replace a fresh STEP whose only difference is its generated timestamp.
            os.replace(staged_model, model)
            if polygon_selection:
                os.replace(staged_recipe, recipe)
            if not same_step_except_timestamp:
                os.replace(staged_step, step)
            os.replace(staged_audit, audit)

        for key in ("model", "source", "step", "audit"):
            path = entry[key].split("?", 1)[0]
            entry[key] = f"{path}?v={sha256(site_path(path).read_bytes())[:12]}"

    refresh_patch_representations(entries)


def refresh_patch_representations(entries: list[dict]) -> None:
    """Publish exact disjoint NURBS covers for the curved and singleton cases."""
    from gallery_curved_patches import build_catalog

    catalog = build_catalog(MANIFEST)
    curved = catalog["entries"]
    linear = catalog["linear_patch_assemblies"]
    if len(curved) != 22 or len(linear) != 7:
        raise ValueError(f"Expected 22 curved and 7 singleton representations; got {len(curved)} and {len(linear)}")
    by_id = {row["id"]: (row, "curved") for row in curved}
    by_id.update({row["id"]: (row, "singleton_linear") for row in linear})
    if sum(row["source_block_count"] for row in curved) != 330:
        raise ValueError("Curved exact extraction source block count changed")

    for entry in entries:
        case_id = entry["id"]
        if case_id not in by_id:
            entry.pop("representation", None)
            continue
        row, family = by_id[case_id]
        model_path, recipe_path, audit_path = (site_path(entry[key]) for key in ("model", "source", "audit"))
        model_doc, audit_doc = read_json(model_path), read_json(audit_path)
        source_blocks = set(block["id"] for block in model_doc.get("blocks", []))
        covered = [cell for patch in row["patches"] for cell in patch.get("cells", [])]
        if (row["exact_extraction_count"] != row["source_block_count"]
                or len(covered) != len(source_blocks) or set(covered) != source_blocks
                or len(covered) != len(set(covered))):
            raise ValueError(f"Exact patch map does not cover every source block once: {case_id}")
        if any(interface.get("status") != "exact_match" for interface in row["interfaces"]):
            raise ValueError(f"Unverified patch-map interface: {case_id}")
        native_patches = int(audit_doc["model_audit"].get("nurbs_patch_count", row["source_block_count"]))
        if family == "singleton_linear" and (native_patches != row["source_block_count"]
                                                or row["patch_count"] >= native_patches):
            raise ValueError(f"Singleton cover does not reduce native NURBS patches: {case_id}")
        if row["model_sha256"] != sha256(model_path.read_bytes()) or row["recipe_sha256"] != sha256(recipe_path.read_bytes()):
            raise ValueError(f"Patch-map source digests are stale: {case_id}")

        payload = {
            "schema": "gallery_exact_patch_map_v1",
            "case_id": case_id,
            "representation_kind": family,
            "method": catalog["method"],
            "scope_note": catalog["scope_note"],
            "source_model_sha256": row["model_sha256"],
            "source_recipe_sha256": row["recipe_sha256"],
            "source_schema": row["source_schema"],
            "source_block_count": row["source_block_count"],
            "native_nurbs_patch_count": native_patches,
            "patch_count": row["patch_count"],
            "exact_extraction_count": row["exact_extraction_count"],
            "exact_interface_count": row["exact_interface_count"] if "exact_interface_count" in row else len(row["interfaces"]),
            "patches": row["patches"],
            "interfaces": row["interfaces"],
            "boundary_faces": row["boundary_faces"],
            "merged_source_blocks": row.get("merged_source_blocks", []),
            "errors": row["errors"],
        }
        path = model_path.with_name(model_path.stem + ".representation.json")
        write_json(path, payload)
        entry["representation"] = f"{path.relative_to(SITE)}?v={sha256(path.read_bytes())[:12]}"
        audit_doc["exact_patch_representation"] = {
            "schema": payload["schema"],
            "kind": family,
            "path": entry["representation"],
            "sha256": sha256(path.read_bytes()),
            "model_sha256": payload["source_model_sha256"],
            "recipe_sha256": payload["source_recipe_sha256"],
            "patch_count": payload["patch_count"],
            "source_block_count": payload["source_block_count"],
            "exact_extraction_count": payload["exact_extraction_count"],
            "exact_interface_count": payload["exact_interface_count"],
        }
        write_json(audit_path, audit_doc)
        audit_url = entry["audit"].split("?", 1)[0]
        entry["audit"] = f"{audit_url}?v={sha256(site_path(audit_url).read_bytes())[:12]}"
    normalize_cover_metrics(entries, by_id)


def normalize_cover_metrics(entries: list[dict], representation_rows: dict | None = None) -> None:
    """Use one primary patch count, while retaining native counts when assembly reduces them."""
    representation_rows = representation_rows or {}
    for entry in entries:
        audit = read_json(site_path(entry["audit"]))
        native = audit["model_audit"]
        native_patches = int(native.get("nurbs_patch_count", native.get("blocks", 0)))
        native_cells = int(native.get("blocks", 0))
        row_pair = representation_rows.get(entry["id"])
        cover_patches = int(row_pair[0]["patch_count"]) if row_pair else native_patches
        source_cells = int(row_pair[0]["source_block_count"]) if row_pair else native_cells
        other_metrics = []
        for metric in entry["metrics"]:
            label = metric["label"]
            if label in {"Exact merged NURBS cover patches", "Source Bézier cells"}:
                continue
            if label in {"Rust model blocks", "Bézier cells"}:
                continue
            if label == "NURBS volumes":
                continue
            if label == "Native NURBS volumes before assembly":
                continue
            other_metrics.append(metric)
        normalized = [
            {"label": "NURBS volumes", "value": cover_patches},
            {"label": "Bézier cells", "value": source_cells},
        ]
        if native_patches != cover_patches:
            normalized.append({"label": "Native NURBS volumes before assembly", "value": native_patches})
        selection = audit.get("partition_selection")
        if selection and not any(metric["label"] == "NURBS partition choice" for metric in other_metrics):
            other_metrics.insert(0, {"label": "NURBS partition choice", "value": selection["selected_mode"]})
        entry["metrics"] = normalized + other_metrics
        summary = (f"The displayed NURBS representation groups {source_cells} Bézier cells into "
                   f"{cover_patches} volumes." if cover_patches < source_cells else
                   f"The displayed NURBS representation uses {cover_patches} volumes across {source_cells} Bézier cells.")
        if entry["id"].startswith("bench_poly_"):
            thickness = next((metric["value"] for metric in entry["metrics"] if metric["label"] == "Stock thickness"), None)
            entry["description"] = (f"{thickness} extrusion. {summary}" if thickness else summary)
        else:
            old = entry.get("description", "")
            if not old.startswith("The displayed NURBS representation"):
                entry["description"] = summary + (" " + old if old else "")


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
    parser.add_argument("--kernel-revision", help="Override the source revision recorded in regeneration audits")
    args = parser.parse_args()
    entries = read_json(MANIFEST)
    if len(entries) != EXPECTED_IDS or len({entry["id"] for entry in entries}) != EXPECTED_IDS:
        raise ValueError(f"Expected {EXPECTED_IDS} unique published gallery entries")
    rebuild(entries, args.cli.resolve(), args.kernel_revision)
    removed = remove_unreferenced_data(entries)
    write_json(MANIFEST, entries)
    version_modules("app.js", ["renderer.js"])
    version_modules("compare.js", ["renderer.js", "comparison-data.json"])
    version_modules("compare.html", ["compare.js", "style.css", "compare.css", "experiments.json"])
    version_modules("index.html", ["app.js", "experiments.json", "style.css"])
    print(f"Regenerated {len(entries)} gallery bundles from tracked recipes with {args.cli}; removed {removed} unreferenced assets")


if __name__ == "__main__":
    main()
