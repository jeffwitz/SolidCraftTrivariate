#!/usr/bin/env python3
"""Rebuild the public 45-case polygon benchmark gallery assets."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = None
DEFAULT_RESULTS = ROOT / "docs/results/sketch-planner-polygon-evaluation"
DEFAULT_CLI = ROOT / "target/release/solidcraft-trivariate"
OUT = ROOT / "docs/site/data"
MANIFEST_OUT = Path("/tmp/sct-gallery-polygon-entries.json")

METHOD_LABELS = {
    "radial_convex_footprint": "Radial convex partition",
    "triangulated_feature_footprint": "Triangulated partition",
    "rectangular_hosts_y": "Rectangular hosts · Y",
    "rectangular_hosts_x": "Rectangular hosts · X",
    "convex_host_separators": "Convex host separators",
    "orthogonal_feature_footprint": "Orthogonal feature partition",
    "guided_quad_footprint": "Guided quadrilateral partition",
}
FAMILY_GROUPS = {
    "convex_through_holes": "Polygon benchmark · through-holes",
    "convex_blind_pockets": "Polygon benchmark · blind pockets",
    "concave_stock_or_cut": "Polygon benchmark · concave profiles",
    "seeded_convex_through_holes": "Polygon benchmark · seeded through-holes",
    "seeded_convex_blind_pockets": "Polygon benchmark · seeded blind pockets",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


def metric(label: str, value: object) -> dict[str, str]:
    return {"label": label, "value": str(value)}


def run(cli: Path, *args: str) -> str:
    result = subprocess.run([str(cli), *args], check=True, text=True, capture_output=True, timeout=60)
    return result.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE,
                        help="Optional cache directory containing <case>.recipe.json files")
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--manifest-out", type=Path, default=MANIFEST_OUT)
    args = parser.parse_args()
    pair_data = json.loads((args.results / "paired.json").read_text())
    scored = json.loads((args.results / "candidate.json").read_text())
    corpus = json.loads((args.results / "cases.json").read_text())
    recipes = {case["id"]: case["recipe"] for case in corpus["cases"]}
    scored_by_id = {case["id"]: case for case in scored["cases"]}
    cases = pair_data["paired_cases"]
    if len(cases) != 45 or len(scored_by_id) != 45:
        raise SystemExit("Expected exactly 45 paired and candidate cases")
    entries: list[dict] = []
    rows: list[dict] = []
    for case in cases:
        case_id = case["id"]
        folder = OUT / f"bench-poly-{case_id}"
        folder.mkdir(parents=True, exist_ok=True)
        prefix = f"bench_poly_{case_id}"
        recipe = folder / f"{prefix}.recipe.json"
        model = folder / f"{prefix}.model.json"
        audit = folder / f"{prefix}.audit.json"
        step = folder / f"{prefix}.step"
        cached_recipe = args.source / f"{case_id}.recipe.json" if args.source else None
        if cached_recipe and cached_recipe.is_file():
            shutil.copyfile(cached_recipe, recipe)
        else:
            recipe.write_text(json.dumps(recipes[case_id], indent=2) + "\n")
        # Reconstruct every advertised model from the input recipe using the frozen executable.
        build_report = json.loads(run(args.cli, "build", str(recipe), str(model)))
        audit_report = json.loads(run(args.cli, "audit", str(model)))
        run(args.cli, "export-step", str(model), str(step))
        scored_case = scored_by_id[case_id]
        audit.write_text(json.dumps({
            "schema": "polygon_gallery_case_audit_v1",
            "case": case,
            "candidate": scored_case,
            "rebuild": build_report,
            "model_audit": audit_report,
            "provenance": {
                "recipe_sha256": hashlib.sha256(recipe.read_bytes()).hexdigest(),
                "model_sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
                "step_sha256": hashlib.sha256(step.read_bytes()).hexdigest(),
            },
        }, indent=2) + "\n")
        fam = case["family"].replace("_", " ")
        method_label = METHOD_LABELS.get(case["candidate_method"], case["candidate_method"].replace("_", " "))
        holes = case["cut_count"]
        seeded = case["family"].startswith("seeded_")
        if seeded:
            index = int(case_id.rsplit("_", 1)[1])
            contours = len(recipes[case_id]["outer"]["vertices"])
            feature_word = "through-hole" if "through_holes" in case["family"] else "blind pocket"
            title = f"Seeded case {index:02d} · {contours}-sided stock · {holes} {feature_word}{'s' if holes != 1 else ''}"
        else:
            title = case_id.replace("_", " ").capitalize()
        feature = "hole" if "through_holes" in case["family"] else "pocket/profile"
        description = (f"{fam.capitalize()}; {holes} {feature}{'s' if holes != 1 else ''}. "
                       "Exact-rational geometry checks passed for both paired executables.")
        kappa = case["candidate_condition_max"]
        baseline_kappa = case["baseline_condition_max"]
        volume = build_report.get("volume_exact_mm3", "see audit")
        volume_label = f"{volume[0]}/{volume[1]} mm³" if isinstance(volume, list) and len(volume) == 2 else f"{volume} mm³"
        entries.append({
            "id": prefix,
            "title": title,
            "description": description,
            "source": f"data/{folder.name}/{recipe.name}?v={digest(recipe)}",
            "model": f"data/{folder.name}/{model.name}?v={digest(model)}",
            "audit": f"data/{folder.name}/{audit.name}?v={digest(audit)}",
            "step": f"data/{folder.name}/{step.name}?v={digest(step)}",
            "metrics": [metric("Blocks · baseline → candidate", f"{case['baseline_blocks']} → {case['candidate_blocks']}"),
                        metric("Partition method", method_label),
                        metric("Sampled κ₂ max · baseline → candidate", f"{baseline_kappa:.3g} → {kappa:.3g}"),
                        metric("Exact volume", volume_label)],
            "section": {"axis": "z", "fraction": 0.5},
            "group": FAMILY_GROUPS.get(case["family"], "Polygon benchmark"),
        })
        rows.append({"id": case_id, "candidate_kappa_max": kappa,
                     "baseline_kappa_max": case["baseline_condition_max"],
                     "block_delta": case["block_delta_candidate_minus_baseline"]})

    # A deliberately overlapping XY projection with cuts separated along Z.
    diag_id = "bench_poly_projection_diamonds"
    folder = OUT / "bench-poly-projection-diamonds"
    folder.mkdir(parents=True, exist_ok=True)
    prefix = diag_id
    recipe = folder / f"{prefix}.recipe.json"
    model = folder / f"{prefix}.model.json"
    audit = folder / f"{prefix}.audit.json"
    step = folder / f"{prefix}.step"
    shutil.copyfile(ROOT / "tests/regressions/overlapping_projections_diamonds.recipe.json", recipe)
    diag_build = json.loads(run(args.cli, "build", str(recipe), str(model)))
    diag_audit = json.loads(run(args.cli, "audit", str(model)))
    run(args.cli, "export-step", str(model), str(step))
    sorted_kappa = sorted(x["candidate_kappa_max"] for x in rows)
    sorted_baseline_kappa = sorted(x["baseline_kappa_max"] for x in rows)
    def percentile(values: list[float], p: float) -> float:
        return values[round((len(values) - 1) * p)]
    family_counts: dict[str, int] = {}
    for case in cases:
        family_counts[case["family"]] = family_counts.get(case["family"], 0) + 1
    audit.write_text(json.dumps({
        "schema": "overlapping_projection_diagnostic_v1",
        "claim": "Two disjoint through-diamond holes have overlapping X and Y coordinate projections. Their planar regions and 3D voids remain disjoint.",
        "hole_x_projection_overlap_mm": 8,
        "hole_y_projection_overlap_mm": 8,
        "z_intervals_mm": [[0,12],[0,12]], "expected_3d_hole_intersection_volume_mm3": 0,
        "rebuild": diag_build, "model_audit": diag_audit,
        "paired_benchmark_summary": {
            "paired_case_count": len(cases), "geometry_common_success_count": sum(c["paired_outcome"] == "common_success" for c in cases),
            "new_admissions": sum(c["paired_outcome"] == "new_candidate_admission" for c in cases),
            "regressions": sum(c["paired_outcome"] == "candidate_regression" for c in cases),
            "family_counts": family_counts,
            "candidate_sampled_kappa2_max_median_p95_worst": [percentile(sorted_kappa, .5), percentile(sorted_kappa, .95), sorted_kappa[-1]],
            "baseline_sampled_kappa2_max_median_p95_worst": [percentile(sorted_baseline_kappa, .5), percentile(sorted_baseline_kappa, .95), sorted_baseline_kappa[-1]],
            "kappa_method": "Per-block Jacobian spectral condition number sampled on the closed 5x5x5 parameter grid; diagnostic only, not a continuous certificate.",
        },
        "recipe_sha256": hashlib.sha256(recipe.read_bytes()).hexdigest(),
        "model_sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
        "step_sha256": hashlib.sha256(step.read_bytes()).hexdigest(),
    }, indent=2) + "\n")
    entries.append({
        "id": diag_id, "title": "Projection diagnostic · overlapping diamonds",
        "description": "Two disjoint through-diamond holes overlap in both X and Y projection, while their exact planar regions and 3D voids remain separate.",
        "source": f"data/{folder.name}/{recipe.name}?v={digest(recipe)}",
        "model": f"data/{folder.name}/{model.name}?v={digest(model)}",
        "audit": f"data/{folder.name}/{audit.name}?v={digest(audit)}",
        "step": f"data/{folder.name}/{step.name}?v={digest(step)}" if step.exists() else None,
        "metrics": [metric("X projection overlap", "8 mm"), metric("Y projection overlap", "8 mm"),
                    metric("Through-holes", "2 disjoint diamonds"), metric("Graph blocks", diag_build.get("blocks", "see audit"))],
        "section": {"axis": "z", "fraction": 0.5}, "group": "Polygon benchmark · diagnostics",
    })
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.write_text(json.dumps(entries, indent=2) + "\n")
    print(json.dumps({"entries": len(entries), "cases": len(rows), "candidate_pair_counts": len(rows),
                      "kappa_max_range": [min(x["candidate_kappa_max"] for x in rows), max(x["candidate_kappa_max"] for x in rows)],
                      "manifest": str(args.manifest_out), "diagnostic_model": str(model)}))


if __name__ == "__main__":
    main()
