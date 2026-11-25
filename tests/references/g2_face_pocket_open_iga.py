#!/usr/bin/env python3
"""Qualify the open-side S3 pocket's C0 space and material maps."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/references"))

import g2_face_pocket_iga as iga


def qualify(model_path: Path, operation_path: Path, source_model_path: Path,
            source_recipe_path: Path, cli_path: Path, package_path: Path) -> dict:
    model = json.loads(model_path.read_text())
    source_model = json.loads(source_model_path.read_text())
    operation_recipe = json.loads(operation_path.read_text())
    geometry = iga.operation_geometry(operation_recipe)
    audit = iga.run_command(cli_path, ["audit", str(model_path)])
    export = iga.run_command(cli_path, ["export-nurbs-space", str(model_path), str(package_path)])
    package = json.loads(package_path.read_text())
    space = iga.load_module("g2_face_open_space", ROOT / "tests/references/nurbs_space_oracle.py")
    patch = iga.load_module("g2_face_open_equilibrium", ROOT / "tests/references/nurbs_patch_equilibrium.py")
    trace = space.package_checks(model, package)
    equilibrium = [patch.integrate(package, order) for order in (4, 8, 16)]
    volume = iga.map_volume_check(model, source_model, source_recipe_path, geometry)
    preservation = iga.source_image_preservation(model, source_model, geometry)

    classifier = iga.load_module("g2_face_open_classifier", ROOT / "tests/references/g2_cross_map_seams.py")
    junctions = iga.load_module("g2_face_open_junctions", ROOT / "tests/references/g2_cross_junctions.py")
    recipe = json.loads(source_recipe_path.read_text())
    profile = next(s for s in recipe["sketches"] if s["direction"] == [0, 0, 1])
    curve = junctions.sketch_curve(profile)
    xlo, xhi = geometry["x_mm"]
    ylo, yhi = geometry["y_mm"]
    ztop, zfloor = geometry["support_z_mm"], geometry["floor_z_mm"]
    ymid = (ylo+yhi)/2
    points = [
        ("open_stock_edge_inside_pocket", [xlo, ymid, (ztop+zfloor)/2], "OUT"),
        ("pocket_void_interior", [(xlo+xhi)/2, ymid, (ztop+zfloor)/2], "OUT"),
        ("pocket_floor", [(xlo+xhi)/2, ymid, zfloor], "ON"),
        ("stock_edge_outside_pocket", [xlo, yhi+1, (ztop+zfloor)/2], "ON"),
        ("new_wall_ymin", [(xlo+xhi)/2, ylo, (ztop+zfloor)/2], "ON"),
        ("new_wall_ymax", [(xlo+xhi)/2, yhi, (ztop+zfloor)/2], "ON"),
    ]
    p29 = junctions.jets.rational_jets(curve, 2.9, "right")["position_mm"]
    p31 = junctions.jets.rational_jets(curve, 3.1, "right")["position_mm"]
    points.extend([
        ("old_profile_edge_above_cut_preserved", [p29[0], p29[1], ztop], "ON"),
        ("old_profile_edge_below_cut_removed", [p31[0], p31[1], ztop], "OUT"),
    ])
    classifications = []
    for name, point, expected in points:
        actual = classifier.model_classify(model, point)["classification"]
        classifications.append({"name": name, "point_mm": point, "expected": expected,
                                "actual": actual, "agreement": actual == expected})

    q16 = equilibrium[-1]["relative_residual_linf"]
    status = "verified" if (
        audit.get("status") == "valid"
        and trace.get("status", "verified") == "verified"
        and volume["status"] == "verified"
        and preservation["status"] == "verified"
        and all(item["agreement"] for item in classifications)
        and q16 <= 1e-8
    ) else "failed"
    return {
        "schema": "g2_face_pocket_open_side_iga_v1",
        "status": status,
        "model": str(model_path), "model_sha256": iga.sha256(model_path),
        "source_model": str(source_model_path), "source_model_sha256": iga.sha256(source_model_path),
        "operation_recipe": str(operation_path),
        "cli": str(cli_path), "cli_sha256": iga.sha256(cli_path),
        "block_count": len(model.get("blocks", [])),
        "build_audit": audit, "space_export": export,
        "space_package_sha256": iga.sha256(package_path),
        "exact_trace_affine_audit": trace,
        "constant_stress_equilibrium_orders": equilibrium,
        "q16_relative_residual_threshold": 1e-8,
        "volume_map_audit": volume,
        "source_map_image_preservation": preservation,
        "open_side_classification_witnesses": {
            "status": "verified" if all(item["agreement"] for item in classifications) else "failed",
            "checks": classifications,
            "method": "inverse Newton over the actual maps with exported exterior labels",
        },
        "limits": ["This is the registered x=stock xmin open-side witness, separate from the interior principal pocket.",
                   "Volume uses Q32 map integration and an independent source-profile Q32/Q64 integral.",
                   "Point classifications and source-map samples are finite diagnostics, not a global image-inclusion proof."],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--operation-recipe", type=Path, required=True)
    parser.add_argument("--source-model", type=Path, required=True)
    parser.add_argument("--source-recipe", type=Path, default=ROOT / "tests/references/g2_cross_sketches.recipe.json")
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--package", type=Path, default=ROOT / "build/g2-face-pocket-open-iga.q3.json")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/results/g2-face-pocket-open-side-iga.json")
    args = parser.parse_args()
    report = qualify(args.model.resolve(), args.operation_recipe.resolve(), args.source_model.resolve(),
                     args.source_recipe.resolve(), args.cli.resolve(), args.package.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "blocks": report["block_count"],
                      "q16_relative_residual": report["constant_stress_equilibrium_orders"][-1]["relative_residual_linf"],
                      "q32_volume_relative_error": report["volume_map_audit"]["relative_total_error"],
                      "output": str(args.output)}))
    return 0 if report["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
