#!/usr/bin/env python3
"""Audit S3 NURBS pocket maps through the C0 Q3 space and verify atomic refusal."""
import argparse
import hashlib
import importlib.util
import json
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load oracle {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pocket_volume_reference(recipe_path: Path, xlo: float, xhi: float,
                            ylo: float, yhi: float, depth: float) -> dict:
    import numpy as np
    junctions = load_module("g2_face_pocket_source_curve", ROOT / "tests/references/g2_cross_junctions.py")
    recipe = json.loads(recipe_path.read_text())
    free = next(s for s in recipe["sketches"] if s["direction"] == [0, 0, 1])
    curve = junctions.sketch_curve(free)
    rows = []
    for order in (32, 64):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        area = 0.0
        for node, weight in zip(nodes, weights):
            y = (ylo+yhi)*0.5 + (yhi-ylo)*0.5*node
            roots = junctions.find_level_roots(curve, float(y), axis=1)
            xs = sorted(junctions.jets.rational_jets(curve, t, "right")["position_mm"][0]
                        for t in roots)
            if len(xs) < 2:
                raise ArithmeticError(f"source profile has fewer than two y={y} crossings")
            material_width = max(0.0, min(xhi, xs[0])-xlo)
            area += float(weight)*(yhi-ylo)*0.5*material_width
        rows.append({"order": order, "removed_area_mm2": area,
                     "removed_volume_mm3": area*depth})
    return {"method": "source periodic NURBS horizontal-ray intersections with Gauss-Legendre in y",
            "x_mm": [xlo,xhi], "y_mm": [ylo,yhi], "depth_mm": depth,
            "quadrature": rows,
            "relative_order_delta": abs(rows[1]["removed_volume_mm3"]-rows[0]["removed_volume_mm3"])
                / max(abs(rows[1]["removed_volume_mm3"]),1.0)}


def operation_geometry(operation_recipe: dict) -> dict:
    operation = operation_recipe["operation"]
    primitive = operation["loop"]["primitive"]
    as_float = lambda pair: float(pair[0])/float(pair[1])
    xmin, ymin = [as_float(v) for v in primitive["min"]]
    xmax, ymax = [as_float(v) for v in primitive["max"]]
    z0 = float(operation["frame"]["origin_mm"][2])
    depth = float(operation["depth_mm"])
    z1 = z0 + float(operation["direction"][2])*depth
    return {"x_mm":[xmin,xmax],"y_mm":[ymin,ymax],"support_z_mm":z0,
            "floor_z_mm":z1,"depth_mm":depth}


def map_volume_check(model: dict, source_model: dict, source_recipe: Path,
                     operation_geometry: dict) -> dict:
    volume_oracle = load_module("g2_face_pocket_volume_oracle", ROOT / "tests/references/g2_cross_model_volume.py")
    source_layers = volume_oracle.integrate(source_model, 32)
    model_layers = volume_oracle.integrate(model, 32)
    xlo,xhi = operation_geometry["x_mm"]
    ylo,yhi = operation_geometry["y_mm"]
    zlo,zhi = sorted((operation_geometry["support_z_mm"],operation_geometry["floor_z_mm"]))
    depth = operation_geometry["depth_mm"]
    pocket = pocket_volume_reference(source_recipe, xlo, xhi, ylo, yhi, depth)
    removal = pocket["quadrature"][-1]["removed_volume_mm3"]
    expected_layers = []
    for (z0,z1), measured in sorted(model_layers.items()):
        expected = 0.0
        for (s0,s1), source_volume in source_layers.items():
            overlap = max(0.0,min(z1,s1)-max(z0,s0))
            expected += source_volume*overlap/(s1-s0)
        cut_depth = max(0.0,min(z1,zhi)-max(z0,zlo))
        expected -= removal*cut_depth/depth
        expected_layers.append({"z_mm":[float(z0),float(z1)],"volume_mm3":measured,
                                "expected_volume_mm3":expected,
                                "relative_error":abs(measured-expected)/max(abs(expected),1.0)})
    total = sum(model_layers.values())
    source_total = sum(source_layers.values())
    expected_total = source_total-removal
    return {"status":"verified" if expected_layers
            and all(row["relative_error"]<1e-8 for row in expected_layers)
            and abs(total-expected_total)/expected_total<1e-8
            and pocket["relative_order_delta"]<1e-8 else "failed",
            "method":"Q32 determinant integration over actual S3 tensor-product Bezier maps",
            "source_model_volume_q32_mm3":source_total,
            "s3_model_volume_q32_mm3":total,
            "expected_removed_volume_mm3":removal,
            "expected_s3_volume_mm3":expected_total,
            "relative_total_error":abs(total-expected_total)/expected_total,
            "max_relative_layer_error":max(row["relative_error"] for row in expected_layers),
            "layers":expected_layers,"independent_source_pocket_integral":pocket,
            "threshold_relative":1e-8}


def pocket_boundary_witnesses(model: dict, source_recipe: Path,
                              operation_geometry: dict) -> dict:
    classifier = load_module("g2_face_pocket_map_classifier", ROOT / "tests/references/g2_cross_map_seams.py")
    junctions = load_module("g2_face_pocket_junctions", ROOT / "tests/references/g2_cross_junctions.py")
    recipe = json.loads(source_recipe.read_text())
    free = next(s for s in recipe["sketches"] if s["direction"] == [0, 0, 1])
    curve = junctions.sketch_curve(free)
    p29 = junctions.jets.rational_jets(curve, 2.9, "right")["position_mm"]
    p31 = junctions.jets.rational_jets(curve, 3.1, "right")["position_mm"]
    xlo,xhi = operation_geometry["x_mm"]
    ylo,yhi = operation_geometry["y_mm"]
    zsupport,zfloor = operation_geometry["support_z_mm"],operation_geometry["floor_z_mm"]
    xleft=xlo-0.1
    xinside=xlo+(xhi-xlo)*0.25
    ymiddle=(ylo+yhi)*0.5
    zmiddle=(zsupport+zfloor)*0.5
    points = [
        ("new_wall_xmin", [xlo, ymiddle, zmiddle], "ON"),
        ("material_left_of_new_wall", [xleft, ymiddle, zmiddle], "IN"),
        ("pocket_void_right_of_new_wall", [xlo+0.1, ymiddle, zmiddle], "OUT"),
        ("new_floor", [xinside, ymiddle, zfloor], "ON"),
        ("old_floor_preserved_x_left_of_pocket", [xleft, ymiddle, zsupport], "ON"),
        ("old_floor_y_above_pocket_preserved", [xinside, yhi+1.0, zsupport], "ON"),
        ("old_profile_edge_y_gt_pocket_preserved", [p29[0], p29[1], zsupport], "ON"),
        ("old_profile_edge_y_lt_pocket_recut", [p31[0], p31[1], zsupport], "OUT"),
        ("old_cap_y_lt_pocket_recut", [xinside, ymiddle, zsupport], "OUT"),
    ]
    rows = []
    for name, point, expected in points:
        actual = classifier.model_classify(model, point)
        rows.append({"name": name, "point_mm": point, "expected": expected,
                     "actual": actual["classification"], "candidate_maps": actual["candidate_maps"],
                     "agreement": actual["classification"] == expected})
    return {"status": "verified" if all(row["agreement"] for row in rows) else "failed",
            "counts": {"samples": len(rows), "agreements": sum(row["agreement"] for row in rows),
                       "disagreements": sum(not row["agreement"] for row in rows)},
            "checks": rows,
            "classification_method": "inverse Newton in every actual S3 volume map; exterior labels distinguish ON from internal shared faces"}


def source_map_preservation(model: dict, source_model: dict) -> dict:
    source = {block["id"]: block for block in source_model.get("blocks", [])}
    target = {block["id"]: block for block in model.get("blocks", [])}
    unchanged, changed, missing = [], [], []
    for block_id, old in source.items():
        new = target.get(block_id)
        if new is None:
            missing.append(block_id)
        elif old.get("degrees") == new.get("degrees") and old.get("controls") == new.get("controls"):
            unchanged.append(block_id)
        else:
            changed.append(block_id)
    unexpected = [block_id for block_id in changed+missing if "_s4_" not in block_id]
    return {"status": "diagnostic_only",
            "source_blocks": len(source), "target_blocks": len(target),
            "unchanged_exact_source_blocks": unchanged,
            "changed_source_blocks": changed, "omitted_source_blocks": missing,
            "unexpected_non_s4_changes": unexpected,
            "contract": "IDs and control nets are diagnostic; geometric image is checked by point samples and layer volume because radial-guide refinement can change every control count"}


def source_image_preservation(model: dict, source_model: dict,
                              operation_geometry: dict) -> dict:
    classifier = load_module("g2_face_pocket_image_classifier", ROOT / "tests/references/g2_cross_map_seams.py")
    space = load_module("g2_face_pocket_image_map", ROOT / "tests/references/nurbs_space_oracle.py")
    xlo,xhi = operation_geometry["x_mm"]
    ylo,yhi = operation_geometry["y_mm"]
    zlo,zhi = sorted((operation_geometry["support_z_mm"],operation_geometry["floor_z_mm"]))
    params = (0.25,0.5,0.75)
    checks = []
    skipped = 0
    for block in source_model.get("blocks", []):
        for u in params:
            for v in params:
                for w in params:
                    p = tuple(float(c) for c in space.map_from_json(block,tuple(
                        Fraction(str(t)) for t in (u,v,w))))
                    eps = 1e-9
                    in_pocket_envelope = (xlo-eps <= p[0] <= xhi+eps
                                          and ylo-eps <= p[1] <= yhi+eps
                                          and zlo-eps <= p[2] <= zhi+eps)
                    if in_pocket_envelope:
                        skipped += 1
                        continue
                    actual = classifier.model_classify(model,p)["classification"]
                    checks.append({"source_block":block["id"],"local_uvw":[u,v,w],
                                   "point_mm":p,"expected":"IN","actual":actual,
                                   "agreement":actual=="IN"})
    failures = [row for row in checks if not row["agreement"]]
    return {"status":"verified" if not failures else "failed",
            "source_map_count":len(source_model.get("blocks",[])),
            "sample_rule":"3x3x3 interior Bernstein nodes per source map; skip only points strictly inside the new pocket envelope",
            "sample_count":len(checks),"expected_removed_envelope_samples":skipped,
            "preserved_material_agreements":len(checks)-len(failures),
            "mismatches":failures[:32],
            "limits":["Samples check image coverage of the old material maps, not a global inclusion proof.",
                      "The independent Q32/Q64 volume balance checks total and per-layer material measure."]}


def run_command(cli: Path, args: list[str]) -> dict:
    process = subprocess.run([str(cli), *args], cwd=ROOT, capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError(f"CLI {' '.join(args)} failed: {process.stderr or process.stdout}")
    return json.loads(process.stdout)


def qualify(model_path: Path, operation_recipe_path: Path,
            source_model_path: Path, source_recipe: Path,
            cli: Path, package_path: Path,
            refusal_recipe: Path | None, refusal_target: Path | None) -> dict:
    space = load_module("g2_face_pocket_space_oracle", ROOT / "tests/references/nurbs_space_oracle.py")
    patch = load_module("g2_face_pocket_patch_oracle", ROOT / "tests/references/nurbs_patch_equilibrium.py")
    audit = run_command(cli, ["audit", str(model_path)])
    model = json.loads(model_path.read_text())
    operation_geometry_data = operation_geometry(json.loads(operation_recipe_path.read_text()))
    source_model = json.loads(source_model_path.read_text())
    export = run_command(cli, ["export-nurbs-space", str(model_path), str(package_path)])
    package = json.loads(package_path.read_text())
    trace = space.package_checks(model, package)
    equilibrium = [patch.integrate(package, order) for order in (4, 8, 16)]
    volume = map_volume_check(model, source_model, source_recipe, operation_geometry_data)
    witnesses = pocket_boundary_witnesses(model, source_recipe, operation_geometry_data)
    preservation = source_map_preservation(model, source_model)
    image_preservation = source_image_preservation(model, source_model, operation_geometry_data)
    if equilibrium[-1]["relative_residual_linf"] > 1e-8:
        raise AssertionError("q16 constant-stress equilibrium exceeded 1e-8")
    faces = Counter(face["feature"] for face in model.get("design_faces", []))

    refusal = None
    if refusal_recipe is not None:
        target = refusal_target or package_path.with_name("preserved-model.json")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(model_path.read_bytes())
        before = sha256(target)
        process = subprocess.run([str(cli), "build", str(refusal_recipe), str(target)],
                                 cwd=ROOT, capture_output=True, text=True)
        after = sha256(target) if target.exists() else None
        refusal = {"recipe": str(refusal_recipe), "target": str(target),
                   "exit_code": process.returncode,
                   "stdout": process.stdout, "stderr": process.stderr,
                   "sha256_before": before, "sha256_after": after,
                   "preserved_atomically": process.returncode != 0 and before == after}

    return {"status": "verified" if audit.get("status") == "valid"
            and trace.get("status", "verified") == "verified"
            and volume["status"] == "verified"
            and witnesses["status"] == "verified"
            and image_preservation["status"] == "verified"
            and equilibrium[-1]["relative_residual_linf"] <= 1e-8
            and (refusal is None or refusal["preserved_atomically"]) else "failed",
            "model": str(model_path), "model_sha256": sha256(model_path),
            "source_model": str(source_model_path), "source_model_sha256": sha256(source_model_path),
            "source_recipe": str(source_recipe),
            "operation_recipe": str(operation_recipe_path),
            "operation_geometry": operation_geometry_data,
            "cli": str(cli), "cli_sha256": sha256(cli),
            "build_audit": audit, "space_export": export,
            "space_package_sha256": sha256(package_path),
            "block_count": len(model.get("blocks", [])),
            "design_face_feature_counts": dict(sorted(faces.items())),
            "exact_trace_affine_audit": trace,
            "volume_map_audit": volume,
            "boundary_and_old_seam_witnesses": witnesses,
            "source_map_preservation": preservation,
            "source_map_image_preservation": image_preservation,
            "constant_stress_equilibrium_orders": equilibrium,
            "q16_relative_residual_threshold": 1e-8,
            "atomic_refusal": refusal,
            "limits": ["The volume maps and Q3 space are audited independently of any external Boolean kernel.",
                       "The refusal check qualifies output preservation for only the supplied counterexample."]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--operation-recipe", type=Path, required=True)
    parser.add_argument("--source-model", type=Path, required=True,
                        help="previous S2 delivered model used as an additive source map oracle")
    parser.add_argument("--source-recipe", type=Path,
                        default=ROOT / "tests/references/g2_cross_sketches.recipe.json")
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--package", type=Path)
    parser.add_argument("--refusal-recipe", type=Path)
    parser.add_argument("--refusal-target", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/results/g2-face-pocket-iga.json")
    args = parser.parse_args()
    package = args.package or ROOT / "build/g2-face-pocket-iga.q3.json"
    report = qualify(args.model.resolve(), args.operation_recipe.resolve(), args.source_model.resolve(), args.source_recipe.resolve(),
                     args.cli.resolve(), package.resolve(),
                     args.refusal_recipe.resolve() if args.refusal_recipe else None,
                     args.refusal_target.resolve() if args.refusal_target else None)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "blocks": report["block_count"],
                      "q16_relative_residual": report["constant_stress_equilibrium_orders"][-1]["relative_residual_linf"],
                      "q32_volume_relative_error": report["volume_map_audit"]["relative_total_error"],
                      "atomic_refusal": None if report["atomic_refusal"] is None
                      else report["atomic_refusal"]["preserved_atomically"],
                      "output": str(args.output)}))
    return 0 if report["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
