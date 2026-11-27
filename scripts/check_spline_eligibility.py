#!/usr/bin/env python3
"""Independent S6 common-parent eligibility checks for two real recipes."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from fractions import Fraction as F
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_cli(cli: str, *args: str) -> str:
    return subprocess.run([cli, *args], check=True, text=True, capture_output=True).stdout


def build(cli: str, recipe_path: Path, output: Path) -> tuple[dict, dict]:
    recipe = json.loads(recipe_path.read_text())
    run_cli(cli, "build", str(recipe_path), str(output))
    model = json.loads(output.read_text())
    audit = json.loads(run_cli(cli, "audit", str(output)))
    assert audit["status"] == "valid"
    return recipe, model


def verify_box_cylinder(recipe: dict, model: dict) -> dict:
    stock = recipe["stock"]
    tool = recipe["tool"]
    assert recipe["operation"] == "difference"
    assert tool["kind"] == "cylinder" and tool["axis"] == "z"
    lo = tuple(map(F, stock["min_mm"]))
    hi = tuple(map(F, stock["max_mm"]))
    cx, cy = map(F, tool["center_mm"])
    radius = F(tool["radius_mm"])
    z0, z1 = map(F, tool["interval_mm"])
    assert radius > 0
    assert lo[0] < cx - radius < cx + radius < hi[0]
    assert lo[1] < cy - radius < cy + radius < hi[1]
    assert z0 < lo[2] and z1 > hi[2]
    assert model["schema"] == "box_cylinder_boolean_v1"
    assert len(model["blocks"]) == 4
    interfaces = model["topology"]["interfaces"]
    assert len(interfaces) == 4
    # The exact positive-radius cylinder crosses both bounded stock end planes
    # and stays strictly inside the lateral stock walls. This is a through
    # tunnel, so the material cannot be an injective image of a tensor cube.
    return {
        "source_schema": model["schema"],
        "source_geometry": "valid_by_native_audit",
        "conversion_status": "unsupported_profile",
        "block_count": len(model["blocks"]),
        "retained_c0_joins": len(interfaces),
        "exact_common_parent_obstruction": "through_tunnel_topology",
        "proof_inputs": {
            "center_mm": [str(cx), str(cy)],
            "radius_mm": str(radius),
            "tool_interval_mm": [str(z0), str(z1)],
            "stock_z_mm": [str(lo[2]), str(hi[2])],
        },
    }


def verify_nurbs_cross(recipe: dict, model: dict) -> dict:
    assert recipe["family"] == "nurbs_cross_bores_v1"
    profiles = recipe["profiles"]
    assert set(profiles) == {"x_min", "x_max", "z_min", "z_max"}
    curves = list(profiles.values())
    assert all(curve["degree"] == 3 for curve in curves)
    assert all(curve["knots"] == curves[0]["knots"] for curve in curves)
    assert model["schema"] == "nurbs_cross_bores_v1"
    assert len(model["blocks"]) == 26
    interfaces = model["topology"]["interfaces"]
    assert len(interfaces) == 36
    # The source recipe supplies four 2D generating curves for two crossing
    # bores; it does not supply one complete 3D material-parent map or a
    # per-block exact restriction operator. The solid remains valid, while
    # conversion to the single-parent profile remains unsupported.
    return {
        "source_schema": model["schema"],
        "source_geometry": "valid_by_native_audit",
        "conversion_status": "unsupported_profile",
        "block_count": len(model["blocks"]),
        "retained_c0_joins": len(interfaces),
        "known_generating_curve_count": len(curves),
        "common_parent_status": "not_proven_from_retained_recipe",
        "source_knot_vector": curves[0]["knots"],
        "source_domain": [curves[0]["knots"][3], curves[0]["knots"][len(curves[0]["controls"])]],
    }


def exact_translation_columns(model: dict, axis: int = 2) -> dict[str, tuple[str, tuple[F, F, F] | None]]:
    result = {}
    for block in model["blocks"]:
        degree = block["degrees"]
        shape = [p + 1 for p in degree]
        if degree[axis] != 1 or shape[0] * shape[1] * shape[2] != len(block["controls"]):
            result[block["id"]] = ("unsupported_column_axis", None)
            continue
        transverse_shape = shape.copy()
        transverse_shape[axis] = 1
        translation = None
        valid = True
        for at in product(*(range(n) for n in transverse_shape)):
            low_index = list(at)
            high_index = list(at)
            high_index[axis] = 1
            flatten = lambda p: (p[0] * shape[1] + p[1]) * shape[2] + p[2]
            low = block["controls"][flatten(low_index)]
            high = block["controls"][flatten(high_index)]
            if low[3] <= 0 or high[3] != low[3]:
                valid = False
                break
            current = tuple(F(high[c] - low[c], low[3]) for c in range(3))
            if translation is not None and translation != current:
                valid = False
                break
            translation = current
        if valid:
            result[block["id"]] = ("certified_exact_prismatic_column", translation)
        else:
            result[block["id"]] = ("not_a_translation_column", None)
    return result


def check_native_report(cli: str, model_path: Path, expected: dict) -> None:
    report = json.loads(run_cli(cli, "spline-eligibility", str(model_path)))
    assert report["source_validity"] == "valid"
    assert report["conversion_status"] == "unsupported_profile"
    assert len(report["retained_joins"]) == expected["retained_c0_joins"]
    if "known_generating_curve_count" in expected:
        assert len(report["known_generating_sources"]) == expected["known_generating_curve_count"]
    model = json.loads(model_path.read_text())
    expected_columns = exact_translation_columns(model)
    native_columns = report["regional_source_columns"]["blocks"]
    assert len(native_columns) == len(expected_columns)
    for item in native_columns:
        expected_status, vector = expected_columns[item["block_id"]]
        assert item["status"] == expected_status
        if vector is not None:
            actual = tuple(
                F(int(entry["numerator"]), int(entry["denominator"]))
                for entry in item["translation_vector_mm"]
            )
            assert actual == vector


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", required=True)
    parser.add_argument("--check-native-report", action="store_true")
    args = parser.parse_args()
    results = []
    with tempfile.TemporaryDirectory(prefix="spline-eligibility-") as temp_name:
        temp = Path(temp_name)
        for recipe_name, verify in (
            ("box_minus_cylinder.recipe.json", verify_box_cylinder),
            ("nurbs_cross_bores.recipe.json", verify_nurbs_cross),
        ):
            recipe_path = ROOT / "tests/references" / recipe_name
            model_path = temp / recipe_name.replace(".recipe.json", ".model.json")
            recipe, model = build(args.cli, recipe_path, model_path)
            result = verify(recipe, model)
            if args.check_native_report:
                check_native_report(args.cli, model_path, result)
            columns = exact_translation_columns(model)
            result["exact_prismatic_columns"] = sum(status == "certified_exact_prismatic_column" for status, _ in columns.values())
            result["nonprismatic_or_unavailable_blocks"] = [
                block_id for block_id, (status, _) in columns.items()
                if status != "certified_exact_prismatic_column"
            ]
            results.append(result)
    print(json.dumps({"status": "independent_eligibility_checks_passed", "models": results}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
