#!/usr/bin/env python3
"""Independent IGA checks for audited multi-loop NURBS models.

The production CLI supplies the audited geometry and its C0 Q3 package. The
space oracle checks exact degree elevation, oriented shared traces, and affine
field values and gradients. The equilibrium oracle independently integrates
constant-stress volume and exterior-traction terms at q4/q8/q16.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_reference(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tests/references" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load reference module {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


space_oracle = load_reference("multiloop_space_oracle", "nurbs_space_oracle.py")
equilibrium_oracle = load_reference("multiloop_equilibrium_oracle", "nurbs_patch_equilibrium.py")


def run_cli(cli: Path, command: list[str]) -> dict:
    process = subprocess.run([str(cli), *command], cwd=ROOT, capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError(f"CLI {' '.join(command)} failed: {process.stderr or process.stdout}")
    try:
        return json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"invalid CLI JSON for {' '.join(command)}: {process.stdout}") from error


def run_case(cli: Path, name: str, recipe_path: Path, work: Path) -> dict:
    model_path = work / f"{name}.model.json"
    package_path = work / f"{name}.space.json"
    build_audit = run_cli(cli, ["build", str(recipe_path), str(model_path)])
    model = json.loads(model_path.read_text())
    package_report = run_cli(cli, ["export-nurbs-space", str(model_path), str(package_path)])
    package = json.loads(package_path.read_text())

    space = space_oracle.package_checks(model, package)
    equilibrium = []
    for order in (4, 8, 16):
        try:
            equilibrium.append(equilibrium_oracle.integrate(package, order))
        except ArithmeticError as error:
            ndof = package["scalar_ndof"]
            bad = next(
                (block["id"] for block in package["blocks"]
                 if _fails_volume_block(block, order, ndof)),
                "unknown block",
            )
            raise ArithmeticError(f"{name} q{order}: {bad}: {error}") from error
    threshold = 1e-8
    if equilibrium[-1]["relative_residual_linf"] > threshold:
        raise AssertionError(
            f"{name}: q16 constant-stress equilibrium residual "
            f"{equilibrium[-1]['relative_residual_linf']:.6g} exceeds {threshold:g}"
        )
    if space["target_tensor_degree"] != 3:
        raise AssertionError(f"{name}: expected Q3 package, got {space['target_tensor_degree']}")

    return {
        "recipe": str(recipe_path.relative_to(ROOT)),
        "model_schema": model.get("schema"),
        "build_audit": build_audit,
        "design_face_feature_counts": dict(sorted(Counter(
            face["feature"] for face in model["design_faces"]
        ).items())),
        "space_export": package_report,
        "space_checks": space,
        "constant_stress_patch": {
            "stress_mpa": equilibrium_oracle.SIGMA.tolist(),
            "traction_convention": "sigma*n on every distinct design face; oriented outward by block side",
            "orders": equilibrium,
            "q16_relative_linf_threshold": threshold,
        },
        "artifacts": {
            "model": str(model_path.relative_to(ROOT)),
            "q3_package": str(package_path.relative_to(ROOT)),
        },
    }


def _fails_volume_block(block: dict, order: int, ndof: int) -> bool:
    try:
        equilibrium_oracle.volume_contribution(block, order, ndof)
    except ArithmeticError:
        return True
    return False


def check_homogeneous_kernel() -> dict:
    controls = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                controls.append([2 * i / 3, 2 * j / 3, 2 * k / 3, 2.0])
    block = {"controls": controls, "dofs": list(range(64))}
    q = (0.17, 0.43, 0.82)
    position, jacobian, det, *_ = equilibrium_oracle.local_data(block, q)
    if max(abs(a - b) for a, b in zip(position, q)) > 2e-14:
        raise AssertionError(f"homogeneous kernel affine position mismatch: {position} != {q}")
    identity_error = max(abs(jacobian[i, j] - (1.0 if i == j else 0.0))
                         for i in range(3) for j in range(3))
    if identity_error > 2e-14 or abs(det - 1.0) > 2e-14:
        raise AssertionError(f"homogeneous kernel affine Jacobian mismatch: det={det}, error={identity_error}")

    scaled = {**block, "controls": [[7 * x for x in p] for p in controls]}
    p2, j2, d2, *_ = equilibrium_oracle.local_data(scaled, q)
    gauge_error = max(
        max(abs(a - b) for a, b in zip(position, p2)),
        max(abs(jacobian[i, j] - j2[i, j]) for i in range(3) for j in range(3)),
        abs(det - d2),
    )
    if gauge_error > 2e-14:
        raise AssertionError(f"homogeneous gauge changed physical geometry: {gauge_error}")
    return {"position_error": max(abs(a - b) for a, b in zip(position, q)),
            "jacobian_linf_error": identity_error, "determinant": det,
            "uniform_homogeneous_gauge_factor": 7, "gauge_invariance_linf_error": gauge_error}


def prospective_recipes(source: dict) -> dict[str, dict]:
    loops = {loop["id"]: loop for loop in source["loops"]}
    first_id = next(iter(loops))
    first = dict(loops[first_id])
    first["z_interval_mm"] = [0, source["height_mm"]]
    through = {**source, "loops": [first]}
    pocket_loop = dict(loops[first_id])
    pocket_loop["z_interval_mm"] = [0, source["height_mm"] * 2 // 5]
    pocket = {**source, "loops": [pocket_loop]}
    two = {**source, "loops": [dict(loop) for loop in source["loops"]]}
    renamed = []
    for index, loop in enumerate(reversed(source["loops"])):
        copy = dict(loop)
        copy["id"] = f"loop-{index + 1}"
        renamed.append(copy)
    permuted = {**source, "loops": renamed}
    return {
        "single-through": through,
        "single-blind-pocket": pocket,
        "independent-two-loop-mixed-depth": two,
        "permuted-renamed-two-loop": permuted,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", required=True, type=Path)
    parser.add_argument("--recipe", type=Path,
                        default=ROOT / "tests/references/nurbs_multiloop.recipe.json",
                        help="canonical multi-loop fixture used to derive prospective variants")
    parser.add_argument("--work", type=Path, default=ROOT / "build/nurbs/multiloop-iga")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs/results/nurbs-multiloop-iga.json")
    args = parser.parse_args()
    start = time.perf_counter()
    kernel_check = check_homogeneous_kernel()
    cli, work = args.cli.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    canonical = json.loads(args.recipe.resolve().read_text())
    cases = {}
    for name, recipe in prospective_recipes(canonical).items():
        recipe_path = work / f"{name}.recipe.json"
        recipe_path.write_text(json.dumps(recipe, separators=(",", ":")) + "\n")
        cases[name] = run_case(cli, name, recipe_path, work)

    report = {
        "schema": "nurbs_multiloop_iga_report_v1",
        "status": "verified",
        "protocol": (
            "Independent reuse of Fraction tensor-Bernstein degree-elevation, D4-oriented "
            "trace identity, and affine value/gradient checks from nurbs_space_oracle; "
            "independent NumPy rational-Q3 constant-stress volume and exterior-traction "
            "integration from nurbs_patch_equilibrium."
        ),
        "environment": {
            "cli": str(cli),
            "cli_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
            "quadrature_orders": [4, 8, 16],
        },
        "homogeneous_kernel_check": kernel_check,
        "cases": cases,
        "elapsed_seconds": time.perf_counter() - start,
        "limits": [
            "The checks apply to the listed audited recipes and their exported Q3 packages.",
            "Quadrature equilibrium is a floating-point diagnostic, not an exact integral proof.",
            "A passing patch test does not prove global geometry injectivity or approximation quality.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "cases": list(cases), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
