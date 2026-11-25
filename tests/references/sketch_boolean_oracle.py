#!/usr/bin/env python3
"""Independent exact oracle for orthogonal sketch Boolean histories (PA8).

The oracle never inspects the native partition algorithm. It decomposes the
input into the exact Cartesian grid induced by polygon vertices and axial
events, classifies each open cell by ordered feature history, then measures
cell volumes and flood-fills material/void adjacency.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from collections import deque
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "target/release/solidcraft-trivariate"
OUTER = {"id": "body", "vertices": [[0, 0], [20, 0], [20, 20], [0, 20]]}
HEIGHT = 10


def operation(name: str, rect: tuple[int, int, int, int], z: tuple[int, int], op: str) -> dict:
    x0, y0, x1, y1 = rect
    return {"id": name, "vertices": [[x0, y0], [x1, y0], [x1, y1], [x0, y1]],
            "z_interval_mm": list(z), "operation": op}


def recipe(ops: list[dict]) -> dict:
    return {"family": "sketch_boolean_layers_v1", "outer": OUTER,
            "operations": ops, "construction": [], "height_mm": HEIGHT}


CASES = {
    "blind_pocket": recipe([operation("pocket", (6, 6, 14, 14), (6, 10), "remove")]),
    "attached_island": recipe([
        operation("pocket", (5, 5, 15, 15), (4, 10), "remove"),
        operation("island", (8, 8, 12, 12), (4, 7), "add"),
    ]),
    "through_hole": recipe([operation("hole", (8, 8, 12, 12), (0, 10), "remove")]),
    "independent_ab": recipe([
        operation("left", (3, 3, 7, 7), (2, 8), "remove"),
        operation("right", (13, 13, 17, 17), (2, 8), "remove"),
    ]),
    "independent_ba": recipe([
        operation("right", (13, 13, 17, 17), (2, 8), "remove"),
        operation("left", (3, 3, 7, 7), (2, 8), "remove"),
    ]),
    "dependent_remove_then_add": recipe([
        operation("pocket", (5, 5, 15, 15), (4, 10), "remove"),
        operation("island", (8, 8, 12, 12), (4, 7), "add"),
    ]),
    "dependent_add_then_remove": recipe([
        operation("island", (8, 8, 12, 12), (4, 7), "add"),
        operation("pocket", (5, 5, 15, 15), (4, 10), "remove"),
    ]),
    "touch_outer_refused": recipe([operation("touch", (18, 6, 20, 14), (3, 8), "remove")]),
    "floating_island_refused": recipe([
        operation("through_void", (4, 4, 16, 16), (0, 10), "remove"),
        operation("floating", (8, 8, 12, 12), (4, 6), "add"),
    ]),
}

EXPECTED_TAGS = {
    "blind_pocket": {"end", "start", "outer/body", "feature/pocket/lower", "feature/pocket/wall"},
    "attached_island": {"end", "start", "outer/body", "feature/pocket/lower", "feature/pocket/wall",
                        "feature/island/upper", "feature/island/wall"},
    "through_hole": {"end", "start", "outer/body", "feature/hole/wall"},
    "independent_ab": {"end", "start", "outer/body", "feature/left/lower", "feature/left/upper",
                       "feature/left/wall", "feature/right/lower", "feature/right/upper", "feature/right/wall"},
    "independent_ba": {"end", "start", "outer/body", "feature/left/lower", "feature/left/upper",
                       "feature/left/wall", "feature/right/lower", "feature/right/upper", "feature/right/wall"},
    "dependent_remove_then_add": {"end", "start", "outer/body", "feature/pocket/lower",
                                   "feature/pocket/wall", "feature/island/upper", "feature/island/wall"},
    "dependent_add_then_remove": {"end", "start", "outer/body", "feature/pocket/lower", "feature/pocket/wall"},
}


def inside_rect_or_polygon(poly: list[list[int]], point: tuple[Fraction, Fraction]) -> bool:
    """Exact even/odd test for a point not on the polygon boundary."""
    x, y = point
    result = False
    for a, b in zip(poly, poly[1:] + poly[:1]):
        ax, ay = map(Fraction, a)
        bx, by = map(Fraction, b)
        if (ay > y) != (by > y):
            xcross = ax + (y - ay) * (bx - ax) / (by - ay)
            if x < xcross:
                result = not result
    return result


def classify(ops: list[dict], point: tuple[Fraction, Fraction], z: Fraction) -> tuple[bool, str | None]:
    state, owner = True, None
    for item in ops:
        lo, hi = item["z_interval_mm"]
        if Fraction(lo) < z < Fraction(hi) and inside_rect_or_polygon(item["vertices"], point):
            state = item["operation"] == "add"
            owner = item["id"]
    return state, owner


def exact_oracle(spec: dict) -> dict:
    outer = spec["outer"]["vertices"]
    ops = spec["operations"]
    xs = sorted({p[0] for poly in [outer] + [o["vertices"] for o in ops] for p in poly})
    ys = sorted({p[1] for poly in [outer] + [o["vertices"] for o in ops] for p in poly})
    zs = sorted({0, spec["height_mm"], *(z for o in ops for z in o["z_interval_mm"])})
    nx, ny, nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
    solid = [[[False for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]
    owner = [[[None for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]
    layer_volume = [0] * nz
    for i, j, k in product(range(nx), range(ny), range(nz)):
        p = (Fraction(xs[i] + xs[i + 1], 2), Fraction(ys[j] + ys[j + 1], 2))
        z = Fraction(zs[k] + zs[k + 1], 2)
        if inside_rect_or_polygon(outer, p):
            solid[i][j][k], owner[i][j][k] = classify(ops, p, z)
        if solid[i][j][k]:
            layer_volume[k] += (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j]) * (zs[k + 1] - zs[k])

    # Extend the grid by one exterior cell on each side so cavity status is
    # computed by adjacency, not inferred from feature names or volume.
    ex, ey, ez = [xs[0] - 1] + xs + [xs[-1] + 1], [ys[0] - 1] + ys + [ys[-1] + 1], [-1] + zs + [zs[-1] + 1]
    shape = (len(ex) - 1, len(ey) - 1, len(ez) - 1)
    def is_solid(i: int, j: int, k: int) -> bool:
        return (1 <= i <= nx and 1 <= j <= ny and 1 <= k <= nz
                and solid[i - 1][j - 1][k - 1])

    visited: set[tuple[int, int, int]] = set()
    void_components = []
    material_components = 0
    for start in product(*(range(n) for n in shape)):
        if start in visited:
            continue
        state = is_solid(*start)
        queue = deque([start])
        visited.add(start)
        touches_exterior = False
        while queue:
            cell = queue.popleft()
            i, j, k = cell
            touches_exterior |= (i in (0, shape[0] - 1) or j in (0, shape[1] - 1)
                                 or k in (0, shape[2] - 1))
            for axis in range(3):
                for delta in (-1, 1):
                    neighbor = list(cell)
                    neighbor[axis] += delta
                    if not 0 <= neighbor[axis] < shape[axis]:
                        continue
                    t = tuple(neighbor)
                    if t not in visited and is_solid(*t) == state:
                        visited.add(t)
                        queue.append(t)
        if state:
            material_components += 1
        elif not touches_exterior:
            void_components.append(len(void_components) + 1)

    volume = sum(layer_volume)
    return {"volume_mm3": volume, "layer_volume_mm3": layer_volume,
            "material_components": material_components,
            "closed_void_components": len(void_components),
            "grid_cells": [nx, ny, nz]}


def run_native(name: str, spec: dict, work: Path) -> dict:
    recipe_path, model_path = work / f"{name}.recipe.json", work / f"{name}.model.json"
    recipe_path.write_text(json.dumps(spec, separators=(",", ":")) + "\n")
    proc = subprocess.run([str(BIN), "build", str(recipe_path), str(model_path)],
                          cwd=ROOT, text=True, capture_output=True)
    report = None
    if proc.stdout.strip():
        try:
            report = json.loads(proc.stdout)
        except json.JSONDecodeError:
            pass
    if proc.returncode:
        return {"accepted": False, "stderr": proc.stderr.strip(), "stdout": proc.stdout.strip()}
    model = json.loads(model_path.read_text())
    tags: dict[str, int] = {}
    for face in model.get("design_faces", []):
        tag = face.get("feature", "<missing>")
        tags[tag] = tags.get(tag, 0) + 1
    return {"accepted": True, "audit": report, "blocks": len(model.get("blocks", [])),
            "design_face_tags": tags,
            "layer_block_counts": model.get("partition_plan", {}).get("layers", []),
            "model": str(model_path.relative_to(ROOT))}


def run_occ(specs: dict[str, dict], work: Path) -> dict[str, dict]:
    """Measure the same CSG histories with an independent FreeCAD/OCCT process."""
    executable = shutil.which("freecadcmd")
    if executable is None:
        return {}
    result_path = work / "occ-volumes.json"
    macro_path = work / "occ-volume-oracle.FCMacro"
    encoded = json.dumps(specs, separators=(",", ":"))
    macro_path.write_text(f'''import json, Part, FreeCAD as App
specs = json.loads(r"""{encoded}""")
out = {{}}
def prism(loop, z0, z1):
    points = [App.Vector(p[0], p[1], 0) for p in loop]
    wire = Part.makePolygon(points + [points[0]])
    shape = Part.Face(wire).extrude(App.Vector(0, 0, z1-z0))
    shape.translate(App.Vector(0, 0, z0))
    return shape
for name, spec in specs.items():
    shape = prism(spec["outer"]["vertices"], 0, spec["height_mm"])
    for item in spec["operations"]:
        tool = prism(item["vertices"], *item["z_interval_mm"])
        shape = shape.fuse(tool) if item["operation"] == "add" else shape.cut(tool)
    out[name] = {{"volume_mm3": shape.Volume, "valid": shape.isValid(), "solids": len(shape.Solids)}}
with open(r"{result_path}", "w") as stream:
    json.dump(out, stream, sort_keys=True)
''')
    proc = subprocess.run([executable, str(macro_path)], cwd=ROOT, text=True,
                          capture_output=True, timeout=180)
    if proc.returncode or not result_path.exists():
        raise RuntimeError(f"FreeCADCmd/OCCT oracle failed: {proc.stdout[-1000:]} {proc.stderr[-1000:]}")
    return json.loads(result_path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "docs/results/sketch-boolean-20261125.json")
    args = parser.parse_args()
    if not BIN.exists():
        raise SystemExit(f"missing release executable: {BIN}; build it first")
    work = ROOT / "build/sketch-boolean-pa8-oracle"
    work.mkdir(parents=True, exist_ok=True)
    results = {}
    occ = run_occ({name: spec for name, spec in CASES.items() if not name.endswith("_refused")}, work)
    for name, spec in CASES.items():
        expected = exact_oracle(spec)
        native = run_native(name, spec, work)
        accepted = native["accepted"]
        if name.endswith("_refused"):
            if accepted:
                raise AssertionError(f"{name}: native unexpectedly accepted invalid case")
            results[name] = {"oracle": "refuse", "native": native}
            continue
        if not accepted:
            raise AssertionError(f"{name}: native refused mandatory case: {native}")
        observed_volume = native["audit"]["volume_exact_mm3"]
        if Fraction(*observed_volume) != expected["volume_mm3"]:
            raise AssertionError(f"{name}: exact volume mismatch {observed_volume} != {expected['volume_mm3']}")
        if native["audit"]["cavities"] != expected["closed_void_components"]:
            raise AssertionError(f"{name}: cavity count mismatch")
        if native["audit"]["determinant_lower_bound_mm3"] <= 0:
            raise AssertionError(f"{name}: nonpositive certified determinant")
        if set(native["design_face_tags"]) != EXPECTED_TAGS[name]:
            raise AssertionError(f"{name}: operation surface attribution mismatch: {native['design_face_tags']}")
        if occ:
            external = occ[name]
            if not external["valid"] or external["solids"] != 1:
                raise AssertionError(f"{name}: OCCT produced an invalid or disconnected solid: {external}")
            if abs(external["volume_mm3"] - float(expected["volume_mm3"])) > 1e-8:
                raise AssertionError(f"{name}: OCCT volume mismatch: {external}")
        native["occ"] = occ.get(name)
        if expected["material_components"] != 1:
            if name != "floating_island_refused":
                raise AssertionError(f"{name}: oracle expects disconnected material")
        results[name] = {"oracle": expected, "native": native}

    ab, ba = results["independent_ab"]["oracle"], results["independent_ba"]["oracle"]
    assert ab == ba, "independent operation permutations changed exact material"
    first = results["dependent_remove_then_add"]["oracle"]["volume_mm3"]
    second = results["dependent_add_then_remove"]["oracle"]["volume_mm3"]
    assert first != second, "dependent operation order did not affect material"

    sentinel = work / "atomic-target.json"
    refusal_targets = {}
    for name in ("touch_outer_refused", "floating_island_refused"):
        sentinel.write_text("preserve-this-model\n")
        recipe_path = work / f"{name}.atomic.recipe.json"
        recipe_path.write_text(json.dumps(CASES[name]))
        proc = subprocess.run([str(BIN), "build", str(recipe_path), str(sentinel)],
                              cwd=ROOT, text=True, capture_output=True)
        unchanged = sentinel.read_text() == "preserve-this-model\n"
        if proc.returncode == 0 or not unchanged:
            raise AssertionError(f"refusal {name} did not preserve the previous target")
        refusal_targets[name] = unchanged
    report = {
        "schema": "sketch_boolean_oracle_report_v1",
        "protocol": "Exact rectilinear cell decomposition from all polygon vertices and z events; ordered last-active operation classification; 3D face-adjacency flood fill.",
        "environment": {"binary": str(BIN.relative_to(ROOT)), "python": "standard library", "exact_arithmetic": "fractions.Fraction and integer cell measures"},
        "cases": results,
        "comparisons": {
            "independent_orders_equal": ab == ba,
            "dependent_order_volumes_mm3": [int(first), int(second)],
            "refusal_preserved_target": refusal_targets,
        },
        "external_oracle": "FreeCADCmd / OpenCascade sequential fuse/cut of exact polygonal prisms" if occ else None,
        "limits": ["Corpus is intentionally orthogonal; arbitrary slanted boundaries are outside this oracle's exact Cartesian decomposition.",
                   "OCCT validates CSG validity and volume; this report does not export/reopen STEP."],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed", "cases": len(results), "report": str(args.output.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
