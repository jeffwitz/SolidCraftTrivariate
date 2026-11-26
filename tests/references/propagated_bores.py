"""Paired public bracket controls. Run under a Python environment with FreeCAD and Part."""
import FreeCAD as App
import Part

import argparse
import hashlib
import json
import math
import random
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from fieldcraft_compact_model_oracle import exact_interface_audit


def cases():
    north = ("x", 4, 32, 20, 4)
    east = ("y", 40, 4, 20, 4)
    low = ("y", 38, 4, 10, 3)
    high = ("y", 38, 4, 30, 3)
    north_low = ("x", 4, 32, 9, 3)
    north_high = ("x", 4, 32, 31, 3)
    layouts = {
        "original": [east, north],
        "in_plane": [east, north, ("y", 54, 4, 20, 2)],
        "east_stack": [low, high, north],
        "north_stack": [north_low, north_high, east],
        "east_offset_stack": [("y", 36, 4, 9, 3), ("y", 40, 4, 25, 4), north],
        "diagonal": [("y", 28, 4, 10, 3), ("y", 48, 4, 30, 3), north],
        "three_levels": [("y", 38, 4, z, 2) for z in (7, 20, 33)] + [north],
        "both_stacked": [low, high, north_low, north_high],
        "axial_z": [("z", 25, 4, 20, 2), east, north],
        "unequal_stack": [low, high, north],
    }
    for name, tools in layouts.items():
        tx, ty = (12, 8) if name == "unequal_stack" else (8, 8)
        recipe = {"family": "orthogonal_sketch_features_v1", "height_mm": 40,
                  "outer": {"id": "profile", "vertices": [[0, 0], [60, 0], [60, ty],
                            [tx, ty], [tx, 50], [0, 50]]},
                  "fillets": [{"id": "bend", "vertex": [tx, ty], "radius_mm": 6}],
                  "bores": [{"id": f"b{i}", "axis": axis, "center_mm": [x, y, z], "radius_mm": r}
                            for i, (axis, x, y, z, r) in enumerate(tools)]}
        yield name, recipe


def expected(recipe):
    tx, ty = recipe["outer"]["vertices"][3]
    radius = recipe["fillets"][0]["radius_mm"]
    volume = (60 * ty + 50 * tx - tx * ty + radius**2 * (1 - math.pi / 4)) * 40
    tools = []
    for bore in recipe["bores"]:
        axis = "xyz".index(bore["axis"])
        tools.append((axis, bore["center_mm"], bore["radius_mm"]))
        volume -= math.pi * bore["radius_mm"]**2 * (tx, ty, 40)[axis]

    def material(p):
        x, y, z = p
        if not (0 < x < 60 and 0 < y < 50 and 0 < z < 40):
            return False
        if x > tx and y > ty:
            if x >= tx + radius or y >= ty + radius:
                return False
            if math.hypot(x-tx-radius, y-ty-radius) < radius:
                return False
        for axis, center, r in tools:
            if sum((p[i]-center[i])**2 for i in range(3) if i != axis) < r*r:
                return False
        return True

    rng = random.Random(170819)
    points = [(rng.uniform(-1, 61), rng.uniform(-1, 51), rng.uniform(-1, 41)) for _ in range(160)]
    for axis, center, r in tools:
        a, b = (axis+1) % 3, (axis+2) % 3
        for angle in range(8):
            for distance in (r-.02, r+.02):
                p = center.copy()
                p[a] += distance * math.cos(angle * math.pi/4)
                p[b] += distance * math.sin(angle * math.pi/4)
                points.append(p)
    return volume, [(p, material(p)) for p in points]


def measure(cli, helper, recipe, path):
    path.mkdir(parents=True, exist_ok=True)
    source, model, step, brep, vtu = (path / f"case.{s}" for s in ("recipe.json", "json", "step", "brep", "vtu"))
    source.write_text(json.dumps(recipe))
    def run(*args):
        return subprocess.run(args, capture_output=True, text=True, timeout=120, check=True).stdout

    try:
        run(cli, "build", source, model)
    except subprocess.CalledProcessError as exc:
        return {"status": "refused", "reason": (exc.stderr or exc.stdout).strip()}
    audit = json.loads(run(cli, "audit", model))
    data = json.loads(model.read_text())
    plan = data["partition_plan"]
    out = {"status": "qualified", "method": plan["method"], "blocks": len(data["blocks"]),
           "reason": plan.get("fallback_reason"), "det_lower_mm3": audit["determinant_lower_bound_mm3"]}
    trace = exact_interface_audit(data)
    labels = [(f["block"], f["side"]) for f in data["design_faces"]]
    exterior = [tuple(pair) for pair in data["topology"]["exterior"]]
    label_ok = sorted(labels) == sorted(exterior) and len(set(labels)) == len(labels)
    bore_tags = {f["feature"] for f in data["design_faces"] if f["feature"].startswith("bore/")}
    label_ok &= bore_tags == {f"bore/{b['id']}/wall" for b in recipe["bores"]}
    if "cut_comparison" in plan:
        out["cut_comparison"] = plan["cut_comparison"]
    if "cut_rounds" in plan:
        out["cut_rounds"] = plan["cut_rounds"]
        out["cut_source_records"] = len(plan["cut_sources"])
    out["interfaces"] = trace["exact_trace_interfaces"]
    run(cli, "export-geometry-vtu", model, "4", vtu)
    arrays = ET.parse(vtu).findall('.//Piece/PointData/DataArray[@Name="condition_number"]')
    samples = [list(map(float, array.text.split())) for array in arrays]
    if len(samples) != out["blocks"] or any(len(v) != 125 for v in samples):
        raise ValueError("incomplete block quality sample")
    values = [value for block in samples for value in block]
    out["max_sampled_kappa"] = max(values)
    run(cli, "export-step", model, step)
    shape = Part.read(str(step))
    shape.exportBrep(str(brep))
    mass = json.loads(run(helper, brep))[-1]["volume_mm3"]
    volume, probes = expected(recipe)
    out.update(volume_relative_error=abs(mass-volume)/volume, probes=len(probes),
               probe_mismatches=sum(bool(shape.isInside(App.Vector(*p), 1e-7, True)) != inside
                                    for p, inside in probes),
               closed_valid=bool(shape.isValid() and len(shape.Solids)==1 and len(shape.Shells)==1
                                 and shape.Solids[0].isClosed()), labels_match=label_ok)
    if not (out["closed_valid"] and label_ok and out["probe_mismatches"] == 0
            and out["volume_relative_error"] <= 2e-9 and out["det_lower_mm3"] > 0
            and len(values) == 125 * out["blocks"] and all(map(math.isfinite,values))):
        out["status"] = "invalid"
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("occ_volume", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    rows = {}
    for name, recipe in cases():
        try:
            rows[name] = measure(args.cli.resolve(), args.occ_volume.resolve(), recipe, work / name)
        except subprocess.TimeoutExpired as exc:
            rows[name] = {"status": "timeout", "stage": str(exc.cmd[1])}
        except (RuntimeError, ValueError, OSError, AssertionError, subprocess.CalledProcessError) as exc:
            rows[name] = {"status": "invalid", "reason": str(exc)}
        print(name, rows[name]["status"], rows[name].get("blocks"), flush=True)
    (work / "scores.json").write_text(json.dumps({
        "binary_sha256": hashlib.sha256(args.cli.read_bytes()).hexdigest(), "cases": rows},
        indent=2) + "\n")
    if any(r["status"] == "invalid" for r in rows.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
