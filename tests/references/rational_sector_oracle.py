#!/usr/bin/env python3
"""Independent analytic and Open CASCADE oracle for rational annular sectors.

Run inside the FreeCAD Python runtime so the STEP round trip uses the installed
OCCT kernel. The reference solid is made by revolving an r-z rectangle.
"""
import FreeCAD as App
import Part

import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coaxial_program import compile_mass_oracle


ROOT = Path(__file__).resolve().parents[2]


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def recipe(n: int, d: int, ri: int, ro: int, height: int) -> dict:
    return {"family": "rational_sector_prism_v1",
            "half_angle_tangent": [n, d], "radii_mm": [ri, ro],
            "height_mm": height}


def split_expected(n: int, d: int) -> bool:
    # 1 + sqrt(2) is irrational; compare without evaluating it in floating point.
    return n > d and (n - d) ** 2 > 2 * d * d


def exact_volume(rec: dict) -> dict:
    n, d = rec["half_angle_tangent"]
    ri, ro = rec["radii_mm"]
    return {"atan_coefficient_mm3": rec["height_mm"] * (ro * ro - ri * ri),
            "argument": [n, d]}


def volume_and_center(rec: dict) -> tuple[float, list[float]]:
    n, d = rec["half_angle_tangent"]
    ri, ro = rec["radii_mm"]
    h = rec["height_mm"]
    t = n / d
    theta = 2.0 * math.atan(t)
    volume = h * (ro * ro - ri * ri) * math.atan(t)
    radial_centroid = (2.0 / 3.0) * (ro * ro + ro * ri + ri * ri) / (ro + ri)
    cx = radial_centroid * math.sin(theta) / theta
    cy = radial_centroid * (1.0 - math.cos(theta)) / theta
    return volume, [cx, cy, h / 2.0]


def reference_solid(rec: dict):
    ri, ro = rec["radii_mm"]
    h = rec["height_mm"]
    n, d = rec["half_angle_tangent"]
    theta = math.degrees(2.0 * math.atan(n / d))
    points = [App.Vector(ri, 0, 0), App.Vector(ro, 0, 0),
              App.Vector(ro, 0, h), App.Vector(ri, 0, h), App.Vector(ri, 0, 0)]
    wire = Part.makePolygon(points)
    face = Part.Face(wire)
    return face.revolve(App.Vector(0, 0, 0), App.Vector(0, 0, 1), theta)


def inside_ideal(rec: dict, point: tuple[float, float, float], eps=1e-10) -> bool:
    x, y, z = point
    n, d = rec["half_angle_tangent"]
    ri, ro = rec["radii_mm"]
    h = rec["height_mm"]
    radius = math.hypot(x, y)
    angle = math.atan2(y, x)
    if angle < -eps:
        angle += 2 * math.pi
    theta = 2 * math.atan(n / d)
    return ri + eps < radius < ro - eps and eps < angle < theta - eps and eps < z < h - eps


def classify_points(rec: dict, shape, seed: int) -> dict:
    n, d = rec["half_angle_tangent"]
    ri, ro = rec["radii_mm"]
    h = rec["height_mm"]
    theta = 2 * math.atan(n / d)
    rng = random.Random(seed)
    points = []
    # Sixteen strict interior points per STEP case.
    for af in (0.07, 0.31, 0.64, 0.93):
        for rf in (0.2, 0.8):
            for zf in (0.2, 0.8):
                radius = ri + rf * (ro - ri)
                angle = af * theta
                points.append((radius * math.cos(angle), radius * math.sin(angle), zf * h))
    # Paired samples on both sides of each analytic radial and angular boundary.
    dr = 1e-6 * min(float(ri), float(ro - ri))
    da = 1e-6 / max(1.0, float(ro))
    z_mid = 0.5 * h
    for af in (0.2, 0.8):
        angle = af * theta
        for boundary, offsets in ((ri, (-dr, dr)), (ro, (-dr, dr))):
            for offset in offsets:
                radius = boundary + offset
                points.append((radius * math.cos(angle), radius * math.sin(angle), z_mid))
    for rf in (0.3, 0.7):
        radius = ri + rf * (ro - ri)
        for boundary, offsets in ((0.0, (-da, da)), (theta, (-da, da))):
            for offset in offsets:
                angle = boundary + offset
                points.append((radius * math.cos(angle), radius * math.sin(angle), z_mid))
    # Deterministic points outside each of the four physical boundaries.
    for _ in range(8):
        angle = rng.uniform(0, theta)
        radius = rng.uniform(ri + 0.1 * (ro - ri), ro - 0.1 * (ro - ri))
        z = rng.uniform(0.15 * h, 0.85 * h)
        points.extend([
            (ri - 0.1 * max(1, ri), 0.0, z),
            (ro + 0.1 * max(1, ro), 0.0, z),
            (radius * math.cos(-0.05 * theta), radius * math.sin(-0.05 * theta), z),
            (radius * math.cos(1.05 * theta), radius * math.sin(1.05 * theta), z),
            (radius * math.cos(angle), radius * math.sin(angle), -0.1 * max(1, h)),
            (radius * math.cos(angle), radius * math.sin(angle), h + 0.1 * max(1, h)),
        ])
    mismatches = []
    records = []
    for point in points:
        expected = inside_ideal(rec, point)
        actual = bool(shape.isInside(App.Vector(*point), 1e-8, True))
        records.append({"point_mm": list(point), "ideal_inside": expected, "step_inside": actual})
        if expected != actual:
            mismatches.append(records[-1])
    return {"decisions": len(records), "mismatches": len(mismatches),
            "witnesses": mismatches[:8], "samples": records}


def control_digest(model: dict) -> str:
    """Hash nets after removing one positive integer gauge per chart."""
    normalized = []
    for block in sorted(model["blocks"], key=lambda b: b["id"]):
        flat = [int(c) for point in block["controls"] for c in point]
        g = 0
        for c in flat:
            g = math.gcd(g, abs(c))
        if not g:
            g = 1
        normalized.append({"id": block["id"], "degrees": block["degrees"],
                           "controls": [[int(c) // g for c in point]
                                        for point in block["controls"]]})
    raw = json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def cli_json(cli: Path, *args: object, expected: bool = True):
    proc = subprocess.run([str(cli), *map(str, args)], cwd=ROOT, capture_output=True, text=True)
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError:
        payload = {"raw_stdout": proc.stdout, "stderr": proc.stderr}
    if expected and proc.returncode:
        raise RuntimeError(f"CLI {args[0]} failed ({proc.returncode}): {payload}")
    return proc, payload


def make_corpus() -> list[dict]:
    named = [
        {"id": "minimum_t_minimum_r_h", "recipe": recipe(1, 4, 1, 2, 1)},
        {"id": "minimum_t_maximum_r_h", "recipe": recipe(1, 4, 1, 64, 64)},
        {"id": "quarter_turn", "recipe": recipe(1, 1, 1, 64, 64)},
        {"id": "below_135_convergent", "recipe": recipe(70, 29, 1, 63, 64)},
        {"id": "above_135_convergent", "recipe": recipe(169, 70, 63, 64, 1)},
        {"id": "intermediate_low", "recipe": recipe(1, 2, 2, 3, 32)},
        {"id": "near_split_below", "recipe": recipe(2, 1, 1, 2, 64)},
        {"id": "maximum_t_minimum_gap", "recipe": recipe(4, 1, 63, 64, 64)},
        {"id": "maximum_t_minimum_h", "recipe": recipe(4, 1, 1, 64, 1)},
        {"id": "equiv_half_angle_a", "recipe": recipe(1, 2, 3, 11, 17)},
        {"id": "equiv_half_angle_b", "recipe": recipe(2, 4, 3, 11, 17)},
        {"id": "oblique_rational", "recipe": recipe(7, 5, 4, 13, 29)},
    ]
    rng = random.Random(2026112524)
    generated = []
    for i in range(72):
        n, d = rng.randint(1, 256), rng.randint(1, 256)
        while not (4 * n >= d and n <= 4 * d):
            n, d = rng.randint(1, 256), rng.randint(1, 256)
        ri = rng.randint(1, 63)
        ro = rng.randint(ri + 1, 64)
        h = rng.randint(1, 64)
        generated.append({"id": f"seed-valid-{i:03d}", "expected_valid": True,
                          "recipe": recipe(n, d, ri, ro, h)})
    invalids = [
        recipe(0, 1, 1, 2, 1), recipe(1, 0, 1, 2, 1),
        recipe(1, 5, 1, 2, 1), recipe(5, 1, 1, 2, 1),
        recipe(1, 1, 2, 2, 4), recipe(1, 1, 1, 2, 0),
        recipe(1, 1, 1, 65, 4), recipe(257, 100, 1, 2, 4),
        recipe(1, 257, 1, 2, 4), recipe(2, 3, 0, 2, 4),
        recipe(2, 3, 1, 2, 65), recipe(3, 4, 64, 64, 4),
    ]
    # Also exercise JSON values that are numerically tempting but have the wrong type.
    invalids.extend([
        {"family": "rational_sector_prism_v1", "half_angle_tangent": [True, 1],
         "radii_mm": [1, 2], "height_mm": 4},
        {"family": "rational_sector_prism_v1", "half_angle_tangent": [1, 1],
         "radii_mm": [1, 2], "height_mm": True},
        {"family": "rational_sector_prism_v1", "half_angle_tangent": [1.0, 2],
         "radii_mm": [1, 2], "height_mm": 4},
        {"family": "rational_sector_prism_v1", "half_angle_tangent": [1, 2],
         "radii_mm": [1, 2], "height_mm": 4.0},
    ])
    while len(invalids) < 24:
        n, d = rng.randint(1, 256), rng.randint(1, 256)
        if 4 * n >= d and n <= 4 * d:
            invalids.append(recipe(n, d, 1, 1, rng.randint(1, 64)))
    generated.extend({"id": f"seed-invalid-{i:03d}", "expected_valid": False,
                      "recipe": rec} for i, rec in enumerate(invalids[:24]))
    for row in named:
        row["expected_valid"] = True
    return named + generated


def audit_mutations(cli: Path, model: dict, folder: Path) -> dict:
    source = folder / "mutation-source.json"
    write_json(source, model)
    results = {}

    common = copy.deepcopy(model)
    for block in common["blocks"]:
        for point in block["controls"]:
            for c in range(4):
                point[c] *= 7
    path = folder / "mutation-common-gauge.json"
    write_json(path, common)
    p, payload = cli_json(cli, "audit", path, expected=False)
    results["common_positive_gauge"] = {"accepted": p.returncode == 0, "report": payload,
                                        "digest_equal": control_digest(common) == control_digest(model)}

    per_chart = copy.deepcopy(model)
    for k, block in enumerate(per_chart["blocks"]):
        factor = 3 + 2 * k
        for point in block["controls"]:
            for c in range(4):
                point[c] *= factor
    path = folder / "mutation-per-chart-gauge.json"
    write_json(path, per_chart)
    p, payload = cli_json(cli, "audit", path, expected=False)
    results["per_chart_positive_gauge"] = {"accepted": p.returncode == 0, "report": payload,
                                           "digest_equal": control_digest(per_chart) == control_digest(model)}

    bad_weight = copy.deepcopy(model)
    bad_weight["blocks"][0]["controls"][0][3] += 1
    path = folder / "mutation-weight-plus-one.json"
    write_json(path, bad_weight)
    p, payload = cli_json(cli, "audit", path, expected=False)
    results["weight_plus_one_rejected"] = {"rejected": p.returncode != 0, "report": payload}

    bad_net = copy.deepcopy(model)
    point = bad_net["blocks"][0]["controls"][4]
    for c in range(4):
        point[c] *= 2
    path = folder / "mutation-one-homogeneous-control-gauge.json"
    write_json(path, bad_net)
    p, payload = cli_json(cli, "audit", path, expected=False)
    results["one_control_gauge_rejected"] = {"rejected": p.returncode != 0, "report": payload}
    return results


def refusal_cases(cli: Path, work: Path) -> list[dict]:
    cases = [
        ("zero_numerator", recipe(0, 1, 1, 2, 1)),
        ("zero_denominator", recipe(1, 0, 1, 2, 1)),
        ("t_below_domain", recipe(1, 5, 1, 2, 1)),
        ("t_above_domain", recipe(5, 1, 1, 2, 1)),
        ("equal_radii", recipe(1, 1, 2, 2, 4)),
        ("height_zero", recipe(1, 1, 1, 2, 0)),
        ("radius_over_limit", recipe(1, 1, 1, 65, 4)),
    ]
    target = work / "preserve-model.json"
    sentinel = b"old valid output sentinel\n"
    rows = []
    for name, rec in cases:
        rp = work / f"invalid-{name}.json"
        write_json(rp, rec)
        target.write_bytes(sentinel)
        p, payload = cli_json(cli, "build", rp, target, expected=False)
        rows.append({"name": name, "refused": p.returncode != 0,
                     "preserved": target.read_bytes() == sentinel, "response": payload})
    return rows


def run(cli: Path, work: Path) -> dict:
    started = time.monotonic()
    work.mkdir(parents=True, exist_ok=True)
    binary_hash = hashlib.sha256(cli.read_bytes()).hexdigest()
    mass_helper = compile_mass_oracle(work)
    cases = make_corpus()
    case_rows = []
    fingerprints = {}
    for item in cases:
        name, rec = item["id"], item["recipe"]
        folder = work / name
        folder.mkdir(parents=True, exist_ok=True)
        rp, mp = folder / "recipe.json", folder / "model.json"
        write_json(rp, rec)
        old = mp.read_bytes() if mp.exists() else None
        if not item.get("expected_valid", True):
            sentinel = b"previous model must survive rejected recipe\n"
            mp.write_bytes(sentinel)
            old = sentinel
        build_proc, build_report = cli_json(cli, "build", rp, mp, expected=False)
        row = {"id": name, "recipe": rec, "expected_split": split_expected(*rec["half_angle_tangent"]),
               "expected_valid": item.get("expected_valid", True),
               "build_returncode": build_proc.returncode, "build_report": build_report}
        if not item.get("expected_valid", True):
            row["preserved_on_refusal"] = mp.read_bytes() == old
            row["qualified"] = build_proc.returncode != 0 and row["preserved_on_refusal"]
            case_rows.append(row)
            print(name, row["qualified"], flush=True)
            continue
        if build_proc.returncode:
            row["qualified"] = False
            case_rows.append(row)
            continue
        model = json.loads(mp.read_text())
        audit_proc, audit_report = cli_json(cli, "audit", mp, expected=False)
        volume, center = volume_and_center(rec)
        n, d = rec["half_angle_tangent"]
        coefficient = exact_volume(rec)
        rate = coefficient["atan_coefficient_mm3"] / (1.0 + (n / d) ** 2)
        angular_sweep = 2 * math.atan(n / d)
        expected_blocks = 2 if row["expected_split"] else 1
        expected_tags = {"angle_start", "angle_end", "inner", "outer", "start", "end"}
        actual_tags = {face["feature"] for face in model["design_faces"]}
        row.update({"audit_returncode": audit_proc.returncode, "audit": audit_report,
                    "expected_volume_exact": coefficient, "analytic_volume_mm3": volume,
                    "analytic_centroid_mm": center, "analytic_theta_rad": angular_sweep,
                    "analytic_volume_rate_dt_mm3": rate,
                    "block_count": len(model["blocks"]), "design_tags": sorted(actual_tags),
                    "partition_fingerprint": control_digest(model)})
        fingerprints[name] = row["partition_fingerprint"]
        row["qualified"] = bool(
            audit_proc.returncode == 0
            and audit_report.get("status") == "valid"
            and audit_report.get("volume_exact") == coefficient
            and len(model["blocks"]) == expected_blocks
            and actual_tags == expected_tags
            and abs(audit_report.get("volume_rate_dt_mm3", math.nan) / rate - 1.0) < 1e-14
        )
        if name == "above_135_convergent":
            row["mutations"] = audit_mutations(cli, model, folder)
            row["qualified"] &= all([
                row["mutations"]["common_positive_gauge"]["accepted"],
                row["mutations"]["common_positive_gauge"]["digest_equal"],
                row["mutations"]["per_chart_positive_gauge"]["accepted"],
                row["mutations"]["per_chart_positive_gauge"]["digest_equal"],
                row["mutations"]["weight_plus_one_rejected"]["rejected"],
                row["mutations"]["one_control_gauge_rejected"]["rejected"],
            ])
        case_rows.append(row)
        print(name, row["qualified"], flush=True)

    # Rationally equivalent parameter pairs must reduce to the same stored net.
    fingerprints_equal = fingerprints.get("equiv_half_angle_a") == fingerprints.get("equiv_half_angle_b")
    for row in case_rows:
        if row["id"] in ("equiv_half_angle_a", "equiv_half_angle_b"):
            row["equivalent_fingerprint"] = fingerprints_equal
            row["qualified"] &= fingerprints_equal

    # STEP round trips: bounds include both parameter limits, both sides of 135 deg,
    # the minimum and maximum radii/gap, and both height limits.
    step_ids = ["minimum_t_minimum_r_h", "minimum_t_maximum_r_h", "quarter_turn",
                "below_135_convergent", "above_135_convergent", "near_split_below",
                "maximum_t_minimum_gap", "maximum_t_minimum_h"]
    steps = []
    all_decisions = 0
    for seed, name in enumerate(step_ids, 1201):
        rec = next(row["recipe"] for row in cases if row["id"] == name)
        folder = work / name
        mp = folder / "model.json"
        step = folder / "roundtrip.step"
        before = time.monotonic()
        p, payload = cli_json(cli, "export-step", mp, step, expected=False)
        row = {"id": name, "export_returncode": p.returncode, "export_report": payload,
               "export_seconds": time.monotonic() - before}
        if p.returncode:
            row["qualified"] = False
            steps.append(row)
            continue
        imported = Part.read(str(step))
        reference = reference_solid(rec)
        solids = list(imported.Solids)
        solid = solids[0] if solids else imported
        expected_volume, expected_center = volume_and_center(rec)
        volume_error = abs(float(solid.Volume) - expected_volume) / expected_volume
        center_error = math.dist(tuple(solid.CenterOfMass), expected_center)
        expected_split = split_expected(*rec["half_angle_tangent"])
        charts = 2 if expected_split else 1
        expected_topology = {"faces": 4 * charts + 2, "edges": 8 * charts + 4,
                             "vertices": 4 * charts + 4}
        actual_topology = {"faces": len(imported.Faces), "edges": len(imported.Edges),
                           "vertices": len(imported.Vertexes)}
        step_brep, reference_brep = folder / "step-import.brep", folder / "ideal-revolution.brep"
        solid.exportBrep(str(step_brep))
        reference.exportBrep(str(reference_brep))
        step_adaptive = json.loads(subprocess.check_output(
            [str(mass_helper), str(step_brep), "--centroid"], text=True))
        reference_adaptive = json.loads(subprocess.check_output(
            [str(mass_helper), str(reference_brep), "--centroid"], text=True))
        step_mass = step_adaptive[-1]
        reference_mass = reference_adaptive[-1]
        adaptive_volume_error = abs(step_mass["volume_mm3"] - expected_volume) / expected_volume
        adaptive_center_error = math.dist(step_mass["centroid_mm"], expected_center)
        brep_pair_volume_error = abs(step_mass["volume_mm3"] - reference_mass["volume_mm3"]) / expected_volume
        brep_pair_center_error = math.dist(step_mass["centroid_mm"], reference_mass["centroid_mm"])
        classifications = classify_points(rec, solid, seed)
        reference_volume_error = abs(float(reference.Volume) - expected_volume) / expected_volume
        reference_center_error = math.dist(tuple(reference.CenterOfMass), expected_center)
        row.update({"solid_count": len(solids), "valid": bool(imported.isValid()),
                    "closed": bool(len(solids) == 1 and solids[0].isClosed()),
                    "step_volume_mm3": float(solid.Volume), "analytic_volume_mm3": expected_volume,
                    "volume_relative_error": volume_error,
                    "center_mm": list(solid.CenterOfMass), "analytic_center_mm": expected_center,
                    "center_error_mm": center_error,
                    "reference_revolution_volume_mm3": float(reference.Volume),
                    "reference_revolution_center_mm": list(reference.CenterOfMass),
                    "reference_volume_relative_error": reference_volume_error,
                    "reference_center_error_mm": reference_center_error,
                    "step_adaptive_mass": step_adaptive,
                    "reference_adaptive_mass": reference_adaptive,
                    "adaptive_volume_relative_error": adaptive_volume_error,
                    "adaptive_center_error_mm": adaptive_center_error,
                    "adaptive_brep_pair_volume_relative_error": brep_pair_volume_error,
                    "adaptive_brep_pair_center_error_mm": brep_pair_center_error,
                    "expected_topology": expected_topology, "step_topology": actual_topology,
                    "topology_matches": actual_topology == expected_topology,
                    "classifications": classifications})
        row["qualified"] = bool(
            row["solid_count"] == 1 and row["valid"] and row["closed"]
            and adaptive_volume_error < 1e-8 and adaptive_center_error < 1e-7
            and brep_pair_volume_error < 1e-8 and brep_pair_center_error < 1e-7
            and reference_volume_error < 1e-12 and reference_center_error < 1e-10
            and actual_topology == expected_topology
            and classifications["decisions"] >= 12 and classifications["mismatches"] == 0
        )
        all_decisions += classifications["decisions"]
        steps.append(row)
        print("STEP", name, row["qualified"], classifications["decisions"], flush=True)

    refusals = refusal_cases(cli, work)
    qualified = all(r.get("qualified", False) for r in case_rows)
    qualified &= len(steps) == 8 and all(s.get("qualified", False) for s in steps)
    qualified &= all_decisions >= 96
    qualified &= all(r["refused"] and r["preserved"] for r in refusals)
    if binary_hash != hashlib.sha256(cli.read_bytes()).hexdigest():
        raise RuntimeError("CLI executable changed during G24 oracle run")
    return {"schema": "rational_sector_independent_oracle_v1",
            "provenance": {"cli": str(cli), "cli_sha256": binary_hash,
                           "checkout_revision": subprocess.check_output(
                               ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                           "python": sys.version.split()[0], "freecad": App.Version()[0],
                           "occ": Part.OCC_VERSION,
                           "command": "rtk env PYTHONPATH=/home/jeff/Code/FreeCAD/.pixi/envs/default/lib LD_LIBRARY_PATH=/home/jeff/Code/FreeCAD/.pixi/envs/default/lib OPENBLAS_NUM_THREADS=1 /home/jeff/Code/FreeCAD/.pixi/envs/default/bin/python tests/references/rational_sector_oracle.py --cli target/release/solidcraft-trivariate --work build/rational-sector-g24-final3 --output docs/results/rational-sector-20261125.json",
                           "elapsed_s": time.monotonic() - started},
            "protocol": {"rational_domain": "1/4 <= n/d <= 4, 1<=n,d<=256",
                         "fixed_seed": 2026112524, "generated_positive_cases": 72,
                         "generated_negative_cases": 24,
                         "radii_mm": "1<=ri<ro<=64", "height_mm": "1<=H<=64",
                         "split_test": "n>d and (n-d)^2>2*d^2",
                         "volume_formula": "H*(ro^2-ri^2)*2*atan(n/d)",
                         "volume_derivative": "H*(ro^2-ri^2)/(1+(n/d)^2)",
                         "step_cases": step_ids, "minimum_step_classifications": 96,
                         "step_volume_relative_tolerance": 1e-8,
                         "step_center_tolerance_mm": 1e-7,
                         "adaptive_occ_integrator": "BRepGProp::VolumePropertiesGK and SurfaceProperties, eps=1e-12",
                         "default_FreeCAD_mass_properties": "diagnostic only; qualification uses both adaptive BREP measurements"},
            "cases": case_rows, "step_roundtrips": steps,
            "classification_decisions": all_decisions, "refusals": refusals,
            "qualified": bool(qualified)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=os.environ.get("SC_CLI"))
    parser.add_argument("--work", type=Path, default=os.environ.get("SC_WORK"))
    parser.add_argument("--output", type=Path, default=os.environ.get("SC_OUTPUT"))
    args = parser.parse_args()
    if not args.cli or not args.work or not args.output:
        parser.error("provide --cli/--work/--output or SC_CLI/SC_WORK/SC_OUTPUT")
    report = run(args.cli.resolve(), args.work.resolve())
    write_json(args.output, report)
    print("G24 rational sector oracle", report["qualified"])
    raise SystemExit(int(not report["qualified"]))
