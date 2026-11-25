#!/usr/bin/env python3
"""Independent SciPy volume oracle and optional FreeCAD STEP check.

The recipe stores the two positive-side meridians in homogeneous form
``(w * half_width, w * ordinate, w)``. Each is reflected about the bore's
mid-plane to obtain a closed symmetric opening. This script does not call the
SolidCraft builder or audit.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import BSpline
from scipy.optimize import brentq


def half_width_at(spec: dict, ordinate: float) -> float:
    """Evaluate the exact rational meridian numerically at |y| in [0, 5]."""
    degree = spec["degree"]
    knots = np.asarray(spec["knots"], dtype=float)
    controls = np.asarray(spec["controls"], dtype=float)
    homogeneous = BSpline(knots, controls, degree, extrapolate=False)
    if not 0.0 <= ordinate <= 5.0:
        raise ValueError("ordinate lies outside the meridian")

    def y_at(t: float) -> float:
        point = np.asarray(homogeneous(t), dtype=float)
        return float(point[1] / point[2])

    def x_at(t: float) -> float:
        point = np.asarray(homogeneous(t), dtype=float)
        return float(point[0] / point[2])

    if ordinate == 0.0:
        return x_at(0.0)
    if ordinate == 5.0:
        return x_at(3.0)
    parameter = brentq(lambda t: y_at(t) - ordinate, 0.0, 3.0, xtol=1e-14)
    return x_at(parameter)


def independent_volume(recipe: dict) -> dict:
    bounds = recipe["bounds_mm"]
    center = recipe["center_mm"]
    mx = recipe["meridians"]["x"]
    mz = recipe["meridians"]["z"]
    for spec in (mx, mz):
        if spec["degree"] != 3 or spec["knots"] != [0, 0, 0, 0, 1, 3, 3, 3, 3]:
            raise ValueError("unexpected reference meridian parameterization")
        weights = [p[2] for p in spec["controls"]]
        ordinate = [p[1] / p[2] for p in spec["controls"]]
        widths = [p[0] / p[2] for p in spec["controls"]]
        if weights != [1, 2, 1, 2, 1] or ordinate != [0, 1, 3, 5, 5]:
            raise ValueError("reference weights or ordinate poles were changed")
        if widths[0] != widths[1] or widths[-1] != 0 or any(w < 0 for w in widths):
            raise ValueError("meridian endpoint tangencies or nonnegative widths changed")
    if [p[0] / p[2] for p in mx["controls"]] != [6, 6, 4, 2, 0]:
        raise ValueError("x-meridian half-width poles changed")
    if [p[0] / p[2] for p in mz["controls"]] != [8, 8, 6, 3, 0]:
        raise ValueError("z-meridian half-width poles changed")
    if bounds != [[0, 60], [0, 40], [0, 30]] or center != [30, 20, 15]:
        raise ValueError("unexpected reference stock or center")

    # At ordinate y relative to the common mid-plane, the two transverse
    # openings remove strips of widths 2*x(y) and 2*z(y); their overlap is
    # counted once by the product of the remaining in-plane dimensions.
    def remaining_area(y: float) -> float:
        x = half_width_at(mx, y)
        z = half_width_at(mz, y)
        if x < 0 or z < 0:
            raise ValueError("negative meridian half-width")
        return (60.0 - 2.0 * x) * (30.0 - 2.0 * z)

    # Split at the knot's physical ordinate, where the rational basis changes.
    y_knots = []
    for spec in (mx, mz):
        c = np.asarray(spec["controls"], dtype=float)
        spline = BSpline(np.asarray(spec["knots"], dtype=float), c, 3, extrapolate=False)
        h = np.asarray(spline(1.0), dtype=float)
        y_knots.append(float(h[1] / h[2]))
    interior = sorted(set(y_knots))
    integral, error = quad(
        remaining_area,
        0.0,
        5.0,
        points=interior,
        epsabs=2e-10,
        epsrel=2e-12,
        limit=400,
    )
    stock_area = (bounds[0][1] - bounds[0][0]) * (bounds[2][1] - bounds[2][0])
    volume = stock_area * ((bounds[1][1] - bounds[1][0]) - 10.0) + 2.0 * integral
    return {
        "oracle": "scipy_bspline_inverse_and_quad",
        "volume_mm3": volume,
        "quad_abs_error_mm3": 2.0 * error,
        "symmetric_half_height_mm": 5.0,
        "meridian_knot_ordinates_mm": interior,
        "cross_section_area_at_center_mm2": remaining_area(0.0),
        "cross_section_area_at_pole_mm2": remaining_area(5.0),
    }


def freecad_step_check(step_path: Path, recipe: dict, reference_volume: float) -> dict:
    """Read STEP in FreeCAD/OCC and classify material/void near both poles."""
    center = recipe["center_mm"]
    points = {
        "north_void": [center[0], center[1] + 4.999, center[2]],
        "north_material": [center[0] + 0.5, center[1] + 4.999, center[2] + 0.5],
        "south_void": [center[0], center[1] - 4.999, center[2]],
        "south_material": [center[0] + 0.5, center[1] - 4.999, center[2] + 0.5],
    }
    expected = {
        "north_void": False,
        "north_material": True,
        "south_void": False,
        "south_material": True,
    }
    freecad = os.environ.get("FREECAD_BIN") or shutil.which("freecadcmd") or shutil.which("freecad") or "/home/jeff/.local/bin/freecadcmd"
    with tempfile.TemporaryDirectory(prefix="solidcraft-nurbs-oval-") as work:
        work_path = Path(work)
        report_path = work_path / "step-report.json"
        brep_path = work_path / "oval.brep"
        macro_path = work_path / "inspect.FCMacro"
        payload = json.dumps({"step": str(step_path.resolve()), "report": str(report_path),
                              "brep": str(brep_path), "points": points})
        macro = f'''import json
import sys
import Part
from FreeCAD import Vector
cfg = json.loads({payload!r})
shape = Part.Shape()
shape.read(cfg["step"])
shape.exportBrep(cfg["brep"])
solids = list(shape.Solids)
classes = {{name: bool(shape.isInside(Vector(*p), 1e-7, True)) for name, p in cfg["points"].items()}}
data = {{"valid": bool(shape.isValid()), "solid_count": len(solids),
        "closed": all(bool(s.isClosed()) for s in solids),
        "volume_mm3": float(shape.Volume), "classifications": classes,
        "python_prefix": sys.prefix}}
with open(cfg["report"], "w", encoding="utf-8") as out:
    json.dump(data, out, sort_keys=True)
'''
        macro_path.write_text(macro, encoding="utf-8")
        completed = subprocess.run(
            [freecad, str(macro_path)],
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        if completed.returncode != 0 or not report_path.exists():
            raise RuntimeError(
                "FreeCAD STEP check failed: "
                + (completed.stderr[-2000:] or completed.stdout[-2000:])
            )
        report = json.loads(report_path.read_text(encoding="utf-8"))
        prefix = Path(report["python_prefix"])
        helper = work_path / "occ-volume"
        source = Path(__file__).with_name("occ_volume.cpp")
        subprocess.run(
            [
                "g++", "-std=c++17", "-O2",
                "-I" + str(prefix / "include/opencascade"), str(source),
                "-L" + str(prefix / "lib"), "-Wl,-rpath," + str(prefix / "lib"),
                "-lTKTopAlgo", "-lTKBRep", "-lTKG3d", "-lTKMath", "-lTKernel",
                "-o", str(helper),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        adaptive = json.loads(subprocess.check_output([str(helper), str(brep_path), "--centroid"], text=True))
        adaptive = adaptive[-1] if isinstance(adaptive, list) else adaptive
        if adaptive.get("tolerance") != 1e-12:
            raise RuntimeError(f"expected the 1e-12 adaptive OCCT measure, got {adaptive}")
    report["volume_default_mm3"] = report["volume_mm3"]
    report["volume_default_relative_error"] = abs(report["volume_default_mm3"] - reference_volume) / reference_volume
    report["adaptive_measure"] = adaptive
    report["adaptive_volume_mm3"] = adaptive["volume_mm3"]
    report["adaptive_volume_relative_error"] = abs(report["adaptive_volume_mm3"] - reference_volume) / reference_volume
    report["adaptive_tolerance"] = adaptive["tolerance"]
    report["volume_default_method"] = "FreeCAD Part.Shape.Volume with default OCCT settings"
    report["volume_adaptive_method"] = "OCCT BRepGProp::VolumePropertiesGK on the FreeCAD-imported STEP BREP"
    report["classifications_match"] = report["classifications"] == expected
    report["volume_relative_threshold"] = 2e-8
    report["volume_threshold_applies_to"] = "adaptive_volume_relative_error"
    report["volume_match"] = report["adaptive_volume_relative_error"] <= report["volume_relative_threshold"]
    if not report["valid"] or report["solid_count"] != 1 or not report["closed"]:
        raise AssertionError(f"STEP is not one valid closed solid: {report}")
    if not report["classifications_match"]:
        raise AssertionError(f"STEP pole-neighborhood classifications differ: {report}")
    if not report["volume_match"]:
        raise AssertionError(f"STEP volume differs from independent SciPy oracle: {report}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "recipe",
        nargs="?",
        type=Path,
        default=Path(__file__).with_name("nurbs_cross_oval.recipe.json"),
    )
    parser.add_argument("--step", type=Path, help="optional SolidCraft STEP export to reread in FreeCAD/OCC")
    parser.add_argument("--freecad", default=os.environ.get("FREECAD_BIN") or shutil.which("freecadcmd") or shutil.which("freecad") or "/home/jeff/.local/bin/freecadcmd",
                        help="FreeCAD or FreeCADCmd executable (default: FREECAD_BIN, PATH, then local installation)")
    args = parser.parse_args()
    recipe = json.loads(args.recipe.read_text(encoding="utf-8"))
    result = independent_volume(recipe)
    if args.step:
        os.environ["FREECAD_BIN"] = args.freecad
        result["step_check"] = freecad_step_check(args.step, recipe, result["volume_mm3"])
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
