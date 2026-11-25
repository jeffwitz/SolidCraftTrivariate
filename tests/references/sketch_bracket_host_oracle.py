"""Compare the automatically partitioned bracket with SolidCraft's saved STEP."""
import FreeCAD as App
import Part

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path


from coaxial_program import compile_mass_oracle
from solidcraft_bridge import samples, shape_member


def replay(cli, host_step, work):
    work.mkdir(parents=True, exist_ok=True)
    recipe = Path(__file__).with_name("orthogonal_sketch_bracket.recipe.json")
    model = work / "automatic.model.json"
    step = work / "automatic.step"
    build = json.loads(subprocess.check_output([cli, "build", recipe, model]))
    subprocess.run([cli, "export-step", model, step], check=True)
    helper = compile_mass_oracle(work)
    shapes = [Part.read(str(path)) for path in (step, host_step)]
    records = []
    for name, shape in zip(("automatic", "solidcraft"), shapes):
        brep = work / (name + ".brep")
        shape.exportBrep(str(brep))
        mass = json.loads(subprocess.check_output([helper, brep, "--centroid"]))[-1]
        records.append(dict(source=name, valid=shape.isValid(), solids=len(shape.Solids),
                            shells=len(shape.Shells), faces=len(shape.Faces), **mass))
    bad = []
    count = 0
    for point in samples(4, 4):
        expected = shape_member(*point, 4, 4)
        observed = [s.isInside(App.Vector(*point), 1e-8, False) for s in shapes]
        if any(v != expected for v in observed):
            bad.append(dict(point_mm=point, expected=expected, observed=observed))
        count += 1
    target = 34080 - 616 * math.pi
    volume_error = max(abs(r["volume_mm3"] / target - 1) for r in records)
    centroid_error = math.dist(*(r["centroid_mm"] for r in records))
    area_error = abs(records[0]["surface_mm2"] / records[1]["surface_mm2"] - 1)
    valid = all(r["valid"] and r["solids"] == r["shells"] == 1 for r in records)
    return dict(cli_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),
                solidcraft_step_sha256=hashlib.sha256(host_step.read_bytes()).hexdigest(),
                cad_engine=Part.OCC_VERSION, build=build, shapes=records,
                exact_volume_mm3=target, max_relative_volume_error=volume_error,
                relative_surface_area_error=area_error, centroid_difference_mm=centroid_error,
                classification_count=count, mismatches=bad, qualified=valid and not bad
                and volume_error <= 1e-8 and area_error <= 1e-8 and centroid_error <= 1e-7)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("cli", "host-step", "work", "output"):
        p.add_argument("--" + name, required=True, type=Path)
    args = p.parse_args()
    result = replay(args.cli.resolve(), args.host_step, args.work)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    raise SystemExit(int(not result["qualified"]))
