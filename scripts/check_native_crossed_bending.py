"""Replay each published partial-interface elasticity case with the native solver."""
import argparse
import json
import math
import subprocess
import tempfile
from pathlib import Path


def replay(binary):
    root = Path(__file__).resolve().parents[1]
    reference = root / "docs/results/crossed-interface-iga/partial_bending_study.json"
    rows = json.loads(reference.read_text())
    recipe = json.loads((root / "tests/jobs/crossed_partial_rise3.recipe.json").read_text())
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        result = work / "solution.json"
        for row in rows:
            rise = row["rise_mm"]
            model = work / f"union-{rise}.json"
            if not model.exists():
                recipe["transverse_rise_mm"] = rise
                source = work / "recipe.json"
                source.write_text(json.dumps(recipe))
                subprocess.run([binary, "crossed-boolean", source, "union", model], check=True,
                               stdout=subprocess.DEVNULL)
            options = [str(row[k]) for k in ("refine_v", "refine_w", "quadrature_order")]
            call = subprocess.run([binary, "crossed-bending", model, *options, result],
                                  capture_output=True, text=True)
            if call.returncode:
                raise RuntimeError(f"Rise {rise}, options {options}: {call.stdout} {call.stderr}")
            native = json.loads(result.read_text())
            for key in ("free_vector_dofs", "partial_trace_rank", "partial_trace_constraint_rows"):
                if native[key] != row[key]:
                    raise ValueError((rise, options, key, native[key], row[key]))
            for key in ("loaded_end_mean_ux_mm", "strain_energy_Nmm", "volume_mm3", "force_N"):
                if not math.isclose(native[key], row[key], rel_tol=3e-8, abs_tol=0):
                    raise ValueError((rise, options, key, native[key], row[key]))
            print(f"rise={rise:+} v,w,q={options} rank={native['partial_trace_rank']} "
                  f"DOFs={native['free_vector_dofs']} ux={native['loaded_end_mean_ux_mm']:.10f} mm",
                  flush=True)
    print(f"Native displacement, energy, volume, force, trace ranks and DOFs match all {len(rows)} independent rows.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    replay(parser.parse_args().binary.resolve())
