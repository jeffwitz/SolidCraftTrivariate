"""Compare native partitions at fixed physical geometry with one frozen CLI."""
import FreeCAD
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import statistics
import subprocess
import time
import numpy as np
from sketch_partition_oracle import jacobian_at

HERE = Path(__file__).resolve().parent


def measure(cli, work, recipe, label):
    path = work / (label + ".recipe.json")
    model = work / (label + ".model.json")
    path.write_text(json.dumps(recipe))
    durations = []
    for i in range(6):
        start = time.perf_counter()
        report = json.loads(subprocess.check_output([cli, "build", path, model]))
        if i:
            durations.append(time.perf_counter()-start)
    design = json.loads(model.read_text())
    points = list(itertools.product((0., .5, 1.), repeat=3))
    worst_condition = 0.
    min_scaled = 1.
    for block in design["blocks"]:
        for uvw in points:
            j = jacobian_at(block, uvw)
            sigmas = np.linalg.svd(j, compute_uv=False)
            worst_condition = max(worst_condition, float(sigmas[0]/sigmas[-1]))
            scaled = np.linalg.det(j)/np.prod(np.linalg.norm(j,axis=1))
            min_scaled = min(min_scaled,float(scaled))
    return dict(recipe=recipe, blocks=len(design["blocks"]), interfaces=report["shared_faces"],
                model_bytes=model.stat().st_size, build_and_audit_median_s=statistics.median(durations),
                repeats_s=durations, min_scaled_jacobian_sample=min_scaled,
                max_condition_sample=worst_condition, audited_min_det_mm3=report.get("determinant_lower_bound_mm3"),
                volume_exact=report.get("volume_exact_mm3",report.get("volume_exact")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cli","work","output"):
        parser.add_argument("--"+name,type=Path,required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    corpus = {entry["id"]:entry["recipe"] for entry in json.loads((HERE/"sketch_partition_corpus.json").read_text())}
    cases = {k:corpus[k] for k in ("plate_blind_pocket","plate_blind_pocket_triangulate",
                                   "c_section_two_holes","c_section_two_holes_triangulate")}
    cases["automatic_bracket"] = json.loads((HERE/"orthogonal_sketch_bracket.recipe.json").read_text())
    cases["specialized_bracket"] = json.loads((HERE/"filleted_bracket.json").read_text())
    cli=args.cli.resolve()
    rows={name:measure(cli,args.work,recipe,name) for name,recipe in cases.items()}
    result=dict(protocol="One warmup then five sequential CLI build+audit wall times, concurrent campaign traffic; 27 Jacobian samples per block.",
                python=platform.python_version(),platform=platform.platform(),numpy=np.__version__,
                cli_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),cases=rows)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
