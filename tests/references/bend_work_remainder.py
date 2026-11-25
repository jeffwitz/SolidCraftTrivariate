"""Resolve the high-order compliance remainder of a sheared, clamped circular bend."""
import argparse
import copy
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from affine_job_derivative import cases, image
from mechanical_job import assemble, independent_load, independently_solve
from sector_job_sensitivity import call, write


def measure(cli, folder, source, specification, matrix, s, order):
    transformation = image(source, 1, matrix, [0,0,0], s)
    write(folder/"recipe.json", transformation)
    write(folder/"job.json", dict(specification, order=order))
    call(cli, "build", folder/"recipe.json", folder/"model.json")
    call(cli, "solve-job", folder/"model.json", folder/"job.json", folder/"result.json")
    model, result = [json.loads((folder/name).read_text()) for name in ["model.json","result.json"]]
    blocks = result["analysis_blocks"]
    mat = specification["material"]
    k, n = assemble(blocks, mat["young_mpa"], mat["poisson"], order)
    f = independent_load(model, result, specification, order, n)[0]
    u, _ = independently_solve(k, f, result, n)
    # With nu=0 the unmounted free stress state has strain sigma/E.
    stress = np.asarray(specification["loads"][0]["stress_mpa"])
    x = np.empty((n,3))
    for b in blocks:
        h = np.asarray(b["controls"], dtype=float)
        x[b["dofs"]] = h[:, :3]/h[:,3,None]
    free = (x @ stress.T / mat["young_mpa"]).ravel()
    delta = free - u
    w = float(f @ u)
    free_work = float(f @ free)
    norm2 = float(delta @ (k @ delta))
    return dict(step=float(s), order=order, cold_work=w, native_work=result["external_work_n_mm"],
                free_work=free_work, energy_defect_squared=norm2,
                identity_relative=abs(free_work-w-norm2)/free_work,
                free_residual=float(np.linalg.norm(k@free-f)/np.linalg.norm(f)),
                field_relative=float(np.linalg.norm(u-result["displacement_mm_xyz"])/np.linalg.norm(u)))


def run(cli, work):
    _, source, job, _, b, _, _ = next(c for c in cases() if c[0]=="curved_stress")
    work.mkdir(exist_ok=True, parents=True)
    rows = []
    for name in ["uniform_yy", "mixed_xy"]:
        spec = copy.deepcopy(job)
        if name=="mixed_xy":
            spec["loads"][0]["stress_mpa"] = [[1,0.3,0],[0.3,2,0],[0,0,0]]
        fit = []
        for q in ([10, 14] if name=="uniform_yy" else [10]):
            folder = work/f"{name}-q{q}"
            folder.mkdir(exist_ok=True)
            zero = measure(cli,folder,source,spec,b,Fraction(0),q)
            direction = dict(schema="solidcraft_affine_velocity_v1", linear=np.asarray(b).tolist(), offset_mm=[0,0,0])
            write(folder/"dir.json",direction)
            call(cli,"differentiate-job",folder/"model.json",folder/"job.json",folder/"dir.json",folder/"jet.json")
            jet=json.loads((folder/"jet.json").read_text())
            tangent=jet["compliance_rate_tangent_n_mm"]
            obs=[]
            for den in [4,8,16,32,64]:
                for sign in [1,-1]:
                    v=measure(cli,folder,source,spec,b,Fraction(sign,den),q)
                    v["absolute_remainder"]=abs(v["cold_work"]-zero["cold_work"]-v["step"]*tangent)
                    obs.append(v)
            powers=[]
            for sign in [1,-1]:
                signed=[r for r in obs if r["step"]*sign>0]
                powers.append([math.log2(a["absolute_remainder"]/c["absolute_remainder"]) for a,c in zip(signed,signed[1:])])
            norms=[[math.log2(a["energy_defect_squared"]/c["energy_defect_squared"])/2
                    for a,c in zip(obs[i::2],obs[i::2][1:])] for i in [0,1]]
            checks=dict(direct_field=max(v["field_relative"] for v in [zero,*obs])<1e-8,
                        equilibrium=max(v["free_residual"] for v in [zero,*obs])<1e-10,
                        energy_identity=max(v["identity_relative"] for v in [zero,*obs])<1e-10)
            if name=="uniform_yy":
                checks.update(quartic=all(3.8<=p<=4.2 for seq in powers for p in seq[-2:]),
                              quadratic_energy_norm=all(1.9<=p<=2.1 for seq in norms for p in seq[-2:]))
            fit.append(dict(order=q, zero=zero, observations=obs, tangent=tangent,
                            signed_remainder_orders=powers, signed_energy_norm_orders=norms, checks=checks))
            print(name,q,checks,powers,flush=True)
        if len(fit)==2:
            difference=max(abs(a["cold_work"]/b["cold_work"]-1)
                           for a,b in zip(fit[0]["observations"],fit[1]["observations"]))
        else:
            difference=None
        rows.append(dict(id=name,job=spec, fits=fit,cross_order_relative_work=difference,
                         passed=all(all(f["checks"].values()) for f in fit)
                                and (difference is None or difference<1e-8)))
    return dict(schema="bend_taylor_independent_oracle_v1",source=source,cases=rows,
                executable_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),qualified=all(v["passed"] for v in rows))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args();data=run(a.cli.resolve(),a.work);write(a.output,data)
    raise SystemExit(int(not data["qualified"]))
