"""Closed-form bar and material-derivative oracle for affine body loads."""
import argparse
import copy
from fractions import Fraction
import json
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial import Polynomial

from affine_job_derivative import image
from equivalent_prism_mechanics import plan


def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def solve_command(cli, *args):
    cmd = subprocess.run([str(cli), *map(str, args)], capture_output=True, text=True)
    if cmd.returncode:
        raise RuntimeError(f"{args[0]}: {cmd.stdout} {cmd.stderr}")


def run(cli, folder):
    folder.mkdir(exist_ok=True, parents=True)
    corners = [[-1, 2], [5, 2], [5, 7], [-1, 7]]
    height, area, young = 8, 30, 2000
    source = dict(family="polygon_prism_v1", height_mm=height, top_shift_mm=[0, 0],
                  vertices=[dict(id=f"p{i}", xy_mm=v) for i, v in enumerate(corners)])
    recipe = folder / "source.json"
    model = folder / "bar.json"
    write(recipe, source)
    solve_command(cli, "build", recipe, model)
    roots = [b["id"] for b in json.loads(model.read_text())["blocks"]]
    rows = []
    for label, a, b, levels in [("quadratic", .00025, 0., [0, 1]), ("cubic", -.0001, .00003, [0, 1, 2])]:
        solution = Polynomial([0., 0., a, b])
        stress = young * solution.deriv()
        density = -stress.deriv()
        traction = stress(height)
        force_z = area * density.integ()(height)
        exact_energy = area * (stress * solution.deriv()).integ()(height) / 2
        for h in levels:
            dz = height / (1 << h)
            job = dict(schema="solidcraft_mechanical_job_v4",
                       material=dict(young_mpa=young, poisson=0), formulation="displacement", order=5,
                       plan=plan(roots, h), supports=[dict(select=dict(feature="start"), components=[0, 1, 2])],
                       loads=[dict(select=dict(feature="end"), vector_mpa=[0, 0, traction])],
                       body_load=dict(constant_n_per_mm3=[0, 0, -2*young*a],
                                      gradient_n_per_mm4=[[0, 0, 0], [0, 0, 0], [0, 0, -6*young*b]]))
            jp, out = folder / "job.json", folder / "solution.json"
            write(jp, job)
            solve_command(cli, "solve-job", model, jp, out)
            res = json.loads(out.read_text())
            disp = np.array(res["displacement_mm_xyz"]).reshape(-1, 3)
            max_err = 0.
            # In each axial interval the derivative of the Ritz solution is the
            # L2 projection of the exact quadratic strain onto linear functions.
            for block in res["analysis_blocks"]:
                dofs = np.array(block["dofs"]).reshape(3, 3, 3)
                geom = np.array(block["controls"], dtype=float).reshape(3, 3, 3, 4)
                z = geom[0, 0, :, 2] / geom[0, 0, :, 3]
                l, r = z[0], z[-1]
                c1 = solution(l) + dz/2 * (solution.deriv()(l) - b*dz**2/2)
                coeff = np.array([solution(l), c1, solution(r)])
                max_err = max(max_err, abs(disp[dofs, 2]-coeff).max(), abs(disp[dofs, :2]).max())
            expected_energy = exact_energy - young*area*b*b*height*dz**4/40
            top_work = area*traction*solution(height)
            expected_body_work = 2*expected_energy-top_work
            check = dict(coefficient=max_err/max(np.linalg.norm(solution(height)), 1e-12),
                         energy=abs(res["strain_energy_n_mm"]/expected_energy-1),
                         body_work=abs(res["body_load"]["work_n_mm"]-expected_body_work)/abs(top_work),
                         surface_work=abs(res["traction_work_n_mm"]/top_work-1),
                         force=np.linalg.norm(np.array(res["body_load"]["resultant_n"])-[0, 0, force_z])/abs(force_z),
                         moment=np.linalg.norm(np.array(res["body_load"]["moment_n_mm"])-np.array([4.5, -2, 0])*force_z)/(5*abs(force_z)),
                         equilibrium=max(res["relative_moment_balance"], res["relative_force_balance"], res["relative_work_energy_error"]))
            rows.append(dict(case=label, h=h, cells=len(res["analysis_blocks"]), errors=check,
                             exact_strain_energy=exact_energy, discrete_strain_energy=res["strain_energy_n_mm"],
                             strain_energy_deficit=exact_energy-res["strain_energy_n_mm"], passed=bool(max(check.values())<=1e-8)))
    defects = [r["strain_energy_deficit"] for r in rows if r["case"]=="cubic"]
    ratios = [x/y for x,y in zip(defects, defects[1:])]
    job["plan"] = plan(roots, 0)
    job["body_load"]["gradient_n_per_mm4"] = [[.01, -.03, .04],[-.02, .01, 0],[.01, -.02, -.12]]
    B, offset = [[0,0,1],[0,1,0],[0,0,0]], [1,0,2]
    jp = folder/"derivative-job.json"
    write(jp, job)
    direction = folder/"direction.json"
    write(direction, dict(schema="solidcraft_affine_velocity_v1", linear=B, offset_mm=offset))
    def perturbed(s):
        rp, mp, op = [folder/f"perturbed.{suffix}" for suffix in ("recipe","model","result")]
        write(rp, image(source, 1, B, offset, s))
        solve_command(cli, "build", rp, mp)
        solve_command(cli, "solve-job", mp, jp, op)
        return json.loads(op.read_text()), mp
    base, mp = perturbed(Fraction(0))
    op = folder/"derivative.json"
    solve_command(cli, "differentiate-job", mp, jp, direction, op)
    dr = json.loads(op.read_text())
    fd = []
    previous = None
    previous_fourth = None
    for denom in (8, 16, 32):
        step = Fraction(1, denom)
        plus, _ = perturbed(step)
        minus, _ = perturbed(-step)
        c = (plus["external_work_n_mm"]-minus["external_work_n_mm"])/(2*float(step))
        f = (np.array(plus["force_n_xyz"])-minus["force_n_xyz"])/(2*float(step))
        grad = dr["compliance_rate_tangent_n_mm"]
        fd.append(dict(step=float(step), compliance=abs(c-grad)/abs(grad),
                       force=np.linalg.norm(f-dr["nodal_force_rate_n_xyz"])/np.linalg.norm(f)))
        if previous is not None:
            c4 = (4*c-previous[0])/3
            f4 = (4*f-previous[1])/3
            fd[-1].update(fourth_order_compliance=abs(c4-grad)/abs(grad),
                          fourth_order_force=np.linalg.norm(f4-dr["nodal_force_rate_n_xyz"])/np.linalg.norm(f4))
            if previous_fourth is not None:
                c6 = (16*c4-previous_fourth[0])/15
                f6 = (16*f4-previous_fourth[1])/15
                fd[-1].update(sixth_order_compliance=abs(c6-grad)/abs(grad),
                              sixth_order_force=np.linalg.norm(f6-dr["nodal_force_rate_n_xyz"])/np.linalg.norm(f6))
            previous_fourth = c4, f4
        previous = c, f
    adjoint = abs(dr["compliance_rate_adjoint_n_mm"]/dr["compliance_rate_tangent_n_mm"]-1)
    rejected = []
    for case in ["old_version", "missing_gradient", "two_densities", "short_vector"]:
        bad = copy.deepcopy(job)
        if case == "old_version": bad["schema"] = "solidcraft_mechanical_job_v3"
        if case == "missing_gradient": del bad["body_load"]["gradient_n_per_mm4"]
        if case == "two_densities": bad["body_load"] = [bad["body_load"]]*2
        if case == "short_vector": bad["body_load"]["constant_n_per_mm3"] = [0, 1]
        write(jp, bad)
        op.write_bytes(b"previous accepted result")
        cmd = subprocess.run([str(cli), "solve-job", str(model), str(jp), str(op)], capture_output=True, text=True)
        rejected.append(dict(case=case, refused=cmd.returncode!=0, preserved=op.read_bytes()==b"previous accepted result", diagnostic=cmd.stdout.strip()))
    return dict(schema="affine_volume_force_oracle_v1", bar=rows, energy_deficit_ratios=ratios,
                centered_derivative=fd, tangent_adjoint=adjoint, refusals=rejected,
                passed=bool(all(r["passed"] for r in rows) and all(15.9<x<16.1 for x in ratios)
                            and max(fd[-1]["sixth_order_compliance"],fd[-1]["sixth_order_force"],adjoint)<1e-5
                            and all(r["refused"] and r["preserved"] for r in rejected)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cli","work","output"): parser.add_argument(name, type=Path)
    args = parser.parse_args()
    report = run(args.cli.resolve(), args.work)
    write(args.output, report)
    print("Body force oracle", report["passed"])
    raise SystemExit(int(not report["passed"]))
