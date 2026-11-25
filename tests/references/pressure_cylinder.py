"""Independent plane-strain boundary-value solution and refinement campaign."""
import argparse
import json
import math
import os
from pathlib import Path
import numpy as np

from refine_tetra_fields import measure


def reference(nu, pi, po):
    # Axisymmetric equilibrium admits u=C*r+D/r. Determine both constants
    # from the prescribed normal stresses using the full isotropic 3D law.
    mu=210000/(2*(1+nu)); lam=2*mu*nu/(1-2*nu)
    def stress(c,d,r):
        strain=np.diag([c-d/r**2,c+d/r**2,0.])
        return 2*mu*strain+lam*strain.trace()*np.eye(3)
    system=[[stress(c,d,r)[0,0] for c,d in ((1,0),(0,1))] for r in (2,4)]
    c,d=np.linalg.solve(system,[-pi,-po])
    def disp(r):return c*r+d/r
    # Independent volume integral of half the strain-stress inner product.
    x,w=np.polynomial.legendre.leggauss(40)
    energy=0.
    for t,weight in zip(x,w):
        r=3+t
        strain=np.diag([c-d/r**2,c+d/r**2,0.])
        energy+=weight*math.pi*8*r*np.sum(strain*stress(c,d,r))
    return dict(inner_displacement_mm=disp(2),outer_displacement_mm=disp(4),
                inner_hoop_mpa=stress(c,d,2)[1,1],outer_hoop_mpa=stress(c,d,4)[1,1],
                strain_energy_n_mm=energy,top_reaction_n=12*math.pi*stress(c,d,3)[2,2])


def checks(r,ref,nu,pi,po):
    same=pi==po
    def relative(a,b):return abs(a-b)/max(abs(b),1e-16)
    exact_diff=max(relative(r["exact_energy_n_mm"],ref["strain_energy_n_mm"]),
                   relative(r["exact_top_reaction_n"],ref["top_reaction_n"]),
                   *(relative(r[s]["exact_radial_displacement_mm"],ref[s+"_displacement_mm"]) for s in ("inner","outer")),
                   *(relative(r[s]["exact_hoop_mpa"],ref[s+"_hoop_mpa"]) for s in ("inner","outer")))
    pressure_area=2*math.pi*8*(2*pi+4*po)
    reactions=r["axial_reactions_n"]
    return dict(analytic_agreement=exact_diff,pressure_force_relative=math.dist(r["pressure_resultant_n"],[0]*3)/pressure_area,
                work_energy_relative=relative(r["work_n_mm"],2*r["strain_energy_n_mm"]),
                axial_reaction_relative=max(relative(reactions[0],-ref["top_reaction_n"]),relative(reactions[1],ref["top_reaction_n"])),
                affine=same,field_quadrature_delta=max(abs(r[key]-r["check_q8_"+key]) for key in ("relative_l2","relative_energy_norm")),
                threshold_satisfied=r["relative_energy_norm"]<=(1e-7 if same else .002)
                  and (r["relative_l2"]<=1e-7 if same else True)
                  and max(r[s]["max_sampled_hoop_relative"] for s in ("inner","outer"))<=.01)


def run(binary,work,out,projected):
    work.mkdir(exist_ok=True,parents=True)
    rows=[]
    for nu,pi,po in ((.3,100,0),(.3,100,30),(.3,40,40),(.49,100,0),(.4999,100,0)):
        result={"nu":nu,"pressures_mpa":[pi,po],"oracle":reference(nu,pi,po),"runs":[]}
        rows.append(result)
        specs=[(h,8,False) for h in range(5)]+[(1,8,True),(4,12,False)]
        for h,q,alternate in specs:
            stem=f"nu{nu}-i{pi}-o{po}-h{h}-q{q}-a{int(alternate)}"
            path=work/(stem+".json")
            if path.exists():
                result["runs"].append(json.loads(path.read_text()));continue
            log=work/(stem+".log")
            command=[str(binary),str(h),str(q),str(nu),str(pi),str(po),str(alternate).lower()]
            if projected:command.append("--projected")
            record=measure(command,os.environ.copy(),log)
            record.update(level=h,order=q,alternate_gauge=alternate,command=command)
            if record["status"]=="completed":
                solved=json.loads(log.read_text())
                record.update(result=solved,checks=checks(solved,result["oracle"],nu,pi,po))
            else:
                record["diagnostic"]=log.read_text()[-1000:]
            path.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
            result["runs"].append(record)
            print(stem,record["status"],record.get("checks",{}).get("threshold_satisfied"),flush=True)
        successful=[r for r in result["runs"] if r["status"]=="completed"]
        level={r["level"]:r["result"] for r in successful if r["order"]==8 and not r["alternate_gauge"]}
        result["convergence"]=[{"level":h,"l2_slope":math.log2(level[h-1]["relative_l2"]/r["relative_l2"]),
            "energy_slope":math.log2(level[h-1]["relative_energy_norm"]/r["relative_energy_norm"])}
            for h,r in level.items() if h-1 in level and pi!=po]
        extra=[r["result"] for r in successful if r["alternate_gauge"]]
        if extra and 1 in level:
            result["gauge_comparison"]={k:abs(extra[0][k]-level[1][k]) for k in ("relative_l2","relative_energy_norm","strain_energy_n_mm")}
        fine=[r for r in successful if r["level"]==4 and r["order"]==12]
        result["qualified"]=bool(fine and len(successful)==len(specs) and fine[0]["checks"]["threshold_satisfied"]
            and all(r["checks"]["analytic_agreement"]<1e-10 and r["checks"]["work_energy_relative"]<1e-8
            and r["checks"]["pressure_force_relative"]<1e-9 and r["checks"]["axial_reaction_relative"]<1e-6
            and r["result"]["relative_residual"]<=1e-10 for r in successful))
        out.write_text(json.dumps({"schema":"lame_cylinder_pressure_v1","projected_bulk":projected,"cases":rows},indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("binary","work","output"):p.add_argument(key,type=Path)
    p.add_argument("--projected",action="store_true")
    a=p.parse_args();run(a.binary.resolve(),a.work,a.output,a.projected)
