"""Measure reduced stiffness nullity on the locally refined sector at q2/q3/q4."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import hstack, kron

from mechanical_job import assemble, projection_from_lifts
from sector_job_sensitivity import call, job, write


def run(cli, work):
    work.mkdir(exist_ok=True, parents=True)
    write(work/"recipe.json", dict(family="rational_sector_prism_v1",radii_mm=[2,6],
                                  half_angle_tangent=[1,2],height_mm=8))
    call(cli,"build",work/"recipe.json",work/"model.json")
    rows=[]
    for q in (2,3,4):
        j=job(dict(select=dict(feature="inner"),pressure_mpa=.3),q)
        j["material"]=dict(young_mpa=210000,poisson=.3)
        j["plan"].update(max_depth=2,max_cells=128,steps=[
            dict(refine=[dict(root="whole",level=0,index=[0,0,0])]),
            dict(refine=[dict(root="whole",level=1,index=[0,0,0])])])
        write(work/"job.json",j)
        call(cli,"solve-job",work/"model.json",work/"job.json",work/"result.json")
        sol=json.loads((work/"result.json").read_text())
        k,n=assemble(sol["analysis_blocks"],210000,.3,q)
        lifts=projection_from_lifts(sol,n)
        c=hstack([kron(lift,np.eye(3)[:,i,None],format="csr") for i,lift in enumerate(lifts)],format="csr")
        reduced=(c.T@k@c).toarray()
        vals,modes=eigh(reduced,driver="evd")
        cutoff=1e-10*vals[-1]
        null=np.flatnonzero(np.abs(vals)<=cutoff)
        rhs=c.T@np.array(sol["force_n_xyz"])
        force_fraction=float(np.linalg.norm(modes[:,null].T@rhs)/np.linalg.norm(rhs))
        # The unreduced local matrices admit rigid fields on every physical control.
        xyz=np.empty((n,3))
        for b in sol["analysis_blocks"]:
            h=np.array(b["controls"],dtype=float)
            xyz[b["dofs"]]=h[:,:3]/h[:,3,None]
        rigid=[np.tile(v,(n,1)) for v in np.eye(3)]+[np.cross(v,xyz-xyz.mean(axis=0)) for v in np.eye(3)]
        residuals=[float(np.linalg.norm(k@v.ravel())/(vals[-1]*np.linalg.norm(v))) for v in rigid]
        rows.append(dict(order=q,scalar_dofs=n,mounted_dimension=len(vals),nullity=len(null),
            lowest_eigenvalues=vals[:16].tolist(),largest_eigenvalue=vals[-1],threshold=cutoff,
            force_fraction_in_kernel=force_fraction,free_rigid_relative_residuals=residuals,
            native_relative_residual=sol["relative_residual"],negative_below_tolerance=int(sum(vals < -cutoff))))
        print(q, len(null), vals[0], force_fraction, flush=True)
    return dict(schema="sector_quadrature_modes_v1",recipe=json.loads((work/"recipe.json").read_text()),
                job=j,cases=rows,executable_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),
                qualified=rows[0]["nullity"]>0 and all(r["nullity"]==0 for r in rows[1:])
                          and max(max(r["free_rigid_relative_residuals"]) for r in rows)<1e-10)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("cli","work","output"):p.add_argument(name,type=Path)
    a=p.parse_args();data=run(a.cli.resolve(),a.work);write(a.output,data)
    raise SystemExit(int(not data["qualified"]))
