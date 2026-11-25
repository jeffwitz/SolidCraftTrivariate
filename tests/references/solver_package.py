"""Independently assemble and solve an exported conforming bracket package."""
import argparse
import json
import os
from pathlib import Path
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy import linalg, optimize, sparse
from scipy.sparse.linalg import spsolve
from immersed_bracket import local_stiffness


def shape(points, controls):
    """Compute rational basis, parametric derivatives and physical Jacobians."""
    t = np.asarray(points, dtype=float).reshape(-1, 3)
    values = np.stack(((1-t)**2, 2*t*(1-t), t*t), axis=2)
    derivatives = np.stack((2*(t-1), 2-4*t, 2*t), axis=2)
    def tensor(a,b,c):
        return np.einsum("qi,qj,qk->qijk", a,b,c).reshape(-1,27)
    pol = tensor(*values.transpose(1,0,2))
    dpol = []
    for axis in range(3):
        dpol.append(tensor(*(derivatives[:,d] if d == axis else values[:,d] for d in range(3))))
    dpol = np.stack(dpol, axis=2)
    h = np.asarray(controls,dtype=float)
    h = h / np.max(h[:,3])
    weights = h[:,3]
    denominator = pol @ weights
    rat = pol*weights/denominator[:,None]
    dw = np.einsum("qia,i->qa",dpol,weights)
    drat = (dpol*weights[None,:,None]-rat[:,:,None]*dw[:,None,:])/denominator[:,None,None]
    cartesian = h[:,:3]/weights[:,None]
    position = rat @ cartesian
    jac = np.einsum("ik,qia->qka",cartesian,drat)
    grad = np.einsum("qia,qak->qik",drat,np.linalg.inv(jac))
    return rat, position, jac, grad


def face_rule(side, uv):
    axis = "uvw".index(side[0])
    axes = [d for d in range(3) if d != axis]
    uvw = np.empty((len(uv),3))
    uvw[:,axis] = float(side[1])
    uvw[:,axes] = uv
    return uvw,axes


def assemble_operators(package):
    """Reassemble the full unconstrained operators from the exported block data."""
    space = package["analysis_space"]
    rows=[];cols=[];entries=[]
    ndof=space["scalar_ndof"]
    coordinates=np.full((ndof,3),np.nan)
    force=np.zeros(3*ndof)
    order=package["quadrature_order"]
    points,weights=np.polynomial.legendre.leggauss(order)
    points=(points+1)/2
    weights=weights/2
    rule=np.asarray(package["quadrature"]["nodes_weights"])
    assert np.max(np.abs(rule-np.column_stack((points,weights)))) < 1e-13
    grid = np.array(np.meshgrid(points,points,points,indexing="ij")).reshape(3,-1).T
    w3 = np.einsum("i,j,k->ijk",weights,weights,weights).ravel()
    facegrid = np.array(np.meshgrid(points,points,indexing="ij")).reshape(2,-1).T
    w2 = np.outer(weights,weights).ravel()
    young,nu=package["material"]["young_mpa"],package["material"]["poisson"]
    mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu))
    load_faces={}
    for face in package["design_faces"]:
        if face["feature"]=="x60":
            load_faces.setdefault(face["block"],[]).append(face["side"])
    volume=0.;area=0.;affine_error=0.;tip_shapes=[]
    for block in space["blocks"]:
        controls=np.asarray(block["controls"],dtype=float)
        ids=np.asarray(block["dofs"])
        local_nodes=controls[:,:3]/controls[:,3,None]
        mask=np.isfinite(coordinates[ids,0])
        if np.any(mask):
            assert np.max(np.abs(coordinates[ids[mask]]-local_nodes[mask]))<1e-11
        coordinates[ids]=local_nodes
        r,x,j,grad=shape(grid,controls)
        det=np.linalg.det(j)
        assert min(det)>0.
        measured=det*w3
        volume+=sum(measured)
        affine_error=max(affine_error,float(np.max(np.abs(np.einsum("ik,qij->qkj",local_nodes,grad)-np.eye(3)))))
        k=local_stiffness(grad,measured,mu,lam)
        dofs=(3*ids[:,None]+np.arange(3)).ravel()
        rows.append(np.repeat(dofs,81));cols.append(np.tile(dofs,81));entries.append(k.ravel())
        for side in load_faces.get(block["id"],[]):
            q,axes=face_rule(side,facegrid)
            rr,_,jj,_=shape(q,controls)
            dw=w2*np.linalg.norm(np.cross(jj[:,:,axes[0]],jj[:,:,axes[1]]),axis=1)
            area+=sum(dw)
            # Reconstruct total load as a constant traction on the design face.
            force[3*ids+1]+=(-1000./320.)*(rr.T@dw)
            target=np.array([60.,4.,20.])
            def residual(t):
                uvw,_=face_rule(side,[t])
                return shape(uvw,controls)[1][0]-target
            solution=optimize.least_squares(residual,[.5,.5],xtol=1e-14,ftol=1e-14,gtol=1e-14)
            if np.linalg.norm(solution.fun)<1e-9 and np.all(solution.x>=-1e-12) and np.all(solution.x<=1+1e-12):
                uvw,_=face_rule(side,[solution.x])
                tip_shapes.append((ids,shape(uvw,controls)[0][0]))
    assert tip_shapes and np.all(np.isfinite(coordinates))
    k=sparse.coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(cols))),
                        shape=(3*ndof,3*ndof)).tocsr()
    rhs_error=np.max(np.abs(force-np.asarray(package["load"]["nodal_rhs_n"])))
    assert affine_error < 1e-10 and rhs_error < 1e-8 and abs(area-320)<1e-9
    return k, force, coordinates, tip_shapes, volume, area, affine_error, rhs_error


def assemble(package):
    k,force,coordinates,tip_shapes,volume,area,affine_error,rhs_error=assemble_operators(package)
    ndof=len(coordinates)
    order=package["quadrature_order"]
    fixed=np.asarray(package["reduction"]["fixed_vector_dofs"])
    free=np.asarray(package["reduction"]["free_vector_dofs"])
    assert np.array_equal(np.sort(np.r_[fixed,free]),np.arange(3*ndof))
    y50=np.flatnonzero(np.abs(coordinates[:,1]-50)<1e-9)
    assert np.array_equal(np.sort(fixed),(3*y50[:,None]+np.arange(3)).ravel())
    rotation=[np.cross(axis,coordinates).ravel() for axis in np.eye(3)]
    translation=[np.tile(axis,ndof) for axis in np.eye(3)]
    rigid=np.column_stack(translation+rotation)
    residual=float(np.linalg.norm(k@rigid)/(np.linalg.norm(k.data)*np.linalg.norm(rigid)))
    assert residual < 1e-12 and np.linalg.matrix_rank(rigid)==6
    eig=linalg.eigvalsh(k.toarray(),driver="evd",check_finite=False)
    zero=int(np.sum(np.abs(eig)<eig[-1]*1e-10))
    assert zero==6 and eig[0] > -eig[-1]*1e-10
    u=np.zeros(len(force));u[free]=spsolve(k[free][:,free],force[free])
    tips=[rr @ u[3*ids+1] for ids,rr in tip_shapes]
    assert np.ptp(tips)<1e-10
    reaction=(k@u-force)[fixed].reshape(-1,3).sum(axis=0)
    return {"scalar_dofs":ndof,"solver":"SciPy sparse direct","volume_mm3":volume,"area_mm2":area,
            "affine_gradient_max_error":affine_error,"rigid_action_relative":residual,
            "free_body_zero_modes":zero,"free_body_eigenvalue_7":float(eig[6]),"rhs_max_error_n":float(rhs_error),
            "point_tip_y_mm":float(tips[0]),"mean_tip_y_mm":float(force@u/-1000.),
            "strain_energy_n_mm":float(u@(k@u)/2),"reaction_n":reaction.tolist(),"q":order}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("package",type=Path)
    p.add_argument("--native",type=Path,required=True)
    args=p.parse_args()
    start=time.monotonic()
    result=assemble(json.loads(args.package.read_text()))
    native=json.loads(args.native.read_text())
    keys=("volume_mm3","point_tip_y_mm","mean_tip_y_mm","strain_energy_n_mm")
    comparison={k:abs(result[k]/native[k]-1.) for k in keys}
    assert max(comparison.values())<=1e-7
    result.update(native_relative_error=comparison,reader_seconds=time.monotonic()-start)
    print(json.dumps(result,indent=2,allow_nan=False))
