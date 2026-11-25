"""Galerkin multigrid preconditioner for the Cartesian quadratic reference."""
import math
import time
import numpy as np
from scipy.interpolate import BSpline
from scipy.sparse import csr_matrix, diags, eye, kron
from scipy.sparse.linalg import LinearOperator, cg, splu


def axis_transfer(fine,coarse,length):
    if fine%coarse or coarse<2:
        raise ValueError("non-nested knot vectors")
    coarse_nodes=np.linspace(0,length,coarse+1)
    knots=np.r_[[0.,0.],coarse_nodes,[length,length]]
    spline=BSpline(knots,np.eye(coarse+2),2,extrapolate=False)
    refinement=fine//coarse
    for j in range(1,fine):
        if j%refinement:
            spline=spline.insert_knot(j*length/fine)
    p=csr_matrix(spline.c)
    if np.max(np.abs(np.squeeze(np.asarray(p.sum(axis=1)))-1.))>1e-13:
        raise ArithmeticError("B-spline insertion lost partition of unity")
    return p


def coarsen(ncell):
    # Unequal integer factors are allowed when each axis remains nested.
    return tuple(next((n//f for f in (2,3,5) if n%f==0 and n//f>=2),n) for n in ncell)


def galerkin_preconditioner(matrix,scaling,free,ncell):
    started=time.perf_counter()
    levels=[]
    audit=[]
    a=matrix.tocsr()
    while len(free)>3500:
        coarse=coarsen(ncell)
        if coarse==ncell:
            raise ValueError("no nested coarser background grid")
        # Keep P scalar until the active rows and columns are known.
        axes=[axis_transfer(n,m,l) for n,m,l in zip(ncell,coarse,(60.,50.,40.))]
        scalar=kron(kron(axes[0],axes[1],format="csr"),axes[2],format="csr")
        if np.any(free.reshape(-1,3)%3!=np.arange(3)):
            raise ValueError("expected complete Cartesian displacement triples")
        small=scalar[free[::3]//3]
        active=np.unique(small.indices)
        active=active[(active//(coarse[2]+2))%(coarse[1]+2)!=coarse[1]+1]
        coarse_free=(3*active[:,None]+np.arange(3)).ravel()
        p=kron(small[:,active],eye(3,format="csr"),format="csr")
        unscaled=diags(1./scaling)@p
        coarse_matrix=(unscaled.T@(a@unscaled)).tocsr()
        coarse_d=coarse_matrix.diagonal()
        if np.any(coarse_d<=0):
            raise ArithmeticError("nonpositive Galerkin coarse diagonal")
        coarse_scaling=1./np.sqrt(coarse_d)
        prolong=(unscaled@diags(coarse_scaling)).tocsr()
        ac=(diags(coarse_scaling)@coarse_matrix@diags(coarse_scaling)).tocsr()
        test=np.random.default_rng(20261125).normal(size=ac.shape[0])
        fine_test=prolong@test
        energy=float(fine_test@(a@fine_test))
        energy_error=abs(float(test@(ac@test))-energy)/abs(energy)
        if energy_error>1e-11:
            raise ArithmeticError("inconsistent Galerkin energy")
        # a has unit diagonal. This Gershgorin bound makes each damped
        # Richardson factor contractive and each symmetric cycle positive.
        bound=float(np.asarray(abs(a).sum(axis=1)).max())*(1+1e-12)
        levels.append((a,prolong,prolong.T.tocsr(),1.8/bound))
        audit.append({"fine_ncell":ncell,"coarse_ncell":coarse,"fine_dofs":a.shape[0],
                      "coarse_dofs":ac.shape[0],"energy_relative":energy_error,
                      "spectral_upper":bound,"P_nnz":prolong.nnz})
        a=ac
        scaling=coarse_scaling
        ncell=coarse
        free=coarse_free
    lu=splu(a.tocsc(),permc_spec="MMD_AT_PLUS_A")
    def apply(level,b):
        if level==len(levels):
            return lu.solve(b)
        system,p,pt,damping=levels[level]
        x=damping*b
        x+=damping*(b-system@x)
        x+=p@apply(level+1,pt@(b-system@x))
        x+=damping*(b-system@x)
        x+=damping*(b-system@x)
        return x
    pre=LinearOperator(matrix.shape,matvec=lambda b:apply(0,b),dtype=float)
    rng=np.random.default_rng(20261125)
    x,y=rng.normal(size=(2,matrix.shape[0]))
    mx,my=pre@x,pre@y
    symmetry=abs(float(x@my-y@mx))/(np.linalg.norm(x)*np.linalg.norm(my)+np.linalg.norm(y)*np.linalg.norm(mx))
    if symmetry>1e-10 or min(x@mx,y@my)<=0.:
        raise ArithmeticError("preconditioner symmetry or energy check failed")
    return pre,{"levels":audit,"coarse_lu_size":a.shape[0],"symmetry_defect":symmetry,
                "coarse_lu_nnz":int(lu.nnz),"setup_s":time.perf_counter()-started}


def multigrid_solve(matrix,scaling,rhs,free,ncell):
    pre,info=galerkin_preconditioner(matrix,scaling,free,ncell)
    start=time.perf_counter()
    count=0
    def observe(_):
        nonlocal count
        count+=1
    y,status=cg(matrix,rhs*scaling,M=pre,rtol=1e-11,atol=0.,maxiter=3000,callback=observe)
    info.update({"iterations":count,"status":int(status),"iterations_s":time.perf_counter()-start})
    if status:
        raise ArithmeticError(f"multigrid PCG did not converge: {info}")
    return scaling*y,info
