"""Independent Cartesian C1 quadratic B-spline elasticity of the filleted bracket.

The material domain is integrated by analytic sections through a Cartesian
background grid. No fitted geometry, mesh or operator is imported.
"""
import argparse
import json
import math
import resource
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu, LinearOperator, onenormest


def splits(lo, hi, values):
    points=sorted({lo,hi,*(x for x in values if lo<x<hi)})
    return zip(points,points[1:])


def interval(lo,hi,nodes,weights):
    return lo+(hi-lo)*nodes, (hi-lo)*weights


def disk_cut_rule(ab,bc,center,radius,nodes,weights,below_only=False):
    """Integrate a rectangle outside a disk, or just below its lower semicircle."""
    a0,a1=ab
    b0,b1=bc
    ca,cb=center
    cuts=[ca-radius,ca,ca+radius]
    for b in (b0,b1):
        h=radius*radius-(b-cb)**2
        if h>=0:
            cuts.extend([ca-math.sqrt(h),ca+math.sqrt(h)])
    arrays=[]
    for l,r in splits(a0,a1,cuts):
        if l==ca-radius:
            x=l+(r-l)*nodes**2
            wx=2*(r-l)*nodes*weights
        elif r==ca+radius:
            x=r-(r-l)*(1-nodes)**2
            wx=2*(r-l)*(1-nodes)*weights
        else:
            x,wx=interval(l,r,nodes,weights)
        for a,weight in zip(x,wx):
            h=radius*radius-(a-ca)**2
            if h<=0:
                runs=[] if below_only else [(b0,b1)]
            else:
                t=math.sqrt(h)
                runs=[(b0,min(b1,cb-t))]
                if not below_only:
                    runs.append((max(b0,cb+t),b1))
            for low,high in runs:
                if low<high:
                    b,wb=interval(low,high,nodes,weights)
                    arrays.append((np.full(len(nodes),a),b,weight*wb))
    if not arrays:
        return (np.array([]),)*3
    return tuple(np.concatenate(v) for v in zip(*arrays))


def tensor_three(ranges,nodes,weights):
    t=[interval(lo,hi,nodes,weights) for lo,hi in ranges]
    xyz=np.meshgrid(*(v[0] for v in t),indexing="ij")
    w=t[0][1][:,None,None]*t[1][1][None,:,None]*t[2][1][None,None,:]
    return tuple(v.ravel() for v in xyz)+(w.ravel(),)


def product_with_axis(rule,t,axis):
    a,b,wa=rule
    c,wc=t
    if not len(a):
        return None
    one=np.repeat(a,len(c))
    two=np.repeat(b,len(c))
    three=np.tile(c,len(a))
    # axis is the independent direction; (a,b) is ordered y,z, x,z, or x,y.
    coords={0:(three,one,two),1:(one,three,two),2:(one,two,three)}[axis]
    return (*coords,(wa[:,None]*wc).ravel())


def material(bounds,nodes,weights,radius_x=4.,radius_y=4.):
    (x0,x1),(y0,y1),(z0,z1)=bounds
    parts=[]
    for xl,xh in splits(x0,x1,[8.,14.]):
        if xh<=8.:
            rule=disk_cut_rule((y0,y1),(z0,z1),(32.,20.),radius_y,nodes,weights)
            part=product_with_axis(rule,interval(xl,xh,nodes,weights),0)
        elif xl>=14.:
            if y0>=8.:
                continue
            upper=min(y1,8.)
            rule=disk_cut_rule((xl,xh),(z0,z1),(40.,20.),radius_x,nodes,weights)
            part=product_with_axis(rule,interval(y0,upper,nodes,weights),1)
        else:
            if y0>=14.:
                continue
            rule=disk_cut_rule((xl,xh),(y0,min(14.,y1)),(14.,14.),6.,nodes,weights,True)
            part=product_with_axis(rule,interval(z0,z1,nodes,weights),2)
        if part is not None:
            parts.append(part)
    if not parts:
        return (np.array([]),)*4
    return tuple(np.concatenate(a) for a in zip(*parts))

class Space:
    def __init__(self,h):
        self.breaks=[]
        self.splines=[]
        for length in (60,50,40):
            count=round(length/h)
            if count<2 or abs(count*h-length)>1e-12:
                raise ValueError("h must divide all dimensions")
            values=np.linspace(0,length,count+1)
            knots=np.r_[[0.,0.],values,[length,length]]
            self.breaks.append(values)
            self.splines.append(BSpline(knots,np.eye(count+2),2,extrapolate=False))
        self.shape=tuple(len(b)+1 for b in self.breaks)
        self.total=int(np.prod(self.shape))

    def basis(self,cell,xyz):
        values=[]
        gradients=[]
        for i,x,spline in zip(cell,xyz,self.splines):
            values.append(spline(x)[:,i:i+3])
            gradients.append(spline(x,nu=1)[:,i:i+3])
        def tensor(a,b,c):
            return (a[:,:,None,None]*b[:,None,:,None]*c[:,None,None,:]).reshape(-1,27)
        n=tensor(*values)
        g=np.empty((len(xyz[0]),27,3))
        for k in range(3):
            n1=list(values)
            n1[k]=gradients[k]
            g[:,:,k]=tensor(*n1)
        return n,g

    def scalar_ids(self,cell):
        i,j,k=cell
        n1,n2=self.shape[1:]
        return (((np.arange(3)+i)[:,None,None]*n1+(np.arange(3)+j)[None,:,None])*n2
                +(np.arange(3)+k)[None,None,:]).ravel()

    def point(self,xyz):
        cell=tuple(min(len(b)-2,max(0,int(np.searchsorted(b,x,side="right")-1)))
                   for b,x in zip(self.breaks,xyz))
        values,grad=self.basis(cell,[np.array([x]) for x in xyz])
        return self.scalar_ids(cell),values[0],grad[0]


def local_stiffness(gradient,weights,mu,lam):
    grad=gradient.reshape(-1,81)
    gram=(grad.T@(weights[:,None]*grad)).reshape(27,3,27,3)
    matrix=lam*gram+mu*gram.transpose(0,3,2,1)
    trace=mu*np.trace(gram,axis1=1,axis2=3)
    for d in range(3):
        matrix[:,d,:,d]+=trace
    return matrix.reshape(81,81)


def full_cell(bounds,radius_x=4.,radius_y=4.):
    """Identify canonical full cells; remaining cells use exact analytic clipping."""
    (xl,xh),(yl,yh),(zl,zh)=bounds
    def outside_disk(l,h,center):
        return max(l-center,0.,center-h)
    dz=outside_disk(zl,zh,20.)
    if xh<=8.:
        return outside_disk(yl,yh,32.)**2+dz**2>=radius_y**2
    if xl>=14.:
        return yh<=8. and outside_disk(xl,xh,40.)**2+dz**2>=radius_x**2
    if xl>=8. and xh<=14. and yh<=14.:
        return (14.-xh)**2+(14.-yh)**2>=36.
    return False


def run(h,order,geometry_only=False,solver="lu",radius_x=4.,radius_y=4.,field_output=None):
    start=time.perf_counter()
    young,nu=210000.,.3
    mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu))
    space=Space(h)
    nodes,weights=leggauss(order)
    nodes,weights=(nodes+1)/2,weights/2
    ncell=tuple(len(b)-1 for b in space.breaks)
    empty=cut=active=qpoints=0
    volume=area=0.
    min_fraction=1.
    integrate_s=assemble_s=0.
    shape=(3*space.total,3*space.total)
    stiffness=None
    rows=[];cols=[];entries=[]
    rhs=np.zeros(shape[0])
    full_matrices={}
    def flush():
        nonlocal stiffness
        if not entries:
            return
        partial=coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(cols))),shape=shape).tocsr()
        stiffness=partial if stiffness is None else stiffness+partial
        rows.clear();cols.clear();entries.clear()

    for cell in np.ndindex(*ncell):
        bounds=[(b[i],b[i+1]) for b,i in zip(space.breaks,cell)]
        (xl,xh),(yl,yh),(zl,zh)=bounds
        if (xl>=14. and yl>=8.) or (xl>=8. and yl>=14.):
            empty+=1
            continue
        begin=time.perf_counter()
        full=full_cell(bounds,radius_x,radius_y)
        if full:
            vol=h**3
            key=tuple(-1 if c==0 else 1 if c==count-1 else 0 for c,count in zip(cell,ncell))
            # The open quadratic basis differs from its translated cardinal
            # form on the first/last span only.
            need=not geometry_only and key not in full_matrices
            if need:
                x,y,z,w=tensor_three(bounds,nodes,weights)
            qpoints+=order**3
        else:
            x,y,z,w=material(bounds,nodes,weights,radius_x,radius_y)
            vol=float(w.sum())
            if vol<=0:
                empty+=1
                integrate_s+=time.perf_counter()-begin
                continue
            cut+=int(vol < h**3*(1-1e-12))
            qpoints+=len(w)
        volume+=vol
        min_fraction=min(min_fraction,vol/h**3)
        active+=1
        integrate_s+=time.perf_counter()-begin
        if geometry_only:
            continue
        begin=time.perf_counter()
        if full and key in full_matrices:
            ke=full_matrices[key]
        else:
            _,grad=space.basis(cell,[x,y,z])
            ke=local_stiffness(grad,w,mu,lam)
            if full:
                full_matrices[key]=ke
        ids=space.scalar_ids(cell)
        ids3=(3*ids[:,None]+np.arange(3)).ravel().astype(np.int32)
        rows.append(np.repeat(ids3,81));cols.append(np.tile(ids3,81));entries.append(ke.ravel())
        if len(entries)>=256:
            flush()
        if xh==60. and yl<8.:
            yp,wy=interval(yl,min(yh,8.),nodes,weights)
            zp,wz=interval(zl,zh,nodes,weights)
            yy,zz=np.meshgrid(yp,zp,indexing="ij")
            force_weights=(wy[:,None]*wz).ravel()
            n,_=space.basis(cell,[np.full(yy.size,60.),yy.ravel(),zz.ravel()])
            rhs[3*ids+1]-=1000./320.*(n.T@force_weights)
            area+=float(force_weights.sum())
        assemble_s+=time.perf_counter()-begin
    result={"h_mm":h,"gauss_order":order,"volume_mm3":volume,
            "radius_x_mm":radius_x,"radius_y_mm":radius_y,
            "volume_target_mm3":34080.-math.pi*(360+8*(radius_x**2+radius_y**2)),"active_cells":active,"cut_cells":cut,
            "empty_cells":empty,"minimum_volume_fraction":min_fraction,"quadrature_points":qpoints,
            "quadrature_s":integrate_s}
    if geometry_only:
        return result
    begin=time.perf_counter()
    flush()
    assemble_s+=time.perf_counter()-begin
    diag=stiffness.diagonal()
    if np.any(diag<0):
        raise ArithmeticError("negative stiffness diagonal")
    active_dof=diag>0
    if np.any(rhs[~active_dof]!=0):
        raise ArithmeticError("load on empty support")
    on_support=np.flatnonzero(active_dof)
    clamp=(on_support//3//space.shape[2])%space.shape[1]==space.shape[1]-1
    free=on_support[~clamp]
    prep=time.perf_counter()
    k=stiffness[free][:,free].tocsc()
    scaling=1/np.sqrt(diag[free])
    scaled=(diags(scaling)@k@diags(scaling)).tocsc()
    prep_s=time.perf_counter()-prep
    solve=time.perf_counter()
    u=np.zeros(shape[0])
    if solver=="multigrid":
        from immersed_multigrid import multigrid_solve
        u[free],solve_info=multigrid_solve(scaled,scaling,rhs[free],free,ncell)
    else:
        lu=splu(scaled,permc_spec="MMD_AT_PLUS_A")
        u[free]=scaling*lu.solve(rhs[free]*scaling)
        # Estimate the inverse 1-norm using the existing direct factor.
        estimate_start=time.perf_counter()
        inverse=LinearOperator(scaled.shape,matvec=lu.solve,rmatvec=lambda v:lu.solve(v,trans="T"))
        np.random.seed(20261125)
        condition=float(abs(scaled).sum(axis=0).max())*float(onenormest(inverse))
        solve_info={"estimated_scaled_condition_1":condition,"condition_estimate_s":time.perf_counter()-estimate_start,
                    "lu_factor_nnz":int(lu.L.nnz+lu.U.nnz),
                    "lu_abs_pivot_ratio_diagnostic":float(np.min(np.abs(lu.U.diagonal()))/np.max(np.abs(lu.U.diagonal())))}
    solve_s=time.perf_counter()-solve
    if field_output is not None:
        np.savez(field_output,h_mm=h,radii_mm=np.array([radius_x,radius_y]),displacement_mm=u.reshape((*space.shape,3)),
                 young_mpa=young,poisson=nu,gauss_order=order)
    residual=stiffness@u-rhs
    residual_norm=float(np.linalg.norm(residual[free])/np.linalg.norm(rhs[free]))
    force=rhs.reshape(-1,3).sum(axis=0)
    reaction=np.zeros_like(u)
    reaction[on_support[clamp]]=residual[on_support[clamp]]
    reaction=reaction.reshape(-1,3).sum(axis=0)
    ids,n,_=space.point([60.,4.,20.])
    tip=float(n@u.reshape(-1,3)[ids,1])
    coordinate=14.-6./math.sqrt(2)
    ids,_,grad=space.point([coordinate,coordinate,20.])
    du=u.reshape(-1,3)[ids].T@grad
    strain=(du+du.T)/2
    stress=2*mu*strain+lam*np.trace(strain)*np.eye(3)
    deviator=stress-np.trace(stress)*np.eye(3)/3.
    result.update({"scalar_active":int(len(on_support)//3),"free_components":int(len(free)),
                   "loaded_area_mm2":area,"force_n":force.tolist(),"reaction_n":reaction.tolist(),
                   "relative_residual":residual_norm,"point_tip_y_mm":tip,"mean_tip_y_mm":-float(rhs@u)/1000.,
                   "von_mises_mpa":float(np.sqrt(1.5*np.sum(deviator**2))),"stress_mpa":stress.tolist(),
                   "strain_energy_n_mm":float(u@(stiffness@u))/2.,"external_work_n_mm":float(rhs@u),
                   "assembly_s":assemble_s,"diagonal_scaling_and_restriction_s":prep_s,
                   "linear_solver":solver,"solver_detail":solve_info,"solver_s":solve_s,"elapsed_s":time.perf_counter()-start,
                   "reduced_nnz":int(k.nnz),
                   "scaling_dynamic_range":float(scaling.max()/scaling.min()),
                   "max_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h",type=float,required=True)
    parser.add_argument("--order",type=int,required=True)
    parser.add_argument("--geometry-only",action="store_true")
    parser.add_argument("--solver",choices=("lu","multigrid"),default="lu")
    parser.add_argument("--radius-x",type=float,default=4.)
    parser.add_argument("--radius-y",type=float,default=4.)
    parser.add_argument("--field-output",help="write the solved B-spline displacement coefficients to NPZ")
    args=parser.parse_args()
    if not 2<=args.order<=16 or not 1<=args.h<=5:
        parser.error("expected order 2..16, h in [1,5] mm")
    print(json.dumps(run(args.h,args.order,args.geometry_only,args.solver,args.radius_x,args.radius_y,args.field_output),indent=2))
