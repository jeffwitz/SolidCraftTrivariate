#!/usr/bin/env python3
"""Independent immersed C1 quadratic B-spline reference for the J4 Kirsch pilot.

The Cartesian background covers the symmetry quadrant. The circular void cuts
the integration domain but does not introduce knots or geometry-fitted basis
functions. The two symmetry conditions are the restrictions of even/odd
displacement fields on the full square. The bore carries zero traction.
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
from scipy.sparse.linalg import splu


def kirsch_tensor(x, y, a, tension):
    """Cartesian in-plane stresses from the traction-free infinite-hole solution."""
    r2 = x * x + y * y
    c, s = x / np.sqrt(r2), y / np.sqrt(r2)
    cos2, sin2 = c*c - s*s, 2*c*s
    d = a*a / r2
    rr = tension/2 * ((1-d) + (1-4*d+3*d*d)*cos2)
    tt = tension/2 * ((1+d) - (1+3*d*d)*cos2)
    rt = -tension/2 * (1+2*d-3*d*d)*sin2
    return (rr*c*c + tt*s*s - 2*rt*c*s,
            rr*s*s + tt*c*c + 2*rt*c*s,
            (rr-tt)*s*c + rt*cos2)


def material_rule(x0, x1, y0, y1, radius, nodes, weights):
    """Tensor integration on smooth portions above the analytic quarter circle."""
    breaks = {x0, x1}
    for y in (0., y0, y1):
        if y <= radius:
            x = math.sqrt(max(0., radius*radius-y*y))
            if x0 < x < x1:
                breaks.add(x)
    xs, ys, measures = [], [], []
    cut = x0*x0 + y0*y0 < radius*radius
    for lo, hi in zip(sorted(breaks), sorted(breaks)[1:]):
        if cut and hi == radius:
            # Smooth the square-root endpoint where the circle meets y=0.
            x = hi - (hi-lo)*(1.-nodes)**2
            wx = weights * (2*(hi-lo)*(1.-nodes))
        else:
            x = lo + nodes*(hi-lo)
            wx = weights*(hi-lo)
        low = np.maximum(y0, np.sqrt(np.maximum(0., radius*radius-x*x)))
        for xx, yy, ww in zip(x, low, wx):
            if yy < y1:
                xs.extend(np.repeat(xx, len(nodes)))
                ys.extend(yy + nodes*(y1-yy))
                measures.extend(ww*(y1-yy)*weights)
    return np.asarray(xs), np.asarray(ys), np.asarray(measures)


def local_basis(spline, points, cell, derivative=0):
    return spline(points, nu=derivative)[:, cell:cell+3]


def gradients(spline, x, y, i, j):
    nx, ny = local_basis(spline,x,i), local_basis(spline,y,j)
    dx, dy = local_basis(spline,x,i,1), local_basis(spline,y,j,1)
    return ((dx[:,:,None]*ny[:,None,:]).reshape(-1,9),
            (nx[:,:,None]*dy[:,None,:]).reshape(-1,9))


def strain_matrix(dx, dy):
    result = np.zeros((len(dx),3,18))
    result[:,0,::2] = dx
    result[:,1,1::2] = dy
    result[:,2,::2] = dy
    result[:,2,1::2] = dx
    return result


def local_dofs(i, j, width):
    nodes = ((i+np.arange(3))[:,None]*width+j+np.arange(3)).ravel()
    return (2*nodes[:,None]+np.arange(2)).ravel()


def run(n, order, radius):
    if n < 2 or not 2 <= order <= 32 or not 0 < radius < 10:
        raise ValueError("expected n>=2, quadrature 2..32 and a bore inside the square")
    start = time.perf_counter()
    side, height, tension, young, nu = 10., 8., 100., 210000., .3
    mu = young / (2*(1+nu))
    lam = young*nu / ((1+nu)*(1-2*nu))
    material = np.array([[lam+2*mu,lam,0.], [lam,lam+2*mu,0.], [0.,0.,mu]])
    breaks = np.linspace(0.,side,n+1)
    knots = np.r_[[0.,0.],breaks,[side,side]]
    width = n+2
    spline = BSpline(knots,np.eye(width),2,extrapolate=False)
    q, wq = leggauss(order)
    q, wq = (1+q)/2, wq/2
    dim = 2*width*width
    rhs = np.zeros(dim)
    rows, cols, entries = [], [], []
    fraction, material_area, cut_count, active_count = 1., 0., 0, 0
    quad_time = 0.
    for i in range(n):
        x0,x1 = breaks[i:i+2]
        for j in range(n):
            y0,y1 = breaks[j:j+2]
            if x1*x1+y1*y1 <= radius*radius:
                continue
            tic = time.perf_counter()
            x,y,measure = material_rule(x0,x1,y0,y1,radius,q,wq)
            quad_time += time.perf_counter()-tic
            area = float(np.sum(measure))
            if area <= 0.:
                continue
            active_count += 1
            if x0*x0 + y0*y0 < radius*radius:
                cut_count += 1
                fraction = min(fraction, area/((x1-x0)*(y1-y0)))
            material_area += area
            b = strain_matrix(*gradients(spline,x,y,i,j))
            k = np.einsum("qai,qaj->ij", b, (material @ b)*(height*measure[:,None,None]), optimize=True)
            dofs = local_dofs(i,j,width)
            rows.extend(np.repeat(dofs,18))
            cols.extend(np.tile(dofs,18))
            entries.extend(k.ravel())
    # Exact external traction from the infinite-domain solution restricted to
    # the finite square. Values on the two symmetry planes generate reactions.
    for t in range(n):
        lo,hi = breaks[t:t+2]
        p = lo+(hi-lo)*q
        basis = local_basis(spline,p,t)
        for axis in (0,1):
            indices = (width-1)*width+t+np.arange(3) if axis == 0 else (t+np.arange(3))*width+width-1
            x,y = (np.full_like(p,side),p) if axis == 0 else (p,np.full_like(p,side))
            xx,yy,xy = kirsch_tensor(x,y,radius,tension)
            traction = np.column_stack((xx,xy) if axis == 0 else (xy,yy))
            force = basis.T @ (traction*(height*(hi-lo)*wq[:,None]))
            for c in (0,1):
                rhs[2*indices+c] += force[:,c]
    k = coo_matrix((entries,(rows,cols)),shape=(dim,dim)).tocsr()
    assembly_end = time.perf_counter()
    active = k.diagonal() > 0
    essential = np.zeros(dim, dtype=bool)
    essential[2*np.arange(width)] = True
    essential[2*(width*np.arange(width))+1] = True
    free = active & ~essential
    reduced = k[free][:,free].tocsc()
    scale = 1./np.sqrt(reduced.diagonal())
    pre = diags(scale)
    scaled = (pre @ reduced @ pre).tocsc()
    precondition_end = time.perf_counter()
    u = np.zeros(dim)
    u[free] = scale*splu(scaled).solve(scale*rhs[free])
    solve_end = time.perf_counter()
    residual = k @ u - rhs
    mask = essential & active
    reaction = np.where(mask,residual,0.).reshape(-1,2).sum(axis=0)
    load = rhs.reshape(-1,2).sum(axis=0)
    force_error = np.linalg.norm(reaction+load)/np.linalg.norm(load)
    # On the upper bore axis the hoop direction is x. Evaluate the material
    # limit from y>=a. A C1 spline has a unique first derivative at knots.
    j = min(n-1, int(np.searchsorted(breaks,radius,side="right")-1))
    dx,dy = gradients(spline,np.array([0.]),np.array([radius]),0,j)
    hoop = float((material @ (strain_matrix(dx,dy)[0] @ u[local_dofs(0,j,width)]))[0])
    if not np.isfinite(u).all() or not math.isfinite(hoop):
        raise ValueError("nonfinite immersed solution")
    return dict(method="immersed_quadratic_bspline", n=n, quadrature=order,
                radius_mm=radius, background_splines=width*width,
                active_vector_dofs=int(np.sum(active)), free_vector_dofs=int(np.sum(free)),
                active_cells=active_count, cut_cells=cut_count, min_material_fraction=fraction,
                nnz=int(reduced.nnz), material_area_mm2=material_area,
                area_relative_error=abs(material_area/(side**2-math.pi*radius**2/4)-1),
                hoop_mpa=hoop, hoop_relative_error=abs(hoop/(3*tension)-1),
                force_balance_relative=force_error,
                free_residual_relative=float(np.linalg.norm(residual[free])/np.linalg.norm(rhs[free])),
                resultant_n=load.tolist(), reaction_n=reaction.tolist(),
                quadrature_geometry_s=quad_time, assembly_including_quadrature_s=assembly_end-start,
                preconditioning_s=precondition_end-assembly_end, solve_s=solve_end-precondition_end,
                total_s=time.perf_counter()-start,
                peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells",type=int)
    parser.add_argument("order",type=int)
    parser.add_argument("radius",type=float)
    args = parser.parse_args()
    print(json.dumps(run(args.cells,args.order,args.radius), allow_nan=False))
