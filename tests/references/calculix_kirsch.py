"""OCC, quadratic Netgen tetrahedra and CalculiX reference for the Kirsch pilot."""

import FreeCAD as App
import Part

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import resource
import subprocess
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

from calculix_bracket import block_ids, node_report
from fillet_probe import basis, read_mesh
from immersed_kirsch import kirsch_tensor


def load(mesh, nodes, radius, order):
    q,w=leggauss(order)
    q,w=(1+q)/2,w/2
    l0=(1-q[:,None])*(1-q[None,:])
    l1=np.broadcast_to(q[:,None],l0.shape)
    l2=(1-q[:,None])*q[None,:]
    bary=np.column_stack((l0.ravel(),l1.ravel(),l2.ravel()))
    weights=((1-q[:,None])*w[:,None]*w[None,:]).ravel()
    a,b,c=bary.T
    shape=np.column_stack((a*(2*a-1),b*(2*b-1),c*(2*c-1),4*a*b,4*b*c,4*c*a))
    nodal=defaultdict(lambda: np.zeros(2))
    area=0.
    for face in mesh.Faces:
        conn=mesh.getElementNodes(face)
        coords=np.array([tuple(nodes[i]) for i in conn])
        planes=np.flatnonzero(np.max(abs(coords[:,:2]-10),axis=0)<1e-8)
        if not len(planes):
            continue
        assert len(planes)==1 and len(coords)==6
        for idx,(i,j) in enumerate(((0,1),(1,2),(2,0)),3):
            assert np.linalg.norm(coords[idx]-(coords[i]+coords[j])/2)<1e-7
        jac=np.linalg.norm(np.cross(coords[1]-coords[0],coords[2]-coords[0]))
        points=bary @ coords[:3]
        xx,yy,xy=kirsch_tensor(points[:,0],points[:,1],radius,100.)
        traction=np.column_stack((xx,xy) if planes[0]==0 else (xy,yy))
        local=shape.T @ (traction*(jac*weights[:,None]))
        for i,value in zip(conn,local):
            nodal[i]+=value
        area+=jac/2
    if abs(area-160.)>1e-6:
        raise ValueError(f"incorrect external area {area}")
    return nodal,area


def stresses(nodes, elements, disp, radius):
    target=np.array([0.,radius,4.])
    result=[]
    for eid,conn in elements.items():
        coords=np.array([nodes[i] for i in conn])
        if np.any(coords.min(axis=0)-1e-5>target) or np.any(coords.max(axis=0)+1e-5<target):
            continue
        # Parametric inversion of the curved isoparametric tetrahedron.
        s=np.linalg.solve((coords[1:4]-coords[0]).T,target-coords[0])
        for _ in range(12):
            bary=np.r_[1-s.sum(),s]
            values,grad=basis(bary)
            residual=values @ coords - target
            if np.linalg.norm(residual)<1e-11: break
            s-=np.linalg.solve(coords.T @ grad,residual)
        if min(bary)<-2e-6 or np.linalg.norm(residual)>1e-9:
            continue
        jac=coords.T @ grad
        if np.linalg.det(jac)<=0:
            raise ValueError("nonpositive tetrahedron at probe")
        gradient=np.array([disp[i] for i in conn]).T @ grad @ np.linalg.inv(jac)
        mu=210000./2.6
        lam=210000.*.3/(1.3*.4)
        result.append({"element":eid,"barycentric":bary.tolist(),
                       "point_error_mm":float(np.linalg.norm(residual)),
                       "hoop_mpa":float(lam*np.trace(gradient)+2*mu*gradient[0,0])})
    if not result:
        raise ValueError("no incident tetrahedron at the material-side probe")
    return result


def run(args):
    root=args.directory.resolve()
    root.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    doc=App.newDocument("KirschReference")
    body=Part.makeBox(10,10,8).cut(Part.makeCylinder(args.radius,8))
    if not body.isValid() or len(body.Solids)!=1 or abs(body.Volume-8*(100-math.pi*args.radius**2/4))>1e-7:
        raise ValueError("incorrect reference solid")
    feature=doc.addObject("PartDesign::Feature","QuarterPlate")
    feature.Shape=body
    obj=doc.addObject("Fem::FemMeshShapeNetgenObject","Tet10")
    obj.Shape=feature
    obj.SecondOrder=True
    obj.MaxSize=args.size
    doc.recompute()
    mesh=obj.FemMesh
    nodes=mesh.Nodes
    if not mesh.VolumeCount or len(mesh.getElementNodes(mesh.Volumes[0]))!=10:
        raise ValueError("C3D10 elements required")
    mesh_done=time.perf_counter()
    sx=[i for i,p in nodes.items() if abs(p.x)<1e-8]
    sy=[i for i,p in nodes.items() if abs(p.y)<1e-8]
    force,area=load(mesh,nodes,args.radius,args.quadrature)
    deck=root/"kirsch.inp"
    mesh.writeABAQUS(str(deck),1,False)
    with deck.open("a") as f:
        f.write("\n*NSET,NSET=SYMX\n"+block_ids(sx)+"\n*NSET,NSET=SYMY\n"+block_ids(sy)+"\n")
        f.write("*MATERIAL,NAME=STEEL\n*ELASTIC\n210000,0.3\n*SOLID SECTION,ELSET=Eall,MATERIAL=STEEL\n")
        f.write("*STEP\n*STATIC\n*BOUNDARY\nSYMX,1,1\nSYMY,2,2\nNall,3,3\n*CLOAD\n")
        for i,val in sorted(force.items()):
            for c in range(2):
                f.write(f"{i},{c+1},{val[c]:.12e}\n")
        f.write("*NODE PRINT,NSET=SYMX\nRF\n*NODE PRINT,NSET=SYMY\nRF\n*NODE PRINT,NSET=Nall\nU\n*END STEP\n")
    deck_done=time.perf_counter()
    with (root/"ccx.log").open("w") as log:
        subprocess.run([args.ccx,"-i","kirsch"],cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
    solve_done=time.perf_counter()
    dat=(root/"kirsch.dat").read_text()
    disp=node_report(dat,"displacements","NALL")
    # CalculiX RF includes nodal external loads at restrained coordinates.
    # Subtract those loads when deriving the support reaction from RF.
    rx=node_report(dat,"forces","SYMX")
    ry=node_report(dat,"forces","SYMY")
    reaction=np.array([sum(rx[i][0]-force.get(i,(0,0))[0] for i in sx),
                       sum(ry[i][1]-force.get(i,(0,0))[1] for i in sy)])
    load_sum=np.sum(list(force.values()),axis=0)
    error=float(np.linalg.norm(reaction+load_sum)/np.linalg.norm(load_sum))
    mesh_nodes,elements=read_mesh(deck.read_text())
    samples=stresses(mesh_nodes,elements,disp,args.radius)
    result=dict(method="curved_tetrahedral_p2",radius_mm=args.radius,size_mm=args.size,
                quadrature=args.quadrature,nodes=len(nodes),tet10=mesh.VolumeCount,
                loaded_area_mm2=area,occ_quarter_volume_mm3=body.Volume,
                hoop_min_mpa=min(s["hoop_mpa"] for s in samples),hoop_max_mpa=max(s["hoop_mpa"] for s in samples),
                probes=samples,reaction_n=reaction.tolist(),resultant_n=load_sum.tolist(),force_balance_relative=error,
                mesh_s=mesh_done-start,load_deck_s=deck_done-mesh_done,solve_s=solve_done-deck_done,
                total_s=time.perf_counter()-start,solver_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                mesh_process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (root/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))
    if error>1e-6: raise ValueError(f"force imbalance {error}")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ccx",required=True)
    p.add_argument("--radius",type=float,required=True)
    p.add_argument("--size",type=float,required=True)
    p.add_argument("--quadrature",type=int,default=6)
    p.add_argument("directory",type=Path)
    run(p.parse_args())
