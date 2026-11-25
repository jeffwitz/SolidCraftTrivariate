"""Recover C3D10 stress at the closest quadratic fillet surface point."""

import argparse
import itertools
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import minimize


EDGES = ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
DL = np.array([[-1.,-1.,-1.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])


def basis(l):
    values = [x*(2*x-1) for x in l]
    gradients = [(4*x-1)*d for x,d in zip(l,DL)]
    for i,j in EDGES:
        values.append(4*l[i]*l[j])
        gradients.append(4*(l[i]*DL[j]+l[j]*DL[i]))
    return np.array(values),np.array(gradients)


def read_mesh(text):
    nodes, elements = {}, {}
    kind=None
    for line in text.splitlines():
        if line.startswith("*"):
            kind = line.upper().split(",")[0]
            continue
        if not line.strip():
            continue
        if kind=="*NODE":
            i,x,y,z=line.split(",")
            nodes[int(i)]=np.array([float(x),float(y),float(z)])
        if kind=="*ELEMENT":
            i,*v=map(int,line.split(","))
            assert len(v)==10
            elements[i]=v
    return nodes,elements


def closest(elements,nodes):
    p=np.array([14-6/math.sqrt(2),14-6/math.sqrt(2),20.])
    best=None
    for eid,conn in elements.items():
        points=np.array([nodes[i] for i in conn])
        if np.any(points.min(axis=0)-1>p) or np.any(points.max(axis=0)+1<p):
            continue
        for face in itertools.combinations(range(4),3):
            face_nodes=list(face)+[4+i for i,e in enumerate(EDGES) if set(e)<=set(face)]
            surface=points[face_nodes]
            if np.max(np.abs(np.linalg.norm(surface[:,:2]-14,axis=1)-6))>1e-4:
                continue
            def bary(uv):
                l=np.zeros(4);l[list(face)]=[1-uv.sum(),*uv];return l
            def residual(uv):
                return basis(bary(uv))[0]@points-p
            def obj(uv):
                err=residual(uv)
                return err@err
            # Start from the projection onto the flat triangle.
            initial=np.linalg.lstsq((surface[1:3]-surface[0]).T,p-surface[0],rcond=None)[0]
            opt=minimize(obj,initial,method="SLSQP",bounds=((0,1),(0,1)),
                         constraints=[{"type":"ineq","fun":lambda x:1-x.sum()}],
                         options={"ftol":1e-16,"maxiter":100})
            if opt.success and (best is None or opt.fun<best[0]):
                best=(opt.fun,eid,bary(opt.x))
    if best is None:
        raise ValueError("no fillet boundary triangle found")
    return best


def run(a):
    original=(a.source/"bracket.inp").read_text()
    nodes,elems=read_mesh(original)
    squared,eid,l=closest(elems,nodes)
    conn=elems[eid]
    deck=original.replace("*STEP\n","*NSET,NSET=PROBE\n"+",".join(map(str,conn))+"\n*STEP\n")
    deck=deck.replace("*END STEP","*NODE PRINT,NSET=PROBE\nU\n*END STEP")
    a.output.mkdir(exist_ok=True,parents=True)
    (a.output/"probe.inp").write_text(deck)
    env=dict(os.environ,LD_LIBRARY_PATH=str(a.ccx.parent.parent/"lib"),OMP_NUM_THREADS="2")
    source_data=(a.source/"bracket.dat").read_text()
    if "for set NALL and time" in source_data:
        lines=source_data.split("for set NALL and time")[-1].splitlines()[1:]
    else:
        with (a.output/"ccx.log").open("w") as log:
            subprocess.run([str(a.ccx),"-i","probe"],cwd=a.output,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        lines=(a.output/"probe.dat").read_text().split("for set PROBE and time")[-1].splitlines()[1:]
    u={}
    for line in lines:
        fields=line.strip().split()
        if len(fields)==4 and fields[0].isdigit():
            u[int(fields[0])]=list(map(float,fields[1:]))
        elif u and line.strip():
            break
    assert set(conn)<=u.keys()
    x=np.array([nodes[i] for i in conn]); disp=np.array([u[i] for i in conn])
    weights,dl=basis(l)
    jac=x.T@dl
    assert np.linalg.det(jac)>0
    gradient=disp.T@dl@np.linalg.inv(jac)
    eps=0.5*(gradient+gradient.T)
    mu=210000/(2*1.3);lam=210000*.3/(1.3*.4)
    sigma=lam*np.trace(eps)*np.eye(3)+2*mu*eps
    deviator=sigma-np.eye(3)*np.trace(sigma)/3
    report=dict(element=eid,barycentric=l.tolist(),point_mm=(weights@x).tolist(),
                distance_mm=math.sqrt(squared),stress_mpa=sigma.tolist(),
                von_mises_mpa=math.sqrt(1.5*np.sum(deviator**2)))
    (a.output/"result.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("source",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--ccx",type=Path,required=True)
    run(p.parse_args())
