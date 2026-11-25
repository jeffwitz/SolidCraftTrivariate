"""Physical L2 and energy distances on common stratified, scrambled Sobol points."""
import argparse
import json
import os
from pathlib import Path
import struct
import sys
import time

os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy.stats import qmc, t as student
from fillet_probe import DL, EDGES, read_mesh
from immersed_bracket import Space
from solver_package import shape

MU=210000./2.6
LAM=210000.*.3/(1.3*.4)


AFFINE=np.array([[.002,.001,-.003],[-.002,.004,.001],[.001,.002,-.001]])
SHIFT=np.array([.3,-.2,.07])


def verify_affine_readers(blocks):
    """Exact nested rational restrictions and spline Greville affine controls."""
    rng=np.random.default_rng(730138)
    uvw=rng.random((109,3))
    block=next(b for b in blocks if b["id"].startswith("hole_x_sector"))
    controls=np.asarray(block["controls"],dtype=float).reshape(3,3,3,4)
    _,x,j,_=shape(uvw,controls.reshape(27,4))
    def child(net,axis,bit):
        c=np.moveaxis(net,axis,0)
        mid=.25*c[0]+.5*c[1]+.25*c[2]
        out=np.stack((mid,.5*(c[1]+c[2]),c[2])) if bit else np.stack((c[0],.5*(c[0]+c[1]),mid))
        return np.moveaxis(out,0,axis)
    native=NativeField.__new__(NativeField)
    native.n=2; native.blocks={"test":np.empty((2,2,2,27,4))}
    for bits in np.ndindex(2,2,2):
        h=controls
        for axis,bit in enumerate(bits): h=child(h,axis,bit)
        h=h.reshape(27,4)
        w=h[:,3];v=(h[:,:3]/w[:,None])@AFFINE.T+SHIFT
        native.blocks["test"][bits]=np.column_stack((w/w.max(),v))
    immersed=ImmersedField.__new__(ImmersedField)
    immersed.h=2.5;immersed.space=Space(2.5)
    greville=[]
    for spl in immersed.space.splines:
        greville.append((spl.t[1:-2]+spl.t[2:-1])*.5)
    xyz=np.stack(np.meshgrid(*greville,indexing="ij"),axis=-1)
    immersed.u=xyz@AFFINE.T+SHIFT
    checks={}
    for name,(u,g) in [("native_q2_h1",native.evaluate("test",uvw,np.linalg.inv(j))),
                       ("immersed_q2",immersed.evaluate(x))]:
        checks[name]={"u_max_mm":float(np.abs(u-(x@AFFINE.T+SHIFT)).max()),
                      "gradient_max":float(np.abs(g-AFFINE).max())}
        assert max(checks[name].values())<1e-11
    return checks


def q2(t):
    pol=np.stack(((1-t)**2,2*t*(1-t),t**2),axis=2)
    dp=np.stack((2*(t-1),2-4*t,2*t),axis=2)
    def tensor(a,b,c):
        return np.einsum("qi,qj,qk->qijk",a,b,c).reshape(-1,27)
    return (tensor(*pol.transpose(1,0,2)),
            np.stack([tensor(*(dp[:,i] if i==a else pol[:,i] for i in range(3))) for a in range(3)],axis=2))


class NativeField:
    def __init__(self,path):
        with path.open("rb") as f:
            assert f.read(8)==b"SCQ2FLD1"
            self.level,self.rx,self.ry,self.q,count=struct.unpack("<BBBBI",f.read(8))
            self.n=2**self.level
            self.blocks={}
            seen=set()
            for _ in range(count):
                length,=struct.unpack("<H",f.read(2))
                name=f.read(length).decode()
                base,*children=name.split("/")
                assert len(children)==self.level and name not in seen
                seen.add(name)
                cell=np.zeros(3,dtype=int)
                for child in children:
                    assert len(child)==3 and set(child)<=set("01")
                    cell=cell*2+np.array(list(child),dtype=int)
                if base not in self.blocks:
                    self.blocks[base]=np.empty((self.n,self.n,self.n,27,4))
                net=np.frombuffer(f.read(27*4*8),dtype="<f8").reshape(27,4).copy()
                assert np.all(net[:,0]>0) and np.isfinite(net).all()
                net[:,0]/=net[:,0].max()
                self.blocks[base][tuple(cell)]=net
            assert count==len(self.blocks)*self.n**3 and f.read()==b""

    def evaluate(self,base,t,invj):
        scaled=t*self.n
        cell=np.minimum(scaled.astype(int),self.n-1)
        net=self.blocks[base][tuple(cell.T)]
        b,db=q2(scaled-cell)
        weights=net[:,:,0]
        denom=np.sum(b*weights,axis=1)
        u=np.einsum("qi,qi,qij->qj",b,weights,net[:,:,1:])/denom[:,None]
        du=np.einsum("qia,qi,qij->qja",db,weights,net[:,:,1:]-u[:,None,:])/denom[:,None,None]
        return u, (self.n*du)@invj


class ImmersedField:
    def __init__(self,path):
        with np.load(path) as data:
            self.h=float(data["h_mm"]);self.radii=data["radii_mm"]
            self.space=Space(self.h);self.u=data["displacement_mm"]
            assert self.u.shape==(*self.space.shape,3)

    def evaluate(self,x):
        cell=np.minimum((x/self.h).astype(int),np.array(self.space.shape)-3)
        assert np.all(cell>=0)
        v=[];d=[]
        for i,s in enumerate(self.space.splines):
            columns=cell[:,i,None]+np.arange(3)
            v.append(np.take_along_axis(s(x[:,i]),columns,axis=1))
            d.append(np.take_along_axis(s(x[:,i],nu=1),columns,axis=1))
        offset=np.array(np.unravel_index(np.arange(27),(3,3,3))).T
        idx=cell[:,None,:]+offset[None,:,:]
        net=self.u[idx[:,:,0],idx[:,:,1],idx[:,:,2]]
        def contract(a,b,c):
            b3=np.einsum("qi,qj,qk->qijk",a,b,c).reshape(-1,27)
            return np.einsum("qi,qij->qj",b3,net)
        return contract(*v),np.stack([contract(*(d[i] if i==j else v[i] for i in range(3))) for j in range(3)],axis=2)


def tetra_basis(l):
    n=np.empty((len(l),10)); dn=np.empty((len(l),10,3))
    n[:,:4]=l*(2*l-1)
    dn[:,:4]=(4*l-1)[:,:,None]*DL
    for k,(i,j) in enumerate(EDGES,4):
        n[:,k]=4*l[:,i]*l[:,j]
        dn[:,k]=4*(l[:,i,None]*DL[j]+l[:,j,None]*DL[i])
    return n,dn


class TetraField:
    def __init__(self,directory):
        nodes,elements=read_mesh((directory/"bracket.inp").read_text())
        max_id=max(nodes); xyz=np.empty((max_id+1,3)); u=np.full_like(xyz,np.nan)
        for i,p in nodes.items(): xyz[i]=p
        lines=(directory/"bracket.dat").read_text().split("for set NALL and time")[-1].splitlines()[1:]
        for line in lines:
            parts=line.split()
            if len(parts)==4 and parts[0].isdigit():
                u[int(parts[0])]=[float(v) for v in parts[1:]]
        conn=np.asarray(list(elements.values()))
        self.x=xyz[conn];self.u=u[conn]
        assert np.all(np.isfinite(self.u))
        self.inv=np.linalg.inv(self.x[:,1:4]-self.x[:,0,None,:])
        # Quadratic Lagrange points do not bound a curved tetrahedron. Bound
        # its entire image using the equivalent degree-two Bernstein hull.
        hull=self.x.copy()
        for k,(i,j) in enumerate(EDGES,4):
            hull[:,k]=2*self.x[:,k]-.5*(self.x[:,i]+self.x[:,j])
        self.lo=hull.min(axis=1)-1e-9;self.hi=hull.max(axis=1)+1e-9
        self.width=2.
        lower=np.floor(self.lo/self.width).astype(int)
        upper=np.floor(self.hi/self.width).astype(int)
        self.origin=lower.min(axis=0);self.dims=upper.max(axis=0)-self.origin+1
        buckets={}
        stride=np.array([self.dims[1]*self.dims[2],self.dims[2],1])
        for eid,(a,b) in enumerate(zip(lower,upper)):
            for i in range(a[0],b[0]+1):
                for j in range(a[1],b[1]+1):
                    base=int((np.array([i,j,a[2]])-self.origin)@stride)
                    for dz in range(b[2]-a[2]+1):
                        buckets.setdefault(base+dz,[]).append(eid)
        self.table=np.full((int(np.prod(self.dims)),max(map(len,buckets.values()))),-1,dtype=np.int32)
        for ix,es in buckets.items():self.table[ix,:len(es)]=es

    def evaluate(self,points):
        cell=np.floor(points/self.width).astype(int)-self.origin
        bucket=np.ravel_multi_index(cell.T,self.dims)
        pp,cc=np.nonzero(self.table[bucket]>=0)
        ee=self.table[bucket[pp],cc]
        keep=np.all((points[pp]>=self.lo[ee])&(points[pp]<=self.hi[ee]),axis=1)
        ii=pp[keep];ee=ee[keep]
        geo=self.x[ee]
        l=np.empty((len(ee),4))
        l[:,1:]=np.einsum("qi,qij->qj",points[ii]-geo[:,0],self.inv[ee])
        for _ in range(9):
            l[:,0]=1.-l[:,1:].sum(axis=1)
            n,dn=tetra_basis(l)
            residual=np.einsum("qi,qic->qc",n,geo)-points[ii]
            jac=np.einsum("qic,qia->qca",geo,dn)
            active=np.linalg.norm(residual,axis=1)>2e-10
            if not active.any():break
            # Polynomial extensions outside candidate tetrahedra can contain
            # singularities; eliminate only these unsuccessful candidate roots.
            active &= np.abs(np.linalg.det(jac))>1e-12
            if not active.any():break
            delta=np.linalg.solve(jac[active],residual[active,:,None])[:,:,0]
            l[active,1:]-=np.clip(delta,-2.,2.)
        l[:,0]=1.-l[:,1:].sum(axis=1)
        n,dn=tetra_basis(l)
        residual=np.linalg.norm(np.einsum("qi,qic->qc",n,geo)-points[ii],axis=1)
        good=(l.min(axis=1)>=-1e-10)&(residual<=1e-8)
        hits=np.bincount(ii[good],minlength=len(points))
        duplicate=int(np.count_nonzero(hits>1))
        inside=hits>0
        # Deterministic selection among boundary-equivalent candidates.
        chosen=np.flatnonzero(good)
        _,index=np.unique(ii[chosen],return_index=True)
        chosen=chosen[index]
        field=np.zeros((len(points),3));gradient=np.zeros((len(points),3,3))
        selected=ee[chosen]; x=self.x[selected]; u=self.u[selected]
        jac=np.einsum("qic,qia->qca",x,dn[chosen])
        assert np.all(np.linalg.det(jac)>0)
        ij=np.linalg.inv(jac)
        field[inside]=np.einsum("qi,qic->qc",n[chosen],u)
        gradient[inside]=np.einsum("qic,qia->qca",u,dn[chosen])@ij
        manufactured=x@AFFINE.T+SHIFT
        value=np.einsum("qi,qic->qc",n[chosen],manufactured)
        affine=np.einsum("qic,qia->qca",manufactured,dn[chosen])@ij
        affine_error=max(float(np.abs(affine-AFFINE).max()),float(np.abs(value-(points[inside]@AFFINE.T+SHIFT)).max())) if len(chosen) else 0.
        return field,gradient,inside,duplicate,affine_error


def density(grad):
    eps=(grad+grad.transpose(0,2,1))*.5
    return 2*MU*np.einsum("qij,qij->q",eps,eps)+LAM*np.trace(eps,axis1=1,axis2=2)**2


def evaluate(args):
    base=Path(args.directory)
    qx=args.radius_ticks
    root=base.parent
    schema="fixed" if qx==64 else "edited"
    package=json.loads((root/f"solver-package/{schema}-q8.json").read_text())
    original=package["analysis_space"]["blocks"]
    efdir=root/"ef1-h1.25" if qx==64 else root/"r14-tet/h1.25"
    tetra_dirs={"tetra_h1.25":efdir}
    for size_text, path in args.additional_tetra:
        size=float(size_text)
        if not np.isfinite(size) or size<=0:
            raise ValueError("Tetrahedral reference size must be finite and positive")
        key=f"tetra_h{size:g}"
        if key in tetra_dirs:
            raise ValueError("Repeated tetrahedral reference size")
        tetra_dirs[key]=Path(path)
    tetra={key:TetraField(path) for key,path in tetra_dirs.items()}
    native={f"native_h{h}":NativeField(base/f"native-x{qx}-h{h}.bin") for h in range(4)}
    immersed={f"immersed_h{h}":ImmersedField(base/f"immersed-r{qx/16:g}-h{h}.npz") for h in (1.25,1)}
    scales={}
    for h in range(4):scales[f"native_h{h}"]=json.loads((base/f"native-x{qx}-h{h}.json").read_text())["strain_energy_n_mm"]
    for h in (1.25,1):scales[f"immersed_h{h}"]=json.loads((base/f"immersed-r{qx/16:g}-h{h}.json").read_text())["strain_energy_n_mm"]
    for key,path in tetra_dirs.items():
        tet_result=json.loads((path/"result.json").read_text())
        scales[key]=-500.*tet_result["mean_tip_y_mm"]
    counts=sorted(set([128,512,2048,2**args.log2_points]))
    counts=[c for c in counts if c<=2**args.log2_points]
    samples={n:[{} for _ in range(8)] for n in counts}
    duplicates=0;affine_error=0.;missing_points=[]
    inverse_counts={k:{"missing":0,"multiple":0} for k in tetra}
    for bi,block in enumerate(sorted(original,key=lambda b:b["id"])):
        for rep in range(8):
            t=qmc.Sobol(3,scramble=True,seed=20261125+1000*bi+rep).random_base2(args.log2_points)
            _,x,j,_=shape(t,block["controls"])
            weights=np.linalg.det(j)
            tetra_fields={}
            ok=np.ones(len(x),dtype=bool)
            for key,ref in tetra.items():
                ut,gt,located,nmulti,ae=ref.evaluate(x)
                duplicates+=nmulti;affine_error=max(affine_error,ae)
                inverse_counts[key]["missing"]+=int(np.count_nonzero(~located))
                inverse_counts[key]["multiple"]+=nmulti
                tetra_fields[key]=(ut,gt)
                ok &= located
            if not ok.all() and len(missing_points)<16:
                missing_points.extend(x[~ok][:16-len(missing_points)].tolist())
            invj=np.linalg.inv(j)
            fields={k:m.evaluate(block["id"],t,invj) for k,m in native.items()}
            fields.update({k:m.evaluate(x) for k,m in immersed.items()})
            fields.update(tetra_fields)
            data={"full_volume":np.ones(len(t)),"missing_volume":~ok,"missing_points":(~ok)/weights}
            for k,(u,g) in fields.items():
                data[k+"_energy"]=density(g)*ok/2.
                data[k+"_l2"]=np.einsum("qi,qi->q",u,u)*ok
            pairs=[(k,b) for b in tetra for k in (*native,*immersed)]
            ordered_tetra=sorted(tetra,key=lambda k:float(k.removeprefix("tetra_h")),reverse=True)
            pairs.extend((a,b) for i,a in enumerate(ordered_tetra) for b in ordered_tetra[i+1:])
            pairs+=[("native_h3","immersed_h1"),("immersed_h1.25","immersed_h1"),("native_h2","native_h3")]
            for a,b in pairs:
                u,g=fields[a];v,h=fields[b];key=f"{a}__{b}"
                data[key+"_energy"]=density(g-h)*ok/2.
                data[key+"_l2"]=np.sum((u-v)**2,axis=1)*ok
            for n in counts:
                row=samples[n][rep]
                for k,v in data.items():
                    row[k]=row.get(k,0.)+float(np.mean((weights*v)[:n]))
        print(f"field {bi+1}/{len(original)} {block['id']}",file=sys.stderr,flush=True)
    output={"radius_x_mm":qx/16,"radius_y_mm":4.,"scrambles":8,"blocks":len(original),
            "affine_tetra_max_error":affine_error,"affine_reader_checks":verify_affine_readers(original),"multiple_tetra_candidates":int(duplicates),
            "tetra_inverse_counts":inverse_counts,"example_unlocated_mm":missing_points,"reference_energy_n_mm":scales,"levels":[]}
    for n in counts:
        rows=samples[n]
        arr=lambda k:np.array([r[k] for r in rows])
        energy={k:float(np.mean(arr(k+"_energy"))) for k in scales}
        distances={}
        for a,b in pairs:
            key=f"{a}__{b}"
            for norm in ("energy","l2"):
                ratios=np.sqrt(arr(key+"_"+norm)/arr(b+"_"+norm))
                mean=float(ratios.mean())
                half=float(student.ppf(.975,7)*ratios.std(ddof=1)/np.sqrt(8))
                distances[key+"_"+norm]={"mean":mean,"empirical_half_width_95":half,
                                        "dispersion_within_budget":half<=max(.001,.05*mean)}
        output["levels"].append({"points_per_block_per_scramble":n,"total_points":len(original)*8*n,
            "unlocated_points":int(round(arr("missing_points").sum()*n)),
            "unlocated_volume_fraction":float(arr("missing_volume").sum()/arr("full_volume").sum()),
            "integrated_volume_mm3":float(arr("full_volume").mean()),"energy_n_mm":energy,
            "relative_energy_recovery":{k:abs(v/scales[k]-1.) for k,v in energy.items()},
            "relative_distances":distances})
    assert affine_error<1e-10
    current=output["levels"][-1]
    previous=output["levels"][-2] if len(counts)>1 else None
    stability={}
    for k,detail in current["relative_distances"].items():
        delta=abs(detail["mean"]-previous["relative_distances"][k]["mean"]) if previous else None
        stability[k]={"last_level_change":delta,"stable":delta is not None and delta<=max(.001,.05*detail["mean"])}
    output["qualification"]={"stability":stability,"energy_recovery":max(current["relative_energy_recovery"].values())<=.015,
            "all_points_located_once":current["unlocated_points"]==0 and duplicates==0,
            "stable_distances":all(r["stable"] for r in stability.values()) and all(r["dispersion_within_budget"] for r in current["relative_distances"].values())}
    return output


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("--radius-ticks",type=int,choices=(56,64),required=True)
    p.add_argument("--log2-points",type=int,choices=(7,9,11,13),default=11)
    p.add_argument("--additional-tetra",nargs=2,metavar=("SIZE_MM","DIRECTORY"),action="append",default=[])
    a=p.parse_args();start=time.monotonic();report=evaluate(a)
    report["elapsed_s"]=time.monotonic()-start
    print(json.dumps(report,indent=2,allow_nan=False))
