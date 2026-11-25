"""Compare adaptive Q2 fields on a common physical Sobol sample with uniform h3."""
import argparse
import json
import os
from pathlib import Path
import struct
import sys
import time

os.environ["OPENBLAS_NUM_THREADS"]="1"
import numpy as np
from scipy.stats import qmc, t as student

from bracket_field_distance import NativeField, density, q2
from solver_package import shape


class AdaptiveField:
    """Locate each parameter point by exact integer dyadic voxel ownership."""
    def __init__(self,path):
        self.voxels={}
        controls=[];scales=[];indices=[]
        with path.open("rb") as f:
            if f.read(8)!=b"SCQ2ADP1":
                raise ValueError("Invalid adaptive field header")
            depth,self.rx,self.ry,self.q,count=struct.unpack("<BBBBI",f.read(8))
            if depth>5:
                raise ValueError("Adaptive reference reader supports levels 0..5")
            self.grid=1<<depth
            for tile in range(count):
                length,=struct.unpack("<H",f.read(2))
                name=f.read(length).decode()
                root,*bits=name.split("/")
                if len(bits)>depth:
                    raise ValueError("Leaf exceeds declared depth")
                ijk=np.zeros(3,dtype=int)
                for child in bits:
                    if len(child)!=3 or not set(child)<={"0","1"}:
                        raise ValueError("Malformed child path")
                    ijk=2*ijk+np.array(list(child),dtype=int)
                stride=self.grid>>len(bits)
                if root not in self.voxels:
                    self.voxels[root]=np.full((self.grid,)*3,-1,dtype=np.int32)
                region=self.voxels[root][tuple(slice(i*stride,(i+1)*stride) for i in ijk)]
                if (region!=-1).any():
                    raise ValueError("Overlapping leaves")
                region[:]=tile
                net=np.frombuffer(f.read(27*4*8),dtype="<f8").reshape(27,4).copy()
                if not np.isfinite(net).all() or net[:,0].min()<=0:
                    raise ValueError("Invalid weighted field coefficients")
                net[:,0]/=net[:,0].max()
                controls.append(net);indices.append(ijk);scales.append(1<<len(bits))
            if f.read() or any(np.any(v<0) for v in self.voxels.values()):
                raise ValueError("Incomplete parameter cover or extraneous coefficients")
        self.nets=np.array(controls);self.indices=np.array(indices);self.scales=np.array(scales)

    def evaluate(self,root,t,invj):
        vox=np.minimum((t*self.grid).astype(int),self.grid-1)
        tile=self.voxels[root][tuple(vox.T)]
        scale=self.scales[tile]
        local=t*scale[:,None]-self.indices[tile]
        b,db=q2(local);net=self.nets[tile];w=net[:,:,0];v=net[:,:,1:]
        denom=(b*w).sum(axis=1)
        u=np.einsum("qi,qi,qij->qj",b,w,v)/denom[:,None]
        grad=np.einsum("qia,qi,qij->qja",db,w,v-u[:,None,:])/denom[:,None,None]
        return u,(grad*scale[:,None,None])@invj


def compare(field_dir,adaptive_dir,ticks,log2=11):
    schema="fixed" if ticks==64 else "edited"
    package=json.loads((field_dir.parent/f"solver-package/{schema}-q8.json").read_text())
    blocks=sorted(package["analysis_space"]["blocks"],key=lambda b:b["id"])
    ref=NativeField(field_dir/f"native-x{ticks}-h3.bin")
    fields={f"uniform_h{h}":NativeField(field_dir/f"native-x{ticks}-h{h}.bin") for h in (0,1,2)}
    for target in sorted(adaptive_dir.glob("round-*.bin")):
        field=AdaptiveField(target)
        assert field.rx==ticks and field.ry==64 and field.q==8 and set(field.voxels)==set(ref.blocks)
        fields[target.stem]=field
    if all(not key.startswith("round-") for key in fields):
        raise ValueError("No completed adaptive checkpoint")
    counts=[n for n in (128,512,2048,8192) if n<=2**log2]
    totals={n:{key:np.zeros((8,2)) for key in (*fields,"reference")} for n in counts}
    for bi,b in enumerate(blocks):
        for s in range(8):
            t=qmc.Sobol(3,scramble=True,seed=20261125+1000*bi+s).random_base2(log2)
            _,_,j,_=shape(t,b["controls"]);ij=np.linalg.inv(j);measure=np.linalg.det(j)
            u,g=ref.evaluate(b["id"],t,ij)
            def add(name,v,h):
                a=measure*np.einsum("ij,ij->i",v,v)
                e=measure*density(h)
                for n in counts:
                    totals[n][name][s]+=[a[:n].mean(),e[:n].mean()]
            add("reference",u,g)
            for name,model in fields.items():
                v,h=model.evaluate(b["id"],t,ij)
                add(name,v-u,h-g)
        print(f"distance {bi+1}/{len(blocks)}",file=sys.stderr,flush=True)
    rows=[]
    for n in counts:
        row={}
        for name in fields:
            rel=np.sqrt(totals[n][name]/totals[n]["reference"])
            mean=rel.mean(axis=0)
            half=student.ppf(.975,7)*rel.std(axis=0,ddof=1)/np.sqrt(8)
            row[name]={"mean_l2_relative":float(mean[0]),"mean_energy_relative":float(mean[1]),
                       "half_width_95_l2":float(half[0]),"half_width_95_energy":float(half[1])}
        rows.append({"points_total":8*n*len(blocks),"points_per_scramble_block":n,"methods":row})
    return {"radius_x_ticks":ticks,"reference":"uniform native h3/q8","levels":rows}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("field_directory",type=Path);p.add_argument("adaptive_directory",type=Path)
    p.add_argument("--radius-ticks",type=int,choices=(56,64),required=True)
    p.add_argument("--log2",type=int,choices=(7,9,11,13),default=11)
    a=p.parse_args();start=time.monotonic()
    result=compare(a.field_directory,a.adaptive_directory,a.radius_ticks,a.log2)
    result["elapsed_s"]=time.monotonic()-start
    print(json.dumps(result,indent=2,allow_nan=False))
