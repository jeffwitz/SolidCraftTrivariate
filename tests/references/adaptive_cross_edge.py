"""Certify the mounted scalar dimension from independent physical box-face equations."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess
import time

PRIME=2147483647


def combine(row,key,delta):
    value=row.get(key,0)+delta
    if value:row[key]=value
    else:row.pop(key,None)


class Rank:
    """A modular pivot minor gives a lower bound on rank over the rationals."""
    def __init__(self):self.pivots={}
    def add(self,entries):
        r={k:(v.numerator*pow(v.denominator,-1,PRIME))%PRIME for k,v in entries.items()}
        r={k:w for k,w in r.items() if w}
        while r:
            i=min(r);coeff=r[i]
            if i not in self.pivots:
                inv=pow(coeff,-1,PRIME)
                self.pivots[i]={k:(v*inv)%PRIME for k,v in r.items()}
                return
            for k,v in self.pivots[i].items():
                new=(r.get(k,0)-coeff*v)%PRIME
                if new:r[k]=new
                else:r.pop(k,None)


def basis(block,point):
    lo,hi=block
    pol=[]
    for k in range(3):
        t=(point[k]-lo[k])/(hi[k]-lo[k])
        assert 0<=t<=1
        pol.append([(1-t)**2,2*t*(1-t),t*t])
    return [(9*i+3*j+k,pol[0][i]*pol[1][j]*pol[2][k]) for i,j,k in product(range(3),repeat=3)
            if pol[0][i] and pol[1][j] and pol[2][k]]


def verify(d):
    blocks=d["analysis_blocks"]
    lift=d["component_lifts"][0]
    assert all(x==lift for x in d["component_lifts"])
    roots=lift["free_roots"]
    lines=[{j:F(n,k) for j,n,k in line["terms"]} for line in lift["lines"]]
    assert all(F(*r["offset"])==0 for r in lift["lines"])
    assert all(lines[id]=={i:F(1)} for i,id in enumerate(roots))
    ids=[z for b in blocks for z in b["dofs"]]
    assert set(roots)<=set(ids)
    rank=Rank();rows=0
    def check(row):
        nonlocal rows
        residual={}
        for k,v in row.items():
            for root,coef in lines[ids[k]].items():combine(residual,root,coef*v)
        if residual:raise AssertionError("Candidate does not satisfy a physical continuity equation")
        rank.add(row);rows+=1
    boxes=[];planes=defaultdict(lambda:[[],[]])
    for i,b in enumerate(blocks):
        pts=[tuple(F(v,p[3]) for v in p[:3]) for p in b["controls"]]
        lo,hi=pts[0],pts[-1]
        assert all(hi[k]>lo[k] for k in range(3))
        assert all(p==tuple(lo[k]+(hi[k]-lo[k])*idx[k]/2 for k in range(3)) for p,idx in zip(pts,product(range(3),repeat=3)))
        assert len({p[3] for p in b["controls"]})==1
        boxes.append((lo,hi))
        for axis in range(3):
            planes[(axis,lo[axis])][0].append(i);planes[(axis,hi[axis])][1].append(i)
    patches=0
    for (axis,coord),(upper_side,lower_side) in planes.items():
        other=[k for k in range(3) if k!=axis]
        for a,b in product(upper_side,lower_side):
            start=[max(boxes[a][0][k],boxes[b][0][k]) for k in other]
            end=[min(boxes[a][1][k],boxes[b][1][k]) for k in other]
            if any(x>=y for x,y in zip(start,end)):continue
            patches+=1
            for i,j in product(range(3),repeat=2):
                pt=[coord]*3
                for k,ix in enumerate((i,j)):pt[other[k]]=start[k]+(end[k]-start[k])*ix/2
                row={}
                for side,sign in ((a,1),(b,-1)):
                    for loc,weight in basis(boxes[side],pt):combine(row,27*side+loc,sign*weight)
                check(row)
    xmin=min(lo[0] for lo,_ in boxes)
    for i,(lo,_) in enumerate(boxes):
        if lo[0]==xmin:
            for j in range(9):check({27*i+j:F(1)})
    rank_mod=len(rank.pivots)
    if len(ids)-rank_mod!=len(roots):
        raise AssertionError(f"Physical kernel dimension bound {len(ids)-rank_mod}, candidate {len(roots)}")
    return dict(blocks=len(blocks),dup_local_dofs=len(ids),mounted_scalar_dofs=len(roots),physical_face_patches=patches,
                tested_equations=rows,prime=PRIME,modular_rank_lower_bound=rank_mod,
                rational_kernel_dimension=len(roots),all_rational_lift_products_zero=True,
                argument="Independent identity rows give full rank of the lift. Exact A*C=0 and the nonzero modular pivot minor give matching rational rank bounds.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli",type=Path);p.add_argument("regression",type=Path);p.add_argument("work",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();start=time.monotonic();data=json.loads(a.regression.read_text());a.work.mkdir(exist_ok=True,parents=True)
    recipe,model,job,result=[a.work/(k+".json") for k in ("recipe","model","job","result")]
    recipe.write_text(json.dumps(data["recipe"]));job.write_text(json.dumps(data["job"]))
    for cmd in (["build",recipe,model],["solve-job",model,job,result]):
        subprocess.run([str(a.cli.resolve()),*cmd],check=True)
    info=verify(json.loads(result.read_text()))
    info.update(schema="cross_edge_exact_physical_rank_v1",elapsed_s=time.monotonic()-start,old_diagnostic=data["old_diagnostic"])
    a.output.write_text(json.dumps(info,indent=2)+"\n");print(info)
