"""Check signed cellular chains and compute Betti numbers by exact GF(2) elimination."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse
from j1_geometry import evaluate


def matrix(cells,nrows):
    rows,cols,values=[],[],[]
    for column,cell in enumerate(cells):
        assert cell["id"]==column
        terms=cell["boundary"]
        assert len({i for i,_ in terms})==len(terms)
        assert all(0<=i<nrows and abs(sign)==1 for i,sign in terms)
        for i,sign in terms:
            rows.append(i);cols.append(column);values.append(sign)
    return sparse.csc_matrix((values,(rows,cols)),shape=(nrows,len(cells)),dtype=np.int64)


def rank_mod_two(boundary):
    pivots={}
    for col in range(boundary.shape[1]):
        ids=boundary.indices[boundary.indptr[col]:boundary.indptr[col+1]]
        vector=sum(1<<int(i) for i in ids)
        while vector:
            i=vector.bit_length()-1
            if i not in pivots:
                pivots[i]=vector
                break
            vector^=pivots[i]
    return len(pivots)


def homology(matrices):
    for first,second in zip(matrices,matrices[1:]):
        assert (first@second).count_nonzero()==0
    dims=[matrices[0].shape[0],*(m.shape[1] for m in matrices)]
    ranks=[rank_mod_two(m.tocsc()) for m in matrices]
    betti=[a-b-c for a,b,c in zip(dims,[0,*ranks],[*ranks,0])]
    euler=sum((-1)**i*dim for i,dim in enumerate(dims))
    assert all(x>=0 for x in betti) and sum((-1)**i*b for i,b in enumerate(betti))==euler
    return dict(cell_counts=dims,boundary_map_ranks_mod_2=ranks,betti_mod_2=betti,euler=euler)


def trace_geometry(c,roots):
    denom=c["grid_denominator"]
    positions=[]
    max_vertex=max_edge=0.
    def pos(root,at):
        return np.array(evaluate(roots[root],np.asarray(at)/denom)[0])
    for i,v in enumerate(c["vertices"]):
        assert i==v["id"] and v["charts"]
        points=np.array([pos(a["root"],a["coordinates"]) for a in v["charts"]])
        positions.append(points[0]);max_vertex=max(max_vertex,float(np.max(np.abs(points-points[0]))))
    for edge in c["edges"]:
        target=np.array([positions[i] for i in edge["vertices"]])
        probes=[]
        for line in edge["provenance"]:
            a,b=np.array(line["start"]),np.array(line["end"])
            assert np.count_nonzero(a!=b)==1
            curve=np.array([pos(line["root"],(1-t)*a+t*b) for t in (0.,.25,.5,.75,1.)])
            if np.max(np.abs(curve[[0,-1]]-target))>1e-9:curve=curve[::-1]
            assert np.max(np.abs(curve[[0,-1]]-target))<1e-9
            probes.append(curve)
        max_edge=max(max_edge,float(np.max(np.abs(np.array(probes)-probes[0]))))
    assert max_vertex<1e-9 and max_edge<1e-9
    return dict(max_chart_vertex_difference_mm=max_vertex,max_chart_edge_difference_mm=max_edge)


def analyze(c,roots,genus):
    start=time.monotonic()
    v,e,f,k=[c[name] for name in ("vertices","edges","faces","cells")]
    de=matrix(e,len(v));df=matrix(f,len(e));dc=matrix(k,len(f))
    interior=set();boundary=set()
    for face,i in zip(f,range(len(f))):
        uses=dc.getrow(i)
        assert sorted(zip(uses.indices,uses.data))==sorted(map(tuple,face["incident_cells"]))
        if len(uses.indices)==1:
            assert face["kind"]=="exterior"
            boundary.add(i)
        else:
            assert face["kind"]=="interface" and set(uses.data)=={-1,1} and len(uses.indices)==2
            interior.add(i)
        # Each exported oriented face cycle must join head to tail.
        order=[e[j]["vertices"][::sign] for j,sign in face["boundary"]]
        assert all(a[1]==b[0] for a,b in zip(order,order[1:]+order[:1]))
    assert boundary==set(c["boundary_faces"])
    volume=homology([de,df,dc])
    boundary_faces=sorted(boundary)
    face_matrix=df[:,boundary_faces]
    bedges=np.flatnonzero(face_matrix.getnnz(axis=1))
    assert np.all(face_matrix.getnnz(axis=1)[bedges]==2)
    edge_matrix=de[:,bedges]
    bvertices=np.flatnonzero(edge_matrix.getnnz(axis=1))
    shell=homology([edge_matrix[bvertices,:],face_matrix[bedges,:]])
    assert volume["betti_mod_2"]==[1,genus,0,0] and shell["betti_mod_2"]==[1,2*genus,1]
    result=dict(material=volume,boundary=shell,maximum_edges_per_face=max(len(ff["boundary"]) for ff in f),
                maximum_faces_per_cell=max(len(cc["boundary"]) for cc in k),elapsed_s=time.monotonic()-start)
    if roots:
        result["chart_agreement"]=trace_geometry(c,roots)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("packages",type=Path,nargs="+")
    args=parser.parse_args();cases=[]
    for file in args.packages:
        data=json.loads(file.read_text())
        if data.get("schema")=="adaptive_oriented_cell_complex_v1":
            result=analyze(data,None,0)
        else:
            roots={b["id"]:b for b in data["source_model"]["blocks"]}
            result=analyze(data["adaptive_mesh"]["oriented_complex"],roots,2)
        cases.append(dict(case=file.stem,**result))
    args.output.write_text(json.dumps({"coefficient_field":"GF(2)","cases":cases},indent=2,allow_nan=False)+"\n")
