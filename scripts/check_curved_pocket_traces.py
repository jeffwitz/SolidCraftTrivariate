#!/usr/bin/env python3
"""Independent exact trace-rank audit for the curved blind-pocket section.

The nine published volume maps reduce to five distinct biquadratic XY section
charts and eight oriented interior edges. This script builds the p=3,4
endpoint-Hermite field with C1 interior knots (r=1), identifies oriented value
and tangent traces, and evaluates physical transverse-gradient matching on each
edge. The rational gradient identity has numerator degree at most p+5, so
p+6 distinct rational samples per span, including endpoints, determine the identity.
The reported ranks use exact sparse Gaussian elimination over Q; two prime-field
calculations provide independent rank crosschecks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

PRIMES = (1_000_003, 1_000_033)
SIDES = ("v0", "u1", "v1", "u0")


def add(row: dict[int, F], col: int, value: F) -> None:
    if value:
        row[col] = row.get(col, F(0)) + value
        if row[col] == 0:
            del row[col]


def mod_rank(rows: list[dict[int, F]], prime: int) -> int:
    basis: dict[int, dict[int, int]] = {}
    for source in rows:
        row = {k: (v.numerator * pow(v.denominator, prime - 2, prime)) % prime
               for k, v in source.items() if v}
        row = {k: v for k, v in row.items() if v}
        while row:
            pivot = min(row)
            old = basis.get(pivot)
            if old is None:
                scale = pow(row[pivot], prime - 2, prime)
                basis[pivot] = {k: v * scale % prime for k, v in row.items() if v}
                break
            factor = row[pivot]
            for k, v in old.items():
                nv = (row.get(k, 0) - factor * v) % prime
                if nv:
                    row[k] = nv
                else:
                    row.pop(k, None)
    return len(basis)


def rational_rank(rows: list[dict[int, F]]) -> int:
    """Exact sparse Gaussian elimination over Q."""
    basis: dict[int, dict[int, F]] = {}
    for source in rows:
        row = {k: v for k, v in source.items() if v}
        while row:
            pivot = min(row)
            old = basis.get(pivot)
            if old is None:
                scale = row[pivot]
                basis[pivot] = {k: v / scale for k, v in row.items() if v}
                break
            factor = row[pivot]
            for k, v in old.items():
                value = row.get(k, F(0)) - factor * v
                if value:
                    row[k] = value
                else:
                    row.pop(k, None)
    return len(basis)


class RationalRankAccumulator:
    """Incrementally report exact row-rank increases for ordered seam groups."""
    def __init__(self):
        self.basis: dict[int, dict[int, F]] = {}

    def add(self, source: dict[int, F]) -> bool:
        row = {k: v for k, v in source.items() if v}
        while row:
            pivot = min(row)
            old = self.basis.get(pivot)
            if old is None:
                scale = row[pivot]
                self.basis[pivot] = {k: v / scale for k, v in row.items() if v}
                return True
            factor = row[pivot]
            for k, v in old.items():
                value = row.get(k, F(0)) - factor * v
                if value:
                    row[k] = value
                else:
                    row.pop(k, None)
        return False


def choose(n: int, k: int) -> int:
    return math.comb(n, k)


def bernstein(degree: int, t: F) -> list[F]:
    return [F(choose(degree, i)) * t**i * (1-t)**(degree-i)
            for i in range(degree + 1)]


def bernstein_derivative(degree: int, t: F) -> list[F]:
    low = bernstein(degree - 1, t)
    return [degree * ((low[i-1] if i else F(0)) -
                      (low[i] if i < degree else F(0)))
            for i in range(degree + 1)]


def descriptors(p: int, n: int) -> list[tuple[str, int]]:
    out = []
    for k in range(n + 1):
        out.extend((("value", k), ("slope", k)))
        if p == 4 and k < n:
            out.append(("middle", k))
    return out


def hermite_spans(p: int, n: int) -> list[list[dict[int, F]]]:
    ds = descriptors(p, n); lookup = {d: i for i, d in enumerate(ds)}; out = []
    for k in range(n):
        rows = [dict() for _ in range(p + 1)]
        rows[0][lookup[("value", k)]] = F(1)
        rows[1][lookup[("value", k)]] = F(1)
        rows[1][lookup[("slope", k)]] = F(1, p)
        if p == 3:
            rows[2][lookup[("value", k+1)]] = F(1)
            rows[2][lookup[("slope", k+1)]] = F(-1, p)
        else:
            rows[2][lookup[("middle", k)]] = F(1)
            rows[3][lookup[("value", k+1)]] = F(1)
            rows[3][lookup[("slope", k+1)]] = F(-1, p)
        rows[p][lookup[("value", k+1)]] = F(1)
        out.append(rows)
    return out


def load_charts(model: dict) -> tuple[list[dict], list[dict]]:
    charts = []
    for block in model["blocks"]:
        if not block["id"].endswith("_below"):
            continue
        if block["degrees"] != [2, 2, 1]:
            raise ValueError("expected the published (2,2,1) curved pocket maps")
        raw = block["controls"]; weights = {q[3] for q in raw}
        if len(weights) != 1 or not next(iter(weights)):
            raise ValueError(f"{block['id']}: nonconstant or zero weights")
        w = next(iter(weights))
        net = [[(F(raw[(i*3+j)*2][0], w), F(raw[(i*3+j)*2][1], w))
                for j in range(3)] for i in range(3)]
        charts.append({"id": block["id"], "net": net})
    if len(charts) != 5:
        raise ValueError(f"expected five section maps, found {len(charts)}")
    by_id = {b["id"]: b for b in model["blocks"]}
    for chart in charts:
        above_id = chart["id"].removesuffix("_below") + "_beside"
        above = by_id.get(above_id)
        if above is None:
            if chart["id"].endswith("footprint_below"):
                continue
            raise ValueError(f"missing thickness continuation {above_id}")
        lower_xy = [(q[0], q[1]) for q in by_id[chart["id"]]["controls"][::2]]
        upper_xy = [(q[0], q[1]) for q in above["controls"][::2]]
        if lower_xy != upper_xy:
            raise ValueError(f"{chart['id']}: XY section map changes with height")

    def side_controls(chart: dict, side: str) -> list[tuple[F, F]]:
        net = chart["net"]
        if side == "v0": return [net[i][0] for i in range(3)]
        if side == "u1": return [net[2][j] for j in range(3)]
        if side == "v1": return [net[i][2] for i in range(2, -1, -1)]
        return [net[0][j] for j in range(2, -1, -1)]

    edge_map: dict[tuple, list[dict]] = {}
    for ci, chart in enumerate(charts):
        for side in SIDES:
            controls = side_controls(chart, side)
            key = tuple(sorted((controls[0], controls[-1])))
            reverse = controls[0] != key[0]
            edge_map.setdefault(key, []).append({"chart": ci, "side": side,
                                                "reverse": reverse,
                                                "controls": controls})
    edges = []
    for key, incidences in edge_map.items():
        if len(incidences) == 2:
            canonical = [inc["controls"] if not inc["reverse"] else list(reversed(inc["controls"]))
                         for inc in incidences]
            if canonical[0] != canonical[1]:
                raise ValueError(f"nonmatching curved section seam {key}")
            edges.append({"key": key, "sides": incidences})
        elif len(incidences) != 1:
            raise ValueError(f"nonmanifold section edge {key}")
    if len(edges) != 8:
        raise ValueError(f"expected eight interior section seams, found {len(edges)}")
    return charts, edges


def edge_axes(side: str) -> tuple[int, int, int]:
    if side == "v0": return 1, 0, 0
    if side == "u1": return 0, 1, 1
    if side == "v1": return 1, 1, 0
    return 0, 0, 1


def map_geometry(chart: dict, u: F, v: F):
    bu, bv = bernstein(2, u), bernstein(2, v)
    du, dv = bernstein_derivative(2, u), bernstein_derivative(2, v)
    x = y = xu = yu = xv = yv = F(0)
    for i in range(3):
        for j in range(3):
            px, py = chart["net"][i][j]
            x += bu[i]*bv[j]*px; y += bu[i]*bv[j]*py
            xu += du[i]*bv[j]*px; yu += du[i]*bv[j]*py
            xv += bu[i]*dv[j]*px; yv += bu[i]*dv[j]*py
    return (x, y), (xu, yu), (xv, yv)


def edge_endpoints(chart: dict, axis: int, endpoint: int):
    net = chart["net"]
    return [net[0][2*endpoint], net[2][2*endpoint]] if axis == 0 else [net[2*endpoint][0], net[2*endpoint][2]]


def node_key(chart: dict, ci: int, i: int, j: int, n: int, ds):
    du, dv = ds[i], ds[j]
    ub = du[1]//n if du[0] == "value" and du[1] in (0,n) else None
    vb = dv[1]//n if dv[0] == "value" and dv[1] in (0,n) else None
    if ub is not None and vb is not None:
        return ("vertex", chart["net"][ub*2][vb*2]), F(1)
    if ub is not None or vb is not None:
        axis = 1 if ub is not None else 0
        endpoint = ub if ub is not None else vb
        tangential = dv if ub is not None else du
        points = edge_endpoints(chart, axis, endpoint)
        if points[0] <= points[1]:
            return ("edge", tuple(points), tangential), F(1)
        kind, knot = tangential
        mapped = (kind, n-knot if kind in ("value", "slope") else n-1-knot)
        return ("edge", tuple(reversed(points)), mapped), F(-1) if kind == "slope" else F(1)
    return ("interior", ci, i, j), F(1)


def full_ids(charts: list[dict], p: int, n: int):
    ds=descriptors(p,n); dim=len(ds); keys=[]
    for ci,c in enumerate(charts):
        keys.append([node_key(c,ci,i,j,n,ds) for i in range(dim) for j in range(dim)])
    unique=sorted({key for row in keys for key,_ in row},key=repr); index={key:i for i,key in enumerate(unique)}
    return [[(index[key],sign) for key,sign in row] for row in keys],len(unique),dim


def trace_local_index(side: str, reverse: bool, d: tuple[str,int], p: int, n: int, dim: int) -> int:
    ds=descriptors(p,n); kind,k=d
    axis_reverse=reverse ^ (side in ("v1","u0"))
    local=n-k if axis_reverse and kind in ("value","slope") else n-1-k if axis_reverse and kind=="middle" else k
    ti=ds.index((kind,local)); normal,fixed,tangent=edge_axes(side)
    ni=ds.index(("value",0 if fixed==0 else n))
    return ni*dim+ti if normal==0 else ti*dim+ni


def trace_sign(reverse: bool, d: tuple[str,int]) -> F:
    return F(-1) if reverse and d[0]=="slope" else F(1)


def physical_gradient_row(charts, edge, p,n,k,tau,local_ids,dim,local_offsets=None):
    ds=descriptors(p,n); spans=hermite_spans(p,n); x0,y0=edge["key"][0]; x1,y1=edge["key"][1]; comp=0 if y0!=y1 else 1
    dets=[]; terms=[]
    for pos,inc in enumerate(edge["sides"]):
        ci,side,reverse=inc["chart"],inc["side"],inc["reverse"]; normal,fixed,tangent=edge_axes(side)
        axis_reverse=reverse ^ (side in ("v1","u0"))
        lk=n-1-k if axis_reverse else k; st=1-tau if axis_reverse else tau
        cell=[0,0]; cell[normal]=0 if fixed==0 else n-1; cell[tangent]=lk
        s=[F(0),F(0)]; s[normal]=F(0) if fixed==0 else F(1); s[tangent]=st
        uv=[(F(cell[q])+s[q])/n for q in range(2)]
        _,(xu,yu),(xv,yv)=map_geometry(charts[ci],*uv); det=xu*yv-xv*yu
        if det==0: raise ValueError("singular chart Jacobian at edge unisolvence point")
        bu,bv=bernstein(p,s[0]),bernstein(p,s[1]); du,dv=bernstein_derivative(p,s[0]),bernstein_derivative(p,s[1])
        hu,hv=spans[cell[0]],spans[cell[1]]; term={}
        for ii in range(p+1):
            for jj in range(p+1):
                coeff=((yv*du[ii]*bv[jj]-yu*bu[ii]*dv[jj]) if comp==0 else
                       (-xv*du[ii]*bv[jj]+xu*bu[ii]*dv[jj]))*n
                for qi,ci_coeff in hu[ii].items():
                    for qj,cj_coeff in hv[jj].items():
                        local_index=qi*dim+qj; gid,orient=local_ids[ci][local_index]
                        col=local_offsets[pos]+local_index if local_offsets else gid
                        add(term,col,coeff*ci_coeff*cj_coeff*orient)
        dets.append(det); terms.append(term)
    out={}
    for col,c in terms[0].items(): add(out,col,c*dets[1])
    for col,c in terms[1].items(): add(out,col,-c*dets[0])
    return out


def transverse_numerator_row(charts, incidence, component, p, n, k, tau, local_ids, dim):
    """One-side physical-gradient numerator row with det(J) left uncleared."""
    ds=descriptors(p,n); spans=hermite_spans(p,n)
    ci,side,reverse=incidence["chart"],incidence["side"],incidence["reverse"]
    normal,fixed,tangent=edge_axes(side); axis_reverse=reverse ^ (side in ("v1","u0"))
    cell=[0,0]; cell[normal]=0 if fixed==0 else n-1; cell[tangent]=n-1-k if axis_reverse else k
    st=1-tau if axis_reverse else tau; s=[F(0),F(0)]; s[normal]=F(0) if fixed==0 else F(1); s[tangent]=st
    uv=[(F(cell[q])+s[q])/n for q in range(2)]
    _,(xu,yu),(xv,yv)=map_geometry(charts[ci],*uv)
    bu,bv=bernstein(p,s[0]),bernstein(p,s[1]); du,dv=bernstein_derivative(p,s[0]),bernstein_derivative(p,s[1])
    hu,hv=spans[cell[0]],spans[cell[1]]; row={}
    for ii in range(p+1):
        for jj in range(p+1):
            coeff=((yv*du[ii]*bv[jj]-yu*bu[ii]*dv[jj]) if component==0 else
                   (-xv*du[ii]*bv[jj]+xu*bu[ii]*dv[jj]))*n
            for qi,ci_coeff in hu[ii].items():
                for qj,cj_coeff in hv[jj].items():
                    gid,orient=local_ids[ci][qi*dim+qj]
                    add(row,gid,coeff*ci_coeff*cj_coeff*orient)
    return row


def value_projection_rows(edge, p, n, local_ids):
    """Bernstein trace coefficients, ordered by canonical edge span."""
    dim=len(descriptors(p,n)); spans=hermite_spans(p,n); incidence=edge["sides"][0]
    ci,side,reverse=incidence["chart"],incidence["side"],incidence["reverse"]
    normal,fixed,tangent=edge_axes(side); axis_reverse=reverse ^ (side in ("v1","u0")); out=[]
    for k in range(n):
        cell=[0,0]; cell[normal]=0 if fixed==0 else n-1; cell[tangent]=n-1-k if axis_reverse else k
        hn,ht=spans[cell[normal]],spans[cell[tangent]]; row=[]
        for q in range(p+1):
            tq=p-q if axis_reverse else q; terms={}
            for a,ca in hn[0 if fixed==0 else p].items():
                for b,cb in ht[tq].items():
                    i,j=(a,b) if normal==0 else (b,a)
                    gid,orient=local_ids[ci][i*dim+j]
                    add(terms,gid,ca*cb*orient)
            row.append(terms)
        out.extend(row)
    return out


def projection_ranks(A_rows, projections):
    base=RationalRankAccumulator()
    for row in A_rows: base.add(row)
    base_rank=len(base.basis); out={}
    for name,rows in projections.items():
        joined=RationalRankAccumulator()
        joined.basis={pivot:dict(row) for pivot,row in base.basis.items()}
        for row in rows: joined.add(row)
        out[name]=len(joined.basis)-base_rank
    return out


def rank_summary(rows):
    ranks={str(prime):mod_rank(rows,prime) for prime in PRIMES}
    exact=rational_rank(rows)
    if exact < max(ranks.values()): raise AssertionError("modular rank exceeds rational rank")
    return {"modular_ranks":ranks,"rational_rank_exact":exact,
            "rows":len(rows)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model",type=Path,default=Path("docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"))
    ap.add_argument("--output",type=Path); args=ap.parse_args()
    raw=args.model.read_bytes(); model=json.loads(raw); charts,edges=load_charts(model); results=[]
    for p in (3,4):
        for n in (1,2,4,8):
            ds=descriptors(p,n); dim=len(ds); ids,whole_dofs,_=full_ids(charts,p,n)
            edge_results=[]; edge_data=[]
            for ei,e in enumerate(edges):
                c0=[]
                a,b=e["sides"]
                for d in ds:
                    ia=trace_local_index(a["side"],a["reverse"],d,p,n,dim)
                    ib=trace_local_index(b["side"],b["reverse"],d,p,n,dim)
                    ar=a["reverse"] ^ (a["side"] in ("v1","u0")); br=b["reverse"] ^ (b["side"] in ("v1","u0"))
                    c0.append({ia:trace_sign(ar,d),dim*dim+ib:-trace_sign(br,d)})
                pair_ids=[[] for _ in charts]
                pair_ids[a["chart"]]=[(q,F(1)) for q in range(dim*dim)]
                pair_ids[b["chart"]]=[(dim*dim+q,F(1)) for q in range(dim*dim)]
                pair_jets=[]; full_jets=[]
                for k in range(n):
                    for r in range(p+6):
                        tau=F(r,p+5)
                        pair_jets.append(physical_gradient_row(charts,e,p,n,k,tau,pair_ids,dim,[0,dim*dim]))
                        full_jets.append(physical_gradient_row(charts,e,p,n,k,tau,ids,dim))
                pair_c0_rank=rational_rank(c0)
                pair_rank=rational_rank(c0+pair_jets)-pair_c0_rank
                edge_results.append({"seam_index":ei,"sides":[[charts[s["chart"]]["id"],s["side"],s["reverse"]] for s in e["sides"]],
                    "value_and_tangent_trace_dimension_per_side":dim,
                    "two_patch_value_and_tangent_constraint_rank_exact":pair_c0_rank,
                    "transverse_rows_including_endpoints":len(pair_jets),
                    "two_patch_transverse_rank_after_value_tangent_exact":pair_rank,
                    "two_patch_scalar_dofs_after_value_tangent":2*dim*dim-pair_c0_rank,
                    "two_patch_scalar_dofs_after_physical_C1":2*dim*dim-pair_c0_rank-pair_rank})
                edge_data.append({"edge":e,"pair_c0":c0,"pair_jets":pair_jets,
                                  "pair_ids":pair_ids,"full_jets":full_jets})
            all_rows=[row for dta in edge_data for row in dta["full_jets"]]
            for ei,dta in enumerate(edge_data):
                e=dta["edge"]; a=e["sides"][0]
                edge_rows=dta["full_jets"]
                other_rows=[row for j,dj in enumerate(edge_data) if j!=ei for row in dj["full_jets"]]
                component=0 if e["key"][0][1]!=e["key"][1][1] else 1
                pval=value_projection_rows(e,p,n,ids); ptrans=[]
                pair_val=value_projection_rows(e,p,n,dta["pair_ids"]); pair_trans=[]
                for k in range(n):
                    for r in range(p+6):
                        tau=F(r,p+5)
                        ptrans.append(transverse_numerator_row(charts,a,component,p,n,k,tau,ids,dim))
                        pair_trans.append(transverse_numerator_row(charts,a,component,p,n,k,tau,dta["pair_ids"],dim))
                projections={"value_bernstein_coefficients":pval,
                             "transverse_numerator_samples":ptrans,
                             "joint_value_and_transverse":pval+ptrans}
                pair_projections={"value_bernstein_coefficients":pair_val,
                                  "transverse_numerator_samples":pair_trans,
                                  "joint_value_and_transverse":pair_val+pair_trans}
                edge_results[ei]["two_patch_observable_projection_ranks"]={
                    "after_value_tangent_gluing":projection_ranks(dta["pair_c0"],pair_projections),
                    "after_physical_C1_gluing":projection_ranks(dta["pair_c0"]+dta["pair_jets"],pair_projections)}
                edge_results[ei]["five_patch_observable_projection_ranks"]={
                    "after_value_tangent_before_C1":projection_ranks([],projections),
                    "after_local_seam_C1_before_other_seams":projection_ranks(edge_rows,projections),
                    "after_other_seven_seams_before_local_C1":projection_ranks(other_rows,projections),
                    "after_all_eight_seams_C1":projection_ranks(all_rows,projections)}
            cumulative=RationalRankAccumulator()
            for ei,dta in enumerate(edge_data):
                inc=sum(cumulative.add(row) for row in dta["full_jets"])
                edge_results[ei]["five_patch_incremental_transverse_rank_exact"]=inc
                edge_results[ei]["five_patch_cumulative_transverse_rank_exact"]=len(cumulative.basis)
            results.append({"p":p,"spans_per_chart_axis":n,"internal_knot_regularity":"C1 (r=1)",
                "section_charts":5,"oriented_interior_seams":8,"dofs_per_chart":dim*dim,
                "whole_five_patch_dofs_before_gluing":5*dim*dim,
                "whole_five_patch_dofs_after_value_tangent_gluing":whole_dofs,
                "value_and_tangent_trace_dimension_per_seam":dim,
                "transverse_unisolvence_points_per_span":p+6,
                "whole_five_patch_transverse_rows":len(all_rows),
                "whole_five_patch_transverse_rank":rank_summary(all_rows),
                "whole_five_patch_scalar_dofs_after_physical_C1":whole_dofs-len(cumulative.basis),"seams":edge_results})
    result={"schema":"curved_pocket_physical_trace_ranks_v1","source_model_sha256":hashlib.sha256(raw).hexdigest(),
        "source_model":str(args.model),"source_volume_maps":9,"section_charts":5,"oriented_section_seams":8,
        "finite_field_primes":list(PRIMES),"method":"Exact Fraction biquadratic geometry and Hermite field rows; physical Cartesian gradient difference after clearing exact Jacobian determinants; p+6 rational unisolvence nodes per span; exact sparse Gaussian elimination over Q with prime-field crosschecks.",
        "scope":"Section-only scalar field, with all value/tangent vertices included. Projection ranks use rank([A;P])-rank(A) over Q and include endpoints and vertices. Transverse observables are first-side physical-gradient numerators with the known nonzero geometry determinant omitted. Axial factor, boundary clamp, 3D solve, and convergence are outside this audit.","results":results}
    data=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output: args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(data)
    else: print(data,end="")

if __name__=="__main__": main()
