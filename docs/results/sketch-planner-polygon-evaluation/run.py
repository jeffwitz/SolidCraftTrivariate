#!/usr/bin/env python3
"""Run paired polygon cases and independently audit serialized model sections."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import time

CORPUS = Path(__file__).with_name("cases.json")
HERE = Path(__file__).resolve().parent
ORACLE_PATH = next((root / "tests/references/pocket_partition_comparison.py"
                    for root in HERE.parents
                    if (root / "tests/references/pocket_partition_comparison.py").is_file()), None)
if ORACLE_PATH is None:
    raise RuntimeError("run.py must live under the repository tree to locate its independent oracle")
spec = importlib.util.spec_from_file_location("section_oracle", ORACLE_PATH)
oracle = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = oracle
spec.loader.exec_module(oracle)


def cond5(model):
    """Per-solid max 2-norm Jacobian condition on the closed 5^3 grid."""
    import numpy as np
    grid = (0.0, 0.25, 0.5, 0.75, 1.0)
    rows = []
    for block in model["blocks"]:
        worst = 0.0
        for u in grid:
            for v in grid:
                for w in grid:
                    sv = np.linalg.svd(oracle.jacobian_at(block, (u, v, w)), compute_uv=False)
                    worst = max(worst, float("inf") if sv[-1] <= 0 else float(sv[0] / sv[-1]))
        rows.append({"block_id": block["id"], "max_condition_2": worst})
    return rows


def checked_run(command, timeout):
    t0 = time.monotonic()
    try:
        p = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        return {"returncode": p.returncode, "stdout": p.stdout, "stderr": p.stderr,
                "elapsed_s": time.monotonic() - t0}
    except subprocess.TimeoutExpired as exc:
        return {"returncode": None, "stdout": "", "stderr": str(exc),
                "elapsed_s": time.monotonic() - t0, "timed_out": True}


def exact_recipe_check(recipe, family=""):
    """Exact integer predicates validate the preregistered input geometry."""
    loops = [recipe["outer"], *recipe.get("cuts", [])]
    def orient(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    def on(p, a, b):
        return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
            and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
    def cross(a, b, c, d):
        return orient(a,b,c)*orient(a,b,d)<0 and orient(c,d,a)*orient(c,d,b)<0
    def inside(p, poly):
        yes = False
        for a,b in zip(poly, poly[1:]+poly[:1]):
            if on(p,a,b): return 0
            if (a[1]>p[1]) != (b[1]>p[1]):
                if (a[0] + F(p[1]-a[1],b[1]-a[1])*(b[0]-a[0])) > p[0]: yes = not yes
        return 1 if yes else -1
    parsed=[]
    for loop in loops:
        pts=[tuple(map(int,p)) for p in loop["vertices"]]
        if len(set(pts)) != len(pts) or len(pts)<3:
            return False, "duplicate or underspecified loop"
        area=sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(pts,pts[1:]+pts[:1]))
        if area == 0: return False, "zero-area loop"
        edges=list(zip(pts,pts[1:]+pts[:1]))
        for i,(a,b) in enumerate(edges):
            for j,(c,d) in enumerate(edges):
                if j<=i or j==i+1 or (i==0 and j==len(edges)-1): continue
                if cross(a,b,c,d) or on(a,c,d) or on(b,c,d) or on(c,a,b) or on(d,a,b):
                    return False, "self-intersecting loop"
        parsed.append((pts,edges))
    def strictly_convex(pts):
        turns=[orient(pts[i],pts[(i+1)%len(pts)],pts[(i+2)%len(pts)])
               for i in range(len(pts))]
        return all(t>0 for t in turns) or all(t<0 for t in turns)
    if "convex_" in family or family.startswith("seeded_convex_"):
        if not all(strictly_convex(points) for points, _ in parsed):
            return False, "case labeled convex contains a non-convex contour"
    if family == "concave_stock_or_cut" and not any(
            not strictly_convex(points) for points, _ in parsed):
        return False, "case labeled concave has no concave contour"
    outer=parsed[0][0]
    for i,(pts,edges) in enumerate(parsed[1:],1):
        if inside(pts[0],outer)!=1: return False, "cut not strictly inside stock"
        if any(cross(a,b,c,d) or on(a,c,d) or on(b,c,d) or on(c,a,b) or on(d,a,b)
               for a,b in edges for c,d in parsed[0][1]): return False, "cut touches stock"
        for j in range(1,i):
            other=parsed[j][0]
            if inside(pts[0],other)>=0 or inside(other[0],pts)>=0: return False, "nested or overlapping cuts"
            if any(cross(a,b,c,d) or on(a,c,d) or on(b,c,d) or on(c,a,b) or on(d,a,b)
                   for a,b in edges for c,d in parsed[j][1]): return False, "cuts touch"
        z0,z1=loops[i]["z_interval_mm"]
        if not 0 <= z0 < z1 <= recipe["height_mm"]: return False, "invalid z interval"
    return True, None


def eval_one(cli, item, work, timeout):
    recipe=item["recipe"]
    ok, why=exact_recipe_check(recipe,item["family"])
    if not ok: raise RuntimeError(f"frozen admissible input invalid {item['id']}: {why}")
    rp=work/f"{item['id']}.recipe.json"; mp=work/f"{item['id']}.model.json"
    rp.write_text(json.dumps(recipe,separators=(",",":"),allow_nan=False)+"\n")
    build=checked_run([str(cli),"build",str(rp),str(mp)],timeout)
    row={"id":item["id"],"family":item["family"],"build_elapsed_s":build["elapsed_s"],
         "build_stdout":build["stdout"].strip(),"build_stderr":build["stderr"].strip(),
         "admitted":False,"geometry_pass":False,"construction_method":None,
         "block_count":None,"condition_grid_5cube_per_block":None}
    if build.get("timed_out"):
        row["outcome"]="timeout"; return row
    if build["returncode"] != 0 or not mp.exists():
        row["outcome"]="rejected"; return row
    try:
        report=json.loads(build["stdout"]); model=json.loads(mp.read_text())
    except Exception as exc:
        row.update(outcome="invalid_model", parse_error=str(exc)); return row
    row["admitted"]=True
    row["block_count"]=len(model.get("blocks",[]))
    row["construction_method"]=(model.get("partition_plan") or {}).get("method")
    row["reported_method"]=report.get("method")
    audit=checked_run([str(cli),"audit",str(mp)],timeout)
    row["audit_elapsed_s"]=audit["elapsed_s"]
    if audit.get("timed_out"):
        row["outcome"]="audit_timeout"; return row
    if audit["returncode"] != 0:
        row["outcome"]="audit_failed"; return row
    try: ar=json.loads(audit["stdout"])
    except Exception as exc:
        row.update(outcome="audit_invalid",audit_parse_error=str(exc)); return row
    row["closed_domain_certificate"]={k:ar.get(k) for k in
        ("status","method","determinant_lower_bound_mm3","hausdorff_error_upper_bound_mm")}
    if ar.get("status") != "valid":
        row["outcome"]="certificate_not_valid"; return row
    levels=sorted({F(0),F(recipe["height_mm"]),
                   *(F(z) for c in recipe.get("cuts",[]) for z in c["z_interval_mm"])})
    layer_rows=[]; all_pass=True; lower_total=0
    for z0,z1 in zip(levels,levels[1:]):
        layer=oracle.audit_layer(model,recipe,z0,z1)
        zmid=(F(z0)+F(z1))/2
        active=[c for c in recipe.get("cuts",[]) if c["z_interval_mm"][0] < zmid < c["z_interval_mm"][1]]
        expected_chi=1-len(active)
        actual_points=set()
        for _,pts,_ in oracle.layer_cells(model,recipe,z0,z1): actual_points.update(pts)
        loops=[recipe["outer"]["vertices"]]+[c["vertices"] for c in active]
        boundary_count=sum(1 for p in actual_points if any(
            oracle.on_segment(p,tuple(a),tuple(b)) for loop in loops
            for a,b in zip(loop,loop[1:]+loop[:1])))
        interior_count=len(actual_points)-boundary_count
        conditional=(F(boundary_count,2)+interior_count-expected_chi)
        mandatory=sum(len(loop) for loop in loops)
        input_bound=math.ceil(mandatory/2)-expected_chi
        exact_volume=oracle.exact_volume(recipe)
        layer["active_hole_count"]=len(active)
        layer["quadrilateral_lower_bounds"]={
            "assumption":"one connected conforming quadrilateral section; all input contour corners mandatory",
            "input_vertex_count":mandatory,
            "input_only_bound_blocks":input_bound,
            "actual_boundary_vertex_count":boundary_count,
            "actual_interior_vertex_count":interior_count,
            "with_inserted_host_vertices_bound_blocks":str(conditional),
            "observed_quads":layer["cells"],
            "attains_conditional_bound":conditional==layer["cells"]}
        lower_total+=input_bound
        all_pass &= (layer["euler"]==expected_chi and layer["area_gap_exact_mm2"]==[0,1]
                     and layer["contained_in_expected_material"]
                     and layer["pairwise_overlap_pairs"]==0
                     and conditional==layer["cells"])
        layer_rows.append(layer)
    volume=oracle.exact_volume(recipe)
    exact_report=report.get("volume_exact_mm3")
    volume_ok=exact_report is not None and F(*map(int,exact_report))==volume
    all_pass &= volume_ok
    row["volume_exact_expected_mm3"]=[volume.numerator,volume.denominator]
    row["volume_exact_cli_mm3"]=exact_report
    row["volume_match"]=volume_ok
    row["input_vertex_lower_bound_total_blocks"]=lower_total
    row["layers"]=layer_rows
    row["geometry_pass"]=bool(all_pass)
    row["outcome"]="valid_model" if all_pass else "invalid_geometry"
    if all_pass:
        row["condition_grid_5cube_per_block"]=cond5(model)
    return row


def summarize(rows):
    out={}
    for family in sorted({r["family"] for r in rows}):
        subset=[r for r in rows if r["family"]==family]
        out[family]={"denominator":len(subset),
            "admitted":sum(r["admitted"] for r in subset),
            "geometric_pass":sum(r["geometry_pass"] for r in subset),
            "rejected":sum(r["outcome"]=="rejected" for r in subset),
            "timeouts":sum("timeout" in r["outcome"] for r in subset),
            "invalid_models":sum(r["outcome"]=="invalid_geometry" for r in subset),
            "methods":dict(Counter(r["construction_method"] for r in subset if r["admitted"]))}
    return out


def condition_summary(rows):
    values=[]
    for r in rows:
        if r["geometry_pass"]:
            values.append(max(x["max_condition_2"] for x in r["condition_grid_5cube_per_block"]))
    return ({"n":len(values),"median":statistics.median(values),
             "p95_nearest_rank":sorted(values)[math.ceil(.95*len(values))-1],"worst":max(values)}
            if values else {"n":0,"median":None,"p95_nearest_rank":None,"worst":None})


def run(cli, tag, work, timeout):
    cases=json.loads(CORPUS.read_text())["cases"]
    work.mkdir(parents=True,exist_ok=True)
    sha=hashlib.sha256(cli.read_bytes()).hexdigest()
    start=time.monotonic(); rows=[]
    for item in cases:
        rows.append(eval_one(cli,item,work,timeout))
        print(tag,item["id"],rows[-1]["outcome"],rows[-1]["block_count"],flush=True)
    if hashlib.sha256(cli.read_bytes()).hexdigest()!=sha:
        raise RuntimeError("executable changed during evaluation")
    return {"tag":tag,"cli_sha256":sha,"elapsed_s":time.monotonic()-start,
            "cases":rows,"summary":summarize(rows),"condition_summary":condition_summary(rows)}


def paired(base, cand):
    b={r["id"]:r for r in base["cases"]}; c={r["id"]:r for r in cand["cases"]}
    rows=[]
    for ident in b.keys() & c.keys():
        x,y=b[ident],c[ident]
        bx=x["geometry_pass"]; cy=y["geometry_pass"]
        rows.append({"id":ident,"family":x["family"],"baseline_pass":bx,"candidate_pass":cy,
            "paired_outcome":"common_success" if bx and cy else
                "newly_admitted" if not bx and cy else
                "regression" if bx and not cy else "common_failure",
            "baseline_blocks":x["block_count"],"candidate_blocks":y["block_count"],
            "block_delta_candidate_minus_baseline":(y["block_count"]-x["block_count"]
                if x["block_count"] is not None and y["block_count"] is not None else None),
            "baseline_condition_max":max((q["max_condition_2"] for q in (x["condition_grid_5cube_per_block"] or [])),default=None),
            "candidate_condition_max":max((q["max_condition_2"] for q in (y["condition_grid_5cube_per_block"] or [])),default=None)})
    return rows


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--cli",type=Path,required=True); ap.add_argument("--tag",required=True)
    ap.add_argument("--work",type=Path,required=True); ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--timeout",type=float,default=60)
    a=ap.parse_args(); cli=a.cli.resolve()
    result=run(cli,a.tag,a.work.resolve(),a.timeout)
    result["corpus_sha256"]=hashlib.sha256(CORPUS.read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
