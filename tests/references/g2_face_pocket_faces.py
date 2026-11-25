#!/usr/bin/env python3
"""Integrate S3 V-Rep boundary faces and compare with the source oracle."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MODELS = {
    "interior_principal": ROOT / "build/s3_face_pocket.delivery.json",
    "open_stock_edge_witness": ROOT / "build/s3_face_pocket_open_side.delivery.json",
}
SOURCE_REPORT = ROOT / "docs/results/g2-face-pocket-geometry-20261125.json"
REPORT = ROOT / "docs/results/g2-face-pocket-faces.json"
SOURCE = runpy.run_path(str(ROOT / "tests/references/g2_face_pocket_geometry.py"))


def bernstein(p, x):
    return [math.comb(p, i)*x**i*(1-x)**(p-i) for i in range(p+1)]


def eval_tensor(net, degrees, uvw):
    values = [bernstein(p, x) for p, x in zip(degrees, uvw)]
    out = [0.0]*4
    for i, bu in enumerate(values[0]):
        for j, bv in enumerate(values[1]):
            for k, bw in enumerate(values[2]):
                c = net[(i*(degrees[1]+1)+j)*(degrees[2]+1)+k]
                q = bu*bv*bw
                for d in range(4):
                    out[d] += q*c[d]
    return out


def derivative_net(net, degrees, axis):
    dims = [p+1 for p in degrees]
    p = degrees[axis]
    ndims = dims[:]
    ndims[axis] -= 1
    out = []
    for i in range(ndims[0]):
        for j in range(ndims[1]):
            for k in range(ndims[2]):
                ijk = [i, j, k]
                nxt = ijk[:]
                nxt[axis] += 1
                flat = lambda a: (a[0]*dims[1]+a[1])*dims[2]+a[2]
                a, b = net[flat(ijk)], net[flat(nxt)]
                out.append([p*(b[d]-a[d]) for d in range(4)])
    nd = degrees[:]
    nd[axis] -= 1
    return out, nd


def homogeneous_eval_with_derivatives(net, degrees, uvw):
    h = eval_tensor(net, degrees, uvw)
    derivs = []
    for axis in range(3):
        dnet, dd = derivative_net(net, degrees, axis)
        derivs.append(eval_tensor(dnet, dd, uvw))
    w = h[3]
    xyz = np.array(h[:3])/w
    dxyz = [ (np.array(d[:3])*w-np.array(h[:3])*d[3])/(w*w) for d in derivs ]
    return xyz, dxyz


def rational_point_jet2(net, degrees, uvw, axis):
    h = eval_tensor(net, degrees, uvw)
    dnet, dd = derivative_net(net, degrees, axis)
    dh = eval_tensor(dnet, dd, uvw)
    ddnet, ddd = derivative_net(dnet, dd, axis)
    ddh = eval_tensor(ddnet, ddd, uvw)
    w = h[3]
    r = np.array(h[:3])/w
    rp = (np.array(dh[:3])*w-np.array(h[:3])*dh[3])/(w*w)
    rpp = (np.array(ddh[:3])*w-np.array(h[:3])*ddh[3])/(w*w) \
        - 2*(np.array(dh[:3])*dh[3])/(w*w) \
        + 2*np.array(h[:3])*(dh[3]**2)/(w**3)
    return r, rp, rpp


def source_point_jet2(t):
    degree = SOURCE["DEG"]
    knots = SOURCE["KNOTS"]
    controls = SOURCE["CONTROLS"]
    def derivative(cp, u, p):
        out=[]
        for i in range(len(cp)-1):
            den=knots[i+p+1]-knots[i+1]
            out.append([p*(cp[i+1][k]-cp[i][k])/den for k in range(3)])
        return out
    cp1 = derivative(controls,None,degree)
    u1 = knots[1:-1]
    cp2=[]
    for i in range(len(cp1)-1):
        den=u1[i+degree]-u1[i+1]
        cp2.append([(degree-1)*(cp1[i+1][k]-cp1[i][k])/den for k in range(3)])
    u2=u1[1:-1]
    h=SOURCE["homogeneous_de_boor"](degree,knots,controls,t)
    h1=SOURCE["homogeneous_de_boor"](degree-1,u1,cp1,t)
    h2=SOURCE["homogeneous_de_boor"](degree-2,u2,cp2,t)
    H,H1,H2=np.array(h),np.array(h1),np.array(h2)
    r=H[:2]/H[2]
    rp=(H1[:2]*H[2]-H[:2]*H1[2])/(H[2]**2)
    rpp=(H2[:2]*H[2]-H[:2]*H2[2])/(H[2]**2) \
        -2*H1[:2]*H1[2]/(H[2]**2)+2*H[:2]*H1[2]**2/(H[2]**3)
    return r,rp,rpp


def curve_edge_jet(block, fixed, s):
    uvw=[.5,.5,.5]
    for axis,value in fixed.items():
        uvw[axis]=value
    uvw[0]=s
    return rational_point_jet2(block["controls"],block["degrees"],uvw,0)


def curvature_vector(rp,rpp):
    speed2=float(np.dot(rp,rp))
    tangent=rp/math.sqrt(speed2)
    return (rpp-tangent*float(np.dot(tangent,rpp)))/speed2


def sampled_edge_comparison(block, t0, t1, z, fixed):
    max_position=max_tangent=max_curvature=0.0
    samples=[]
    for q in np.linspace(0.0,1.0,17):
        r,rp,rpp=curve_edge_jet(block,fixed,float(q))
        rs,rps,rpps=source_point_jet2(t0+(t1-t0)*float(q))
        rs=np.array([rs[0],rs[1],z])
        rps=np.array([rps[0],rps[1],0.0])
        rpps=np.array([rpps[0],rpps[1],0.0])
        dp=float(np.linalg.norm(r-rs))
        dt=float(np.dot(rp/np.linalg.norm(rp),rps/np.linalg.norm(rps)))
        dk=float(np.linalg.norm(curvature_vector(rp,rpp)-curvature_vector(rps,rpps)))
        max_position=max(max_position,dp)
        max_tangent=max(max_tangent,1.0-dt)
        max_curvature=max(max_curvature,dk)
        samples.append({"fraction":float(q),"position_error_mm":dp,"tangent_dot":dt,"curvature_vector_error_per_mm":dk})
    return {"sample_count":len(samples),"max_position_error_mm":max_position,
            "max_tangent_1_minus_dot":max_tangent,"max_curvature_vector_error_per_mm":max_curvature,
            "samples":samples}


def measure_trace(blocks, faces, name, trace):
    by_id={b["id"]:b for b in blocks}
    entries=[]
    for face in faces:
        block_id=face["block"]
        if face["feature"] != trace["feature"] or face["side"] != trace["side"]:
            continue
        if not trace["block_match"](block_id):
            continue
        b=by_id[block_id]
        edge_len=0.0
        nodes,weights=np.polynomial.legendre.leggauss(32)
        for x,w in zip(nodes,weights):
            _,rp,_=curve_edge_jet(b,trace["fixed"],float((x+1)/2))
            edge_len += (w/2)*np.linalg.norm(rp)
        comparison=sampled_edge_comparison(b,*trace["source_t"],trace["z"],trace["fixed"])
        entries.append({"block":block_id,"feature":face["feature"],"face_side":face["side"],
                        "source_parameter_interval":trace["source_t"],"integrated_edge_length_mm":edge_len,
                        "direct_source_jet_comparison":comparison})
    return entries


def find_trace_from_source(blocks, faces, feature, source_t, z):
    """Find a V-Rep boundary edge by matching its endpoints to source points."""
    by_id={b["id"]:b for b in blocks}
    target_a=np.array([*SOURCE["de_boor"](SOURCE["DEG"],SOURCE["KNOTS"],SOURCE["CONTROLS"],source_t[0]),z])
    target_b=np.array([*SOURCE["de_boor"](SOURCE["DEG"],SOURCE["KNOTS"],SOURCE["CONTROLS"],source_t[1]),z])
    candidates=[]
    for face in faces:
        if face["feature"]!=feature or face["side"]!="w1" or not face["block"].startswith("s3_i_"):
            continue
        block=by_id[face["block"]]
        a=curve_edge_jet(block,{1:1.0,2:1.0},0.0)[0]
        b=curve_edge_jet(block,{1:1.0,2:1.0},1.0)[0]
        direct=np.linalg.norm(a-target_a)+np.linalg.norm(b-target_b)
        reverse=np.linalg.norm(a-target_b)+np.linalg.norm(b-target_a)
        candidates.append((min(direct,reverse),face,block,reverse<direct))
    candidates.sort(key=lambda x:x[0])
    if not candidates or candidates[0][0]>2e-8:
        return []
    _,face,block,reversed_orientation=candidates[0]
    directed_t=source_t[::-1] if reversed_orientation else source_t
    comparison=sampled_edge_comparison(block,*directed_t,z,{1:1.0,2:1.0})
    nodes,weights=np.polynomial.legendre.leggauss(32)
    length=sum((w/2)*np.linalg.norm(curve_edge_jet(block,{1:1.0,2:1.0},(x+1)/2)[1])
                for x,w in zip(nodes,weights))
    return [{"block":face["block"],"feature":face["feature"],"face_side":face["side"],
             "source_parameter_interval":source_t,"edge_orientation_reversed":bool(reversed_orientation),
             "integrated_edge_length_mm":length,"direct_source_jet_comparison":comparison}]


def face_area(block, side, order):
    degrees = block["degrees"]
    net = block["controls"]
    fixed_axis = "uvw".index(side[0])
    fixed_value = 1.0 if side[1] == "1" else 0.0
    free = [a for a in range(3) if a != fixed_axis]
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes = (nodes+1.0)/2.0
    weights = weights/2.0
    area = 0.0
    for i, u in enumerate(nodes):
        for j, v in enumerate(nodes):
            uvw = [0.0, 0.0, 0.0]
            uvw[fixed_axis] = fixed_value
            uvw[free[0]], uvw[free[1]] = u, v
            _, ders = homogeneous_eval_with_derivatives(net, degrees, uvw)
            area += weights[i]*weights[j]*np.linalg.norm(np.cross(ders[free[0]],ders[free[1]]))
    return area


def main():
    source = json.loads(SOURCE_REPORT.read_text())
    cases = {}
    for name, path in MODELS.items():
        if not path.exists():
            cases[name] = {"status": "model_unavailable", "path": str(path.relative_to(ROOT))}
            continue
        model = json.loads(path.read_text())
        by_id = {b["id"]: b for b in model["blocks"]}
        faces = model["design_faces"]
        selected = {}
        for order in (16, 32):
            areas = {}
            for face in faces:
                feature = face["feature"]
                if feature.startswith("pocket/"):
                    key = feature
                    areas[key] = areas.get(key, 0.0)+face_area(by_id[face["block"]],face["side"],order)
            selected[str(order)] = areas
        area16, area32 = selected["16"], selected["32"]
        expected = source["cases"][name]
        exp_floor = expected["new_surfaces"]["floor"]["exposed_area_mm2"]
        side = expected["new_surfaces"]["side_walls"]
        depth = float(source["pocket"]["support_z_mm"]-source["pocket"]["floor_z_mm"])
        exp_wall_parts = {"pocket/wall_ymin": depth*side["y_min_wall_length_mm"],
                          "pocket/wall_ymax": depth*side["y_max_wall_length_mm"],
                          "pocket/wall_xmin": side["x_min_wall_area_mm2"]}
        exp_walls = sum(exp_wall_parts.values())
        trace_results = {}
        source_arc=expected["old_arc_at_z4"]["old_floor_left_edge_arc_lengths_mm"]
        edge_a,edge_b=source_arc["edge_parameter_interval"]
        cut_a,cut_b=source_arc["truncated_portion_interval"]
        retained_left=[edge_a,min(edge_b,cut_a)]
        retained_right=[max(edge_a,cut_b),edge_b]
        if retained_left[1]>retained_left[0]:
            trace_results["retained_z4_left"]=find_trace_from_source(
                model["blocks"],faces,"rect_x/zmin",retained_left,4.0)
        if retained_right[1]>retained_right[0]:
            trace_results["retained_z4_right"]=find_trace_from_source(
                model["blocks"],faces,"rect_x/zmin",retained_right,4.0)
        trace_results["relocated_floor_z2"]=find_trace_from_source(
            model["blocks"],faces,"pocket/floor",[cut_a,cut_b],2.0)
        measured_z4=sum(e["integrated_edge_length_mm"]
                        for key in ("retained_z4_left","retained_z4_right")
                        for e in trace_results.get(key,[]))
        measured_floor=sum(e["integrated_edge_length_mm"]
                           for e in trace_results.get("relocated_floor_z2",[]))
        expected_z4=source_arc["retained_left_portion_length"]+source_arc["retained_right_portion_length"]
        topology=model.get("topology",{})
        cases[name] = {
            "status": "measured",
            "model_path": str(path.relative_to(ROOT)),
            "model_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "topology_counts": {"blocks":len(model["blocks"]),"exterior_sides":len(topology.get("exterior",[])),
                                 "interfaces":len(topology.get("interfaces",[]))},
            "pocket_face_areas_mm2_by_gauss_order": selected,
            "source_oracle": {"floor_area_mm2": exp_floor, "side_wall_area_mm2": exp_walls,
                              "side_wall_components_mm2": exp_wall_parts,
                              "truncated_arc_length_relocated_mm": expected["old_arc_at_z4"]["old_floor_left_edge_arc_lengths_mm"]["truncated_portion_relocated_on_floor_z2_length"]},
            "differences_order32_mm2": {
                "floor": area32.get("pocket/floor", 0.0)-exp_floor,
                "sum_walls": sum(v for k,v in area32.items() if k.startswith("pocket/wall_"))-exp_walls,
                "wall_components": {key:area32.get(key,0.0)-value for key,value in exp_wall_parts.items()}},
            "quadrature_delta_order32_minus16_mm2": {
                key: area32.get(key,0.0)-area16.get(key,0.0) for key in sorted(set(area16)|set(area32))},
            "model_face_counts": {key: sum(1 for f in faces if f["feature"]==key)
                                  for key in sorted({f["feature"] for f in faces if f["feature"].startswith("pocket/")})},
            "source_edge_traces": trace_results,
            "arc_length_comparison_mm": {
                "z4_retained_source_expected": expected_z4,
                "z4_retained_model_measured": measured_z4,
                "z4_retained_model_minus_source": measured_z4-expected_z4,
                "z2_relocated_source_expected": source_arc["truncated_portion_relocated_on_floor_z2_length"],
                "z2_relocated_model_measured": measured_floor,
                "z2_relocated_model_minus_source": measured_floor-source_arc["truncated_portion_relocated_on_floor_z2_length"]},
        }
    report = {"method": "direct tensor-product Bernstein integration of rational V-Rep boundary faces; 16/32 Gauss rules; source expectations from independent Cox/de Boor oracle", "cases": cases}
    REPORT.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
