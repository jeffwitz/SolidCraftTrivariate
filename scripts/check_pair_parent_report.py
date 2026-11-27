#!/usr/bin/env python3
"""Independent exact replay of native pairwise Bézier-parent reports.

This Fraction oracle parses the source nets and report only. It performs its
own axis orientation, Bernstein degree elevation, gauge comparison, normal-jet
test, and de Casteljau reconstruction of every recovered parent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb, factorial
from pathlib import Path


def index(degrees, p):
    dims = [d + 1 for d in degrees]
    return (p[0] * dims[1] + p[1]) * dims[2] + p[2]


def net_from_block(block):
    degrees = [int(x) for x in block["degrees"]]
    net = [tuple(Fraction(int(x)) for x in point) for point in block["controls"]]
    if len(net) != (degrees[0]+1)*(degrees[1]+1)*(degrees[2]+1):
        raise ValueError("source net count mismatch")
    if any(point[3] <= 0 for point in net):
        raise ValueError("source net has a nonpositive weight")
    return degrees, net


def reverse_axis(net, degrees, axis):
    out = [None] * len(net)
    for p in product(*(range(d+1) for d in degrees)):
        source = list(p)
        source[axis] = degrees[axis] - p[axis]
        out[index(degrees,p)] = net[index(degrees,source)]
    return out


def orient_second(net, source_degrees, first_axis, second_axis, second_upper, orientation):
    first_free = [a for a in range(3) if a != first_axis]
    second_free = [a for a in range(3) if a != second_axis]
    source_axis = [None] * 3
    source_axis[first_axis] = second_axis
    flips = [False] * 3
    flips[first_axis] = second_upper
    for out_tangent in range(2):
        in_tangent = 1-out_tangent if orientation["swap"] else out_tangent
        source_axis[first_free[out_tangent]] = second_free[in_tangent]
        if orientation["swap"]:
            flips[first_free[out_tangent]] = (
                orientation["flip_second"] if out_tangent == 0 else orientation["flip_first"]
            )
        else:
            flips[first_free[out_tangent]] = (
                orientation["flip_first"] if out_tangent == 0 else orientation["flip_second"]
            )
    degrees = [source_degrees[a] for a in source_axis]
    out = [None] * len(net)
    for p in product(*(range(d+1) for d in degrees)):
        source = [None] * 3
        for axis, src_axis in enumerate(source_axis):
            source[src_axis] = source_degrees[src_axis]-p[axis] if flips[axis] else p[axis]
        out[index(degrees,p)] = net[index(source_degrees,source)]
    return degrees, out


def elevate_line(line):
    p = len(line)-1
    result = [line[0]]
    for i in range(1,p+1):
        alpha = Fraction(i,p+1)
        result.append(tuple(alpha*a+(1-alpha)*b for a,b in zip(line[i-1],line[i])))
    result.append(line[-1])
    return result


def elevate_axis(net, degrees, axis, target):
    degrees = list(degrees)
    while degrees[axis] < target:
        old_degrees = list(degrees)
        degrees[axis] += 1
        out = [None] * ((degrees[0]+1)*(degrees[1]+1)*(degrees[2]+1))
        other = [a for a in range(3) if a != axis]
        for a in range(old_degrees[other[0]]+1):
            for b in range(old_degrees[other[1]]+1):
                line_indices = []
                for k in range(old_degrees[axis]+1):
                    p = [0,0,0]
                    p[axis], p[other[0]], p[other[1]] = k,a,b
                    line_indices.append(index(old_degrees,p))
                elevated = elevate_line([net[i] for i in line_indices])
                for k,point in enumerate(elevated):
                    p = [0,0,0]
                    p[axis], p[other[0]], p[other[1]] = k,a,b
                    out[index(degrees,p)] = point
        net = out
    return degrees, net


def face(net, degrees, axis, upper):
    free = [a for a in range(3) if a != axis]
    fixed = degrees[axis] if upper else 0
    result = []
    for i in range(degrees[free[0]]+1):
        for j in range(degrees[free[1]]+1):
            p = [0,0,0]
            p[axis],p[free[0]],p[free[1]] = fixed,i,j
            result.append(net[index(degrees,p)])
    return result


def gauge(first, second):
    if len(first) != len(second) or not first:
        raise ValueError("trace dimensions differ")
    scale = first[0][3]/second[0][3]
    if scale <= 0 or any(a != tuple(scale*x for x in b) for a,b in zip(first,second)):
        raise ValueError("traces do not agree up to one positive homogeneous gauge")
    return scale


def normal_jet(net, degrees, axis, upper, order):
    free = [a for a in range(3) if a != axis]
    p = degrees[axis]
    factor = Fraction(factorial(p),factorial(p-order))
    result = []
    for i in range(degrees[free[0]]+1):
        for j in range(degrees[free[1]]+1):
            value = [Fraction(0)]*4
            for h in range(order+1):
                normal_index = p-h if upper else h
                sign = (-1)**h if upper else (-1)**(order-h)
                coefficient = factor*sign*comb(order,h)
                point = [0,0,0]
                point[axis],point[free[0]],point[free[1]] = normal_index,i,j
                for c in range(4):
                    value[c] += coefficient*net[index(degrees,point)][c]
            result.append(tuple(value))
    return result


def full_jet_candidate(first, second, degrees, axis):
    ratio = None
    for left,right in zip(normal_jet(first,degrees,axis,True,1),normal_jet(second,degrees,axis,False,1)):
        for a,b in zip(left,right):
            if b:
                candidate = a/b
                if candidate <= 0 or (ratio is not None and ratio != candidate):
                    return False, None
                ratio = candidate
            elif a:
                return False, None
    if ratio is None:
        return False, None
    for order in range(2,degrees[axis]+1):
        left = normal_jet(first,degrees,axis,True,order)
        right = normal_jet(second,degrees,axis,False,order)
        scale = ratio**order
        if any(a != tuple(scale*x for x in b) for a,b in zip(left,right)):
            return False, ratio
    return True, ratio


def split_line(line, t):
    levels = [list(line)]
    while len(levels[-1]) > 1:
        previous = levels[-1]
        levels.append([
            tuple((1-t)*a+t*b for a,b in zip(previous[i],previous[i+1]))
            for i in range(len(previous)-1)
        ])
    left = [level[0] for level in levels]
    right = [levels[-1-i][i] for i in range(len(levels))]
    return left,right


def split_net(net,degrees,axis,t):
    left=[None]*len(net)
    right=[None]*len(net)
    free=[a for a in range(3) if a!=axis]
    for i in range(degrees[free[0]]+1):
        for j in range(degrees[free[1]]+1):
            p=[0,0,0]
            p[free[0]],p[free[1]]=i,j
            line=[]
            for k in range(degrees[axis]+1):
                p[axis]=k
                line.append(net[index(degrees,p)])
            l,r=split_line(line,t)
            for k in range(degrees[axis]+1):
                p[axis]=k
                left[index(degrees,p)]=l[k]
                right[index(degrees,p)]=r[k]
    return left,right


def q(value):
    return Fraction(int(value[0]),int(value[1]))


def check(model_path, report_path):
    model_bytes=Path(model_path).read_bytes()
    report_bytes=Path(report_path).read_bytes()
    model=json.loads(model_bytes)
    report=json.loads(report_bytes)
    blocks={b["id"]:b for b in model["blocks"]}
    interfaces=model["topology"]["interfaces"]
    rows=report["interfaces"]
    if report.get("schema")!="pair_bezier_parent_report_v1" or len(rows)!=len(interfaces):
        raise ValueError("native pair report does not match the source interface list")
    recovered=refused=indeterminate=0
    jet_candidates=jet_refused=0
    jet_candidates_by_axis={axis:0 for axis in "uvw"}
    max_parent_coefficients=0
    axis_counts={axis:{"interfaces":0,"recovered":0,"refused":0,"indeterminate":0} for axis in "uvw"}
    for ordinal,(source,row) in enumerate(zip(interfaces,rows)):
        if row["interface_index"]!=ordinal or row["first"]!=source["first"] or row["second"]!=source["second"]:
            raise ValueError(f"interface order/side metadata differs at {ordinal}")
        orientation=source.get("orientation",{"swap":False,"flip_first":False,"flip_second":False})
        if row["orientation"]!=orientation:
            raise ValueError(f"interface orientation differs at {ordinal}")
        first_id,first_face=source["first"]
        second_id,second_face=source["second"]
        axis_counts[first_face[0]]["interfaces"]+=1
        first_axis,first_upper="uvw".index(first_face[0]),first_face[1]=="1"
        second_axis,second_upper="uvw".index(second_face[0]),second_face[1]=="1"
        da,a=net_from_block(blocks[first_id])
        db,b=net_from_block(blocks[second_id])
        if not first_upper:
            a=reverse_axis(a,da,first_axis)
        db,b=orient_second(b,db,first_axis,second_axis,second_upper,orientation)
        target=[max(x,y) for x,y in zip(da,db)]
        for axis in range(3):
            da,a=elevate_axis(a,da,axis,target[axis])
            db,b=elevate_axis(b,db,axis,target[axis])
        try:
            scale=gauge(face(a,target,first_axis,True),face(b,target,first_axis,False))
            b=[tuple(scale*x for x in p) for p in b]
            candidate,ratio=full_jet_candidate(a,b,target,first_axis)
        except (ValueError,ZeroDivisionError):
            candidate,ratio=False,None
        jet_candidates+=int(candidate)
        jet_candidates_by_axis[first_face[0]]+=int(candidate)
        status=row["status"]
        if status=="recovered":
            recovered+=1
            axis_counts[first_face[0]]["recovered"]+=1
            rec=row["recovery"]
            if rec.get("source_block_ids") != [first_id,second_id] or rec.get("identity_kind") != "mathematical_parent_recovery_not_feature_provenance":
                raise ValueError(f"recovery lost pair provenance or overclaims feature identity at {ordinal}")
            if not candidate:
                raise ValueError(f"Rust recovered a pair that fails independent normal jets at {ordinal}")
            degrees=[int(x) for x in rec["degrees"]]
            parent=[tuple(q(x) for x in point) for point in rec["controls_homogeneous"]]
            t=q(rec["split_parameter"])
            if degrees!=target or t<=0 or t>=1:
                raise ValueError(f"recovered parent degree or split parameter is invalid at {ordinal}")
            if q(rec["first_to_second_interval_ratio"])!=t/(1-t):
                raise ValueError(f"recovered interval ratio is incorrect at {ordinal}")
            if q(rec["second_to_first_gauge"])!=scale:
                raise ValueError(f"reported projective gauge is incorrect at {ordinal}")
            for p in parent:
                if p[3]<=0:
                    raise ValueError(f"recovered parent has nonpositive weight at {ordinal}")
            left,right=split_net(parent,degrees,first_axis,t)
            if left!=a or right!=b:
                raise ValueError(f"recovered parent does not reproduce full input nets at {ordinal}")
            # The first-jet interval ratio is independently recovered from
            # exact endpoint derivative control nets, not from the report.
            if ratio!=t/(1-t):
                raise ValueError(f"parent interval ratio disagrees with independent jets at {ordinal}")
            max_parent_coefficients=max(max_parent_coefficients,len(parent))
        elif status=="refused":
            refused+=1
            axis_counts[first_face[0]]["refused"]+=1
            jet_refused+=int(candidate)
            if row.get("recovery") is not None:
                raise ValueError(f"refused pair unexpectedly contains a parent at {ordinal}")
        elif status=="indeterminate":
            indeterminate+=1
            axis_counts[first_face[0]]["indeterminate"]+=1
        else:
            raise ValueError(f"unknown pair status at {ordinal}")
    if (recovered,refused,indeterminate)!=(report["recovered_count"],report["refused_count"],report["indeterminate_count"]):
        raise ValueError("native report summary counts disagree with pair records")
    native_axis={entry["axis"]:entry for entry in report["counts_by_first_interface_axis"]}
    for axis,counts in axis_counts.items():
        if (counts["interfaces"],counts["recovered"],counts["refused"],counts["indeterminate"]) != (
            native_axis[axis]["interface_count"],native_axis[axis]["recovered_count"],
            native_axis[axis]["refused_count"],native_axis[axis]["indeterminate_count"]):
            raise ValueError(f"native per-axis counts disagree for {axis}")
    result={
        "schema":"pair_parent_independent_oracle_v1",
        "model":Path(model_path).name,
        "report":Path(report_path).name,
        "model_sha256":sha256(model_bytes).hexdigest(),
        "native_report_sha256":sha256(report_bytes).hexdigest(),
        "interfaces":len(interfaces),
        "counts_by_first_interface_axis":axis_counts,
        "recovered_and_exactly_replayed":recovered,
        "refused":refused,
        "indeterminate":indeterminate,
        "independent_full_normal_jet_candidates":jet_candidates,
        "independent_full_normal_jet_candidates_by_axis":jet_candidates_by_axis,
        "full_jet_candidates_refused_by_exact_parent_reconstruction":jet_refused,
        "max_recovered_parent_control_count":max_parent_coefficients,
        "status":"pass",
    }
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model",type=Path)
    parser.add_argument("report",type=Path)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=check(args.model,args.report)
    output=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output)
    print(output,end="")


if __name__=="__main__":
    main()
