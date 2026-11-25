#!/usr/bin/env python3
"""Independent circle intersection and implicit-jet oracle for P7."""
from __future__ import annotations

import json
import math
import platform
import pathlib
import subprocess
import sys

MAX_CENTER = 1_000_000
MAX_RADIUS = 1_000_000
ROOT = pathlib.Path(__file__).resolve().parents[2]


def max_error(expected, actual):
    if expected is None or actual is None:
        return 0.0 if expected is None and actual is None else math.inf
    if isinstance(expected,(int,float)) and isinstance(actual,(int,float)):
        return abs(expected-actual)
    if isinstance(expected,list) and isinstance(actual,list) and len(expected)==len(actual):
        return max((max_error(a,b) for a,b in zip(expected,actual)),default=0.0)
    return 0.0 if expected==actual else math.inf


def call_rust_probe(payload):
    subprocess.run(["rtk","cargo","build","--release","--locked","--offline","--example","p6b_oracle_probe"],
                   cwd=ROOT,check=True,capture_output=True,text=True)
    binary=ROOT/"target/release/examples/p6b_oracle_probe"
    completed=subprocess.run([str(binary),"p7"],input=json.dumps(payload,allow_nan=False),
                             cwd=ROOT,check=True,capture_output=True,text=True)
    return json.loads(completed.stdout)


def classify(c1: tuple[int, int, int], c2: tuple[int, int, int]) -> dict:
    if any(abs(v) > MAX_CENTER for v in (*c1[:2], *c2[:2])) or any(
        r <= 0 or r > MAX_RADIUS for r in (c1[2], c2[2])
    ):
        return {"accepted": False, "reason": "outside supported circle domain"}
    dx, dy = c2[0]-c1[0], c2[1]-c1[1]
    d2 = dx*dx + dy*dy
    r1, r2 = c1[2], c2[2]
    if d2 == 0:
        if r1 == r2:
            return {"accepted": True, "relation": "Coincident", "cardinality": "Infinite", "count": 0,"rank":None,"points":[],"tangent_point_mm":None}
        return {"accepted": True, "relation": "ConcentricDisjoint", "cardinality": "Zero", "count": 0,"rank":None,"points":[],"tangent_point_mm":None}
    outer=(r1+r2)*(r1+r2)-d2
    inner=d2-(r1-r2)*(r1-r2)
    if outer < 0:
        return {"accepted": True, "relation": "Disjoint", "cardinality": "Zero", "count": 0,"rank":None,"points":[],"tangent_point_mm":None}
    if inner < 0:
        return {"accepted": True, "relation": "StrictContainment", "cardinality": "Zero", "count": 0,"rank":None,"points":[],"tangent_point_mm":None}
    d = math.sqrt(d2)
    a = (r1*r1-r2*r2+d2)/(2.0*d)
    base = (c1[0]/1000.0 + (a/d)*(dx/1000.0),
            c1[1]/1000.0 + (a/d)*(dy/1000.0))
    if outer == 0:
        return {"accepted": True, "relation": "TangentExternal", "cardinality": "One", "count": 1,
                "rank":1,"points":[],"tangent_point_mm": list(base), "derivatives": None}
    if inner == 0:
        return {"accepted": True, "relation": "TangentInternal", "cardinality": "One", "count": 1,
                "rank":1,"points":[],"tangent_point_mm": list(base), "derivatives": None}
    h = math.sqrt(max(0.0, r1*r1-a*a))
    ex, ey = dx/d, dy/d
    points = []
    for sign in (1.0, -1.0):
        p = (base[0]-sign*(h*ey)/1000.0, base[1]+sign*(h*ex)/1000.0)
        # Implicit function theorem: solve the 2x2 constraint Jacobian for
        # each physical parameter, independently of the closed-form point jet.
        x1, y1 = p[0]-c1[0]/1000.0, p[1]-c1[1]/1000.0
        x2, y2 = p[0]-c2[0]/1000.0, p[1]-c2[1]/1000.0
        matrix = ((2*x1, 2*y1), (2*x2, 2*y2))
        det = matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
        jets = []
        jet_residual = 0.0
        param_partials = ((-2*x1, 0.0), (-2*y1, 0.0), (-2*r1/1000.0, 0.0),
                          (0.0, -2*x2), (0.0, -2*y2), (0.0, -2*r2/1000.0))
        for df1, df2 in param_partials:
            b0, b1 = -df1, -df2
            dp=((b0*matrix[1][1]-matrix[0][1]*b1)/det,
                (matrix[0][0]*b1-b0*matrix[1][0])/det)
            jets.append(dp)
            jet_residual=max(jet_residual,abs(matrix[0][0]*dp[0]+matrix[0][1]*dp[1]+df1),
                             abs(matrix[1][0]*dp[0]+matrix[1][1]*dp[1]+df2))
        residual = max(abs(x1*x1+y1*y1-(r1/1000.0)**2),
                       abs(x2*x2+y2*y2-(r2/1000.0)**2))
        relative_residual=residual/max((r1/1000.0)**2,(r2/1000.0)**2)
        gram_trace = sum(v*v for v in (2*x1,2*y1,2*x2,2*y2))
        gram_det = det*det
        disc = math.sqrt(max(0.0, gram_trace*gram_trace-4*gram_det))
        eig_hi = (gram_trace+disc)/2.0
        sigma_max=math.sqrt(eig_hi)
        sigma_min=abs(det)/sigma_max
        condition=sigma_max/sigma_min
        points.append({"point_mm": list(p), "jet_mm_per_mm": [list(v) for v in jets],
                       "constraint_residual_mm2": residual,
                       "constraint_residual_relative_to_radius2":relative_residual,
                       "implicit_condition": condition,
                       "implicit_jet_residual_mm":jet_residual})
    return {"accepted": True, "relation": "Transverse", "cardinality": "Two",
            "count": 2, "rank": 2, "points": points,"tangent_point_mm":None}


def main() -> dict:
    cases = [
        ("transverse_axis", (0,0,5000), (6000,0,5000)),
        ("transverse_oblique", (0,0,5000), (6000,1000,5000)),
        ("external_tangent", (0,0,5000), (10000,0,5000)),
        ("external_1um_inside", (0,0,5000), (9999,0,5000)),
        ("external_1um_outside", (0,0,5000), (10001,0,5000)),
        ("internal_tangent", (0,0,5000), (4000,0,1000)),
        ("internal_1um_inside", (0,0,5000), (3999,0,1000)),
        ("internal_1um_outside", (0,0,5000), (4001,0,1000)),
        ("disjoint", (0,0,5000), (11000,0,5000)),
        ("strict_containment", (0,0,5000), (2000,0,1000)),
        ("concentric_disjoint", (0,0,5000), (0,0,4000)),
        ("coincident", (0,0,5000), (0,0,5000)),
        ("near_coincident_1um", (0,0,1000000), (1,0,1000000)),
        ("at_center_bound", (1000000,0,1000000), (994000,0,1000000)),
        ("at_radius_bound", (0,0,1000000), (1000000,0,1000000)),
        ("beyond_center_bound", (1000001,0,1000), (0,0,1000)),
        ("zero_radius", (0,0,0), (1,0,1)),
        ("negative_radius", (0,0,-1), (1,0,1)),
        ("beyond_radius_bound", (0,0,1000001), (1,0,1)),
        ("i64_min_center", (-9223372036854775808,0,1), (0,0,1)),
        ("i64_max_radius", (0,0,9223372036854775807), (1,0,1)),
    ]
    out=[]
    max_res=0.0
    max_relative_res=0.0
    max_jet_res=0.0
    native_errors={"point_mm":0.0,"jet_mm_per_mm":0.0,"condition":0.0,"tangent_point_mm":0.0}
    native=call_rust_probe([{"case":name,"first":a,"second":b} for name,a,b in cases])
    if len(native)!=len(cases):
        raise AssertionError("Rust P7 probe returned a different case count")
    native_by_name={x["case"]:x for x in native}
    for name,a,b in cases:
        result=classify(a,b)
        row={"case":name,"first":a,"second":b,**result}
        rust=native_by_name[name]
        row["rust"]=rust
        if rust["accepted"] != result["accepted"]:
            raise AssertionError(f"P7 admission mismatch: {name}")
        if result["accepted"]:
            for key in ("relation","cardinality","count","rank"):
                if rust.get(key)!=result.get(key):
                    raise AssertionError(f"P7 {key} mismatch: {name}: {rust.get(key)} vs {result.get(key)}")
            tangent_error=max_error(result.get("tangent_point_mm"),rust.get("tangent_point_mm"))
            native_errors["tangent_point_mm"]=max(native_errors["tangent_point_mm"],tangent_error)
            if tangent_error>2e-12:
                raise AssertionError(f"P7 tangent point mismatch: {name}")
            if len(result.get("points",[]))!=len(rust.get("points",[])):
                raise AssertionError(f"P7 point count mismatch: {name}")
            for refp,rp in zip(result.get("points",[]),rust.get("points",[])):
                point_error=max_error(refp["point_mm"],rp["point_mm"])
                native_errors["point_mm"]=max(native_errors["point_mm"],point_error)
                if point_error>2e-11:
                    raise AssertionError(f"P7 point coordinate mismatch: {name}")
                jet_error=max_error(refp["jet_mm_per_mm"],rp["jet_mm_per_mm"])
                native_errors["jet_mm_per_mm"]=max(native_errors["jet_mm_per_mm"],jet_error)
                if jet_error>3e-6:
                    raise AssertionError(f"P7 implicit jet mismatch: {name}")
                if abs(refp["constraint_residual_mm2"]-rp["constraint_residual_mm2"])>2e-12:
                    raise AssertionError(f"P7 constraint residual mismatch: {name}")
                condition_error=abs(refp["implicit_condition"]-rp["condition"])
                native_errors["condition"]=max(native_errors["condition"],condition_error)
                if condition_error>max(1e-8,1e-10*refp["implicit_condition"]):
                    raise AssertionError(f"P7 condition mismatch: {name}")
        for point in result.get("points",[]):
            max_res=max(max_res,point["constraint_residual_mm2"])
            max_relative_res=max(max_relative_res,point["constraint_residual_relative_to_radius2"])
            max_jet_res=max(max_jet_res,point["implicit_jet_residual_mm"])
            point["fd_jet_error_max_mm_per_mm"]=None
        out.append(row)
    near_coincident=classify((0,0,1000000),(1,0,1000000))
    # Centered integer-micrometre perturbations validate each implicit column.
    base_circles=((0,0,5000),(6000,1000,5000))
    base=classify(*base_circles)
    max_fd=0.0
    for p in range(6):
        ci,field=divmod(p,3)
        plus=[list(v) for v in base_circles]; minus=[list(v) for v in base_circles]
        plus[ci][field]+=1; minus[ci][field]-=1
        rp=classify(tuple(plus[0]),tuple(plus[1]))["points"][0]["point_mm"]
        rm=classify(tuple(minus[0]),tuple(minus[1]))["points"][0]["point_mm"]
        jet=base["points"][0]["jet_mm_per_mm"][p]
        for axis in range(2):
            fd=(rp[axis]-rm[axis])/0.002
            max_fd=max(max_fd,abs(fd-jet[axis]))
    # Rigid symmetries independently map the ordered axis-side intersections.
    square=classify((0,0,5000),(6000,0,5000))
    rotated=classify((0,0,5000),(0,6000,5000))
    rotation_error=max(abs(rotated["points"][i]["point_mm"][0]+square["points"][i]["point_mm"][1])
                       for i in range(2))
    rotation_error=max(rotation_error,max(abs(rotated["points"][i]["point_mm"][1]-square["points"][i]["point_mm"][0])
                                          for i in range(2)))
    swapped=classify((6000,0,5000),(0,0,5000))
    swap_error=max(abs(swapped["points"][i]["point_mm"][1]-square["points"][1-i]["point_mm"][1]) for i in range(2))
    translated=classify((7000,-3000,5000),(13000,-3000,5000))
    translation_error=max(abs(translated["points"][i]["point_mm"][axis]
                              -square["points"][i]["point_mm"][axis]
                              -(7.0 if axis==0 else -3.0))
                          for i in range(2) for axis in range(2))
    oblique=classify((0,0,5000),(6000,1000,5000))
    reflected=classify((0,0,5000),(6000,-1000,5000))
    reflection_error=max(abs(reflected["points"][i]["point_mm"][0]-oblique["points"][1-i]["point_mm"][0])
                         for i in range(2))
    reflection_error=max(reflection_error,max(abs(reflected["points"][i]["point_mm"][1]+oblique["points"][1-i]["point_mm"][1])
                                               for i in range(2)))
    tay=[]
    for step in (100,200):
        perturbed=classify((0,0,5000+step),(6000,1000,5000))["points"][0]["point_mm"]
        h=step/1000.0
        jet=base["points"][0]["jet_mm_per_mm"][2]
        tay.append(math.dist(perturbed,[base["points"][0]["point_mm"][axis]+jet[axis]*h for axis in range(2)]))
    return {"schema":"p7_sketch_intersections_oracle_v1",
            "provenance":{"source_baseline":"7b61094+ff2c36f","checkout_revision":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
                          "python":platform.python_version(),"protocol":"docs/lot-p6b-oracle.md",
                          "commands":["rtk python3 tests/references/sketch_intersections_oracle.py",
                                      "rtk cargo build --release --locked --offline --example p6b_oracle_probe",
                                      "rtk cargo test --locked --offline sketch_intersections"]},
            "native_probe_case_count":len(native),"native_reference_comparison":"passed",
            "max_native_reference_absolute_errors":native_errors,
            "acceptance":{"ordinary_constraint_residual_mm2_max":2e-12,
                          "implicit_jet_residual_mm_max":2e-11,
                          "centered_fd_jet_error_mm_per_mm_max":3e-6,
                          "taylor_ratio_target":4.0,"taylor_ratio_tolerance":0.15},
            "case_count":len(out),"accepted_transverse_count":sum(x.get("relation")=="Transverse" for x in out),
            "max_constraint_residual_mm2":max_res,"max_centered_fd_jet_error_mm_per_mm":max_fd,
            "max_constraint_residual_relative_to_radius2":max_relative_res,
            "max_implicit_jet_residual_mm":max_jet_res,
            "rotation_covariance_error_mm":rotation_error,"translation_covariance_error_mm":translation_error,
            "reflection_covariance_error_mm":reflection_error,"swapped_order_error_mm":swap_error,
            "taylor_remainder_mm":tay,"taylor_ratio":tay[1]/tay[0],
            "near_coincident_condition":near_coincident["points"][0]["implicit_condition"],
            "cases":out}


if __name__ == "__main__":
    print(json.dumps(main(),indent=2,allow_nan=False))
