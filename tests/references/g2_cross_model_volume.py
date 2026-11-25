#!/usr/bin/env python3
"""Integrate the rational trivariate maps and compare each layer to the two sketch cuts."""
import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def basis(n, x):
    return np.array([math.comb(n, i) * x**i * (1 - x)**(n-i) for i in range(n+1)]).T


def gradient(n, x):
    b = basis(n-1, x)
    return n * (np.pad(b, ((0,0),(1,0))) - np.pad(b, ((0,0),(0,1))))


def integrate(model, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    volumes = defaultdict(float)
    for block in model["blocks"]:
        p, q, r = block["degrees"]
        net = np.array(block["controls"], dtype=float).reshape(p+1,q+1,r+1,4)
        assert r == 1 and np.array_equal(net[:,:,0,[0,1,3]], net[:,:,1,[0,1,3]])
        z = net[:,:,:,2] / net[:,:,:,3]
        z0, z1 = z[0,0,:]
        assert np.all(z[:,:,0] == z0) and np.all(z[:,:,1] == z1) and z1 > z0
        plane = net[:,:,0,:]
        def evaluate(left, right):
            return np.einsum("ai,bj,ijc->abc", left, right, plane)
        h = evaluate(basis(p,nodes), basis(q,nodes))
        du = evaluate(gradient(p,nodes), basis(q,nodes))
        dv = evaluate(basis(p,nodes), gradient(q,nodes))
        w = h[:,:,3:4]
        a = (du[:,:,:2]*w - h[:,:,:2]*du[:,:,3:4]) / w**2
        b = (dv[:,:,:2]*w - h[:,:,:2]*dv[:,:,3:4]) / w**2
        det = (a[:,:,0]*b[:,:,1] - a[:,:,1]*b[:,:,0])*(z1-z0)
        assert np.all(det > 0)
        volumes[(z0,z1)] += float(weights @ det @ weights)
    return volumes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model",type=Path)
    parser.add_argument("source_report",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    model=json.loads(args.model.read_text())
    source=json.loads(args.source_report.read_text())["independent_volume_reference"]
    reference=source["remaining_stock_volume_mm3"]
    bounds=model["source_recipe"]["bounds_xyz_mm"]
    stock_area=(bounds[0][1]-bounds[0][0])*(bounds[1][1]-bounds[1][0])
    hole_area=source["periodic_profile_removed_volume_mm3"]/(bounds[2][1]-bounds[2][0])
    lo,hi=[a/b for a,b in model["sketch_plan"]["z_cut_levels_mm"]]
    extra_area=(source["rectangle_removed_volume_mm3"]-source["intersection_volume_q64_mm3"])/(hi-lo)
    rows=[]
    for order in (8,16,32):
        by_layer=integrate(model,order)
        total=sum(by_layer.values())
        layer_errors=[]
        for (a,b), volume in by_layer.items():
            expected=(stock_area-hole_area)*(b-a)-max(0,min(hi,b)-max(lo,a))*extra_area
            layer_errors.append(abs(volume-expected)/expected)
        rows.append({"order":order,"volume_mm3":total,"relative_source_error":abs(total-reference)/reference,"max_relative_layer_error":max(layer_errors),
                     "layers":[{"z_mm":[a,b],"volume_mm3":v} for (a,b),v in sorted(by_layer.items())]})
    assert max(rows[-1]["relative_source_error"],rows[-1]["max_relative_layer_error"]) < 1.e-8
    args.output.write_text(json.dumps({"schema":"cross_plate_map_quadrature_v1","model_sha256":hashlib.sha256(args.model.read_bytes()).hexdigest(),
        "source_reference_mm3":reference,"quadrature":rows},indent=2)+"\n")


if __name__ == "__main__":
    main()
