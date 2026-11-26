#!/usr/bin/env python3
"""Exact torus, rank-stratum, and volume audit for emitted bracket fillet maps."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction as F
from math import gcd
from pathlib import Path

import numpy as np

from bracket_edge_edit_oracle import (
    ROOT, SOURCE, SAMPLES, add, bernstein_power, evaluate_net, model_map, mul,
    scale, verify_full_circle_sections,
)

_SPEC = importlib.util.spec_from_file_location(
    "s3_stock_rim_oracle", ROOT / "tests/references/s3_stock_rim.py")
_CERT = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CERT)


def read_net(span):
    du, dv, dw = span["degrees"]
    flat = span["controls_xyzw"]
    assert len(flat) == (du + 1) * (dv + 1) * (dw + 1)
    return [[[flat[(i * (dv + 1) + j) * (dw + 1) + k] for k in range(dw + 1)]
            for j in range(dv + 1)] for i in range(du + 1)]


def bivariate_control_poly(net, component):
    """Expand a 3x3 bivariate Bernstein control net to exact power form."""
    out = {}
    for i in range(3):
        for j in range(3):
            value = F(net[i][j][component])
            for pu, cu in bernstein_power(2, i).items():
                for pv, cv in bernstein_power(2, j).items():
                    out[(pu, pv)] = out.get((pu, pv), F(0)) + value * cu * cv
    return {k: v for k, v in out.items() if v}


def torus_residual(net, bore, meridional_center_q):
    """Clear W^4 in the implicit torus equation on radial parameter v=1."""
    axis = bore["axis"]
    center = [F(v, 4) for v in bore["center_quarters"]]
    center[axis] = F(meridional_center_q, 4)
    radius = F(bore["radius_quarters"], 4)
    fillet = F(1, 4)
    major = radius + fillet
    # The current native net axes are (angle, radial, meridional). Select its
    # inner bore/fillet surface at radial=1 and use angle/meridional parameters.
    surface = [[net[i][1][j] for j in range(3)] for i in range(3)]
    components = [bivariate_control_poly(surface, c) for c in range(4)]
    weight = components[3]
    offsets = [add(components[c], scale(weight, -center[c])) for c in range(3)]
    transverse = [c for c in range(3) if c != axis]
    radial2 = add(mul(offsets[transverse[0]], offsets[transverse[0]]),
                  mul(offsets[transverse[1]], offsets[transverse[1]]))
    axial2 = mul(offsets[axis], offsets[axis])
    sum_term = add(add(radial2, axial2), scale(mul(weight, weight), major**2 - fillet**2))
    return add(mul(sum_term, sum_term),
               scale(mul(mul(weight, weight), radial2), -4 * major**2))


def verify_meridional_endpoints(source_block, net, bore, span, radius):
    axis = bore["axis"]
    center = [F(v, 4) for v in bore["center_quarters"]]
    transverse = [(axis + 1) % 3, (axis + 2) % 3]
    source_net = model_map(source_block)
    lo_q = F(span["station_quarters"][0], 4)
    hi_q = F(span["station_quarters"][1], 4)
    lower = span["span_index"] == 0
    for u in SAMPLES:
        original = evaluate_net(source_net, (u, F(1), F(0)))
        direction = [original[a] - center[a] for a in transverse]
        assert sum(x*x for x in direction) == radius*radius
        for t, target_r, target_axis in (
            (F(0), radius + F(1,4), lo_q) if lower else (F(0), radius, lo_q),
            (F(1), radius, hi_q) if lower else (F(1), radius + F(1,4), hi_q),
        ):
            point = evaluate_net(net, (u, F(1), t))
            offset = [point[a] - center[a] for a in transverse]
            assert offset[0]*direction[1] - offset[1]*direction[0] == 0
            projected_radius = sum(offset[i]*direction[i] for i in range(2)) / radius
            assert projected_radius == target_r
            assert point[axis] == target_axis


def normalize_controls(controls):
    divisor = 0
    for row in controls:
        for value in row:
            divisor = gcd(divisor, abs(int(value)))
    assert divisor > 0
    return [[int(value) // divisor for value in row] for row in controls]


def verify_cli_equivalence(source, candidate, cli_path):
    cli_raw = cli_path.read_bytes()
    cli = json.loads(cli_raw)
    assert cli["schema"] == "fieldcraft_bore_rim_spans_v1"
    assert cli["state"] == "fillet" and cli["singular_cap_contact"] is True
    assert cli["parent_topology"] == source["topology"]
    assert cli["parent_design_faces"] == source["design_faces"]
    direct = {(e["logical_id"], e["span_index"]): e for e in candidate["elements"]}
    compared = 0
    for volume in cli["volumes"]:
        intervals = volume["parameter_intervals_32nds"]
        spans = volume["spans"]
        assert len(intervals) == len(spans)
        assert len(spans) == len([key for key in direct if key[0] == volume["id"]])
        for index, (interval, raw) in enumerate(zip(intervals, spans)):
            other = direct[(volume["id"], index)]
            assert interval == other["parameter_interval_32nds"]
            assert raw["id"] == f"{volume['id']}#{index}"
            assert raw["degrees"] == other["degrees"]
            assert normalize_controls(raw["controls"]) == normalize_controls(other["controls_xyzw"])
            compared += 1
    assert compared == len(candidate["elements"]) == 30
    return hashlib.sha256(cli_raw).hexdigest(), compared


def run(source_path: Path, output_path: Path, cli_path: Path):
    source_raw = source_path.read_bytes()
    source = json.loads(source_raw)
    candidate = json.loads(output_path.read_text())
    cli_sha, cli_maps_compared = verify_cli_equivalence(source, candidate, cli_path)
    assert candidate["schema"] == "fixed_partition_bore_edge_v1"
    assert candidate["state"] == "fillet"
    assert candidate["source_audit_status"] == "valid"
    assert candidate["singular_boundary_qualified"] is True
    assert candidate["determinant_lower_bound_mm3"] == 0
    assert candidate["edge_parameter_quarters"] == 1
    assert len(candidate["logical_blocks"]) == 10
    assert candidate["interfaces"] == 53
    source_blocks = {b["id"]: b for b in source["blocks"]}
    bore_configs = {b["chart_prefix"]: b for b in candidate["bores"]}
    assert set(bore_configs) == {"wing_east_", "wing_north_"}
    by_id = {}
    for element in candidate["elements"]:
        by_id.setdefault(element["logical_id"], []).append(element)
    assert set(by_id) == set(source_blocks)

    # Exact integer Bernstein determinant check on all 30 maps. The only
    # permitted zero strata are the bore/fillet tangent-contact rails.
    strata_by_id = {}
    for element in candidate["elements"]:
        strata = _CERT.certify(np.array(read_net(element), dtype=object))
        strata_by_id[(element["logical_id"], element["span_index"])] = strata

    bore_summaries = []
    mouth_count = 0
    curved_patch_count = 0
    for prefix, bore in bore_configs.items():
        axis = int(bore["axis"])
        center = [F(x, 4) for x in bore["center_quarters"]]
        radius = F(bore["radius_quarters"], 4)
        assert radius == 4
        chart_ids = sorted(k for k in source_blocks if k.startswith(prefix))
        assert len(chart_ids) == 4
        transverse = [(axis + 1) % 3, (axis + 2) % 3]
        center_2d = [center[a] for a in transverse]
        source_arcs = verify_full_circle_sections(
            source_blocks, chart_ids, center, transverse, radius)
        all_chart_records = []
        bounds = None
        for chart_id in chart_ids:
            source_block = source_blocks[chart_id]
            axial = [F(p[axis] * 4, p[3]) for p in source_block["controls"]]
            assert all(x.denominator == 1 for x in axial)
            chart_bounds = [int(min(axial)), int(max(axial))]
            bounds = chart_bounds if bounds is None else bounds
            assert chart_bounds == bounds
            spans = sorted(by_id[chart_id], key=lambda s: s["span_index"])
            assert len(spans) == 3
            assert [s["span_index"] for s in spans] == [0, 1, 2]
            assert [s["station_quarters"] for s in spans] == [[bounds[0], bounds[0]+1],
                                                               [bounds[0]+1, bounds[1]-1],
                                                               [bounds[1]-1, bounds[1]]]
            assert [s["parameter_interval_32nds"] for s in spans] == [[0, 1], [1, 31], [31, 32]]
            chart_patch_rows = []
            for span in spans:
                net = read_net(span)
                if span["span_index"] in (0, 2):
                    assert span["degrees"] == [2, 1, 2]
                    center_station = span["station_quarters"][0] + 1 if span["span_index"] == 0 else span["station_quarters"][1] - 1
                    residual = torus_residual(net, bore, center_station)
                    assert residual == {}, (chart_id, span["span_index"], residual)
                    verify_meridional_endpoints(source_block, net, bore, span, radius)
                    curved_patch_count += 1
                    expected_axis_rail = (-1 if span["span_index"] == 0 else 1)
                    zero = sorted(tuple(x) for x in strata_by_id[(chart_id, span["span_index"])] if x[1] == 1 and x[2] == expected_axis_rail)
                    assert zero, (chart_id, span["span_index"], strata_by_id[(chart_id, span["span_index"])] )
                    chart_patch_rows.append({
                        "span": span["span_index"],
                        "station_quarters": span["station_quarters"],
                        "parameter_interval_32nds": span["parameter_interval_32nds"],
                        "torus_identity_exact": True,
                        "zero_jacobian_strata": [list(x) for x in strata_by_id[(chart_id, span["span_index"])]],
                    })
                else:
                    assert span["degrees"] == [2, 1, 1]
                    assert not strata_by_id[(chart_id, span["span_index"])], (chart_id, strata_by_id[(chart_id, span["span_index"])])
                    chart_patch_rows.append({"span": 1, "cylindrical_radius_mm": 4.0, "zero_jacobian_strata": []})
                # At the untouched outer radial face v=0, every span must be
                # the exact source map restricted to its fixed source interval.
                for u in SAMPLES:
                    for t in SAMPLES:
                        source_t = F(span["parameter_interval_32nds"][0], 32) + F(span["parameter_interval_32nds"][1]-span["parameter_interval_32nds"][0], 32) * t
                        newp = evaluate_net(net, (u, F(0), t))
                        oldp = evaluate_net(model_map(source_block), (u, F(0), source_t))
                        assert newp == oldp, (chart_id, span["span_index"], u, t)
            all_chart_records.append({"chart": chart_id, "spans": chart_patch_rows})
        mouth_count += 2
        bore_summaries.append({
            "prefix": prefix,
            "axis": axis,
            "center_quarters": bore["center_quarters"],
            "radius_mm": float(radius),
            "fillet_radius_mm": 0.25,
            "axial_bounds_quarters": bounds,
            "fixed_parameter_breaks_32nds": [0, 1, 31, 32],
            "source_circle_arcs": source_arcs,
            "charts": all_chart_records,
        })

    # The two elbow maps are also split at exactly the same declared knots.
    for chart_id, source_block in source_blocks.items():
        if chart_id.startswith("wing_"):
            continue
        spans = sorted(by_id[chart_id], key=lambda s: s["span_index"])
        assert len(spans) == 3
        intervals = [[0, 1], [1, 31], [31, 32]]
        for span, interval in zip(spans, intervals):
            assert span["parameter_interval_32nds"] == interval
            net = read_net(span)
            assert not strata_by_id[(chart_id, span["span_index"])]
            for u in SAMPLES:
                for v in SAMPLES:
                    for w in SAMPLES:
                        local = evaluate_net(net, (u, v, w))
                        source_v = F(interval[0], 32) + F(interval[1]-interval[0], 32)*v
                        original = evaluate_net(model_map(source_block), (u, source_v, w))
                        assert local == original, (chart_id, span["span_index"], u, v, w)

    # Exact four-fillet analytic volume; this is a design reference, not a
    # volume read from a valid serialized B-rep (the maps have singular rails).
    R, r = F(4), F(1, 4)
    coeff_pi = 8 * R * r**2 + F(20, 3) * r**3
    coeff_pi2 = -2 * R * r**2 - 2 * r**3
    expected = float(coeff_pi) * math.pi + float(coeff_pi2) * math.pi**2
    assert curved_patch_count == 16 and mouth_count == 4
    assert candidate["jacobian_scope"] == "zero at certified bore/cap contact rails; strictly positive away from those rails"
    return {
        "schema": "independent_bracket_fillet_oracle_v1",
        "source_sha256": hashlib.sha256(source_raw).hexdigest(),
        "candidate_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "common_cli_sha256": cli_sha,
        "common_cli_maps_exact_up_to_positive_homogeneous_scaling": cli_maps_compared,
        "state": candidate["state"],
        "logical_blocks": len(candidate["logical_blocks"]),
        "elements": len(candidate["elements"]),
        "interfaces_declared": candidate["interfaces"],
        "physical_fillet_edges": mouth_count,
        "quarter_fillet_patches_checked": curved_patch_count,
        "exact_torus_identity_all_patches": True,
        "exact_jacobian_strata_per_map": [
            {"logical_id": key[0], "span_index": key[1], "zero_strata": [list(x) for x in value]}
            for key, value in sorted(strata_by_id.items())
        ],
        "no_zero_strata_on_regular_spans": True,
        "all_singular_strata_confined_to_16_fillet_patches": True,
        "volume_formula": "4 pi(2 R r^2(1-pi/4)+r^3(5/3-pi/2)) = (101/48)pi - (17/32)pi^2 mm^3",
        "expected_four_fillet_removed_volume_mm3": expected,
        "volume_status": "analytic reference only; zero-Jacobian contact rails mean no regular closed-cube chart qualification",
        "bores": bore_summaries,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--candidate", type=Path, default=ROOT / "build/bracket-edge-edit-fillet.json")
    parser.add_argument("--cli-output", type=Path, default=ROOT / "build/fieldcraft-x64-fillet.json")
    parser.add_argument("--output", type=Path, default=ROOT / "build/bracket-fillet-independent.json")
    args = parser.parse_args()
    report = run(args.source, args.candidate, args.cli_output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
