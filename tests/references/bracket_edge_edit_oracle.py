#!/usr/bin/env python3
"""Independent exact-radius and removed-volume oracle for bracket end chamfers.

Input is the JSON emitted by ``examples/bracket_edge_edit.rs``. The check uses
the frozen published source model, exact Fraction polynomial identities for
the two edited bore walls, exact restriction checks on untouched chart walls,
and an analytic cross-section integral. It does not trust the Rust audit's
radius or determinant fields.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction as F
from pathlib import Path
from math import comb

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"
TICKS_PER_MM = 4
SAMPLES = tuple(map(F, (0, F(1, 4), F(1, 2), F(3, 4), 1)))
_CERT_PATH = ROOT / "tests/references/s3_stock_rim.py"
_SPEC = importlib.util.spec_from_file_location("s3_stock_rim_oracle", _CERT_PATH)
_CERT_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CERT_MODULE)


def add(a, b):
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, F(0)) + value
        if out[key] == 0:
            del out[key]
    return out


def scale(a, c):
    return {key: value * c for key, value in a.items() if value * c}


def mul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            out[(i + k, j + l)] = out.get((i + k, j + l), F(0)) + x * y
    return {key: value for key, value in out.items() if value}


def bernstein_power(n, i):
    """Power coefficients of B_i^n(t), as a one-variable sparse polynomial."""
    return {i + k: F(comb(n, i) * comb(n - i, k) * (-1) ** k) for k in range(n - i + 1)}


def evaluate_net(net, q):
    """Evaluate homogeneous tensor Bezier net exactly at three parameters."""
    degrees = [len(net) - 1, len(net[0]) - 1, len(net[0][0]) - 1]
    h = [F(0)] * 4
    for i in range(degrees[0] + 1):
        bi = F(comb(degrees[0], i)) * q[0]**i * (1-q[0])**(degrees[0]-i)
        for j in range(degrees[1] + 1):
            bj = F(comb(degrees[1], j)) * q[1]**j * (1-q[1])**(degrees[1]-j)
            for k in range(degrees[2] + 1):
                bk = F(comb(degrees[2], k)) * q[2]**k * (1-q[2])**(degrees[2]-k)
                for c in range(4):
                    h[c] += F(net[i][j][k][c]) * bi * bj * bk
    if h[3] <= 0:
        raise AssertionError("nonpositive homogeneous denominator")
    return tuple(v / h[3] for v in h[:3])


def read_net(span):
    du, dv, dw = span["degrees"]
    flat = span["controls_xyzw"]
    assert len(flat) == (du + 1) * (dv + 1) * (dw + 1)
    return [[[flat[(i * (dv + 1) + j) * (dw + 1) + k] for k in range(dw + 1)]
            for j in range(dv + 1)] for i in range(du + 1)]


def affine_parameter_poly(a0, a1):
    return {(0, 0): a0, (0, 1): a1 - a0} if a0 != a1 else {(0, 0): a0}


def circle_residual(net, transverse_axes, center_mm, radius_endpoints, radial_index=1):
    # The map order is (angle, radial, axial). Restrict to the requested radial
    # endpoint and expand the angular/axial surface exactly in Q[u,t].
    du, dt = len(net) - 1, len(net[0][0]) - 1
    components = []
    for c in (*transverse_axes, 3):
        plane = [[net[i][radial_index][k][c] for k in range(dt + 1)] for i in range(du + 1)]
        poly = {}
        for i in range(du + 1):
            for j in range(dt + 1):
                value = F(plane[i][j])
                for pu, cu in bernstein_power(du, i).items():
                    for pv, cv in bernstein_power(dt, j).items():
                        poly[(pu, pv)] = poly.get((pu, pv), F(0)) + value * cu * cv
        components.append({k: v for k, v in poly.items() if v})
    X, Y, W = components
    X = add(X, scale(W, -center_mm[0]))
    Y = add(Y, scale(W, -center_mm[1]))
    radius = affine_parameter_poly(*radius_endpoints)
    residual = add(add(mul(X, X), mul(Y, Y)), scale(mul(mul(radius, radius), mul(W, W)), -1))
    return residual


def model_map(block):
    du, dv, dw = block["degrees"]
    flat = block["controls"]
    return [[[flat[(i * (dv + 1) + j) * (dw + 1) + k] for k in range(dw + 1)]
            for j in range(dv + 1)] for i in range(du + 1)]


def eval_model(block, q):
    return evaluate_net(model_map(block), q)


def eval_span(net, q):
    return evaluate_net(net, q)


def verify_full_circle_sections(source_blocks, chart_ids, center, transverse, radius):
    arcs = []
    for chart_id in chart_ids:
        net = model_map(source_blocks[chart_id])
        ends = [eval_model(source_blocks[chart_id], (u, F(1), F(0))) for u in (F(0), F(1))]
        uv = [tuple(point[a] - center[a] for a in transverse) for point in ends]
        for x, y in uv:
            assert x*x + y*y == radius*radius, (chart_id, x, y)
        cross = uv[0][0]*uv[1][1] - uv[0][1]*uv[1][0]
        dot = uv[0][0]*uv[1][0] + uv[0][1]*uv[1][1]
        assert abs(cross) == radius*radius and dot == 0, (chart_id, cross, dot)
        arcs.append((chart_id, uv[0], uv[1], 1 if cross > 0 else -1))
    assert len({arc[1] for arc in arcs}) == 4
    assert len({arc[2] for arc in arcs}) == 4
    assert len({arc[3] for arc in arcs}) == 1
    for _, _, end, _ in arcs:
        assert sum(1 for _, start, _, _ in arcs if start == end) == 1
    cursor = arcs[0][1]
    seen = set()
    for _ in arcs:
        match = next(i for i, arc in enumerate(arcs) if arc[1] == cursor)
        assert match not in seen
        seen.add(match)
        cursor = arcs[match][2]
    assert len(seen) == 4 and cursor == arcs[0][1]
    return [{"chart": a[0], "start_offset_mm": list(map(str, a[1])), "end_offset_mm": list(map(str, a[2])), "orientation": a[3]} for a in arcs]


def run(source_path: Path, result_path: Path):
    source_bytes = source_path.read_bytes()
    source = json.loads(source_bytes)
    result = json.loads(result_path.read_text())
    assert result["schema"] == "fixed_partition_bore_chamfer_v1"
    source_blocks = {b["id"]: b for b in source["blocks"]}
    assert result["source_model"] == str(source_path.relative_to(ROOT))
    assert result["source_audit_status"] == "valid"
    assert len(source_blocks) == 10 and len(result["logical_blocks"]) == 10
    assert len(set(result["logical_blocks"])) == 10
    assert len(result["elements"]) == 30
    assert result["interfaces"] == 53

    bore_configs = {b["chart_prefix"]: b for b in result["chamfers"]}
    assert set(bore_configs) == {"wing_east_", "wing_north_"}
    all_elements = result["elements"]
    by_id = {}
    for span in all_elements:
        by_id.setdefault(span["logical_id"], []).append(span)
    assert set(by_id) == set(source_blocks)
    for span in all_elements:
        zero_strata = _CERT_MODULE.certify(np.array(read_net(span), dtype=object))
        assert not zero_strata, (span["logical_id"], span["span_index"], zero_strata)

    spec_summary = []
    for prefix, cfg in bore_configs.items():
        assert cfg["both_axial_ends"]
        radius = F(cfg["radius_quarters"], TICKS_PER_MM)
        width = F(cfg["width_quarters"], TICKS_PER_MM)
        axis = int(cfg["axis"])
        center = [F(x, TICKS_PER_MM) for x in cfg["center_quarters"]]
        transverse = [(axis + 1) % 3, (axis + 2) % 3]
        center_2d = [center[a] for a in transverse]
        chart_ids = sorted(k for k in source_blocks if k.startswith(prefix))
        assert len(chart_ids) == 4, (prefix, chart_ids)
        circle_arcs = verify_full_circle_sections(source_blocks, chart_ids, center, transverse, radius)

        # Source axial support is derived from homogeneous source controls.
        axial_values = []
        for chart_id in chart_ids:
            block = source_blocks[chart_id]
            for p in block["controls"]:
                axial_values.append(F(p[axis], p[3]))
        lo, hi = min(axial_values), max(axial_values)
        loq, hiq = lo * TICKS_PER_MM, hi * TICKS_PER_MM
        assert loq.denominator == hiq.denominator == 1
        loq, hiq = int(loq), int(hiq)
        expected = [[loq, loq + int(cfg["width_quarters"])],
                    [loq + int(cfg["width_quarters"]), hiq - int(cfg["width_quarters"])],
                    [hiq - int(cfg["width_quarters"]), hiq]]
        assert expected[1][0] < expected[1][1]

        for chart_id in chart_ids:
            spans = sorted(by_id[chart_id], key=lambda x: x["span_index"])
            assert len(spans) == 3 and [s["span_index"] for s in spans] == [0, 1, 2]
            assert [s["station_quarters"] for s in spans] == expected
            src = source_blocks[chart_id]
            srcnet = model_map(src)
            # First prove that the frozen source boundary is the requested
            # full-radius circular bore, across the complete angle/axial patch.
            assert circle_residual(srcnet, transverse, center_2d, (radius, radius), radial_index=1) == {}
            # The outside radial boundary (edited v=0, source radial=0) must
            # remain the same surface after the axial interval restriction.
            for span in spans:
                net = read_net(span)
                assert span["degrees"] == [2, 1, 1]
                station = span["station_quarters"]
                z0, z1 = F(station[0], TICKS_PER_MM), F(station[1], TICKS_PER_MM)
                for u in SAMPLES:
                    for v in SAMPLES:
                        physical_axis = z0 + (z1 - z0) * v
                        source_axis = (physical_axis - lo) / (hi - lo)
                        pnew = eval_span(net, (u, F(0), v))
                        pold = eval_model(src, (u, F(0), source_axis))
                        assert pnew == pold, (chart_id, span["span_index"], u, v, pnew, pold)

            # Exact full-span circle identities prove the intended linear
            # taper at every angular and axial parameter, not only samples.
            for span in (spans[0], spans[2]):
                net = read_net(span)
                s0, s1 = span["station_quarters"]
                r0 = radius + width if s0 == loq else radius
                r1 = radius + width if s1 == hiq else radius
                residual = circle_residual(net, transverse, center_2d, (r0, r1), radial_index=1)
                assert residual == {}, (chart_id, span["span_index"], residual)
            # Untreated middle span retains a constant-radius R bore.
            mid = read_net(spans[1])
            assert circle_residual(mid, transverse, center_2d, (radius, radius), radial_index=1) == {}
        # Both independent mouths of each bore remove this exact section volume.
        added_volume_over_pi = 2 * (radius * width**2 + width**3 / 3)
        spec_summary.append({
            "prefix": prefix,
            "axis": axis,
            "radius_mm": float(radius),
            "width_mm": float(width),
            "axial_bounds_mm": [float(lo), float(hi)],
            "four_chart_sections": len(chart_ids),
            "full_circle_arc_coverage": circle_arcs,
            "mouths_verified": 2,
            "exact_removed_volume_pi_coefficient_mm3": str(added_volume_over_pi),
        })

    # Non-bore elbow charts must be exact parameter restrictions of their
    # originals, including the two extra matched blend stations.
    unchanged = []
    for block_id, block in source_blocks.items():
        if block_id.startswith("wing_"):
            continue
        spans = sorted(by_id[block_id], key=lambda x: x["span_index"])
        assert len(spans) == 3
        src = model_map(block)
        expected_intervals = [[0, 1], [1, 31], [31, 32]]
        if "parameter_interval_32nds" in spans[0]:
            assert [span["parameter_interval_32nds"] for span in spans] == expected_intervals
        for span in spans:
            net = read_net(span)
            assert span["degrees"] == block["degrees"]
            lo_t, hi_t = expected_intervals[span["span_index"]]
            for u in SAMPLES:
                for v in SAMPLES:
                    for w in SAMPLES:
                        local = eval_span(net, (u, v, w))
                        original = evaluate_net(src, (u, F(lo_t, 32) + F(hi_t-lo_t, 32) * v, w))
                        assert local == original, (block_id, span["span_index"], u, v, w)
        unchanged.append(block_id)

    total_pi = sum((F(x["exact_removed_volume_pi_coefficient_mm3"]) for x in spec_summary), F(0))
    report = {
        "schema": "independent_bore_chamfer_oracle_v1",
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "logical_block_count": len(result["logical_blocks"]),
        "element_count": len(all_elements),
        "interfaces_reported": result["interfaces"],
        "physical_bore_mouth_edge_denominator": 4,
        "quarter_arc_sections_checked": 16,
        "independently_checked_bores": spec_summary,
        "untouched_non_bore_maps_exact": unchanged,
        "total_exact_removed_volume_pi_coefficient_mm3": str(total_pi),
        "total_removed_volume_mm3": float(total_pi) * math.pi,
        "radius_method": "exact bivariate polynomial identity (X-cx W)^2+(Y-cy W)^2-r(v)^2 W^2 == 0 on each complete rational Bezier bore-wall span",
        "volume_method": "independent cross-section integral of pi*(r(z)^2-R^2), both axial mouths of each certified full-circle bore",
        "independent_jacobian_method": "integer Bernstein tensor determinant certificate reused from tests/references/s3_stock_rim.py; nonnegative coefficients and no zero stratum required on every edited span",
    }
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--result", type=Path, default=ROOT / "build/bracket-edge-edit.json")
    parser.add_argument("--output", type=Path, default=ROOT / "build/bracket-edge-edit-independent.json")
    args = parser.parse_args()
    report = run(args.source, args.result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
