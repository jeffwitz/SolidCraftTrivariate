#!/usr/bin/env python3
"""Exact local normal-jet audit for the five compact FieldCraft brackets.

The audit tests whether each declared face join can be one split of a single
Bezier parent in the blocks' local coordinates. Homogeneous control nets are
degree-elevated exactly, face orientation and a positive projective gauge are
accounted for, and all normal derivatives through the common polynomial
degree are compared over the complete tangential Bernstein net.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from math import comb, factorial
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = ROOT / "docs/site/data"
IDS = [
    "fieldcraft_x40_y80",
    "fieldcraft_x50_y74",
    "fieldcraft_x51_y74",
    "fieldcraft_x64_y64",
    "fieldcraft_x80_y40",
]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def elevate_line(line, target):
    degree = len(line) - 1
    out = [line[0]]
    for i in range(1, degree + 1):
        a = F(i, degree + 1)
        out.append(tuple(a * x + (1 - a) * y for x, y in zip(line[i - 1], line[i])))
    out.append(line[-1])
    return out


def elevate_axis(net, degrees, axis, target):
    while degrees[axis] < target:
        old_shape = [d + 1 for d in degrees]
        new_degrees = degrees.copy()
        new_degrees[axis] += 1
        new_shape = [d + 1 for d in new_degrees]
        next_net = [None] * (new_shape[0] * new_shape[1] * new_shape[2])
        other = [a for a in range(3) if a != axis]
        for i in range(old_shape[other[0]]):
            for j in range(old_shape[other[1]]):
                old_indices, new_indices = [], []
                for k in range(old_shape[axis]):
                    p = [0, 0, 0]
                    p[axis] = k
                    p[other[0]], p[other[1]] = i, j
                    old_indices.append(
                        (p[0] * old_shape[1] + p[1]) * old_shape[2] + p[2]
                    )
                for k in range(new_shape[axis]):
                    p = [0, 0, 0]
                    p[axis] = k
                    p[other[0]], p[other[1]] = i, j
                    new_indices.append(
                        (p[0] * new_shape[1] + p[1]) * new_shape[2] + p[2]
                    )
                line = [net[index] for index in old_indices]
                elevated = elevate_line(line, target)
                for index, point in zip(new_indices, elevated):
                    next_net[index] = point
        degrees[axis] += 1
        net = next_net
    return net


def common_degree(block_a, block_b, side_a, side_b, orient):
    da, db = [int(x) for x in block_a["degrees"]], [int(x) for x in block_b["degrees"]]
    aa, _, free_a = parse_side(side_a)
    ab, _, free_b = parse_side(side_b)
    target_a, target_b = da.copy(), db.copy()
    target_normal = max(da[aa], db[ab])
    target_a[aa] = target_normal
    target_b[ab] = target_normal
    for k in range(2):
        other_k = 1 - k if orient.get("swap", False) else k
        target = max(da[free_a[k]], db[free_b[other_k]])
        target_a[free_a[k]] = target
        target_b[free_b[other_k]] = target

    # Elevate each axis independently; do not mutate input blocks.
    def converted(block, target):
        net = [tuple(F(x) for x in p) for p in block["controls"]]
        deg = [int(x) for x in block["degrees"]]
        for axis in range(3):
            net = elevate_axis(net, deg, axis, target[axis])
        return net

    return (
        target_a,
        target_b,
        converted(block_a, target_a),
        converted(block_b, target_b),
    )


def dims(degrees):
    return [d + 1 for d in degrees]


def index(degrees, p):
    d = dims(degrees)
    return (p[0] * d[1] + p[1]) * d[2] + p[2]


def parse_side(side):
    axis = "uvw".index(side[0])
    return axis, int(side[1]), [a for a in range(3) if a != axis]


def boundary_and_jet(net, degree, axis, endpoint, q):
    d = dims(degree)
    free = [a for a in range(3) if a != axis]
    result = []
    p = degree[axis]
    factor = factorial(p) // factorial(p - q)
    for a in range(d[free[0]]):
        row = []
        for b in range(d[free[1]]):
            value = [F(0)] * 4
            for j in range(q + 1):
                pos = [0, 0, 0]
                pos[axis] = j if endpoint == 0 else p - j
                pos[free[0]], pos[free[1]] = a, b
                coeff = F(
                    ((-1) ** (q - j) if endpoint == 0 else (-1) ** j)
                    * comb(q, j)
                    * factor
                )
                for c, x in enumerate(net[index(degree, pos)]):
                    value[c] += coeff * x
            row.append(tuple(value))
        result.append(row)
    return result


def orient_rows(rows, orient, shape):
    na, nb = shape
    out = []
    for a in range(na):
        line = []
        for b in range(nb):
            x = na - 1 - a if orient.get("flip_first", False) else a
            y = nb - 1 - b if orient.get("flip_second", False) else b
            sx, sy = (y, x) if orient.get("swap", False) else (x, y)
            line.append(rows[sx][sy])
        out.append(line)
    return out


def flatten(rows):
    return [p for row in rows for p in row]


def qjson(value):
    return [value.numerator, value.denominator]


def gauge(a, b):
    first, second = flatten(a), flatten(b)
    if len(first) != len(second) or not first or second[0][3] <= 0:
        return None
    scale = first[0][3] / second[0][3]
    if scale > 0 and all(
        x == scale * y for va, vb in zip(first, second) for x, y in zip(va, vb)
    ):
        return scale
    return None


def jet_restriction(net_a, net_b, degree_a, degree_b, side_a, side_b, orient):
    aa, ea, _ = parse_side(side_a)
    ab, eb, _ = parse_side(side_b)
    faces_a = boundary_and_jet(net_a, degree_a, aa, ea, 0)
    free_a = [a for a in range(3) if a != aa]
    target_shape = (degree_a[free_a[0]] + 1, degree_a[free_a[1]] + 1)
    faces_b = orient_rows(
        boundary_and_jet(net_b, degree_b, ab, eb, 0), orient, target_shape
    )
    g = gauge(faces_a, faces_b)
    if g is None:
        return {"face_gauge": "mismatch", "full_jet_match": False}
    order = degree_a[aa]
    ratios = []
    details = []
    for q in range(1, order + 1):
        ja = flatten(boundary_and_jet(net_a, degree_a, aa, ea, q))
        jb = flatten(
            orient_rows(
                boundary_and_jet(net_b, degree_b, ab, eb, q), orient, target_shape
            )
        )
        # Align increasing parent coordinates across either endpoint orientation.
        # The unknown positive interval-length ratio lambda must satisfy
        # A_q=lambda^q*g*B_q.
        scalar = None
        valid = True
        witness = None
        for point_index, (va, vb) in enumerate(zip(ja, jb)):
            for component, (x, y) in enumerate(zip(va, vb)):
                rhs = g * y * ((-1) ** q if ea == eb else 1)
                if rhs:
                    r = x / rhs
                    if r <= 0 or (scalar is not None and r != scalar):
                        valid = False
                        witness = {
                            "point_index": point_index,
                            "homogeneous_component": component,
                            "first_jet": qjson(x),
                            "gauge_scaled_second_jet": qjson(rhs),
                            "reason": "nonpositive_or_nonconstant_exact_jet_ratio",
                        }
                        break
                    scalar = r
                elif x:
                    valid = False
                    witness = {
                        "point_index": point_index,
                        "homogeneous_component": component,
                        "first_jet": qjson(x),
                        "gauge_scaled_second_jet": qjson(rhs),
                        "reason": "zero_second_jet_but_nonzero_first_jet",
                    }
                    break
            if not valid:
                break
        both_zero = scalar is None and valid and q > 1
        if scalar is None and not both_zero:
            valid = False
            if witness is None:
                witness = {"reason": "all_compared_second_jet_coefficients_are_zero"}
        ratios.append(scalar if valid else None)
        details.append(
            {
                "jet_order": q,
                "scaled_ratio_power": None
                if scalar is None or not valid
                else qjson(scalar),
                "match": valid,
                "both_jets_zero": both_zero,
                "exact_failure_witness": witness,
            }
        )
    # lambda is recovered from q=1. Check higher orders obey its powers.
    lam = ratios[0] if ratios else None
    cross_order_failure = None
    if lam is not None:
        for q, r in enumerate(ratios, 1):
            if r != lam**q and not details[q - 1]["both_jets_zero"]:
                cross_order_failure = {
                    "jet_order": q,
                    "observed_ratio": None if r is None else qjson(r),
                    "required_ratio_power": qjson(lam**q),
                    "reason": "local interval ratio is not consistent across jet orders",
                }
                break
    full = lam is not None and cross_order_failure is None
    return {
        "face_gauge_second_to_first": [g.numerator, g.denominator],
        "normal_degree_first": degree_a[aa],
        "normal_degree_second": degree_b[ab],
        "compared_jet_orders": order,
        "jet_scaled_ratios": details,
        "candidate_interval_length_ratio_first_to_second": None
        if lam is None
        else qjson(lam),
        "cross_order_failure": cross_order_failure,
        "full_jet_match": full,
    }


def audit_model(path):
    raw = path.read_bytes()
    model = json.loads(raw)
    source = model.get("source_model", {})
    blocks = {b["id"]: b for b in model["blocks"]}
    edges = model["topology"]["interfaces"]
    entries = []
    for i, e in enumerate(edges):
        aid, aside = e["first"]
        bid, bside = e["second"]
        degree_a, degree_b, net_a, net_b = common_degree(
            blocks[aid], blocks[bid], aside, bside, e.get("orientation", {})
        )
        jet = jet_restriction(
            net_a, net_b, degree_a, degree_b, aside, bside, e.get("orientation", {})
        )
        entries.append(
            {
                "interface_index": i,
                "first": e["first"],
                "second": e["second"],
                "orientation": e.get("orientation", {}),
                "common_degree_first_axes": degree_a,
                "common_degree_second_axes": degree_b,
                "parent_support_candidate": jet["full_jet_match"],
                "jet_evidence": jet,
                "local_parameter_provenance": {
                    "first_block": aid,
                    "second_block": bid,
                    "first_source_column": None,
                    "second_source_column": None,
                    "normalized_interval_lengths": "not recorded in model; ratio is only a candidate inferred from jets",
                },
            }
        )
    return {
        "model": path.name,
        "model_sha256": sha(raw),
        "source_schema": source.get("schema"),
        "partition_method": source.get("partition_plan", {}).get("method"),
        "block_count": len(model["blocks"]),
        "interface_count": len(edges),
        "interfaces": entries,
        "candidate_count": sum(x["parent_support_candidate"] for x in entries),
        "no_certified_parent_subdivisions": not any(
            x["parent_support_candidate"] for x in entries
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    reports = []
    for ident in IDS:
        path = MODEL_DIR / ident / f"{ident}.model.json"
        audit = json.loads(
            subprocess.run(
                [str(args.cli), "audit", str(path)],
                check=True,
                text=True,
                capture_output=True,
            ).stdout
        )
        if audit.get("status") != "valid":
            raise RuntimeError(f"invalid source {ident}: {audit}")
        eligibility = json.loads(
            subprocess.run(
                [str(args.cli), "spline-eligibility", str(path)],
                check=True,
                text=True,
                capture_output=True,
            ).stdout
        )
        report = audit_model(path)
        report["source_parent_report"] = {
            "conversion_status": eligibility.get("conversion_status"),
            "common_parent_obstruction": eligibility.get("common_parent_obstruction"),
            "known_source_restriction_status": eligibility.get(
                "known_source_restriction_status"
            ),
            "source_schema": eligibility.get("source_schema"),
            "source_validity": eligibility.get("source_validity"),
            "regional_translation_columns": [
                {
                    "block_id": b["block_id"],
                    "parametric_axis": b["parametric_axis"],
                    "translation_vector_mm": b["translation_vector_mm"],
                    "status": b["status"],
                }
                for b in eligibility.get("regional_source_columns", {}).get(
                    "blocks", []
                )
            ],
            "retained_join_count": len(eligibility.get("retained_joins", [])),
            "retained_join_reason_counts": {
                reason: sum(
                    j.get("reason") == reason
                    for j in eligibility.get("retained_joins", [])
                )
                for reason in sorted(
                    {j.get("reason") for j in eligibility.get("retained_joins", [])}
                )
            },
        }
        reports.append(report)
    summary = {
        "schema": "fieldcraft_same_parent_support_audit_v1",
        "method": "exact oriented homogeneous Bernstein normal jets, positive gauge, through common normal degree",
        "claim_scope": "necessary same-parent split test on every declared full-face join; inferred interval ratios are candidates, not provenance",
        "source_models": reports,
        "cli_sha256": sha(args.cli.read_bytes()),
        "totals": {
            "models": len(reports),
            "interfaces": sum(r["interface_count"] for r in reports),
            "full_jet_candidates": sum(r["candidate_count"] for r in reports),
            "interfaces_without_full_jet_candidate": sum(
                r["interface_count"] - r["candidate_count"] for r in reports
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["totals"], sort_keys=True))


if __name__ == "__main__":
    main()
