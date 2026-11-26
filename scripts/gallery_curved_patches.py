#!/usr/bin/env python3
"""Build exact clamped-NURBS volume patches for the gallery's curved models.

Join compatible pairs across one shared full face, using exact homogeneous
degree elevation and orientation-preserving parameter maps. Keep remaining
cells as single-span clamped volumes. The emitted span maps allow independent
verification against the complete source control nets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from urllib.parse import urlparse

SCHEMA = "gallery_curved_patches_v1"
TARGET_SCHEMAS = {
    "orthogonal_sketch_features_v1",
    "curved_sketch_layers_v1",
    "affine_image_v1",
}


def _path(base: Path, reference: str) -> Path:
    return base / urlparse(reference).path


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _dims(degrees):
    return tuple(int(d) + 1 for d in degrees)


def _grid(block):
    dims = _dims(block["degrees"])
    controls = block.get("controls", block.get("controls_homogeneous"))
    if len(controls) != math.prod(dims):
        raise ValueError(f"{block.get('id')}: control count does not match degrees")
    return [
        [
            [
                tuple(map(int, controls[(i * dims[1] + j) * dims[2] + k]))
                for k in range(dims[2])
            ]
            for j in range(dims[1])
        ]
        for i in range(dims[0])
    ]


def _flatten(grid):
    return [list(p) for plane in grid for row in plane for p in row]


def _knots(degree):
    return [[0] * (degree + 1) + [1] * (degree + 1)]


def _exact_patch(block):
    degrees = list(map(int, block["degrees"]))
    controls = [list(map(int, p)) for p in block["controls"]]
    if not degrees or any(d < 1 for d in degrees):
        raise ValueError(f"{block.get('id')}: degrees must be positive")
    if any(len(p) != 4 or p[3] <= 0 for p in controls):
        raise ValueError(
            f"{block.get('id')}: expected positive-weight homogeneous XYZW controls"
        )
    if len(controls) != math.prod(d + 1 for d in degrees):
        raise ValueError(f"{block.get('id')}: control count does not match degrees")
    # Clamped single-span NURBS extraction is the identity on the Bezier net.
    patch = {
        "id": block["id"],
        "degrees": degrees,
        "knots": [[0] * (d + 1) + [1] * (d + 1) for d in degrees],
        "controls_homogeneous": controls,
    }
    extracted = _flatten(_grid(patch))
    return patch, extracted == controls


def _face(block, side):
    grid = _grid(block)
    degrees = list(map(int, block["degrees"]))
    if len(side) != 2 or side[0] not in "uvw" or side[1] not in "01":
        raise ValueError(f"unsupported side {side!r}")
    axis = "uvw".index(side[0])
    index = 0 if side[1] == "0" else degrees[axis]
    if axis == 0:
        net = [
            [grid[index][j][k] for k in range(degrees[2] + 1)]
            for j in range(degrees[1] + 1)
        ]
    elif axis == 1:
        net = [
            [grid[i][index][k] for k in range(degrees[2] + 1)]
            for i in range(degrees[0] + 1)
        ]
    else:
        net = [
            [grid[i][j][index] for j in range(degrees[1] + 1)]
            for i in range(degrees[0] + 1)
        ]
    return net


def _orient(net, orientation):
    out = [row[:] for row in net]
    if orientation.get("swap"):
        if len({len(row) for row in out}) != 1:
            raise ValueError("cannot swap nonrectangular face net")
        out = [list(row) for row in zip(*out)]
    if orientation.get("flip_first"):
        out.reverse()
    if orientation.get("flip_second"):
        out = [list(reversed(row)) for row in out]
    return out


def _elevate_axis(net, axis, target):
    """Exactly elevate one Bernstein direction in homogeneous coordinates."""
    current = len(net) - 1 if axis == 0 else len(net[0]) - 1
    while current < target:
        if axis == 0:
            nxt = []
            for i in range(current + 2):
                row = []
                for j in range(len(net[0])):
                    if i == 0:
                        p = net[0][j]
                    elif i == current + 1:
                        p = net[-1][j]
                    else:
                        left, right = net[i - 1][j], net[i][j]
                        p = tuple(
                            Fraction(i, current + 1) * left[c]
                            + Fraction(current + 1 - i, current + 1) * right[c]
                            for c in range(4)
                        )
                    row.append(p)
                nxt.append(row)
        else:
            nxt = []
            for row in net:
                newrow = []
                for i in range(current + 2):
                    if i == 0:
                        p = row[0]
                    elif i == current + 1:
                        p = row[-1]
                    else:
                        left, right = row[i - 1], row[i]
                        p = tuple(
                            Fraction(i, current + 1) * left[c]
                            + Fraction(current + 1 - i, current + 1) * right[c]
                            for c in range(4)
                        )
                    newrow.append(p)
                nxt.append(newrow)
        net = nxt
        current += 1
    return net


def _elevate_pair(a, b):
    shape_a = (len(a), len(a[0]))
    shape_b = (len(b), len(b[0]))
    targets = (max(shape_a[0], shape_b[0]) - 1, max(shape_a[1], shape_b[1]) - 1)
    return _elevate_axis(_elevate_axis(a, 0, targets[0]), 1, targets[1]), _elevate_axis(
        _elevate_axis(b, 0, targets[0]), 1, targets[1]
    )


def _projective_equal(a, b):
    """Exact equality of homogeneous nets up to one common nonzero scale."""
    if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
        return False
    flat_a = [v for row in a for p in row for v in p]
    flat_b = [v for row in b for p in row for v in p]
    pivot = next((i for i, (x, y) in enumerate(zip(flat_a, flat_b)) if x or y), None)
    if pivot is None:
        return True
    x, y = flat_a[pivot], flat_b[pivot]
    if not x or not y:
        return False
    ratio = Fraction(x, y)
    return all(Fraction(u) == ratio * v for u, v in zip(flat_a, flat_b))


def _interface_result(model, iface, blocks):
    first_id, first_side = iface["first"]
    second_id, second_side = iface["second"]
    if first_id not in blocks or second_id not in blocks:
        return {"status": "unsupported", "reason": "interface references missing block"}
    orient = iface.get("orientation", {})
    try:
        a, b = (
            _face(blocks[first_id], first_side),
            _face(blocks[second_id], second_side),
        )
        b = _orient(b, orient)
        a, b = _elevate_pair(a, b)
        equal = _projective_equal(a, b)
        return {
            "status": "exact_match" if equal else "mismatch",
            "first": [first_id, first_side],
            "second": [second_id, second_side],
            "orientation": orient,
            "homogeneous_net_equal_up_to_common_scale": equal,
        }
    except (ValueError, IndexError) as exc:
        return {
            "status": "unsupported",
            "first": [first_id, first_side],
            "second": [second_id, second_side],
            "reason": str(exc),
        }


def build_model_patch_assembly(model: dict) -> dict:
    """Return a deterministic disjoint patch cover and exact trace audit for one model."""
    errors = []
    block_map = {}
    for block in model.get("blocks", []):
        try:
            _, exact = _exact_patch(block)
            block_map[block["id"]] = block
            if not exact:
                errors.append(
                    {
                        "block": block["id"],
                        "reason": "Bezier extraction did not reproduce controls",
                    }
                )
        except (ValueError, KeyError, TypeError) as exc:
            errors.append({"block": block.get("id", "?"), "reason": str(exc)})
    interfaces = [
        _interface_result(model, i, block_map)
        for i in model.get("topology", {}).get("interfaces", [])
    ]
    if any(i["status"] != "exact_match" for i in interfaces):
        errors.append(
            {
                "reason": "one or more source interfaces failed exact homogeneous trace verification"
            }
        )
    final, remapped, boundary = _finalize_model(
        model, block_map, model.get("topology", {}).get("interfaces", [])
    )
    owners = [bid for patch in final for bid in patch["source_blocks"]]
    if sorted(owners) != sorted(block_map) or len(owners) != len(set(owners)):
        errors.append(
            {"reason": "final patches do not form a disjoint source-block cover"}
        )
    return {
        "source_block_count": len(block_map),
        "patch_count": len(final),
        "interface_count": len(interfaces),
        "exact_extraction_count": len(block_map)
        - sum(e["reason"].startswith("Bezier extraction") for e in errors),
        "exact_interface_count": sum(i["status"] == "exact_match" for i in interfaces),
        "patches": final,
        "interfaces": remapped,
        "boundary_faces": boundary,
        "errors": errors,
    }


def build_catalog(experiments_path: str | Path) -> dict:
    experiments_path = Path(experiments_path)
    root = experiments_path.parent
    entries = []
    linear = []
    for experiment in json.loads(experiments_path.read_text()):
        model_path = _path(root, experiment["model"])
        model = json.loads(model_path.read_bytes())
        is_curved = model.get("schema") in TARGET_SCHEMAS and any(
            max(b.get("degrees", [1])) > 1 for b in model.get("blocks", [])
        )
        patches = model.get("nurbs_patches", [])
        is_polygon = (
            model.get("schema") == "sketch_layer_partition_v1"
            and bool(patches)
            and all(len(patch["cells"]) == 1 for patch in patches)
        )
        if not (is_curved or is_polygon):
            continue
        assembly = build_model_patch_assembly(model)
        if assembly["errors"]:
            raise ValueError(
                f"Invalid source for {experiment['id']}: {assembly['errors']}"
            )
        final = assembly["patches"]
        used = {
            bid
            for p in final
            if len(p["source_blocks"]) > 1
            for bid in p["source_blocks"]
        }
        item = {
            "id": experiment["id"],
            "model_sha256": _sha(model_path),
            "recipe_sha256": _sha(_path(root, experiment["source"]))
            if experiment.get("source") and _path(root, experiment["source"]).is_file()
            else None,
            "source_schema": model.get("schema"),
            **assembly,
            "merged_source_blocks": sorted(used),
        }
        (entries if is_curved else linear).append(item)
    return {
        "schema": SCHEMA,
        "method": "exact_homogeneous_bezier_to_clamped_nurbs",
        "scope_note": "Final patch sets are disjoint unions of exact single-span Bezier volumes and compatible two-span clamped NURBS volumes. Exact piece extraction and interface traces are checked; no global injectivity, shape, or material coverage certificate is implied.",
        "entries": entries,
        "linear_patch_assemblies": linear,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiments", default="docs/site/experiments.json")
    parser.add_argument("--output", default="build/gallery-curved-patches.json")
    args = parser.parse_args()
    result = build_catalog(args.experiments)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"wrote {len(result['entries'])} models to {output}")
    for item in result["entries"] + result["linear_patch_assemblies"]:
        merged = sum(len(p["source_blocks"]) > 1 for p in item["patches"])
        print(
            f"{item['id']}: source_blocks={item['source_block_count']} patches={item['patch_count']} interfaces={item['interface_count']} extracted={item['exact_extraction_count']} merged={merged}"
        )


def _elevate_volume(net, axis, target):
    # Nets use axis 0/1/2 as nested list dimensions.
    def get_shape(n):
        return (len(n), len(n[0]), len(n[0][0]))

    while get_shape(net)[axis] - 1 < target:
        shape = get_shape(net)
        degree = shape[axis] - 1
        nxt_shape = list(shape)
        nxt_shape[axis] += 1
        nxt = [
            [[None for _ in range(nxt_shape[2])] for _ in range(nxt_shape[1])]
            for _ in range(nxt_shape[0])
        ]
        for i in range(nxt_shape[0]):
            for j in range(nxt_shape[1]):
                for k in range(nxt_shape[2]):
                    idx = [i, j, k]
                    t = idx[axis]
                    if t == 0:
                        lo = hi = 0
                        alpha = Fraction(0)
                    elif t == degree + 1:
                        lo = hi = degree
                        alpha = Fraction(0)
                    else:
                        lo = t - 1
                        hi = t
                        alpha = Fraction(t, degree + 1)
                    a = idx.copy()
                    b = idx.copy()
                    a[axis] = lo
                    b[axis] = hi
                    pa = net[a[0]][a[1]][a[2]]
                    pb = net[b[0]][b[1]][b[2]]
                    nxt[i][j][k] = tuple(
                        alpha * pa[c] + (1 - alpha) * pb[c] for c in range(4)
                    )
        net = nxt
    return net


def _permute_block(block, normal_axis, transverse_order, flips):
    old = _grid(block)
    old_deg = list(map(int, block["degrees"]))
    axes = [normal_axis, *transverse_order]
    shape = [old_deg[a] + 1 for a in axes]
    out = [
        [[None for _ in range(shape[2])] for _ in range(shape[1])]
        for _ in range(shape[0])
    ]
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                coords = [i, j, k]
                for q, flip in enumerate(flips):
                    if flip:
                        coords[q] = shape[q] - 1 - coords[q]
                src = [0, 0, 0]
                for q, ax in enumerate(axes):
                    src[ax] = coords[q]
                out[i][j][k] = old[src[0]][src[1]][src[2]]
    return out, [old_deg[a] for a in axes], axes, list(flips)


def _band_pair(a, b, iface):
    """Return exact two-span NURBS union after signed axis permutation."""
    first, second = a, b
    s1, s2 = iface["first"], iface["second"]
    o = iface.get("orientation", {})
    axis1 = "uvw".index(s1[1][0])
    axis2 = "uvw".index(s2[1][0])
    t1 = [x for x in range(3) if x != axis1]
    t2 = [x for x in range(3) if x != axis2]
    map2 = [t2[1], t2[0]] if o.get("swap") else t2
    ga, da, axa, flipa = _permute_block(first, axis1, t1, [False, False, False])
    gb, db, axb, flipb = _permute_block(
        second,
        axis2,
        map2,
        [False, o.get("flip_first", False), o.get("flip_second", False)],
    )
    # Put the first member's shared face at the end and second member's at the start.
    if s1[1][1] == "0":
        ga.reverse()
        flipa[0] = not flipa[0]
    if s2[1][1] == "1":
        gb.reverse()
        flipb[0] = not flipb[0]

    def parity(axes, flips):
        inversions = sum(axes[i] > axes[j] for i in range(3) for j in range(i + 1, 3))
        return (-1 if inversions % 2 else 1) * (-1 if sum(flips) % 2 else 1)

    # A common transverse reversal preserves the seam while making both source-to-patch maps positive.
    if parity(axa, flipa) < 0:
        ga = [plane[::-1] for plane in ga]
        gb = [plane[::-1] for plane in gb]
        flipa[1] = not flipa[1]
        flipb[1] = not flipb[1]
    if parity(axa, flipa) < 0 or parity(axb, flipb) < 0:
        return None, "cannot preserve positive orientation on both source spans"
    targets = [max(da[i], db[i]) for i in range(3)]
    for q in range(3):
        ga = _elevate_volume(ga, q, targets[q])
        gb = _elevate_volume(gb, q, targets[q])
    sa = ga[-1]
    sb = gb[0]
    if not _projective_equal(sa, sb):
        return None, "elevated homogeneous seam mismatch"
    af = [v for row in sa for p in row for v in p]
    bf = [v for row in sb for p in row for v in p]
    pivot = next((i for i, (x, y) in enumerate(zip(af, bf)) if x or y), None)
    if pivot is None:
        return None, "zero seam net"
    scale = Fraction(af[pivot], bf[pivot])
    gb = [
        [[tuple(Fraction(v) * scale for v in p) for p in row] for row in plane]
        for plane in gb
    ]
    merged = ga + gb[1:]
    deg = targets
    shape = [len(merged), len(merged[0]), len(merged[0][0])]
    join_degree = deg[0]
    knots = [
        [0] * (join_degree + 1) + [1] * join_degree + [2] * (join_degree + 1),
        *[[0] * (d + 1) + [1] * (d + 1) for d in deg[1:]],
    ]
    denominator = math.lcm(
        *(
            Fraction(v).denominator
            for plane in merged
            for row in plane
            for p in row
            for v in p
        )
    )
    merged = [
        [
            [tuple(int(Fraction(v) * denominator) for v in p) for p in row]
            for row in plane
        ]
        for plane in merged
    ]
    patch = {
        "id": f"band:{first['id']}+{second['id']}",
        "source_blocks": [first["id"], second["id"]],
        "cells": [first["id"], second["id"]],
        "patch_kind": "two_span_nurbs",
        "degrees": deg,
        "knots": knots,
        "control_shape": shape,
        "controls_homogeneous": [
            list(p) for plane in merged for row in plane for p in row
        ],
        "homogeneous_common_scale_denominator": denominator,
        "spans": [
            {
                "source_block": first["id"],
                "axis_map": axa,
                "flips": flipa,
                "degree_elevation": [targets[q] - da[q] for q in range(3)],
                "parameter_interval": [0, 1],
            },
            {
                "source_block": second["id"],
                "axis_map": axb,
                "flips": flipb,
                "degree_elevation": [targets[q] - db[q] for q in range(3)],
                "parameter_interval": [1, 2],
            },
        ],
        "exact_piece_extraction_expected": True,
    }
    return patch, None


def _finalize_model(model, blocks, interfaces):
    used = set()
    merged = []
    pair_counts = {}
    for iface in interfaces:
        x = iface["first"][0]
        y = iface["second"][0]
        pair_counts[frozenset((x, y))] = pair_counts.get(frozenset((x, y)), 0) + 1
    for iface in interfaces:
        if iface.get("status") not in (None, "exact_match"):
            continue
        a_id, a_side = iface["first"]
        b_id, b_side = iface["second"]
        if (
            a_id not in blocks
            or b_id not in blocks
            or a_id == b_id
            or a_id in used
            or b_id in used
        ):
            continue
        if pair_counts.get(frozenset((a_id, b_id)), 0) != 1:
            continue
        candidate, reason = _band_pair(blocks[a_id], blocks[b_id], iface)
        if candidate is None:
            continue
        used.update((a_id, b_id))
        merged.append(candidate)
    final = []
    owner = {}
    for p in merged:
        final.append(p)
        for bid in p["source_blocks"]:
            owner[bid] = p["id"]
    for bid, b in blocks.items():
        if bid in used:
            continue
        patch, ok = _exact_patch(b)
        patch.update(
            {
                "source_blocks": [bid],
                "cells": [bid],
                "patch_kind": "single_bezier",
                "control_shape": list(_dims(patch["degrees"])),
                "spans": [
                    {
                        "source_block": bid,
                        "axis_map": [0, 1, 2],
                        "flips": [False] * 3,
                        "degree_elevation": [0] * 3,
                        "parameter_interval": [0, 1],
                    }
                ],
            }
        )
        final.append(patch)
        owner[bid] = bid

    def mapped_side(bid, side):
        patch_id = owner[bid]
        patch = next(p for p in final if p["id"] == patch_id)
        span = next(x for x in patch["spans"] if x["source_block"] == bid)
        source_axis = "uvw".index(side[0])
        patch_axis = span["axis_map"].index(source_axis)
        endpoint = int(side[1])
        endpoint = 1 - endpoint if span["flips"][patch_axis] else endpoint
        return {
            "patch_id": patch_id,
            "source_block": bid,
            "side": side,
            "patch_side": "uvw"[patch_axis] + str(endpoint),
            "axis_map": span["axis_map"],
            "flips": span["flips"],
            "parameter_interval": span["parameter_interval"],
        }

    remapped = []
    for iface in interfaces:
        (a, sa), (b, sb) = iface["first"], iface["second"]
        if owner[a] == owner[b]:
            continue
        remapped.append(
            {
                "first": mapped_side(a, sa),
                "second": mapped_side(b, sb),
                "orientation": iface.get("orientation", {}),
                "status": "exact_match",
            }
        )
    boundary = [
        {**mapped_side(bid, side), "status": "extracted_exactly_from_model"}
        for bid, side in model.get("topology", {}).get("exterior", [])
        if bid in owner
    ]
    return final, remapped, boundary


if __name__ == "__main__":
    main()
