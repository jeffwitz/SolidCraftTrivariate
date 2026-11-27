#!/usr/bin/env python3
"""Publish exact NURBS patch covers for the seven tall-cylinder Boolean states.

Exact maximal chains across the extrusion axis are merged first. Remaining
compatible cell pairs may then be assembled in any orientation; unmergeable
cells remain single-span NURBS volumes. The standard gallery patch-map verifier
independently checks every source coefficient, sampled evaluation, exterior
face, and surviving interface.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import subprocess
import tempfile
from fractions import Fraction
from pathlib import Path

STATES = ("box", "cylinder", "union", "intersection", "a-minus-b", "b-minus-a", "xor")
SCHEMA = "gallery_exact_patch_map_v1"
NOTE = (
    "Exact homogeneous extraction and deterministic axial-chain/pair planning. "
    "The native plan is checked against every source span and topology record. "
    "Candidate merges are not a global optimality claim. This is a representation "
    "and trace check, not a global injectivity or material-coverage certificate."
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_one(root: Path, state: str) -> dict:
    sys.path.insert(0, str(root / "scripts"))
    import gallery_curved_patches as curved
    from check_gallery_patch_maps import verify_patch_samples, verify_patch_topology

    model_path = root / "docs/site/box-cylinder/tall" / f"{state}.model.json"
    recipe_path = root / "docs/site/box-cylinder/tall" / f"{state}.recipe.json"
    model = json.loads(model_path.read_bytes())
    blocks = {block["id"]: block for block in model.get("blocks", [])}
    if not blocks:
        raise ValueError(f"{state}: model has no blocks")
    extracted = {}
    for bid, block in blocks.items():
        patch, exact = curved._exact_patch(block)
        if not exact:
            raise ValueError(f"{state}/{bid}: exact Bezier extraction changed controls")
        extracted[bid] = patch

    interfaces = model.get("topology", {}).get("interfaces", [])
    checked = [curved._interface_result(model, iface, blocks) for iface in interfaces]
    bad = [item for item in checked if item["status"] != "exact_match"]
    if bad:
        raise ValueError(f"{state}: source trace mismatch: {bad}")

    chains = _maximal_w_chains(curved, blocks, interfaces, checked)
    chained = {bid for chain in chains if len(chain) > 1 for bid in chain}
    patches = [
        _chain_patch(curved, blocks, chain) for chain in chains if len(chain) > 1
    ]

    # Let the established exact pair assembler consider all remaining
    # orientations after maximal extrusion chains have been reserved.
    residual_ids = [bid for bid in blocks if bid not in chained]
    residual_id_set = set(residual_ids)
    residual_blocks = {bid: blocks[bid] for bid in residual_ids}
    residual_model = dict(model)
    residual_model["topology"] = dict(model.get("topology", {}))
    residual_model["topology"]["interfaces"] = [
        iface
        for iface in interfaces
        if iface["first"][0] in residual_id_set
        and iface["second"][0] in residual_id_set
    ]
    residual_model["topology"]["exterior"] = [
        face
        for face in model.get("topology", {}).get("exterior", [])
        if face[0] in residual_id_set
    ]
    residual_patches, _, _ = curved._finalize_model(
        residual_model, residual_blocks, residual_model["topology"]["interfaces"]
    )
    patches.extend(residual_patches)

    owner, spans = {}, {}
    for patch in patches:
        for span in patch["spans"]:
            bid = span["source_block"]
            owner[bid], spans[bid] = patch["id"], span

    def mapped_side(bid, side):
        patch = next(p for p in patches if p["id"] == owner[bid])
        span = spans[bid]
        source_axis = "uvw".index(side[0])
        patch_axis = span["axis_map"].index(source_axis)
        endpoint = int(side[1])
        if span["flips"][patch_axis]:
            endpoint = 1 - endpoint
        return {
            "patch_id": patch["id"],
            "source_block": bid,
            "side": side,
            "patch_side": "uvw"[patch_axis] + str(endpoint),
            "axis_map": span["axis_map"],
            "flips": span["flips"],
            "parameter_interval": span["parameter_interval"],
        }

    mapped_interfaces = []
    for iface in interfaces:
        (a, sa), (b, sb) = iface["first"], iface["second"]
        if owner[a] != owner[b]:
            mapped_interfaces.append(
                {
                    "first": mapped_side(a, sa),
                    "second": mapped_side(b, sb),
                    "orientation": iface.get("orientation", {}),
                    "status": "exact_match",
                }
            )
    boundary_faces = [
        {**mapped_side(bid, side), "status": "extracted_exactly_from_model"}
        for bid, side in model.get("topology", {}).get("exterior", [])
    ]

    source_cells = set(blocks)
    listed_cells = [cell for patch in patches for cell in patch["cells"]]
    if set(listed_cells) != source_cells or len(listed_cells) != len(source_cells):
        raise ValueError(f"{state}: patch cover does not own each source block once")
    sampled = verify_patch_samples(patches, model)
    verify_patch_topology(
        {
            "patches": patches,
            "interfaces": mapped_interfaces,
            "boundary_faces": boundary_faces,
        },
        model,
    )
    merged = [
        cell for patch in patches if len(patch["cells"]) > 1 for cell in patch["cells"]
    ]
    output = {
        "schema": SCHEMA,
        "case_id": f"box-cylinder-tall-{state}",
        "representation_kind": "exact_tall_cylinder_axial_nurbs_cover",
        "method": "exact_homogeneous_bezier_to_clamped_nurbs_chain_and_pair_assembly",
        "scope_note": NOTE,
        "source_model": f"{state}.model.json",
        "source_model_sha256": sha(model_path),
        "source_recipe": f"{state}.recipe.json",
        "source_recipe_sha256": sha(recipe_path),
        "source_schema": model.get("schema"),
        "source_block_count": len(blocks),
        "native_nurbs_patch_count": len(blocks),
        "patch_count": len(patches),
        "exact_extraction_count": len(blocks),
        "exact_interface_count": len(checked),
        "source_interface_count": len(interfaces),
        "independent_sample_evaluations": sampled,
        "axial_chain_patch_count": sum(
            patch["id"].startswith("axial:") for patch in patches
        ),
        "axial_chain_source_block_count": sum(
            len(patch["cells"]) for patch in patches if patch["id"].startswith("axial:")
        ),
        "additional_exact_pair_patch_count": sum(
            patch["patch_kind"] == "two_span_nurbs" for patch in patches
        ),
        "patches": patches,
        "interfaces": mapped_interfaces,
        "boundary_faces": boundary_faces,
        "merged_source_blocks": sorted(merged),
        "errors": [],
    }
    return output


def _maximal_w_chains(curved, blocks, interfaces, audits):
    """Find deterministic, nonbranching chains across exact w1/w0 seams."""
    forward, predecessor = {}, {}
    for iface, audit in zip(interfaces, audits):
        if audit["status"] != "exact_match":
            continue
        (a, sa), (b, sb) = iface["first"], iface["second"]
        orientation = iface.get("orientation", {})
        if any(
            orientation.get(k, False) for k in ("swap", "flip_first", "flip_second")
        ):
            continue
        if sa == "w1" and sb == "w0":
            lo, hi = a, b
        elif sa == "w0" and sb == "w1":
            lo, hi = b, a
        else:
            continue
        left, right = blocks[lo], blocks[hi]
        if left["degrees"][2] != right["degrees"][2]:
            continue
        if left["degrees"][:2] != right["degrees"][:2]:
            continue
        if not curved._projective_equal(
            curved._face(left, "w1"), curved._face(right, "w0")
        ):
            continue
        if lo in forward or hi in predecessor:
            continue
        forward[lo], predecessor[hi] = hi, lo
    chains, consumed = [], set()
    starts = sorted(set(forward) - set(predecessor))
    for start in starts:
        chain, current, seen = [start], start, {start}
        while current in forward:
            nxt = forward[current]
            if nxt in seen:
                break
            chain.append(nxt)
            seen.add(nxt)
            current = nxt
        if len(chain) > 1:
            chains.append(chain)
            consumed.update(chain)
    return chains


def _chain_patch(curved, blocks, chain):
    """Make a maximally chained NURBS patch with w as its first axis."""
    nets = []
    previous_seam = None
    for bid in chain:
        block = blocks[bid]
        grid = curved._grid(block)
        # Reorder source (u,v,w) to patch (w,u,v).
        net = [
            [[grid[i][j][k] for j in range(len(grid[0]))] for i in range(len(grid))]
            for k in range(len(grid[0][0]))
        ]
        if previous_seam is not None:
            current_seam = net[0]
            if not curved._projective_equal(previous_seam, current_seam):
                raise ValueError(f"{bid}: chain seam lost exact equality")
            flat_a = [x for row in previous_seam for point in row for x in point]
            flat_b = [x for row in current_seam for point in row for x in point]
            index = next(i for i, (a, b) in enumerate(zip(flat_a, flat_b)) if a or b)
            factor = Fraction(flat_a[index], flat_b[index])
            net = [
                [
                    [tuple(Fraction(x) * factor for x in point) for point in row]
                    for row in plane
                ]
                for plane in net
            ]
            nets.extend(net[1:])
        else:
            nets.extend(net)
        previous_seam = net[-1]
    normal_degree = int(blocks[chain[0]]["degrees"][2])
    transverse_degrees = [int(x) for x in blocks[chain[0]]["degrees"][:2]]
    degree = [normal_degree, *transverse_degrees]
    denominator = math.lcm(
        *(
            Fraction(x).denominator
            for plane in nets
            for row in plane
            for point in row
            for x in point
        )
    )
    integer_nets = [
        [
            [[int(Fraction(x) * denominator) for x in point] for point in row]
            for row in plane
        ]
        for plane in nets
    ]
    spans = [
        {
            "source_block": bid,
            "axis_map": [2, 0, 1],
            "flips": [False] * 3,
            "degree_elevation": [0] * 3,
            "parameter_interval": [i, i + 1],
        }
        for i, bid in enumerate(chain)
    ]
    return {
        "id": "axial:" + "+".join(chain),
        "source_blocks": list(chain),
        "cells": list(chain),
        "patch_kind": "multi_span_nurbs",
        "degrees": degree,
        "knots": [
            [0] * (degree[0] + 1)
            + [k for k in range(1, len(chain)) for _ in range(degree[0])]
            + [len(chain)] * (degree[0] + 1),
            *[[0] * (d + 1) + [1] * (d + 1) for d in degree[1:]],
        ],
        "control_shape": [
            len(integer_nets),
            len(integer_nets[0]),
            len(integer_nets[0][0]),
        ],
        "controls_homogeneous": [
            point for plane in integer_nets for row in plane for point in row
        ],
        "homogeneous_common_scale_denominator": denominator,
        "spans": spans,
        "exact_piece_extraction_expected": True,
    }


def verify_cylinder_seams(cover, model, expected_segments):
    """Classify each original same-patch face as an internal knot or an identified end."""
    from collections import Counter

    locations = {}
    for patch in cover["patches"]:
        for span in patch["spans"]:
            if span["axis_map"] != [0, 1, 2] or any(span["flips"]):
                raise ValueError("Unexpected cylinder tensor orientation")
            locations[span["source_block"]] = (patch, span["parameter_intervals"])

    expected = Counter()
    for iface in model["topology"]["interfaces"]:
        (cell_a, side_a), (cell_b, side_b) = iface["first"], iface["second"]
        pa, ranges_a = locations[cell_a]
        pb, ranges_b = locations[cell_b]
        if pa["id"] != pb["id"]:
            continue
        axis = "uvw".index(side_a[0])
        if (
            side_a[0] != side_b[0]
            or side_a[1] == side_b[1]
            or any(iface["orientation"].values())
        ):
            raise ValueError("Undeclared axis change inside cylinder tensor")
        if any(ranges_a[d] != ranges_b[d] for d in range(3) if d != axis):
            raise ValueError(
                "Seam spans do not coincide in their tangential coordinates"
            )
        ends = [ranges_a[axis][int(side_a[1])], ranges_b[axis][int(side_b[1])]]
        if ends[0] != ends[1]:
            if sorted(ends) != [0, pa["knots"][axis][-1]] or axis != 0:
                raise ValueError("Invalid identification across a cylinder patch")
            expected[tuple(sorted(((cell_a, side_a), (cell_b, side_b))))] += 1

    actual = Counter()
    for seam in cover.get("self_seams", []):
        faces = [seam[key] for key in ("first", "second")]
        if any(
            face["patch_id"] != locations[face["source_block"]][0]["id"]
            or face["patch_side"] != face["side"]
            for face in faces
        ) or any(seam["orientation"].values()):
            raise ValueError("Invalid cylinder self seam")
        actual[
            tuple(sorted((face["source_block"], face["side"]) for face in faces))
        ] += 1
    if expected != actual or sum(actual.values()) != expected_segments:
        raise ValueError("Cylinder cover omits or duplicates a required self seam")


def publish_cylinder(root, outdir, check):
    """Use the native periodic cylinder planner and independently verify it."""
    from check_gallery_patch_maps import verify_patch_samples, verify_patch_topology

    model_path = outdir / "cylinder.model.json"
    model = json.loads(model_path.read_bytes())
    recipe_path = outdir / "cylinder.recipe.json"
    cli = root / "target/release/solidcraft-trivariate"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for count in (2, 3):
            filename = f"cylinder-{count}.patches.json"
            temp_path = tmp / filename
            subprocess.run(
                [cli, "cover-cylinder", str(model_path), str(count), str(temp_path)],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            cover = json.loads(temp_path.read_bytes())
            cover["source_model_sha256"] = sha(model_path)
            cover["source_recipe_sha256"] = sha(recipe_path)
            data = (json.dumps(cover, separators=(",", ":"), sort_keys=True) + "\n").encode()
            if len(cover["patches"]) != count:
                raise ValueError(f"Unexpected cylinder cover size for {filename}")
            verify_patch_samples(cover["patches"], model)
            verify_patch_topology(cover, model)
            verify_cylinder_seams(cover, model, 3 if count == 2 else 0)
            destinations = [outdir / filename]
            if count == 2:
                destinations.append(outdir / "cylinder.patches.json")
            for path in destinations:
                if check:
                    if path.read_bytes() != data:
                        raise ValueError(f"Stale cylinder cover: {path}")
                else:
                    path.write_bytes(data)
            print(
                f"cylinder: 15 source cells extracted exactly from {count} NURBS volumes"
            )


def assemble_native(root: Path, state: str, compare_legacy: bool = False) -> dict:
    """Plan and assemble a tall Boolean cover natively, then verify independently."""
    from check_gallery_patch_maps import verify_patch_samples, verify_patch_topology

    model_path = root / "docs/site/box-cylinder/tall" / f"{state}.model.json"
    recipe_path = root / "docs/site/box-cylinder/tall" / f"{state}.recipe.json"
    model = json.loads(model_path.read_bytes())
    cli = root / "target/release/solidcraft-trivariate"
    if not cli.is_file():
        raise FileNotFoundError(f"release CLI is required for native patch assembly: {cli}")
    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)
        result_path = temporary / f"{state}.patches.json"
        subprocess.run(
            [cli, "plan-patches", model_path, result_path],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        native = json.loads(result_path.read_bytes())

    sample_count = verify_patch_samples(native["patches"], model)
    verify_patch_topology(native, model)
    patches = native["patches"]
    if compare_legacy:
        # The historical Python algorithm is retained as a validation-only
        # comparison. Its nets never contribute to published output.
        legacy = build_one(root, state)
        verify_patch_samples(legacy["patches"], model)
        verify_patch_topology(legacy, model)
        legacy_groups = {frozenset(patch["cells"]) for patch in legacy["patches"]}
        native_groups = {frozenset(patch["cells"]) for patch in patches}
        print(
            f"{state}: legacy partition "
            f"{'matches' if legacy_groups == native_groups else 'differs from'} native"
        )

    # Preserve the versioned public map contract while making the Rust nets,
    # ordering, knots, orientation and source maps the sole generated data.
    merged = sorted(
        cell for patch in patches if len(patch["cells"]) > 1 for cell in patch["cells"]
    )
    result = {
        "schema": SCHEMA,
        "case_id": f"box-cylinder-tall-{state}",
        "representation_kind": "exact_tall_cylinder_axial_nurbs_cover",
        "method": "native_exact_pair_and_axial_chain_planning",
        "scope_note": NOTE,
        "source_model": f"{state}.model.json",
        "source_model_sha256": sha(model_path),
        "source_recipe": f"{state}.recipe.json",
        "source_recipe_sha256": sha(recipe_path),
        "source_schema": model.get("schema"),
        "source_block_count": len(model.get("blocks", [])),
        "native_nurbs_patch_count": len(model.get("blocks", [])),
        "patch_count": len(patches),
        "exact_extraction_count": len(model.get("blocks", [])),
        "exact_interface_count": len(model.get("topology", {}).get("interfaces", [])),
        "source_interface_count": len(model.get("topology", {}).get("interfaces", [])),
        "independent_sample_evaluations": sample_count,
        "axial_chain_patch_count": sum(patch["id"].startswith("axial:") for patch in patches),
        "axial_chain_source_block_count": sum(
            len(patch["cells"]) for patch in patches if patch["id"].startswith("axial:")
        ),
        "additional_exact_pair_patch_count": sum(
            patch["patch_kind"] == "two_span_nurbs" for patch in patches
        ),
        "patches": patches,
        "interfaces": native["interfaces"],
        "boundary_faces": native["boundary_faces"],
        "self_seams": native.get("self_seams", []),
        "merged_source_blocks": merged,
        "errors": [],
    }
    verify_patch_samples(result["patches"], model)
    verify_patch_topology(result, model)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--states", nargs="*", default=STATES, choices=STATES)
    parser.add_argument(
        "--check",
        action="store_true",
        help="recompute every map and fail if a published payload is stale",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    outdir = root / "docs/site/box-cylinder/tall"
    if not args.check:
        outdir.mkdir(parents=True, exist_ok=True)
    for state in args.states:
        if state == "cylinder":
            publish_cylinder(root, outdir, args.check)
            continue
        payload = assemble_native(root, state, compare_legacy=args.check)
        path = outdir / f"{state}.patches.json"
        encoded = json.dumps(payload, separators=(",", ":"), sort_keys=True) + "\n"
        if args.check:
            if not path.is_file() or path.read_text() != encoded:
                raise SystemExit(f"stale or missing patch map: {path}")
        else:
            path.write_text(encoded)
        print(
            f"{state}: cells={payload['source_block_count']} patches={payload['patch_count']} "
            f"interfaces={payload['source_interface_count']} "
            f"axial_chains={payload['axial_chain_patch_count']} "
            f"axial_cells={payload['axial_chain_source_block_count']} "
            f"additional_pairs={payload['additional_exact_pair_patch_count']} "
            f"samples={payload['independent_sample_evaluations']} "
            f"{'checked' if args.check else 'wrote'} {path}"
        )


if __name__ == "__main__":
    main()
