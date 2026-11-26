#!/usr/bin/env python3
"""Independent paired benchmark for polygonal sketch NURBS partitions.

The script treats recipes as immutable inputs, invokes the same CLI in its
default and ``nurbs_partition=optimize`` modes, and inspects serialized models
without calling the Rust audit. Exact homogeneous comparisons use Fractions;
sampled Jacobian metrics are diagnostics only.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frac(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        return Fraction(value)
    if isinstance(value, float):
        return Fraction(str(value))
    raise TypeError(f"unsupported exact scalar {value!r}")


def point_key(control):
    if len(control) != 4:
        raise ValueError("homogeneous point must have four coordinates")
    q = tuple(frac(x) for x in control)
    if q[3] == 0:
        raise ValueError("homogeneous point has zero weight")
    # Euclidean coordinates make this invariant under any nonzero common scale.
    return tuple(x / q[3] for x in q[:3])


def normalized_net(controls):
    """Return an oriented homogeneous net modulo a single common scale."""
    rows = [tuple(frac(x) for x in row) for row in controls]
    if not rows or any(len(row) != 4 for row in rows):
        raise ValueError("expected a nonempty homogeneous control net")
    scale = next((x for row in rows for x in row if x), None)
    if scale is None:
        raise ValueError("zero homogeneous control net")
    return tuple(tuple(x / scale for x in row) for row in rows)


def invoke(cli: Path, recipe_path: Path, mode: str, out: Path, timeout: int):
    recipe = json.loads(recipe_path.read_text())
    if mode == "optimize":
        recipe["nurbs_partition"] = "optimize"
    elif mode == "assemble":
        recipe["nurbs_partition"] = "assemble"
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="sct-nurbs-recipe-") as tmp:
        source = Path(tmp) / "recipe.json"
        source.write_text(json.dumps(recipe, separators=(",", ":")))
        proc = subprocess.run([str(cli), "build", str(source), str(out)],
                              capture_output=True, text=True, timeout=timeout)
    elapsed = time.monotonic() - started
    if proc.returncode:
        return {"status": "rejected", "elapsed_s": elapsed,
                "returncode": proc.returncode,
                "diagnostic": (proc.stderr or proc.stdout).strip()[-1200:]}
    return {"status": "built", "elapsed_s": elapsed, "model_sha256": digest(out)}


def block_metrics(model):
    """Use the established independent spectral and thickness diagnostics."""
    from evaluate_polygon_extrusion_quality import summarize_model
    return summarize_model(model)


def face_net(block, side):
    axis = "uvw".index(side[0])
    fixed = int(side[1])
    varying = [a for a in range(3) if a != axis]
    values = []
    for a in (0, 1):
        for b in (0, 1):
            bits = [0, 0, 0]
            bits[axis] = fixed
            bits[varying[0]], bits[varying[1]] = a, b
            values.append(tuple(frac(x) for x in block["controls"][(bits[0]*2+bits[1])*2+bits[2]]))
    return values


def orient_face(net, orientation, first):
    # Incidence flags map the first face's parameterization into the second's.
    if not first:
        return net
    flip_a = orientation["flip_first"]
    flip_b = orientation["flip_second"]
    swap = orientation["swap"]
    out = []
    for a in (0, 1):
        for b in (0, 1):
            x, y = (b, a) if swap else (a, b)
            if flip_a:
                x = 1 - x
            if flip_b:
                y = 1 - y
            out.append(net[x*2+y])
    return out


def conformity_check(model):
    by_id = {b.get("id"): b for b in model.get("blocks", [])}
    checked = 0
    mismatches = []
    for interface in model.get("interfaces", model.get("topology", {}).get("interfaces", [])):
        try:
            left_id, left_side = interface["first"]
            right_id, right_side = interface["second"]
            left = face_net(by_id[left_id], left_side)
            right = face_net(by_id[right_id], right_side)
            orient = interface["orientation"]
            a = normalized_net(orient_face(left, orient, True))
            b = normalized_net(orient_face(right, orient, False))
            checked += 1
            if a != b:
                mismatches.append([left_id, left_side, right_id, right_side])
        except (KeyError, IndexError, ValueError, TypeError) as exc:
            mismatches.append({"interface": interface, "error": str(exc)})
    return {"interfaces_checked_exactly": checked,
            "interface_control_mismatch_count": len(mismatches),
            "mismatch_examples": mismatches[:8]}


def exact_extraction_check(model):
    """Check each consecutive section pair against all eight block controls.

    Each section is row-major in ``(r,z)``; consecutive sections span ``s``.
    Only the four positive XY rotations are admitted, with z preserved. The
    rotation is derived from control coordinates and then checked against the
    serialized ``source_xy_rotations`` metadata.
    """
    blocks = {b.get("id"): b for b in model.get("blocks", [])}
    patches = model.get("nurbs_patches", [])
    spans = 0
    failures = []
    referenced = []
    rotation_map = (
        lambda s, r: (s, r),
        lambda s, r: (1-r, s),
        lambda s, r: (1-s, 1-r),
        lambda s, r: (r, 1-s),
    )
    rotation_results = []
    for pi, patch in enumerate(patches):
        sections = patch.get("sections", [])
        cells = patch.get("cells", [])
        if len(cells) != max(0, len(sections) - 1):
            failures.append(f"patch {pi}: {len(sections)} sections do not define {len(cells)} consecutive cells")
            continue
        for si, block_id in enumerate(cells):
            spans += 1
            referenced.append(block_id)
            if block_id not in blocks:
                failures.append(f"patch {pi} span {si}: unknown cell {block_id!r}")
                continue
            lo, hi = sections[si], sections[si + 1]
            if len(lo) != 4 or len(hi) != 4:
                failures.append(f"patch {pi} span {si}: each end section must have four controls")
                continue
            source = [tuple(frac(x) for x in q) for q in list(lo) + list(hi)]
            target = [tuple(frac(x) for x in q) for q in blocks[block_id].get("controls", [])]
            if len(target) != 8 or any(q[3] <= 0 for q in source + target):
                failures.append(f"patch {pi} span {si}: expected eight positive-weight controls")
                continue
            try:
                target_norm = normalized_net(target)
            except (ValueError, TypeError, ZeroDivisionError) as exc:
                failures.append(f"patch {pi} span {si}: {exc}")
                continue
            matching_rotations = []
            for rotation, transform in enumerate(rotation_map):
                arranged = [None] * 8
                for s, r, z in product((0, 1), repeat=3):
                    u, v = transform(s, r)
                    src_index = (s * 2 + r) * 2 + z
                    dst_index = (u * 2 + v) * 2 + z
                    arranged[dst_index] = source[src_index]
                try:
                    if normalized_net(arranged) == target_norm:
                        matching_rotations.append(rotation)
                except (ValueError, TypeError, ZeroDivisionError):
                    continue
            if not matching_rotations:
                failures.append(f"patch {pi} span {si} ({block_id}): no exact positive-orientation extraction modulo one common scale")
                rotation_results.append(None)
                continue
            if len(matching_rotations) > 1:
                failures.append(f"patch {pi} span {si} ({block_id}): ambiguous exact XY rotation {matching_rotations}")
            inferred = matching_rotations[0]
            rotation_results.append(inferred)
            declared = patch.get("source_xy_rotations")
            if not isinstance(declared, list) or len(declared) != len(cells):
                failures.append(f"patch {pi}: source_xy_rotations must contain one integer per span")
            elif declared[si] != inferred:
                failures.append(f"patch {pi} span {si}: declared source XY rotation {declared[si]!r} != independently inferred {inferred}")
    counts = Counter(referenced)
    duplicates = [key for key, count in counts.items() if count != 1]
    missing = sorted(set(blocks) - set(referenced)) if patches else []
    if duplicates:
        failures.append(f"extracted spans do not reference cells exactly once: {duplicates[:8]}")
    if missing:
        failures.append(f"extraction omits {len(missing)} cells: {missing[:8]}")
    return {"exact_spans_checked": spans,
            "positive_xy_rotations_tested": 4,
            "inferred_source_xy_rotations": rotation_results,
            "referenced_cell_count": len(set(referenced)),
            "unreferenced_cell_count": len(missing),
            "extraction_failures": failures}


def exact_q1_volume(model):
    """Integrate the trilinear Jacobian determinant exactly with rational Simpson."""
    total = Fraction(0)
    samples = ((Fraction(0), Fraction(1, 6)),
               (Fraction(1, 2), Fraction(4, 6)),
               (Fraction(1), Fraction(1, 6)))
    for block in model.get("blocks", []):
        controls = block.get("controls", [])
        if block.get("degrees") != [1, 1, 1] or len(controls) != 8:
            raise ValueError(f"{block.get('id')}: exact volume requires eight Q1 controls")
        points = []
        for c in controls:
            h = tuple(frac(x) for x in c)
            if h[3] == 0:
                raise ValueError(f"{block.get('id')}: zero control weight")
            points.append(tuple(x / h[3] for x in h[:3]))
        for u, wu in samples:
            for v, wv in samples:
                for w, ww in samples:
                    cols = [[Fraction(0) for _ in range(3)] for _ in range(3)]
                    for i, j, k in product((0, 1), repeat=3):
                        p = points[(i * 2 + j) * 2 + k]
                        f = (v if j else 1-v, w if k else 1-w, u if i else 1-u)
                        signs = (1 if i else -1, 1 if j else -1, 1 if k else -1)
                        factors = (signs[0] * f[0] * f[1],
                                   signs[1] * f[2] * f[1],
                                   signs[2] * f[2] * f[0])
                        for a in range(3):
                            for d in range(3):
                                cols[a][d] += factors[a] * p[d]
                    a, b, c = cols
                    det = (a[0]*(b[1]*c[2]-b[2]*c[1])
                           - a[1]*(b[0]*c[2]-b[2]*c[0])
                           + a[2]*(b[0]*c[1]-b[1]*c[0]))
                    total += wu * wv * ww * det
    return total


def inspect_model(path: Path):
    model = json.loads(path.read_text())
    blocks = model.get("blocks", [])
    ids = {b.get("id") for b in blocks}
    interfaces = model.get("interfaces", model.get("topology", {}).get("interfaces", []))
    patches = model.get("nurbs_patches", [])
    patch_errors = []
    exact_keys = set()
    cell_ids = []
    for pi, patch in enumerate(patches):
        degrees = patch.get("degrees")
        knots = patch.get("knots")
        sections = patch.get("sections", [])
        if degrees != [1, 1, 1]:
            patch_errors.append(f"patch {pi}: degrees are not [1,1,1]")
        if not isinstance(knots, list) or len(knots) != 3:
            patch_errors.append(f"patch {pi}: expected three knot vectors")
        else:
            control_counts = [len(sections), 2, 2]
            for axis, (vector, p, count) in enumerate(zip(knots, degrees or [], control_counts)):
                try:
                    values = [frac(x) for x in vector]
                    if len(values) != count + p + 1:
                        raise ValueError(f"length {len(values)} != {count + p + 1}")
                    if values != sorted(values):
                        raise ValueError("knots are decreasing")
                    if values[0] != 0 or values[-1] != (count - p):
                        raise ValueError(f"domain [{values[0]}, {values[-1]}] is not [0, {count-p}]")
                    if sum(x == values[0] for x in values) != p + 1 or sum(x == values[-1] for x in values) != p + 1:
                        raise ValueError("end knots are not clamped")
                    interior = values[p+1:-(p+1)]
                    if any(interior.count(k) > p for k in set(interior)):
                        raise ValueError("internal knot multiplicity exceeds degree")
                except (ValueError, TypeError, ZeroDivisionError) as exc:
                    patch_errors.append(f"patch {pi} knot axis {axis}: {exc}")
        for section_index, section in enumerate(sections):
            if len(section) != 4:
                patch_errors.append(f"patch {pi} section {section_index}: expected four controls")
                continue
            try:
                net = normalized_net(section)
                exact_keys.add(net)
            except (ValueError, TypeError, ZeroDivisionError) as exc:
                patch_errors.append(f"patch {pi} section {section_index}: {exc}")
        for cell in patch.get("cells", []):
            cell_ids.append(cell)
            if cell not in ids:
                patch_errors.append(f"patch {pi}: unknown cell id {cell!r}")
    if patches:
        duplicate_cells = sorted(k for k, n in Counter(cell_ids).items() if n != 1)
        missing_cells = sorted(ids - set(cell_ids))
        if duplicate_cells:
            patch_errors.append(f"patch cell references are not unique: {duplicate_cells[:8]}")
        if missing_cells:
            patch_errors.append(f"patches omit {len(missing_cells)} model cells; examples: {missing_cells[:8]}")
    # Q1 scalar DOFs are the unique geometric corners of the trilinear blocks.
    dofs = set()
    for block in blocks:
        controls = block.get("controls", [])
        if len(controls) == 8 and block.get("degrees") == [1, 1, 1]:
            try:
                dofs.update(point_key(c) for c in controls)
            except (ValueError, TypeError, ZeroDivisionError) as exc:
                patch_errors.append(f"block {block.get('id')}: {exc}")
    return {
        "model_sha256": digest(path), "block_count": len(blocks),
        "interface_count": len(interfaces), "nurbs_patch_count": len(patches),
        "nurbs_section_count": sum(len(p.get("sections", [])) for p in patches),
        "nurbs_cell_reference_count": len(cell_ids),
        "nurbs_unique_oriented_section_nets_modulo_scale": len(exact_keys),
        "q1_coordinate_deduplicated_corner_count_dof_estimate": len(dofs),
        "exact_q1_volume_mm3": str(exact_q1_volume(model)),
        "exact_conformity": conformity_check(model),
        "exact_extraction": exact_extraction_check(model),
        "nurbs_patch_errors": patch_errors,
        "quality": block_metrics(model),
    }


def corpus_cases(corpus_dir: Path):
    files = sorted(corpus_dir.glob("bench-poly-*/*.recipe.json"))
    if not files:
        raise FileNotFoundError(f"no bench-poly recipe files under {corpus_dir}")
    return files


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", required=True, type=Path)
    parser.add_argument("--corpus", type=Path, default=ROOT / "docs/site/data")
    parser.add_argument("--work", type=Path, default=ROOT / "build/sketch-nurbs-evaluation")
    parser.add_argument("--output", type=Path, default=ROOT / "build/sketch-nurbs-evaluation.json")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--inspect-only", action="store_true",
                        help="replay exact checks over existing assembled model files")
    parser.add_argument("--reuse-existing", action="store_true",
                        help="reuse already written mode models and build only missing files")
    parser.add_argument("--expect-cases", type=int, default=46)
    args = parser.parse_args()
    cli = args.cli.resolve()
    corpus_dir = args.corpus.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    rows = []
    for recipe_path in corpus_cases(corpus_dir):
        case_id = recipe_path.name.removesuffix(".recipe.json")
        row = {"id": case_id, "recipe": str(recipe_path),
               "recipe_sha256": digest(recipe_path)}
        for mode in ("assemble", "optimize"):
            model_path = args.work / f"{case_id}.{mode}.model.json"
            if args.inspect_only or (args.reuse_existing and model_path.is_file()):
                if model_path.is_file():
                    run = {"status": "built", "model_sha256": digest(model_path),
                           "elapsed_s": 0.0, "replayed": True}
                else:
                    run = {"status": "missing", "diagnostic": f"saved model not found: {model_path}"}
            else:
                run = invoke(cli, recipe_path, mode, model_path, args.timeout)
            row[mode] = run
            if run["status"] == "built":
                row[mode]["model"] = inspect_model(model_path)
        rows.append(row)
        print(f"processed {len(rows)}/{len(corpus_cases(corpus_dir))}: {case_id}", flush=True)
    summary = {
        "schema": "sketch_nurbs_partition_benchmark_v1",
        "cli": str(cli), "cli_sha256": digest(cli),
        "corpus_directory": str(corpus_dir), "case_count": len(rows),
        "mode_counts": {mode: Counter(r[mode]["status"] for r in rows)
                        for mode in ("assemble", "optimize")},
        "block_delta_counts_optimize_minus_assemble": {
            "reductions": sum(1 for r in rows if r["assemble"]["status"] == r["optimize"]["status"] == "built"
                            and r["optimize"]["model"]["block_count"] < r["assemble"]["model"]["block_count"]),
            "unchanged": sum(1 for r in rows if r["assemble"]["status"] == r["optimize"]["status"] == "built"
                        and r["optimize"]["model"]["block_count"] == r["assemble"]["model"]["block_count"]),
            "increases": sum(1 for r in rows if r["assemble"]["status"] == r["optimize"]["status"] == "built"
                            and r["optimize"]["model"]["block_count"] > r["assemble"]["model"]["block_count"]),
        },
        "total_block_delta_optimize_minus_assemble": sum(
            r["optimize"]["model"]["block_count"] - r["assemble"]["model"]["block_count"]
            for r in rows if r["assemble"]["status"] == r["optimize"]["status"] == "built"),
        "mode_execution_counts": {
            mode: {"reused_models": sum(bool(r[mode].get("replayed")) for r in rows),
                   "fresh_cli_builds": sum(not bool(r[mode].get("replayed")) and r[mode]["status"] == "built" for r in rows)}
            for mode in ("assemble", "optimize")},
        "validation_history": [
            "An initial checker draft normalized homogeneous nets before rotating them and applied the inverse XY map. That draft was rejected during review and is not a scored result.",
            "The final checker first reorders raw homogeneous controls by the independently enumerated positive XY map, then normalizes and compares the complete eight-control span against its Q1 cell.",
        ],
        "cases": rows,
        "limitations": [
            "Quality values use the independent 5x5x5 spectral and thickness diagnostic and are not certificates.",
            "The corner-count DOF estimate deduplicates exact geometric coordinates; coincident but topologically distinct vertices require separate identity data.",
            "The extraction oracle compares all eight homogeneous controls modulo one common scale and the four positive XY rotations only.",
            "When --reuse-existing is used, cached model elapsed times are zero and timing comparisons must be omitted.",
            "Conformity beyond serialized shared block identities requires interpretation of the patch section indexing contract.",
        ],
    }
    failed = []
    if len(rows) != args.expect_cases:
        failed.append(f"expected {args.expect_cases} corpus cases, found {len(rows)}")
    for row in rows:
        for mode in ("assemble", "optimize"):
            result = row[mode]
            if result["status"] != "built":
                failed.append(f"{row['id']} {mode}: {result.get('diagnostic', result['status'])}")
                continue
            model = result["model"]
            if model["nurbs_patch_errors"]:
                failed.append(f"{row['id']} {mode}: {model['nurbs_patch_errors'][0]}")
            if mode == "optimize" and model["nurbs_patch_count"] == 0:
                failed.append(f"{row['id']} optimize: no NURBS patches emitted")
            extraction = model["exact_extraction"]
            if extraction["extraction_failures"]:
                failed.append(f"{row['id']} {mode}: {extraction['extraction_failures'][0]}")
            conformity = model["exact_conformity"]
            if conformity["interface_control_mismatch_count"]:
                failed.append(f"{row['id']} {mode}: {conformity['interface_control_mismatch_count']} exact interface mismatches")
            if conformity["interfaces_checked_exactly"] != model["interface_count"]:
                failed.append(f"{row['id']} {mode}: checked {conformity['interfaces_checked_exactly']} of {model['interface_count']} interfaces")
        if row["assemble"]["status"] == row["optimize"]["status"] == "built":
            base_volume = row["assemble"]["model"]["exact_q1_volume_mm3"]
            opt_volume = row["optimize"]["model"]["exact_q1_volume_mm3"]
            if base_volume != opt_volume:
                failed.append(f"{row['id']}: exact Q1 volumes differ: {base_volume} vs {opt_volume}")
    summary["verdict"] = "PASS" if not failed else "FAIL"
    summary["failures"] = failed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"{summary['verdict']}: {len(rows)} cases; assemble={dict(summary['mode_counts']['assemble'])}; optimize={dict(summary['mode_counts']['optimize'])}; report={args.output}")
    if failed:
        for failure in failed[:20]:
            print(f"FAIL: {failure}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
