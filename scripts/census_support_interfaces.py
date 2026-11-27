#!/usr/bin/env python3
"""Rebuild and enumerate exact block interfaces in frozen polygon corpora.

The recipe identity is the SHA-256 of canonical JSON. Corpus memberships are
kept on each identity, so deduplication never erases which evaluation supplied
an input. Interface checks compare rational homogeneous face nets after exact
degree elevation and the stored orientation. This certifies a matching face
trace only; it does not infer a common parent support or higher continuity.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
POLY_CASES = ROOT / "docs/results/sketch-planner-polygon-evaluation/cases.json"
POLY_SUPPLEMENTS = ROOT / "docs/results/polygon-extrusion-quality/ring-lookahead/supplemental-cases.json"
HELDOUT_MANIFEST = ROOT / "docs/results/planar-external-review/inputs.sha256.json"
HELDOUT_DIR = ROOT / "docs/results/planar-external-review"


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_polygon_recipes() -> list[dict]:
    rows: list[dict] = []
    historical = json.loads(POLY_CASES.read_text())
    for case in historical["cases"]:
        recipe = case["recipe"]
        rows.append({
            "corpus": "historical_45", "case_id": case["id"],
            "expected_admissible": bool(case.get("admissible", True)),
            "recipe": recipe, "raw_sha256": None,
        })
    supplements = json.loads(POLY_SUPPLEMENTS.read_text())
    for case in supplements["cases"]:
        recipe = case["recipe"]
        rows.append({
            "corpus": "supplemental_12", "case_id": case["id"],
            "expected_admissible": True, "recipe": recipe, "raw_sha256": None,
        })
    manifest = json.loads(HELDOUT_MANIFEST.read_text())
    for item in manifest["inputs"]:
        path = HELDOUT_DIR / item["path"]
        data = path.read_bytes()
        if sha256(data) != item["sha256"]:
            raise ValueError(f"held-out source hash mismatch: {item['path']}")
        rows.append({
            "corpus": "heldout_11", "case_id": Path(item["path"]).name.removesuffix(".recipe.json"),
            "expected_admissible": True, "recipe": json.loads(data),
            "raw_sha256": item["sha256"],
        })
    return rows


def grouped_recipes(rows: list[dict]) -> list[dict]:
    grouped: dict[str, dict] = {}
    for row in rows:
        payload = canonical_bytes(row["recipe"])
        digest = sha256(payload)
        entry = grouped.setdefault(digest, {
            "recipe_sha256": digest, "recipe": row["recipe"],
            "memberships": [], "expected_admissible": row["expected_admissible"],
        })
        if entry["recipe"] != row["recipe"]:
            raise AssertionError("canonical recipe hash collision")
        entry["memberships"].append({
            "corpus": row["corpus"], "case_id": row["case_id"],
            "raw_sha256": row["raw_sha256"],
        })
        if entry["expected_admissible"] != row["expected_admissible"]:
            raise ValueError(f"conflicting expected outcomes for recipe {digest}")
    return [grouped[key] for key in sorted(grouped)]


def _face_axes(side: str) -> tuple[int, tuple[int, int], int]:
    if len(side) != 2 or side[0] not in "uvw" or side[1] not in "01":
        raise ValueError(f"invalid block side {side!r}")
    fixed = "uvw".index(side[0])
    return fixed, tuple(axis for axis in range(3) if axis != fixed), int(side[1])


def _transform_control(control: list[int], transform: dict) -> tuple[Fraction, ...]:
    matrix = transform.get("matrix", [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    translation = transform.get("translation", [0, 0, 0])
    denom = F(transform.get("denominator", 1))
    xyz = [sum(F(matrix[a][j]) * control[j] for j in range(3))
           + denom * translation[a] * control[3] for a in range(3)]
    return tuple(xyz + [denom * control[3]])


def face_control_net(block: dict, side: str, transform: dict | None = None):
    degrees = block["degrees"]
    dims = [int(p) + 1 for p in degrees]
    fixed, free, endpoint = _face_axes(side)
    face_index = endpoint * (dims[fixed] - 1)
    rows = []
    for a in range(dims[free[0]]):
        row = []
        for b in range(dims[free[1]]):
            ix = [0, 0, 0]
            ix[fixed] = face_index
            ix[free[0]], ix[free[1]] = a, b
            flat = (ix[0] * dims[1] + ix[1]) * dims[2] + ix[2]
            row.append(_transform_control(block["controls"][flat], transform or {}))
        rows.append(row)
    return tuple(degrees[i] for i in free), rows


def orient_second(rows: list[list[tuple[Fraction, ...]]], orient: dict,
                  shape: tuple[int, int]) -> list[list[tuple[Fraction, ...]]]:
    out = []
    for a in range(shape[0]):
        row = []
        for b in range(shape[1]):
            x = shape[0] - 1 - a if orient.get("flip_first", False) else a
            y = shape[1] - 1 - b if orient.get("flip_second", False) else b
            sa, sb = (y, x) if orient.get("swap", False) else (x, y)
            row.append(rows[sa][sb])
        out.append(row)
    return out


def elevate_curve(values: list[tuple[Fraction, ...]]) -> list[tuple[Fraction, ...]]:
    degree = len(values) - 1
    elevated = [values[0]]
    for i in range(1, degree + 1):
        alpha = F(i, degree + 1)
        elevated.append(tuple(alpha * left + (1 - alpha) * right
                              for left, right in zip(values[i - 1], values[i])))
    elevated.append(values[-1])
    return elevated


def elevate_face(rows: list[list[tuple[Fraction, ...]]], target: tuple[int, int]):
    result = rows
    while len(result) - 1 < target[0]:
        columns = list(zip(*result))
        result = [list(row) for row in zip(*[elevate_curve(list(c)) for c in columns])]
    while len(result[0]) - 1 < target[1]:
        result = [elevate_curve(row) for row in result]
    return result


def exact_homogeneous_face_match(first, second, orientation: dict) -> tuple[bool, dict]:
    first_degree, a = first
    second_degree, b = second
    a = [[tuple(F(x) for x in point) for point in row] for row in a]
    b = [[tuple(F(x) for x in point) for point in row] for row in b]
    if orientation.get("swap", False):
        second_degree = tuple(reversed(second_degree))
    target = tuple(max(x, y) for x, y in zip(first_degree, second_degree))
    a = elevate_face(a, target)
    second_target = tuple(reversed(target)) if orientation.get("swap", False) else target
    b = elevate_face(b, second_target)
    b = orient_second(b, orientation, (target[0] + 1, target[1] + 1))
    flat_a, flat_b = [x for row in a for x in row], [x for row in b for x in row]
    ref_a = next((v[3] for v in flat_a if v[3]), None)
    ref_b = next((v[3] for v in flat_b if v[3]), None)
    if ref_a is None or ref_b is None:
        return False, {"reason": "face has no nonzero homogeneous weight"}
    scale = ref_a / ref_b
    if scale <= 0:
        return False, {"reason": "nonpositive homogeneous gauge ratio"}
    matches = all(x == tuple(scale * y for y in z) for x, z in zip(flat_a, flat_b))
    return matches, {
        "degree_after_elevation": list(target),
        "positive_gauge_second_to_first": [scale.numerator, scale.denominator],
    }


def census_model_interfaces(model: dict) -> list[dict]:
    topology = model.get("topology", {})
    interfaces = topology.get("interfaces", [])
    blocks = {b["id"]: b for b in model.get("blocks", [])}
    transform = model.get("transform", {})
    rows = []
    for index, edge in enumerate(interfaces):
        a_id, a_side = edge["first"]
        b_id, b_side = edge["second"]
        exact, method = exact_homogeneous_face_match(
            face_control_net(blocks[a_id], a_side, transform),
            face_control_net(blocks[b_id], b_side, transform), edge.get("orientation", {}))
        rows.append({
            "interface_index": index,
            "first": {"block": a_id, "side": a_side},
            "second": {"block": b_id, "side": b_side},
            "orientation": edge.get("orientation", {}),
            "relation": "exact_full_face_homogeneous_trace" if exact else "trace_mismatch",
            "trace_verified": exact,
            "verification_method": "exact_fraction_homogeneous_control_net_after_orientation_and_degree_elevation",
            "details": method,
            "support_continuity_claim": "not_inferred_from_trace_equality",
        })
    return rows


def compile_inherited_spaces(cli: Path, work: Path, model_path: Path, digest: str,
                            timeout: int) -> dict:
    """Compile the requested pz=2 baseline and inherited-C1 packages natively."""
    plan_dir = work / "native-plans"
    output_dir = work / "native-fields"
    plan_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {}
    for label, continuity in (("pz2_c0", 0), ("pz2_inherited_c1", 1)):
        plan_path = plan_dir / f"{label}.plan.json"
        plan_path.write_bytes(canonical_bytes({
            "schema": "inherited_field_plan_v1", "axial_degree": 2,
            "axial_continuity": continuity, "restore_strip": False,
            "transverse_interfaces": [],
        }) + b"\n")
        package_path = output_dir / f"{digest}.{label}.field.json"
        command = [str(cli), "export-inherited-space", str(model_path),
                   str(plan_path), str(package_path)]
        try:
            run = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            result[label] = {"status": "timeout", "timeout_seconds": timeout,
                             "stdout": exc.stdout, "stderr": exc.stderr}
            continue
        if run.returncode != 0:
            try:
                observed = json.loads(run.stdout)
            except (json.JSONDecodeError, TypeError):
                observed = {}
            result[label] = {"status": "refused", "returncode": run.returncode,
                             "stdout": run.stdout, "stderr": run.stderr,
                             "reason": observed.get("reason", run.stderr.strip() or "unspecified refusal"),
                             "command": command}
            continue
        pack = json.loads(package_path.read_text())
        field = pack["field"]
        result[label] = {
            "status": "exported", "package_path": str(package_path.relative_to(ROOT)),
            "command": command, "cli_observation": json.loads(run.stdout),
            "scalar_dofs": field["scalar_dofs"],
            "c0_scalar_dofs": field["c0_scalar_dofs"],
            "inherited_join_count": field["inherited_join_count"],
            "retained_c0_join_count": field["retained_c0_join_count"],
            "transverse_c1_join_count": field.get("transverse_c1_join_count", 0),
        }
    return result


def build_one(cli: Path, work: Path, entry: dict, timeout: int, reuse_existing: bool = False) -> dict:
    digest = entry["recipe_sha256"]
    recipe_path = work / "recipes" / f"{digest}.recipe.json"
    model_path = work / "models" / f"{digest}.model.json"
    recipe_path.parent.mkdir(parents=True, exist_ok=True)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    recipe_path.write_bytes(canonical_bytes(entry["recipe"]) + b"\n")
    build_command = [str(cli), "build", str(recipe_path), str(model_path)]
    try:
        reused = bool(reuse_existing and model_path.is_file())
        built = None if reused else subprocess.run(build_command, text=True, capture_output=True, timeout=timeout)
        build_output = (json.loads(built.stdout) if built.stdout.strip().startswith("{") else None) if built else None
        if built and built.returncode != 0:
            return {"status": "refused" if entry["expected_admissible"] else "expected_refusal",
                    "build_returncode": built.returncode, "build_stdout": built.stdout,
                    "build_stderr": built.stderr, "command": build_command,
                    "interfaces": []}
        if reused and not model_path.is_file():
            raise ValueError(f"cached model missing: {model_path}")
        audited = subprocess.run([str(cli), "audit", str(model_path)], text=True,
                                 capture_output=True, timeout=timeout)
        audit_output = json.loads(audited.stdout) if audited.stdout.strip().startswith("{") else None
        if audited.returncode != 0:
            return {"status": "audit_failed", "build_returncode": 0,
                    "audit_returncode": audited.returncode, "audit_stdout": audited.stdout,
                    "audit_stderr": audited.stderr, "command": build_command,
                    "interfaces": []}
        model = json.loads(model_path.read_text())
        iface_rows = census_model_interfaces(model)
        native = compile_inherited_spaces(cli, work, model_path, digest, timeout)
        return {
            "status": "admitted" if entry["expected_admissible"] else "unexpected_admission",
            "build_returncode": 0, "audit_returncode": audited.returncode,
            "build_reused_from_completed_run": reused,
            "command": build_command, "build_observation": build_output,
            "audit_observation": audit_output,
            "model_path": str(model_path.relative_to(ROOT)) if model_path.is_relative_to(ROOT) else str(model_path),
            "model_sha256": sha256(model_path.read_bytes()),
            "block_count": len(model.get("blocks", [])),
            "declared_interface_count": len(model.get("topology", {}).get("interfaces", [])),
            "exact_trace_count": sum(r["trace_verified"] for r in iface_rows),
            "native_inherited_space": native,
            "interfaces": iface_rows,
        }
    except subprocess.TimeoutExpired as exc:
        return {"status": "timeout", "command": build_command,
                "timeout_seconds": timeout, "stdout": exc.stdout, "stderr": exc.stderr,
                "interfaces": []}


def polygon_census(cli: Path, work: Path, timeout: int, reuse_existing: bool = False) -> dict:
    source_rows = load_polygon_recipes()
    entries = grouped_recipes(source_rows)
    for index, entry in enumerate(entries, start=1):
        result = build_one(cli, work, entry, timeout, reuse_existing)
        entry.update(result)
        entry.pop("recipe", None)
        native = entry.get("native_inherited_space", {})
        print(f"[{index}/{len(entries)}] {entry['recipe_sha256'][:12]} {entry['status']} "
              f"blocks={entry.get('block_count', 0)} interfaces={entry.get('declared_interface_count', 0)} "
              f"native={native.get('pz2_inherited_c1', {}).get('status', 'not-run')}", flush=True)
    return {
        "source_counts": dict(Counter(r["corpus"] for r in source_rows)),
        "source_entries": len(source_rows),
        "unique_canonical_recipe_hashes": len(entries),
        "duplicate_hash_groups": [
            {"recipe_sha256": e["recipe_sha256"], "memberships": e["memberships"]}
            for e in entries if len(e["memberships"]) > 1
        ],
        "entries": entries,
    }


def _state_record(state_id: str, path: Path, status: str, note: str = "") -> dict:
    model = json.loads(path.read_text())
    rows = census_model_interfaces(model)
    return {
        "state_id": state_id, "source": str(path.relative_to(ROOT)), "status": status,
        "model_sha256": sha256(path.read_bytes()),
        "note": note, "block_count": len(model.get("blocks", [])),
        "declared_interface_count": len(model.get("topology", {}).get("interfaces", [])),
        "exact_trace_count": sum(r["trace_verified"] for r in rows),
        "interfaces": rows,
        "other_topology_records": {k: len(v) for k, v in model.get("topology", {}).items()
                                   if k != "interfaces" and isinstance(v, list)},
    }


def named_states_census() -> dict:
    states = []
    tall_dir = ROOT / "docs/site/box-cylinder/tall"
    for name in ("box", "cylinder", "union", "intersection", "a-minus-b", "b-minus-a", "xor"):
        states.append(_state_record(f"tall_cylinder_{name}", tall_dir / f"{name}.model.json",
                                    "qualified_reference", "Tall-cylinder Boolean state."))
    box_dir = ROOT / "docs/site/box-cylinder"
    for name in ("classic", "result"):
        states.append(_state_record(f"box_cylinder_{name}", box_dir / f"{name}.model.json",
                                    "qualified_reference", "Centered through-cylinder box/cylinder model."))
    fieldcraft_dir = ROOT / "docs/site/data"
    for path in sorted(fieldcraft_dir.glob("fieldcraft_*/fieldcraft_*.model.json")):
        states.append(_state_record(path.parent.name, path, "qualified_compact_fieldcraft",
                                    "Compact FieldCraft model; exact trace does not imply common parent support."))
    s2 = s2_census()
    return {"states": states, "s2_crossed_union": s2}


def _raw_face_controls(block: dict, axis: int, side: int):
    dims = [int(p) + 1 for p in block["degrees"]]
    free = tuple(i for i in range(3) if i != axis)
    rows = []
    for a in range(dims[free[0]]):
        row = []
        for b in range(dims[free[1]]):
            ix = [0, 0, 0]
            ix[axis] = side * (dims[axis] - 1)
            ix[free[0]], ix[free[1]] = a, b
            flat = (ix[0] * dims[1] + ix[1]) * dims[2] + ix[2]
            row.append(tuple(F(x) for x in block["controls"][flat]))
        rows.append(row)
    return tuple(int(block["degrees"][i]) for i in free), rows


def s2_census() -> dict:
    """Enumerate seven complete band joins and four wing/core partial joins."""
    result = []
    base = ROOT / "docs/site/crossed-tools"
    full_pairs = []
    for col in (0, 1, 2):
        bands = [2, 3] if col != 1 else list(range(6))
        for band in bands[:-1]:
            full_pairs.append((f"band{band}_col{col}_z" + ("4_8" if col != 1 else "0_12"), "u1",
                               f"band{band + 1}_col{col}_z" + ("4_8" if col != 1 else "0_12"), "u0"))
    partial_pairs = []
    for band in (2, 3):
        for col, wing_side, core_side in ((0, "v0", "v1"), (2, "v1", "v0")):
            wing_z = "z4_8"
            partial_pairs.append((f"band{band}_col{col}_{wing_z}", wing_side,
                                  f"band{band}_col1_z0_12", core_side))
    for rise, model_name in ((-3, "union--3.json"), (0, "union-0.json"), (3, "union-3.json")):
        model = json.loads((base / model_name).read_text())
        blocks = {b["id"]: b for b in model["blocks"]}
        rows = []
        for i, (a, sa, b, sb) in enumerate(full_pairs):
            a_block, b_block = blocks[a], blocks[b]
            # u is axis 0; remaining local axes are v,w.
            a_net, b_net = _raw_face_controls(a_block, 0, 1), _raw_face_controls(b_block, 0, 0)
            exact, detail = exact_homogeneous_face_match(a_net, b_net, {})
            rows.append({"interface_index": i, "first": {"block": a, "side": sa},
                         "second": {"block": b, "side": sb},
                         "relation": "exact_full_face_homogeneous_trace" if exact else "trace_mismatch",
                         "trace_verified": exact, "evidence": detail,
                         "support_continuity_claim": "not_inferred_from_trace_equality"})
        for j, (a, sa, b, sb) in enumerate(partial_pairs, start=len(rows)):
            rows.append({"interface_index": j, "first": {"block": a, "side": sa},
                         "second": {"block": b, "side": sb},
                         "relation": "partial_wing_to_core_trace",
                         "trace_verified": None,
                         "verification_method": "partial_physical_footprint_with_preimage_aware_trace_constraints",
                         "evidence": {
                             "full_face_equality": False,
                             "source": "docs/results/crossed-interface-iga/README.md",
                             "trace_method": "exact rational constraint null-space construction; polynomial degree bound and seven distinct ordinate nodes; held-out probes are numerical diagnostics",
                             "documented_held_out_trace_defect_upper": 8e-16,
                             "documented_geometry_mismatch_mm": 5.7e-14,
                         },
                         "support_continuity_claim": "partial_coupling_only_not_a_shared_full_face_support"})
        result.append({
            "state_id": f"S2_crossed_union_rise_{rise:+d}mm", "source": f"docs/site/crossed-tools/{model_name}",
            "status": "qualified_geometry_with_partial_interfaces",
            "block_count": len(model["blocks"]), "interface_count": len(rows),
            "full_exact_trace_count": sum(r["trace_verified"] is True for r in rows),
            "partial_trace_count": sum(r["relation"] == "partial_wing_to_core_trace" for r in rows),
            "non_face_collapsed_polar_line_constraints": 2,
            "interfaces": rows,
        })
    return {"states": result, "scope_note": "Seven full face joins and four partial wing/core faces per rise; the two collapsed polar-line constraints are separately counted and are not faces."}


def make_report(cli: Path, work: Path, timeout: int, reuse_existing: bool = False) -> dict:
    polygon = polygon_census(cli, work, timeout, reuse_existing)
    named = named_states_census()
    all_statuses = Counter(e["status"] for e in polygon["entries"])
    total_interfaces = sum(e.get("declared_interface_count", 0) for e in polygon["entries"])
    verified = sum(e.get("exact_trace_count", 0) for e in polygon["entries"])
    failed_rows = [(e["recipe_sha256"], row) for e in polygon["entries"] for row in e.get("interfaces", [])
                   if row["trace_verified"] is not True]
    by_corpus = {}
    for corpus, source_n in polygon["source_counts"].items():
        selected = [e for e in polygon["entries"]
                    if any(m["corpus"] == corpus for m in e["memberships"])]
        native = {name: Counter() for name in ("pz2_c0", "pz2_inherited_c1")}
        refusal_reasons = {name: Counter() for name in native}
        common_pairs = []
        for entry in selected:
            for name, record in entry.get("native_inherited_space", {}).items():
                native[name]["exported_unique_models"] += record.get("status") == "exported"
                native[name]["refused_unique_models"] += record.get("status") == "refused"
                native[name]["timeout_unique_models"] += record.get("status") == "timeout"
                if record.get("status") == "exported":
                    for key in ("scalar_dofs", "inherited_join_count", "retained_c0_join_count"):
                        native[name][f"sum_{key}_over_unique_models"] += record[key]
                elif record.get("status") == "refused":
                    refusal_reasons[name][record.get("reason", "unspecified refusal")] += 1
            c0 = entry.get("native_inherited_space", {}).get("pz2_c0", {})
            c1 = entry.get("native_inherited_space", {}).get("pz2_inherited_c1", {})
            if c0.get("status") == "exported" and c1.get("status") == "exported":
                common_pairs.append((c0, c1))
        by_corpus[corpus] = {
            "source_recipe_memberships": source_n,
            "unique_recipe_hashes_in_corpus": len(selected),
            "interface_observations_once_per_hash_in_corpus": sum(
                e.get("declared_interface_count", 0) for e in selected),
            "duplicate_memberships": [m for e in selected if len(e["memberships"]) > 1
                                       for m in e["memberships"] if m["corpus"] == corpus],
            "native_spaces": {k: dict(v) for k, v in native.items()},
            "native_refusal_reasons": {k: dict(v) for k, v in refusal_reasons.items()},
            "native_common_compilation_pairs": {
                "unique_model_count": len(common_pairs),
                "c0_scalar_dofs_sum": sum(a["scalar_dofs"] for a, _ in common_pairs),
                "inherited_c1_scalar_dofs_sum": sum(b["scalar_dofs"] for _, b in common_pairs),
                "scalar_dof_delta_c1_minus_c0": sum(b["scalar_dofs"] - a["scalar_dofs"]
                                                      for a, b in common_pairs),
                "inherited_join_count_sum": sum(b["inherited_join_count"] for _, b in common_pairs),
            },
        }
    named_states = named["states"]
    s2_states = named["s2_crossed_union"]["states"]
    named_summary = {
        "named_model_state_observations": len(named_states),
        "named_model_interface_observations": sum(s["declared_interface_count"] for s in named_states),
        "named_model_exact_trace_observations": sum(s["exact_trace_count"] for s in named_states),
        "s2_state_observations": len(s2_states),
        "s2_interface_observations": sum(s["interface_count"] for s in s2_states),
        "s2_full_face_exact_trace_observations": sum(s["full_exact_trace_count"] for s in s2_states),
        "s2_partial_interface_observations": sum(s["partial_trace_count"] for s in s2_states),
        "s2_collapsed_polar_line_constraints_not_counted_as_faces": sum(
            s["non_face_collapsed_polar_line_constraints"] for s in s2_states),
        "aggregation_note": "These are state observations, not unique interface geometries; same template across parameter states remains separately attributable.",
    }
    manifests = [POLY_CASES, POLY_SUPPLEMENTS, HELDOUT_MANIFEST]
    try:
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
                                  capture_output=True, check=True).stdout.strip()
    except (subprocess.SubprocessError, OSError):
        revision = "unavailable"
    return {
        "schema": "support_interface_census_v1",
        "run_provenance": {
            "source_revision": revision,
            "python_version": sys.version.split()[0],
            "cli_path": str(cli), "cli_sha256": sha256(cli.read_bytes()),
            "model_execution": "built_and_audited" if not reuse_existing else "audited_prebuilt_models_from_completed_build_run",
            "source_manifests": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p.read_bytes())}
                                  for p in manifests],
        },
        "method": {
            "recipe_deduplication": "SHA-256 over canonical JSON; each source corpus membership retained",
            "face_trace": "exact Fraction homogeneous control-net equality after orientation and exact degree elevation",
            "interpretation_limit": "A shared exact face trace does not prove a common parent spline support, G1/G2 continuity, or an analysis-space coupling.",
            "aggregation": "recipe interfaces counted once per unique canonical recipe; named state observations are reported separately and are not pooled as unique geometries",
            "native_candidate": "export-inherited-space with axial_degree=2, axial_continuity=1, restore_strip=false, transverse_interfaces=[]; baseline is the same plan with axial_continuity=0",
        },
        "polygon_recipes": polygon,
        "named_state_summary": named_summary,
        "polygon_summary": {
            "source_recipe_entries": polygon["source_entries"],
            "unique_recipe_hashes": polygon["unique_canonical_recipe_hashes"],
            "build_and_audit_statuses": dict(all_statuses),
            "interfaces_declared_on_successful_unique_recipes": total_interfaces,
            "interfaces_with_exact_homogeneous_trace": verified,
            "trace_mismatch_or_unavailable": len(failed_rows),
            "unresolved_interfaces": [{"recipe_sha256": h, "interface": row} for h, row in failed_rows],
            "no_silent_refusal": True,
            "by_corpus": by_corpus,
            "native_candidate_vs_c0": {
                "deduplication_unit": "unique canonical recipe hash",
                "membership_rule": "each corpus membership is reported separately in by_corpus, even when a model hash is shared",
                "c0_baseline": "pz2_c0",
                "candidate": "pz2_inherited_c1",
            },
        },
        "named_states": named,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--work", type=Path, default=ROOT / "build/support-census")
    parser.add_argument("--output", type=Path, default=ROOT / "build/support-census/census.json")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--reuse-existing", action="store_true",
                        help="audit already built models under --work; useful after a reporting-only fix")
    args = parser.parse_args()
    report = make_report(args.cli.resolve(), args.work.resolve(), args.timeout, args.reuse_existing)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    summary = report["polygon_summary"]
    print(json.dumps({"output": str(args.output), "polygon_summary": summary,
                      "named_state_count": len(report["named_states"]["states"]),
                      "s2_state_count": len(report["named_states"]["s2_crossed_union"]["states"])}))
    bad = summary["trace_mismatch_or_unavailable"] > 0 or any(
        entry["status"] != ("admitted" if entry["expected_admissible"] else "expected_refusal")
        for entry in report["polygon_recipes"]["entries"]
    )
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
