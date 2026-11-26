#!/usr/bin/env python3
"""Compare the latest planner binary with the frozen convex host-ring corpus.

The 45 input recipes, reference scores, exact section oracle, and diagnostic
quality metric are taken from the checked-in polygon evaluation. This script
does not regenerate the corpus or need the prototype archive at runtime.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import random
import os
import platform
from pathlib import Path
import statistics
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "docs/results/sketch-planner-polygon-evaluation"
CORPUS = EVAL_DIR / "cases.json"
OLD_CANDIDATE = EVAL_DIR / "candidate.json"
EXPECTED_CORPUS_SHA256 = "f10964bb0417a702a394f66aa96e148b6ead67685a3a9f9f541cea8d6ab98e7b"
EXPECTED_REFERENCE_CLI_SHA256 = "684c7c6a89423f798648800e28c53d1065c23b8aaa7e5407409714cea4a36f70"
DEFAULT_CLI = Path("/tmp/gallery-minimal-49f6/target/release/solidcraft-trivariate")
DEFAULT_OUTPUT = ROOT / "docs/results/convex-host-rings.json"
DEFAULT_WORK = Path("/tmp/convex-host-rings-latest")

# Prototype observations copied from the reviewed 30-case extract. This table
# is intentionally self-contained: running the evaluation does not open the
# original archive. The archived Python quality values are diagnostic samples,
# not continuous certificates.
PROTO_OBSERVATIONS: dict[str, tuple[int, float, float]] = {
    "seeded_convex_00": (12, 21.099906424888086, 146.89872892919297),
    "seeded_convex_01": (48, 22.354411677204656, 486.3244985738964),
    "seeded_convex_02": (19, 30.21501035369444, 172.72703439966415),
    "seeded_convex_03": (14, 11.434452901602901, 405.959681771679),
    "seeded_convex_04": (12, 31.603415321521602, 215.53333721628002),
    "seeded_convex_05": (57, 14.689085634726743, 496.57614035433653),
    "seeded_convex_06": (9, 16.58278318260224, 93.80511145085347),
    "seeded_convex_07": (28, 45.0323103798767, 596.1511924634355),
    "seeded_convex_08": (33, 44.75674060283433, 116.86064155592015),
    "seeded_convex_09": (14, 27.171961725072705, 892.90247224466),
    "seeded_convex_10": (24, 19.087580767097243, 133.21185972527363),
    "seeded_convex_11": (72, 36.7882525295242, 271.55565535439627),
    "seeded_convex_12": (5, 4.351017061804353, 118.14371963716529),
    "seeded_convex_13": (48, 15.656016924455162, 474.0684949717208),
    "seeded_convex_14": (30, 20.58121473368661, 99.89083507473516),
    "seeded_convex_15": (29, 19.33602487746356, 267.91136845056167),
    "seeded_convex_16": (22, 43.252118657324445, 153.38680522976173),
    "seeded_convex_17": (42, 29.570105683999085, 585.955200329631),
    "seeded_convex_18": (12, 21.099906424888086, 77.77256050241715),
    "seeded_convex_19": (48, 40.68449482150903, 432.98640441588327),
    "seeded_convex_20": (27, 16.85366420108068, 160.5092700969279),
    "seeded_convex_21": (24, 39.50912512278567, 669.3917352376196),
    "seeded_convex_22": (22, 45.689614419051, 210.89525812323842),
    "seeded_convex_23": (74, 36.30508212376244, 378.998303956225),
    "seeded_convex_24": (9, 36.26596724622451, 161.43495357671443),
    "seeded_convex_25": (28, 36.962367847544755, 450.1099261322755),
    "seeded_convex_26": (33, 44.75674060283433, 132.33761468243796),
    "seeded_convex_27": (29, 46.211610416550506, 396.3663097700276),
    "seeded_convex_28": (16, 14.592631527330207, 79.59178304683428),
    "seeded_convex_29": (74, 121.26557860914045, 379.18272347006706),
}
PROTO_ARCHIVE_SHA256 = "d5238579f7478caec0e03e3e2b86b4520354ecffa4ee8116f3525d80fd27d409"
PROTO_REVIEW_SNAPSHOT_SHA256 = "74725bbfa08004cdb0623f61c7b2f957ecf4c1cb8b35bd3b9474180f45f25109"
FRESH_PROBE_SOURCE_SHA256 = "9ccca9f5a55dfd59ea8f9710c1fce5b32067f1d081e7be7c03c6c22827b0cce0"
FRESH_PROTO_BLOCK_COUNTS = (31, 44, 32, 49, 15, 51, 26, 42, 37, 47, 15, 53)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_runner():
    source = EVAL_DIR / "run.py"
    spec = importlib.util.spec_from_file_location("frozen_polygon_runner", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load frozen evaluation runner: {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def condition_max(row: dict) -> float | None:
    samples = row.get("condition_grid_5cube_per_block") or []
    values = [sample.get("max_condition_2") for sample in samples
              if sample.get("max_condition_2") is not None]
    return max(values) if values else None


def quality_snapshot(row: dict) -> dict:
    cert = row.get("closed_domain_certificate") or {}
    return {
        "sampled_jacobian_condition_2_max_5cube": condition_max(row),
        "continuous_determinant_lower_bound_mm3": cert.get("determinant_lower_bound_mm3"),
        "certificate_status": cert.get("status"),
        "exact_volume_match": row.get("volume_match"),
        "geometry_pass": row.get("geometry_pass", False),
        "outcome": row.get("outcome"),
    }


def pair_case(old: dict, new: dict) -> dict:
    old_pass = bool(old.get("geometry_pass"))
    new_pass = bool(new.get("geometry_pass"))
    old_blocks = old.get("block_count")
    new_blocks = new.get("block_count")
    proto = PROTO_OBSERVATIONS.get(old["id"])
    proto_blocks = proto[0] if proto else None
    old_quality = quality_snapshot(old)
    new_quality = quality_snapshot(new)
    return {
        "id": old["id"],
        "family": old["family"],
        "old_candidate": {
            "outcome": old.get("outcome"),
            "admitted": old.get("admitted", False),
            "blocks": old_blocks,
            "quality": old_quality,
        },
        "latest_release": {
            "outcome": new.get("outcome"),
            "admitted": new.get("admitted", False),
            "blocks": new_blocks,
            "quality": new_quality,
        },
        "prototype": {
            "blocks": proto_blocks,
            "sampled_jacobian_condition_2_max_5cube": proto[2] if proto else None,
            "sampled_jacobian_condition_2_max_5cube_inplane_5x5": proto[1] if proto else None,
            "quality_observed": proto is not None,
            "quality_note": "Prototype samples are diagnostics, not continuous certificates." if proto else None,
        },
        "block_delta_latest_minus_old": (new_blocks - old_blocks
                                          if isinstance(new_blocks, int) and isinstance(old_blocks, int)
                                          else None),
        "block_delta_latest_minus_prototype": (new_blocks - proto_blocks
                                                if isinstance(new_blocks, int) and isinstance(proto_blocks, int)
                                                else None),
        "outcome_change": (
            "common_success" if old_pass and new_pass else
            "regression" if old_pass and not new_pass else
            "newly_valid" if not old_pass and new_pass else
            "common_failure"
        ),
        "quality_condition_ratio_latest_over_old": (
            new_quality["sampled_jacobian_condition_2_max_5cube"] /
            old_quality["sampled_jacobian_condition_2_max_5cube"]
            if new_quality["sampled_jacobian_condition_2_max_5cube"] is not None
            and old_quality["sampled_jacobian_condition_2_max_5cube"] not in (None, 0)
            else None
        ),
        "quality_condition_ratio_latest_over_prototype": (
            new_quality["sampled_jacobian_condition_2_max_5cube"] / proto[2]
            if proto and new_quality["sampled_jacobian_condition_2_max_5cube"] is not None
            and proto[2] != 0 else None
        ),
    }


def summarize(pairs: list[dict]) -> dict:
    from collections import Counter
    result = {}
    families = sorted({pair["family"] for pair in pairs})
    for family in ["all", *families]:
        rows = pairs if family == "all" else [p for p in pairs if p["family"] == family]
        latest = [p["latest_release"] for p in rows]
        old = [p["old_candidate"] for p in rows]
        block_deltas = [p["block_delta_latest_minus_old"] for p in rows
                        if p["block_delta_latest_minus_old"] is not None]
        proto_deltas = [p["block_delta_latest_minus_prototype"] for p in rows
                        if p["block_delta_latest_minus_prototype"] is not None]
        q_new = [p["quality_condition_ratio_latest_over_old"] for p in rows
                 if p["quality_condition_ratio_latest_over_old"] is not None]
        prototype_ratios = [p["quality_condition_ratio_latest_over_prototype"] for p in rows
                            if p["quality_condition_ratio_latest_over_prototype"] is not None]
        result[family] = {
            "denominator": len(rows),
            "old_candidate_valid": sum(x["quality"]["geometry_pass"] for x in old),
            "latest_release_valid": sum(x["quality"]["geometry_pass"] for x in latest),
            "old_candidate_failures": dict(Counter(x["quality"]["outcome"] for x in old
                                                    if not x["quality"]["geometry_pass"])),
            "latest_release_failures": dict(Counter(x["quality"]["outcome"] for x in latest
                                                     if not x["quality"]["geometry_pass"])),
            "cases_with_old_and_new_block_counts": len(block_deltas),
            "old_blocks_total": sum(x["blocks"] for x in old if isinstance(x["blocks"], int)),
            "latest_blocks_total": sum(x["blocks"] for x in latest if isinstance(x["blocks"], int)),
            "block_delta_latest_minus_old_total": sum(block_deltas),
            "block_delta_latest_minus_old_median": statistics.median(block_deltas)
                if block_deltas else None,
            "cases_with_proto_counts": sum(p["prototype"]["blocks"] is not None for p in rows),
            "block_delta_latest_minus_proto_total": sum(proto_deltas),
            "block_delta_latest_minus_proto_median": statistics.median(proto_deltas)
                if proto_deltas else None,
            "sampled_condition_ratio_median_latest_over_old": statistics.median(q_new)
                if q_new else None,
            "sampled_condition_ratio_median_latest_over_prototype": statistics.median(prototype_ratios)
                if prototype_ratios else None,
            "latest_release_condition_summary": latest_condition_summary(rows),
        }
    return result


def seeded_breakdown(pairs: list[dict]) -> dict:
    rows = [pair for pair in pairs if pair["id"].startswith("seeded_convex_")]
    return {
        "case_count": len(rows),
        "old_candidate_blocks_total": sum(
            pair["old_candidate"]["blocks"] for pair in rows
            if isinstance(pair["old_candidate"]["blocks"], int)),
        "latest_release_blocks_total": sum(
            pair["latest_release"]["blocks"] for pair in rows
            if isinstance(pair["latest_release"]["blocks"], int)),
        "prototype_blocks_total": sum(
            pair["prototype"]["blocks"] for pair in rows
            if isinstance(pair["prototype"]["blocks"], int)),
        "old_candidate_valid": sum(pair["old_candidate"]["quality"]["geometry_pass"]
                                    for pair in rows),
        "latest_release_valid": sum(pair["latest_release"]["quality"]["geometry_pass"]
                                     for pair in rows),
        "per_case": [{
            "id": pair["id"],
            "old_candidate_blocks": pair["old_candidate"]["blocks"],
            "latest_release_blocks": pair["latest_release"]["blocks"],
            "prototype_blocks": pair["prototype"]["blocks"],
            "old_candidate_sampled_kappa3_max": pair["old_candidate"]["quality"][
                "sampled_jacobian_condition_2_max_5cube"],
            "latest_release_sampled_kappa3_max": pair["latest_release"]["quality"][
                "sampled_jacobian_condition_2_max_5cube"],
            "prototype_sampled_kappa3_max": pair["prototype"][
                "sampled_jacobian_condition_2_max_5cube"],
            "outcome_change": pair["outcome_change"],
        } for pair in rows],
    }


def latest_condition_summary(rows: list[dict]) -> dict:
    vals = [p["latest_release"]["quality"]["sampled_jacobian_condition_2_max_5cube"]
            for p in rows
            if p["latest_release"]["quality"]["geometry_pass"] and
            p["latest_release"]["quality"]["sampled_jacobian_condition_2_max_5cube"] is not None]
    vals.sort()
    if not vals:
        return {"n": 0, "median": None, "p95_nearest_rank": None, "worst": None}
    import math
    return {"n": len(vals), "median": vals[len(vals)//2],
            "p95_nearest_rank": vals[math.ceil(.95 * len(vals)) - 1], "worst": vals[-1]}


def cross(a: tuple[int, int], b: tuple[int, int], c: tuple[int, int]) -> int:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def convex_hull(points: list[tuple[int, int]]) -> list[list[int]]:
    """Exact monotone-chain hull, returned clockwise with a canonical start."""
    ordered = sorted(set(points))
    lower: list[tuple[int, int]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    ccw = lower[:-1] + upper[:-1]
    clockwise = list(reversed(ccw))
    # Rotate the reversed hull to the lexicographically smallest point. The
    # source only fixes a contour cycle; rotation does not alter its geometry.
    start = min(range(len(clockwise)), key=clockwise.__getitem__)
    return [[x, y] for x, y in clockwise[start:] + clockwise[:start]]


def fresh_cases() -> list[dict]:
    """Regenerate the 12 additional convex stock/cut probes, stdlib only."""
    rng = random.Random(560127)
    stock = [[0, 10], [20, 0], [124, 0], [150, 20],
             [150, 85], [129, 100], [14, 100], [0, 82]]
    cases = []
    for i in range(12):
        cuts = []
        for j in range(2 + i % 2):
            cx = 23 + 50 * j + rng.randint(-5, 5)
            cy = rng.randint(25, 75)
            points = [(cx + rng.randint(-12, 12), cy + rng.randint(-12, 12))
                      for _ in range(12)]
            vertices = convex_hull(points)
            z_interval = [4, 10] if (i + j) % 3 == 0 else [0, 10]
            cuts.append({"id": f"cut_{j + 1}", "vertices": vertices,
                         "z_interval_mm": z_interval})
        recipe = {"family": "sketch_layer_partition_v1", "height_mm": 10,
                  "outer": {"id": "stock", "vertices": stock},
                  "cuts": cuts, "construction": []}
        cases.append({"id": f"fresh_convex_{i:02d}",
                      "family": "fresh_convex_host_rings", "recipe": recipe})
    return cases


def signed_area2(vertices: list[list[int]]) -> int:
    return sum(a[0] * b[1] - a[1] * b[0]
               for a, b in zip(vertices, vertices[1:] + vertices[:1]))


def exact_original_geometry(recipe: dict) -> dict:
    """Independent shoelace area and layer integration using integer inputs."""
    from fractions import Fraction
    stock_area = Fraction(abs(signed_area2(recipe["outer"]["vertices"])), 2)
    cut_areas = [Fraction(abs(signed_area2(cut["vertices"])), 2)
                 for cut in recipe["cuts"]]
    levels = sorted({0, recipe["height_mm"],
                     *(int(z) for cut in recipe["cuts"] for z in cut["z_interval_mm"])})
    layers = []
    volume = Fraction(0)
    for z0, z1 in zip(levels, levels[1:]):
        zmid = Fraction(z0 + z1, 2)
        active = [(idx, cut) for idx, cut in enumerate(recipe["cuts"])
                  if cut["z_interval_mm"][0] < zmid < cut["z_interval_mm"][1]]
        area = stock_area - sum((cut_areas[idx] for idx, _ in active), Fraction(0))
        volume += area * (z1 - z0)
        layers.append({"z_interval_mm": [z0, z1],
                       "stock_area_mm2": [stock_area.numerator, stock_area.denominator],
                       "active_cut_ids": [cut["id"] for _, cut in active],
                       "active_cut_area_mm2": [
                           [cut_areas[idx].numerator, cut_areas[idx].denominator]
                           for idx, _ in active],
                       "expected_material_area_mm2": [area.numerator, area.denominator]})
    return {
        "stock_area_mm2": [stock_area.numerator, stock_area.denominator],
        "cut_areas_mm2": {cut["id"]: [cut_areas[idx].numerator, cut_areas[idx].denominator]
                           for idx, cut in enumerate(recipe["cuts"])},
        "layers": layers,
        "volume_exact_mm3": [volume.numerator, volume.denominator],
    }


def evaluate_fresh_cases(cli: Path, work: Path, timeout: float) -> dict:
    runner = load_runner()
    cases = fresh_cases()
    rows = []
    work.mkdir(parents=True, exist_ok=True)
    for index, item in enumerate(cases):
        row = runner.eval_one(cli, item, work, timeout)
        exact = exact_original_geometry(item["recipe"])
        exact_matches_runner = row.get("volume_exact_expected_mm3") == exact["volume_exact_mm3"]
        row.update({
            "recipe": item["recipe"],
            "prototype_block_count": FRESH_PROTO_BLOCK_COUNTS[index],
            "block_delta_latest_minus_prototype": (
                row["block_count"] - FRESH_PROTO_BLOCK_COUNTS[index]
                if isinstance(row.get("block_count"), int) else None),
            "exact_original_geometry": exact,
            "exact_original_volume_matches_independent_runner": exact_matches_runner,
            "recipe_sha256": hashlib.sha256(
                json.dumps(item["recipe"], separators=(",", ":"), allow_nan=False).encode()
            ).hexdigest(),
        })
        if not exact_matches_runner:
            row["geometry_pass"] = False
            row["outcome"] = "exact_original_volume_disagreement"
        rows.append(row)
        print("supplemental", item["id"], row["outcome"], row.get("block_count"), flush=True)
    canonical = json.dumps([item["recipe"] for item in cases], separators=(",", ":"),
                           allow_nan=False).encode()
    valid_rows = [row for row in rows if row.get("geometry_pass")]
    latest_blocks = [row["block_count"] for row in rows if isinstance(row.get("block_count"), int)]
    proto_blocks = [row["prototype_block_count"] for row in rows]
    condition_values = [condition_max(row) for row in valid_rows if condition_max(row) is not None]
    return {
        "description": "Twelve additional seeded convex stock/cut geometries; not part of the frozen 45-case corpus.",
        "generator": "Python random.Random seed 560127; exact integer monotone-chain hull; stock vertices and cut/Z-interval rules reproduced from the independent fresh-probe definition.",
        "source_definition_sha256": FRESH_PROBE_SOURCE_SHA256,
        "runtime_dependency_on_source_or_archive": False,
        "recipe_set_sha256": hashlib.sha256(canonical).hexdigest(),
        "denominator": len(rows),
        "valid": sum(bool(row.get("geometry_pass")) for row in rows),
        "summary": {
            "latest_release_blocks_total": sum(latest_blocks),
            "prototype_blocks_total_context_only": sum(proto_blocks),
            "block_delta_latest_minus_prototype": sum(latest_blocks) - sum(proto_blocks)
                if len(latest_blocks) == len(proto_blocks) else None,
            "latest_release_sampled_kappa3_max_median": statistics.median(condition_values)
                if condition_values else None,
            "latest_release_sampled_kappa3_max_worst": max(condition_values)
                if condition_values else None,
            "exact_direct_volume_matches": sum(bool(row.get(
                "exact_original_volume_matches_independent_runner")) for row in rows),
        },
        "failures": [{"id": row["id"], "outcome": row.get("outcome"),
                      "stderr": row.get("build_stderr", "")}
                     for row in rows if not row.get("geometry_pass")],
        "prototype_counts_are_context_only": True,
        "cases": rows,
    }


def redact_machine_paths(value, replacements: dict[str, str]):
    if isinstance(value, str):
        for source, label in sorted(replacements.items(), key=lambda item: len(item[0]), reverse=True):
            if source:
                value = value.replace(source, label)
        return value
    if isinstance(value, list):
        return [redact_machine_paths(item, replacements) for item in value]
    if isinstance(value, dict):
        return {key: redact_machine_paths(item, replacements) for key, item in value.items()}
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()
    if not args.cli.is_file():
        raise SystemExit(f"release executable not found: {args.cli}")

    corpus_hash = sha256(CORPUS)
    if corpus_hash != EXPECTED_CORPUS_SHA256:
        raise SystemExit(f"frozen corpus SHA-256 mismatch: {corpus_hash}")
    reference_path = OLD_CANDIDATE
    ref_hash = sha256(reference_path)
    frozen = json.loads(reference_path.read_text())
    if len(frozen.get("cases", [])) != 45:
        raise SystemExit("reference dataset must contain exactly 45 frozen cases")
    if frozen.get("cli_sha256") != EXPECTED_REFERENCE_CLI_SHA256:
        raise SystemExit(f"unexpected reference executable SHA-256: {frozen.get('cli_sha256')}")

    runner = load_runner()
    executable_hash_before = sha256(args.cli)
    result = runner.run(args.cli.resolve(), "latest-release",
                        args.work.resolve(), args.timeout)
    executable_hash_after = sha256(args.cli)
    if executable_hash_before != executable_hash_after:
        raise SystemExit("release executable changed during evaluation")
    old_by_id = {row["id"]: row for row in frozen["cases"]}
    new_by_id = {row["id"]: row for row in result["cases"]}
    if old_by_id.keys() != new_by_id.keys():
        raise SystemExit("latest result IDs do not match the frozen reference dataset")
    pairs = [pair_case(old_by_id[case_id], new_by_id[case_id])
             for case_id in old_by_id]
    git_head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              text=True, capture_output=True, check=True).stdout.strip()
    git_status = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain"],
                                text=True, capture_output=True, check=True).stdout
    script_hash = sha256(Path(__file__))
    supplemental = evaluate_fresh_cases(
        args.cli.resolve(), args.work.resolve() / "supplemental", args.timeout)
    output = {
        "schema": "convex_host_rings_release_comparison_v1",
        "claim_scope": "Public deterministic synthetic corpus only. The separately preregistered holdout is no longer blind.",
        "frozen_inputs": {
            "corpus": "docs/results/sketch-planner-polygon-evaluation/cases.json",
            "corpus_sha256": corpus_hash,
            "case_count": len(pairs),
            "reference_dataset": reference_path.relative_to(ROOT).as_posix(),
            "reference_dataset_sha256": ref_hash,
            "reference_tag": frozen.get("tag"),
            "reference_cli_sha256": frozen.get("cli_sha256"),
        },
        "prototype_snapshot": {
            "observed_case_count": len(PROTO_OBSERVATIONS),
            "original_archive_sha256_for_provenance_only": PROTO_ARCHIVE_SHA256,
            "reviewed_extract_sha256_for_provenance_only": PROTO_REVIEW_SNAPSHOT_SHA256,
            "runtime_archive_dependency": False,
        },
        "latest_release": {
            "binary": args.cli.name,
            "binary_sha256_before": executable_hash_before,
            "binary_sha256_after": executable_hash_after,
            "source_checkout_head": git_head,
            "source_worktree_dirty": bool(git_status.strip()),
            "source_worktree_status_sha256": hashlib.sha256(git_status.encode()).hexdigest(),
            "evaluator_script_sha256": script_hash,
            "timeout_s_per_build_or_audit": args.timeout,
            "elapsed_s": result["elapsed_s"],
        },
        "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                        "cpu_count": os.cpu_count()},
        "methodology": {
            "geometry": "Fresh build and audit of each frozen recipe, independent exact Fraction-based section topology/area/containment/volume checks, plus the CLI closed-domain determinant certificate.",
            "quality": "Per-solid maximum 2-norm Jacobian condition number sampled on the closed 5x5x5 parameter grid; candidate and prototype samples are diagnostics only, not continuous certificates. The prototype also records an in-plane 5x5 sample.",
            "failures": "Refusals, build errors, audit failures, invalid geometry and timeouts remain in the fixed 45-case denominator.",
            "prototype": "Thirty archived prototype block-count and sampled-condition observations are embedded as literals. Their source archive is not read at runtime; the prototype did not provide the kernel's continuous certificate.",
            "supplemental": "Twelve additional deterministic geometries are generated with the standard library, checked with the Rust CLI and independent exact section oracle, then volume is cross-checked against direct exact shoelace-area integration.",
        },
        "summary": summarize(pairs),
        "seeded_30_case_breakdown": seeded_breakdown(pairs),
        "cases": pairs,
        "latest_run_details": result,
        "supplemental_fresh_cases": supplemental,
    }
    output = redact_machine_paths(output, {
        str(args.cli.resolve()): "solidcraft-trivariate",
        str(args.work.resolve()): "<evaluation-workdir>",
        str(ROOT): ".",
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "cases": len(pairs),
                      "new_valid": sum(r["geometry_pass"] for r in result["cases"]),
                      "binary_sha256": executable_hash_after}))


if __name__ == "__main__":
    main()
