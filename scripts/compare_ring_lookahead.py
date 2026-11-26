#!/usr/bin/env python3
"""Paired independent geometry and quality evaluation of ring lookahead."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
RUN_PATH = ROOT / "docs/results/sketch-planner-polygon-evaluation/run.py"
METRICS_PATH = ROOT / "scripts/evaluate_polygon_extrusion_quality.py"


def import_from_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


runner = import_from_path("frozen_polygon_runner", RUN_PATH)
planar = import_from_path("independent_planar_metrics", METRICS_PATH)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def nearest_rank(values, p):
    if not values:
        return None
    return sorted(values)[max(0, math.ceil(p * len(values)) - 1)]


def source_revision(binary_path):
    path = Path(binary_path).resolve()
    for directory in (path.parent, *path.parents):
        if (directory / ".git").exists():
            result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=directory,
                                    text=True, capture_output=True)
            if result.returncode == 0:
                return result.stdout.strip()
    return None


def baseline_inputs(ledger, corpus_path):
    frozen = json.loads(Path(corpus_path).read_text())
    cases = frozen["cases"]
    if len(cases) != 45 or len({c["id"] for c in cases}) != 45:
        raise ValueError("expected exactly 45 unique frozen historical cases")
    baseline = ledger["latest_run_details"]["cases"]
    supplemental = ledger["supplemental_fresh_cases"]["cases"]
    if len(baseline) != 45 or len(supplemental) != 12:
        raise ValueError("baseline ledger must contain 45 historical and 12 supplemental cases")
    by_id = {c["id"]: c for c in baseline}
    if set(by_id) != {c["id"] for c in cases}:
        raise ValueError("baseline historical IDs do not match frozen cases.json")
    supplement_by_id = {c["id"]: c for c in supplemental}
    supplements = []
    for case in supplemental:
        recipe = case.get("recipe")
        if not recipe:
            raise ValueError(f"missing frozen supplemental recipe for {case['id']}")
        supplements.append({"id": case["id"], "family": case["family"], "recipe": recipe})
    inputs = [{"id": c["id"], "family": c["family"], "recipe": c["recipe"]} for c in cases]
    inputs += supplements
    if len(inputs) != 57 or len({x["id"] for x in inputs}) != 57:
        raise ValueError("expected 57 unique total case IDs")
    return inputs, by_id, supplement_by_id


def load_model_metrics(path):
    path = Path(path)
    return {**planar.summarize_model(json.loads(path.read_text())),
            "model_path": str(path), "model_sha256": sha256(path)}


def metric_pair(base, candidate):
    ratios = {
        "planar_spectral_kappa": candidate["planar_spectral_kappa"] / base["planar_spectral_kappa"],
        "thickness_direction_aspect": candidate["thickness_direction_aspect"] / base["thickness_direction_aspect"],
        "max_3d_spectral_kappa": candidate["max_3d_spectral_kappa"] / base["max_3d_spectral_kappa"],
        "corner_sine_worsening_ratio": (base["minimum_corner_sine"] / candidate["minimum_corner_sine"]
                                        if candidate["minimum_corner_sine"] else math.inf),
    }
    def worse(ratio):
        return ratio > 1.0 + 1e-10
    return {
        "baseline": base,
        "candidate": candidate,
        "candidate_over_baseline_ratio": ratios,
        "worsened": {"planar_spectral_kappa": worse(ratios["planar_spectral_kappa"]),
                     "thickness_direction_aspect": worse(ratios["thickness_direction_aspect"]),
                     "max_3d_spectral_kappa": worse(ratios["max_3d_spectral_kappa"]),
                     "minimum_corner_sine": worse(ratios["corner_sine_worsening_ratio"])},
    }


def make_case_rows(inputs, baseline_rows, baseline_models, candidate_rows, candidate_models):
    rows = []
    for item in inputs:
        ident = item["id"]
        b, c = baseline_rows[ident], candidate_rows[ident]
        base_metric = baseline_models.get(ident)
        cand_metric = candidate_models.get(ident)
        quality = metric_pair(base_metric, cand_metric) if base_metric and cand_metric else None
        bcount, ccount = b.get("block_count"), c.get("block_count")
        elapsed_b = (b.get("build_elapsed_s") or 0.0) + (b.get("audit_elapsed_s") or 0.0)
        elapsed_c = (c.get("build_elapsed_s") or 0.0) + (c.get("audit_elapsed_s") or 0.0)
        rows.append({
            "id": ident,
            "family": item["family"],
            "baseline_geometry_pass": bool(b.get("geometry_pass")),
            "candidate_geometry_pass": bool(c.get("geometry_pass")),
            "baseline_outcome": b.get("outcome"),
            "candidate_outcome": c.get("outcome"),
            "geometry_regression": bool(b.get("geometry_pass")) and not bool(c.get("geometry_pass")),
            "volume_match_candidate": c.get("volume_match"),
            "baseline_volume_exact_mm3": b.get("volume_exact_expected_mm3"),
            "candidate_volume_exact_expected_mm3": c.get("volume_exact_expected_mm3"),
            "candidate_volume_exact_cli_mm3": c.get("volume_exact_cli_mm3"),
            "candidate_certificate": c.get("closed_domain_certificate"),
            "baseline_blocks": bcount,
            "candidate_blocks": ccount,
            "baseline_construction_method": b.get("construction_method"),
            "candidate_construction_method": c.get("construction_method"),
            "block_delta_candidate_minus_baseline": ccount - bcount if bcount is not None and ccount is not None else None,
            "same_block_count": bcount == ccount if bcount is not None and ccount is not None else None,
            "baseline_model_sha256": base_metric.get("model_sha256") if base_metric else None,
            "candidate_model_sha256": cand_metric.get("model_sha256") if cand_metric else None,
            "original_recipe_sha256": hashlib.sha256(json.dumps(item["recipe"], sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "baseline_build_audit_s": elapsed_b,
            "candidate_build_audit_s": elapsed_c,
            "quality": quality,
        })
    return rows


def summarize(rows):
    outcomes = {}
    ids = {r["id"]: r for r in rows}
    geometry_passes = {"baseline": sum(r["baseline_geometry_pass"] for r in rows),
                       "candidate": sum(r["candidate_geometry_pass"] for r in rows)}
    compared = [r for r in rows if r["quality"]]
    block_deltas = [r["block_delta_candidate_minus_baseline"] for r in rows
                    if r["block_delta_candidate_minus_baseline"] is not None]
    common_success = [r for r in rows if r["baseline_geometry_pass"] and r["candidate_geometry_pass"]]
    timings_b = [r["baseline_build_audit_s"] for r in common_success]
    timings_c = [r["candidate_build_audit_s"] for r in common_success]
    quality_summary = {}
    fields = ("planar_spectral_kappa", "thickness_direction_aspect", "max_3d_spectral_kappa",
              "corner_sine_worsening_ratio")
    for field in fields:
        adverse = [r["quality"]["candidate_over_baseline_ratio"][field] for r in compared]
        if field == "corner_sine_worsening_ratio":
            worsened = [r["id"] for r in compared if r["quality"]["worsened"]["minimum_corner_sine"]]
        else:
            worsened = [r["id"] for r in compared if r["quality"]["worsened"][field]]
        quality_summary[field] = {"n": len(adverse), "median_worsening_ratio": statistics.median(adverse) if adverse else None,
                                  "p95_worsening_ratio": nearest_rank(adverse, .95),
                                  "max_worsening_ratio": max(adverse) if adverse else None,
                                  "worsened_case_ids": worsened}
    for family in sorted({r["family"] for r in rows}):
        group = [r for r in rows if r["family"] == family]
        outcomes[family] = {"n": len(group), "geometry_pass_baseline": sum(r["baseline_geometry_pass"] for r in group),
                            "geometry_pass_candidate": sum(r["candidate_geometry_pass"] for r in group),
                            "geometry_regressions": [r["id"] for r in group if r["geometry_regression"]],
                            "same_block_count": sum(r["same_block_count"] is True for r in group),
                            "baseline_blocks": sum(r["baseline_blocks"] or 0 for r in group),
                            "candidate_blocks": sum(r["candidate_blocks"] or 0 for r in group),
                            "block_delta_total": sum(r["block_delta_candidate_minus_baseline"] or 0 for r in group)}
    return {
        "case_count": len(rows),
        "geometry_passes": geometry_passes,
        "geometry_regressions": [r["id"] for r in rows if r["geometry_regression"]],
        "block_count_comparison": {"n": len(block_deltas), "equal": sum(r["same_block_count"] is True for r in rows),
                                   "increased": sum(d > 0 for d in block_deltas), "decreased": sum(d < 0 for d in block_deltas),
                                   "baseline_total": sum(r["baseline_blocks"] or 0 for r in rows),
                                   "candidate_total": sum(r["candidate_blocks"] or 0 for r in rows),
                                   "delta_total": sum(block_deltas), "delta_median": statistics.median(block_deltas) if block_deltas else None},
        "quality_worsening": quality_summary,
        "cost_common_success_s": {"n": len(common_success), "baseline_median": statistics.median(timings_b) if timings_b else None,
                                   "candidate_median": statistics.median(timings_c) if timings_c else None,
                                   "baseline_total": sum(timings_b), "candidate_total": sum(timings_c)},
        "by_family": outcomes,
    }


def render_report(result):
    summary = result["summary"]
    blocks = summary["block_count_comparison"]
    costs = summary["cost_common_success_s"]
    lines = [
        "# Paired evaluation: planar ring lookahead",
        "",
        f"Candidate source revision: `{result['candidate'].get('source_revision')}`; binary SHA-256: `{result['candidate']['cli_sha256']}`.",
        f"Baseline ledger SHA-256: `{result['baseline_provenance']['ledger_sha256']}`.",
        "",
        "## Cohort and geometric gate",
        "",
        f"The cohort contains {summary['case_count']} fixed recipes: 45 historical cases and 12 supplements. The candidate recipe adds only the explicit `partition_strategy=ring_planar_lookahead` selector. Both sides use the frozen `run.py` build/audit and exact Fraction section checks.",
        "",
        f"Geometry passes: baseline {summary['geometry_passes']['baseline']}/{summary['case_count']}; candidate {summary['geometry_passes']['candidate']}/{summary['case_count']}.",
        f"Block counts: {blocks['equal']} equal, {blocks['decreased']} decreased, {blocks['increased']} increased; total {blocks['baseline_total']} → {blocks['candidate_total']} (delta {blocks['delta_total']}).",
        f"Common-success build+audit cost: median {costs['baseline_median']:.4f} s baseline / {costs['candidate_median']:.4f} s candidate; cohort sums {costs['baseline_total']:.3f} s / {costs['candidate_total']:.3f} s. These are descriptive single-run timings.",
        "",
        "## Independent quality comparison",
        "",
        "All ratios are adverse-direction ratios where values above 1 mean worse: candidate/baseline for condition numbers and thickness aspect, baseline/candidate for minimum corner sine. Metrics use the independent planar evaluator; extrema are sampled, not continuous bounds. A relative change must exceed 1e−10 to be counted as worsening.",
        "",
        "| Metric | Paired n | Median adverse ratio | p95 | Maximum | Cases worsened |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for field, label in (("planar_spectral_kappa", "Planar spectral κ"),
                         ("thickness_direction_aspect", "Normal thickness aspect"),
                         ("max_3d_spectral_kappa", "Full 3D spectral κ"),
                         ("corner_sine_worsening_ratio", "Minimum corner sine")):
        row = summary["quality_worsening"][field]
        lines.append(f"| {label} | {row['n']} | {row['median_worsening_ratio']:.6g} | {row['p95_worsening_ratio']:.6g} | {row['max_worsening_ratio']:.6g} | {len(row['worsened_case_ids'])} |")
    regressions = summary["geometry_regressions"]
    lines += ["", f"Geometry regression cases: {', '.join(regressions) if regressions else 'none'}.", "",
              "## Per-family counts", "", "| Family | n | Geometry baseline/candidate | Equal blocks | Block delta |", "| --- | ---: | ---: | ---: | ---: |"]
    for family, values in summary["by_family"].items():
        lines.append(f"| {family} | {values['n']} | {values['geometry_pass_baseline']}/{values['geometry_pass_candidate']} | {values['same_block_count']} | {values['block_delta_total']:+d} |")
    worsened_rows = [r for r in result["cases"] if r["quality"] and any(r["quality"]["worsened"].values())]
    lines += ["", "## Cases with at least one quality regression", "",
              "| Case | Block Δ | Planar κ ratio | Thickness ratio | 3D κ ratio | Corner-sine adverse ratio |", "| --- | ---: | ---: | ---: | ---: | ---: |"]
    if worsened_rows:
        for r in worsened_rows:
            q = r["quality"]["candidate_over_baseline_ratio"]
            lines.append(f"| {r['id']} | {r['block_delta_candidate_minus_baseline']} | {q['planar_spectral_kappa']:.6g} | {q['thickness_direction_aspect']:.6g} | {q['max_3d_spectral_kappa']:.6g} | {q['corner_sine_worsening_ratio']:.6g} |")
    else:
        lines.append("| none | — | — | — | — | — |")
    lines += ["", "## Reproduction", "", "Run the comparison script from the checkout root with the frozen release binary, baseline review ledger, and baseline model directory:", "",
              "```sh",
              "python3 scripts/compare_ring_lookahead.py \\",
              "  --cli /path/to/solidcraft-trivariate \\",
              "  --baseline-ledger /path/to/clean-ledger.json \\",
              "  --baseline-models /path/to/replay \\",
              "  --work build/ring-lookahead-candidate \\",
              "  --output-dir docs/results/polygon-extrusion-quality/ring-lookahead",
              "```", "",
              "Machine-readable paired results are in `paired-results.json`; frozen supplement recipes are in `supplemental-cases.json`. The candidate command also writes its detailed raw run ledger under the selected work directory.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", required=True, type=Path)
    parser.add_argument("--baseline-ledger", required=True, type=Path)
    parser.add_argument("--baseline-models", required=True, type=Path,
                        help="directory with historical cases and supplemental/ models")
    parser.add_argument("--corpus", default=ROOT / "docs/results/sketch-planner-polygon-evaluation/cases.json", type=Path)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--strategy", default="ring_planar_lookahead")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    baseline_ledger = json.loads(args.baseline_ledger.read_text())
    inputs, baseline_main, baseline_supp = baseline_inputs(baseline_ledger, args.corpus)
    supplemental_input_path = args.output_dir / "supplemental-cases.json"
    supplemental_inputs = [x for x in inputs if x["id"].startswith("fresh_convex_")]
    supplemental_input_path.write_text(json.dumps({"description": "Frozen supplemental case recipes from the baseline review ledger; unchanged for candidate build except the explicit strategy selector.",
                                                   "recipe_set_sha256": baseline_ledger["supplemental_fresh_cases"]["recipe_set_sha256"],
                                                   "cases": supplemental_inputs}, indent=2, sort_keys=True) + "\n")
    candidate_rows = []
    for index, item in enumerate(inputs, start=1):
        candidate_item = {**item, "recipe": dict(item["recipe"])}
        candidate_item["recipe"]["partition_strategy"] = args.strategy
        row = runner.eval_one(args.cli.resolve(), candidate_item, args.work, args.timeout)
        row["strategy"] = args.strategy
        row["original_recipe_sha256"] = hashlib.sha256(json.dumps(item["recipe"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        model_path = args.work / f"{item['id']}.model.json"
        row["model_sha256"] = sha256(model_path) if model_path.is_file() else None
        candidate_rows.append(row)
        print(f"[{index}/57] {item['id']}: {row.get('outcome')} blocks={row.get('block_count')}", flush=True)
    candidate_ledger = {"tag": args.strategy, "cli": str(args.cli.resolve()), "cli_sha256": sha256(args.cli),
                        "source_revision": source_revision(args.cli), "recipe_strategy": args.strategy,
                        "timeout_s_per_build_or_audit": args.timeout, "cases": candidate_rows}
    raw_path = args.work / "candidate-run-ledger.json"
    raw_path.write_text(json.dumps(candidate_ledger, indent=2, sort_keys=True) + "\n")

    base_rows = {**baseline_main, **baseline_supp}
    baseline_metrics = {}
    candidate_metrics = {}
    for item in inputs:
        ident = item["id"]
        subdir = "supplemental" if ident.startswith("fresh_convex_") else ""
        base_path = args.baseline_models / subdir / f"{ident}.model.json"
        candidate_path = args.work / f"{ident}.model.json"
        if base_path.is_file():
            try:
                baseline_metrics[ident] = load_model_metrics(base_path)
            except Exception as exc:
                baseline_metrics[ident] = {"metric_error": str(exc)}
        if candidate_path.is_file():
            try:
                candidate_metrics[ident] = load_model_metrics(candidate_path)
            except Exception as exc:
                candidate_metrics[ident] = {"metric_error": str(exc)}
    # Quality results only compare if both independent summaries parsed successfully.
    valid_base = {k: v for k, v in baseline_metrics.items() if "metric_error" not in v}
    valid_candidate = {k: v for k, v in candidate_metrics.items() if "metric_error" not in v}
    pairs = make_case_rows(inputs, base_rows, valid_base, {r["id"]: r for r in candidate_rows}, valid_candidate)
    summary = summarize(pairs)
    summary["quality_metric_coverage"] = {"baseline": len(valid_base), "candidate": len(valid_candidate),
                                          "paired": sum(r["quality"] is not None for r in pairs),
                                          "baseline_metric_errors": {k: v["metric_error"] for k, v in baseline_metrics.items() if "metric_error" in v},
                                          "candidate_metric_errors": {k: v["metric_error"] for k, v in candidate_metrics.items() if "metric_error" in v}}
    result = {
        "schema": "polygon_ring_lookahead_paired_evaluation_v1",
        "claim_scope": "Same 45 frozen historical and 12 frozen supplemental convex polygon recipes; explicit opt-in strategy key only changes construction selection.",
        "candidate": candidate_ledger,
        "baseline_provenance": {"ledger_sha256": sha256(args.baseline_ledger),
                                "baseline_model_root": str(args.baseline_models),
                                "baseline_source_revision": baseline_ledger.get("latest_release", {}).get("source_checkout_head"),
                                "supplemental_recipe_set_sha256": baseline_ledger.get("supplemental_fresh_cases", {}).get("recipe_set_sha256"),
                                "baseline_run_ledger_sha256": sha256(args.baseline_ledger)},
        "input_provenance": {"historical_corpus_path": str(args.corpus), "historical_corpus_sha256": sha256(args.corpus),
                             "supplemental_recipes_path": str(supplemental_input_path),
                             "supplemental_recipes_sha256": sha256(supplemental_input_path)},
        "methodology": {"geometric_gate": "docs/results/sketch-planner-polygon-evaluation/run.py eval_one; exact independent Fraction section topology/area/containment/volume checks and CLI closed-domain certificate.",
                        "planar_metrics": "scripts/evaluate_polygon_extrusion_quality.py summarize_model; 5x5x5 planar spectral kappa and normal thickness aspect, corner sine at corners; separate full 3D spectral kappa.",
                        "worsening_threshold_relative": 1e-10,
                        "cost": "per-case build plus audit wall time; descriptive measurements"},
        "summary": summary,
        "cases": pairs,
        "raw_candidate_ledger_path": str(raw_path),
        "raw_candidate_ledger_sha256": sha256(raw_path),
    }
    (args.output_dir / "paired-results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (args.output_dir / "README.md").write_text(render_report(result))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
