#!/usr/bin/env python3
"""Rebuild the current FieldCraft bracket and measure its PA5 energy convergence."""

from __future__ import annotations

import xml.etree.ElementTree as ET
import argparse
import itertools
import json
import hashlib
import math
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build" / "fieldcraft-bracket-pa5"
RECIPE = ROOT / "docs" / "site" / "data" / "fieldcraft_x64_y64" / "fieldcraft_x64_y64.recipe.json"
JOB = ROOT / "tests" / "references" / "mechanical_jobs" / "fieldcraft_bracket_hex_pa5.json"


def cli(*args: str) -> None:
    subprocess.run(
        [str(CLI), *args],
        cwd=ROOT,
        check=True,
    )


def uniform_plan(block_ids: list[str], level: int) -> dict:
    steps = []
    for depth in range(level):
        width = 1 << depth
        leaves = [
            {"root": root, "level": depth, "index": list(index)}
            for root in block_ids
            for index in itertools.product(range(width), repeat=3)
        ]
        steps.append({"refine": leaves})
    return {
        "schema": "adaptive_refinement_plan_v1",
        "max_depth": level,
        "max_cells": max(1, len(block_ids) * (8**level) * 2),
        "steps": steps,
    }


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def solve(family: str, model_path: Path, order: int, level: int) -> dict:
    model = load_json(model_path)
    job = load_json(JOB)
    job["order"] = order
    # Select the exterior L-profile face by its persisted design provenance.
    tags = {face["feature"] for face in model["design_faces"]}
    profile_tag = "outer/profile" if "outer/profile" in tags else "outer/outline"
    for selection in [job["supports"][0]["select"], job["loads"][0]["select"]]:
        selection["feature"] = profile_tag
    job["plan"] = uniform_plan([b["id"] for b in model["blocks"]], level)
    job_path = OUT / f"{family}-q{order}-h{level}.job.json"
    result_path = OUT / f"{family}-q{order}-h{level}.result.json"
    job_path.write_text(json.dumps(job, indent=2) + "\n")
    start = time.perf_counter()
    cli("solve-job", str(model_path), str(job_path), str(result_path))
    result = load_json(result_path)
    result["_wall_seconds"] = time.perf_counter() - start
    return result


def number(obj: dict, key: str) -> float:
    return float(obj[key])


def render_existing(records: list[dict]) -> None:
    import os
    mpl_config = OUT / "matplotlib-config"
    mpl_config.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    import matplotlib.pyplot as plt

    q16 = [r for r in records if r["quadrature_order"] == 16 and r["family"] == "10-block"]
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    rows = sorted(q16, key=lambda row: row["scalar_dofs"])
    ax.plot([r["scalar_dofs"] for r in rows], [r["strain_energy_n_mm"] for r in rows],
            color="#147D92", marker="o", linewidth=2.1, label="Current sketch partition")
    ax.set_xscale("log")
    dof_ticks = sorted({r["scalar_dofs"] for r in q16})
    ax.set_xticks(dof_ticks, [f"{x:,}" for x in dof_ticks])
    ax.set_xlabel("Global scalar Q2 degrees of freedom")
    ax.set_ylabel("Strain energy (N·mm)")
    ax.set_title("FieldCraft bracket: energy versus analysis size")
    ax.grid(True, which="both", color="#D9DEE5", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, title="Geometric partition")
    fig.text(0.01, -0.015, "Same geometry, supports, +2000 N traction, steel, and q16 integration; full uniform octree refinement.", fontsize=8.5, color="#475569")
    for extension in ("png", "svg"):
        fig.savefig(OUT / f"fieldcraft-bracket-pa5-energy-dof.{extension}", dpi=220 if extension == "png" else None, bbox_inches="tight")
    plt.close(fig)


def geometry_condition(model: Path) -> float:
    vtu = OUT / f"{model.stem}.quality.vtu"
    cli("export-geometry-vtu", str(model), "4", str(vtu))
    arrays = ET.parse(vtu).findall('.//Piece/PointData/DataArray[@Name="condition_number"]')
    samples = [list(map(float, array.text.split())) for array in arrays]
    if len(samples) != len(load_json(model)["blocks"]) or any(len(v) != 125 for v in samples):
        raise ValueError("incomplete block quality sample")
    return max(max(block) for block in samples)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-level", type=int, default=2, choices=range(0, 4))
    parser.add_argument("--cli", type=Path, default=ROOT / "target" / "release" / "solidcraft-trivariate")
    parser.add_argument("--orders", type=int, nargs="+", default=[8, 16, 24, 32])
    parser.add_argument("--render-only", action="store_true", help="re-render the current build JSON without running mechanics")
    args = parser.parse_args()
    if set(args.orders) != {8, 16, 24, 32}:
        parser.error("the quadrature study requires exactly orders 8, 16, 24, and 32")
    global CLI
    CLI = args.cli.resolve()
    if not CLI.is_file():
        raise FileNotFoundError(f"native CLI binary not found: {CLI}")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.render_only:
        render_existing(load_json(OUT / "fieldcraft-bracket-pa5-energy-dof.json")["records"])
        return

    ten_model = OUT / "fieldcraft_x64_y64_10.model.json"
    cli("build", str(RECIPE), str(ten_model))
    family = "10-block"
    expected = len(load_json(ten_model)["blocks"])
    if expected != 10:
        raise RuntimeError(f"FieldCraft sketch produced {expected} blocks; expected 10")
    cli("audit", str(ten_model))
    records = []
    levels = range(args.max_level + 1)
    for level in levels:
        for order in (args.orders if level == 0 else [16]):
            result = solve(family, ten_model, order, level)
            actual_blocks = len(result["analysis_blocks"])
            expected_blocks = expected * 8**level
            if actual_blocks != expected_blocks:
                raise RuntimeError(f"{family} h{level}: got {actual_blocks} cells, expected {expected_blocks}")
            energy = number(result, "strain_energy_n_mm")
            area = number(result["loads"][0], "area_mm2")
            force = result["applied_resultant_n"]
            if not math.isclose(area, 320.0, rel_tol=2e-10, abs_tol=2e-9):
                raise RuntimeError(f"{family} h{level} q{order}: load area mismatch: {area}")
            if not math.isclose(float(force[0]), 2000.0, rel_tol=2e-10, abs_tol=2e-8):
                raise RuntimeError(f"{family} h{level} q{order}: load resultant mismatch: {force}")
            if order >= 16 and not math.isclose(number(result, "volume_mm3"), 32144.778925, rel_tol=5e-8, abs_tol=1e-7):
                raise RuntimeError(f"{family} h{level} q{order}: matched geometry volume mismatch")
            if not math.isfinite(energy) or energy <= 0:
                raise RuntimeError(f"{family} h{level} q{order}: invalid strain energy")
            records.append({
                "family": family,
                "base_blocks": expected,
                "blocks": len(result["analysis_blocks"]),
                "level": level,
                "quadrature_order": order,
                "scalar_dofs": int(result["scalar_dofs"]),
                "vector_dofs": 3 * int(result["scalar_dofs"]),
                "volume_mm3": number(result, "volume_mm3"),
                "strain_energy_n_mm": energy,
                "external_work_n_mm": number(result, "external_work_n_mm"),
                "relative_work_energy_error": number(result, "relative_work_energy_error"),
                "relative_residual": number(result, "relative_residual"),
                "iterations": int(result["iterations"]),
                "loaded_area_mm2": area,
                "applied_force_n": force,
                "support_reaction_n": result["support_reaction_n"],
                "relative_force_balance": number(result, "relative_force_balance"),
                "relative_moment_balance": number(result, "relative_moment_balance"),
                "wall_seconds": number(result, "_wall_seconds"),
                "source_result": str((OUT / f"{family}-q{order}-h{level}.result.json").relative_to(ROOT)),
            })


    h0 = {r["quadrature_order"]: r for r in records if r["family"] == family and r["level"] == 0}
    q8, q16, q32 = h0[8], h0[16], h0[32]
    q24 = h0[24]
    q8["h0_q8_to_q16_energy_relative_difference"] = abs(q16["strain_energy_n_mm"] / q8["strain_energy_n_mm"] - 1)
    q16["h0_q16_to_q32_energy_relative_difference"] = abs(q32["strain_energy_n_mm"] / q16["strain_energy_n_mm"] - 1)
    q16["h0_q16_to_q32_volume_relative_difference"] = abs(q32["volume_mm3"] / q16["volume_mm3"] - 1)
    q24["h0_q24_to_q32_energy_relative_difference"] = abs(q32["strain_energy_n_mm"] / q24["strain_energy_n_mm"] - 1)
    q24["h0_q24_to_q32_volume_relative_difference"] = abs(q32["volume_mm3"] / q24["volume_mm3"] - 1)

    output = {
        "schema": "fieldcraft_bracket_pa5_energy_dof_v1",
        "protocol": {
            "case": "FieldCraft filleted-bracket-fillet-stress PA5 load transfer",
            "geometry": "same original bracket, 60 x 50 x 40 mm L outline, R6 inside fillet, two diameter-8 mm through-bores",
            "geometries": {
                "10-block": "rebuilt from docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.recipe.json",
            },
            "supports": "all displacement components fixed on y=0 design face selected by outer/outline provenance",
            "traction": "+6.25 MPa in x on y=50 face; 320 mm^2 loaded area; target resultant +2000 N",
            "material": {"young_mpa": 210000.0, "poisson": 0.3},
            "formulation": "native Rust Q2 displacement elasticity solve via solve-job; full octree refinement of every block at each level",
            "quadrature_orders": args.orders,
            "levels": {family: list(levels)},
            "unit_convention": "N, mm, MPa; strain energy N mm",
            "timing_note": "wall_seconds is end-to-end solve-job time; includes model audit, refinement, assembly, and solve",
            "quadrature_qualification": "h0 q8, q16, q24, q32 are compared for the compact space; q16/q32 energy agreement is checked explicitly. q32 provides the volume cross-check. All refined runs use q16.",
            "sampled_geometry_condition_number_max": {
                family: geometry_condition(ten_model)
            },
            "condition_number_method": "Independent export-geometry-vtu sampling on {0, 1/4, 1/2, 3/4, 1}^3 in every block; see docs/results/fieldcraft-compact.md.",
            "cli_binary": CLI.name,
            "cli_sha256": hashlib.sha256(CLI.read_bytes()).hexdigest(),
            "input_sha256": {
                "ten_block_recipe": hashlib.sha256(RECIPE.read_bytes()).hexdigest(),
                "mechanical_job_fixture": hashlib.sha256(JOB.read_bytes()).hexdigest(),
            },
            "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        },
        "records": records,
    }
    json_path = OUT / "fieldcraft-bracket-pa5-energy-dof.json"
    json_path.write_text(json.dumps(output, indent=2) + "\n")
    render_existing(records)
    print(json.dumps({"result": str(json_path), "records": len(records), "plot": str(OUT / "fieldcraft-bracket-pa5-energy-dof.png")}, indent=2))


if __name__ == "__main__":
    main()
