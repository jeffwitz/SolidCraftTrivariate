"""Run a low-cost matched-mesh HEX8/HEX20 study for the tilted S2 union."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))


def load_study(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def model_from_cli(binary: Path, rise: int) -> dict:
    recipe = json.loads(
        (ROOT / "tests/references/g2_cross_sketches.recipe.json").read_text()
    )
    document = json.dumps(
        dict(
            family="inclined_crossed_tools_v1",
            sketch_pair=recipe,
            transverse_rise_mm=rise,
        )
    )
    with tempfile.TemporaryDirectory(prefix="s2-hex20-") as tmp:
        source, output = Path(tmp) / "recipe.json", Path(tmp) / "model.json"
        source.write_text(document)
        subprocess.run(
            [str(binary), "crossed-boolean", str(source), "union", str(output)],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(output.read_text())


def iga_dof_counts(study, model: dict) -> dict:
    blocks = model["blocks"]
    nets = [np.asarray(block["controls"], float) for block in blocks]
    import re

    ids = {
        tuple(
            map(int, re.fullmatch(r"band(\d+)_col(\d+)_z.*", block["id"]).groups())
        ): i
        for i, block in enumerate(blocks)
    }
    full, _ = study.constraints(nets, ids, np.linspace(0, 1, 7))
    basis = study.exact_trace_space(nets, ids, full)
    clamp = []
    for (band, col), block in ids.items():
        if col == 0:
            clamp.extend(16 * block + np.arange(16).reshape(4, 2, 2)[:, 1, :].ravel())
    free_basis = null_space(basis[clamp, :], rcond=1e-11)
    return {
        "charts": len(nets),
        "scalar_total": int(basis.shape[1]),
        "vector_total": int(3 * basis.shape[1]),
        "scalar_free_after_clamp": int(free_basis.shape[1]),
        "vector_free_after_clamp": int(3 * free_basis.shape[1]),
    }


def brief_hex(row: dict, node_count: int) -> dict:
    total_dofs = 3 * node_count
    return {
        "element": row["element"],
        "elements": row["hexes"],
        "nodes": node_count,
        "total_displacement_dofs": total_dofs,
        "fixed_displacement_dofs": total_dofs - row["free_dofs"],
        "free_displacement_dofs": row["free_dofs"],
        "loaded_end_traction_weighted_mean_ux_mm": row["ux_mm"],
        "strain_energy_Nmm": row["energy_Nmm"],
        "volume_mm3": row["volume_mm3"],
        "minimum_gauss_jacobian": row["min_gauss_det"],
        "affine_energy_relative_error": row["affine_energy_relative_error"],
        "traction_moment_max_gap_N": row["max_gap_N"],
        "resultant_force_N": row["force_N"],
    }


def topology_counts(cells: np.ndarray) -> dict:
    face_local = (
        (0, 1, 2, 3),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
    )
    owners = {}
    parent = list(range(len(cells)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for eid, cell in enumerate(cells):
        for local in face_local:
            key = tuple(sorted(int(cell[j]) for j in local))
            owners.setdefault(key, []).append(eid)
    overused = sum(len(owners[f]) > 2 for f in owners)
    for incident in owners.values():
        if len(incident) == 2:
            a, b = map(root, incident)
            parent[a] = b
    return {
        "unique_corner_faces": len(owners),
        "boundary_faces": sum(len(v) == 1 for v in owners.values()),
        "interior_faces": sum(len(v) == 2 for v in owners.values()),
        "overused_faces": overused,
        "cell_components": len({root(i) for i in range(len(cells))}),
    }


def free_dofs(X: np.ndarray, C: np.ndarray, quadratic: bool, hexstudy) -> dict:
    nodes = len(X)
    if quadratic:
        X, _ = hexstudy.elevate_hexes(X, C)
        nodes = len(X)
    total = 3 * nodes
    fixed_nodes = int(np.count_nonzero(X[:, 0] < 1e-9))
    return {
        "nodes": nodes,
        "total_dofs": total,
        "fixed_dofs": 3 * fixed_nodes,
        "free_dofs": total - 3 * fixed_nodes,
    }


def enriched_targets() -> tuple[dict, str]:
    path = HERE / "smooth_bending.json"
    document = json.loads(path.read_text())
    run = next(
        row
        for row in document["runs"]
        if row["rise_mm"] == 3 and row["spans_vw"] == [2, 2]
    )
    targets = {
        field["continuity"]: {
            "free_vector_dofs": field["free_vector_dofs"],
            "loaded_end_traction_weighted_mean_ux_mm": field["mean_ux_mm"],
            "relative_equilibrium_residual": field["relative_residual"],
            "integration_cells": 88,
            "source_volume_mm3": run["volume_mm3"],
        }
        for field in run["fields"]
    }
    return targets, hashlib.sha256(path.read_bytes()).hexdigest()


def find_dof_brackets(model, hexstudy, targets: dict, rise: int) -> tuple[dict, dict]:
    candidates = []
    attempted = 0
    rejected = 0
    for circumferential in (1, 2):
        for wing_a in range(1, 14):
            for core in (1, 2, 3):
                for wing_b in (1, 2, 3, 4):
                    divisions = (wing_a, core, wing_b)
                    attempted += 1
                    options = {
                        "circumferential": circumferential,
                        "nv_by_col": divisions,
                        "height": 1,
                    }
                    try:
                        X, C, _, traces, _ = hexstudy.mesh(model, 1, rise, **options)
                    except (AssertionError, ValueError, IndexError):
                        rejected += 1
                        continue
                    topology = topology_counts(C)
                    if topology["cell_components"] != 1 or topology["overused_faces"]:
                        rejected += 1
                        continue
                    dofs = free_dofs(X, C, True, hexstudy)
                    candidates.append(
                        {
                            "options": options,
                            "free_dofs": dofs["free_dofs"],
                            "nodes": dofs["nodes"],
                            "total_dofs": dofs["total_dofs"],
                            "elements": len(C),
                            "partial_trace_quads": len(traces),
                            "topology": topology,
                        }
                    )
    results = {}
    selected_configs = {}
    for continuity, target in targets.items():
        desired = target["free_vector_dofs"]
        below = [c for c in candidates if c["free_dofs"] < desired]
        exact = [c for c in candidates if c["free_dofs"] == desired]
        above = [c for c in candidates if c["free_dofs"] > desired]
        pick = exact[:1] if exact else []
        if not pick and below:
            pick.append(max(below, key=lambda c: c["free_dofs"]))
        if not exact and above:
            pick.append(min(above, key=lambda c: c["free_dofs"]))
        match_status = (
            "exact" if exact else "bracketed" if below and above else "one_sided"
        )
        for candidate in pick:
            key = tuple(
                candidate["options"][k]
                if k != "nv_by_col"
                else tuple(candidate["options"][k])
                for k in ("circumferential", "nv_by_col", "height")
            )
            if key not in selected_configs:
                selected_configs[key] = candidate
        results[continuity] = {
            "target": target,
            "match_status": match_status,
            "selected": pick,
        }
    return results, {
        "candidates_validated": len(candidates),
        "candidates_attempted": attempted,
        "candidates_rejected": rejected,
        "configs": selected_configs,
        "candidate_dof_values": sorted({c["free_dofs"] for c in candidates}),
        "candidate_records": candidates,
    }


def nearest_dof_scan(candidate_records: list[dict], target: int) -> dict:
    values = sorted({row["free_dofs"] for row in candidate_records})
    lower = max((value for value in values if value < target), default=None)
    upper = min((value for value in values if value > target), default=None)
    exact = target in values
    return {
        "target_free_vector_dofs": target,
        "exact_match_found": exact,
        "nearest_lower": None
        if lower is None
        else {
            "free_dofs": lower,
            "difference": lower - target,
            "difference_percent": 100.0 * (lower - target) / target,
        },
        "nearest_upper": None
        if upper is None
        else {
            "free_dofs": upper,
            "difference": upper - target,
            "difference_percent": 100.0 * (upper - target) / target,
        },
    }


def additional_iga_targets() -> tuple[dict, list[int], str]:
    path = HERE / "conforming_bending_study.json"
    data = json.loads(path.read_text())
    cases = {
        name: {
            "free_vector_dofs": row["vector_dofs_free"],
            "loaded_end_traction_weighted_mean_ux_mm": row["loaded_end_mean_ux_mm"],
            "field_degree": row["degrees_uvw"],
            "subdivisions_vw": row["subdivisions_vw"],
            "integration_cells": row["element_count"],
        }
        for name, row in data["cases"].items()
        if name.startswith("p322")
    }
    return (
        cases,
        [594, 996, 1152, 1926, 2268],
        hashlib.sha256(path.read_bytes()).hexdigest(),
    )


def matched_iga_target() -> tuple[dict, str]:
    path = HERE / "matched_iga.json"
    data = json.loads(path.read_text())
    return {
        "free_vector_dofs": data["vector_dofs_free"],
        "loaded_end_traction_weighted_mean_ux_mm": data["loaded_end_mean_ux_mm"],
        "field_degree": data["degrees_uvw"],
        "subdivisions_vw": data["subdivisions_vw"],
        "integration_cells": data["element_count"],
        "source_volume_mm3": data["geometry_volume_mm3"],
        "quadrature_order": data["gauss_order"],
        "relative_equilibrium_residual": data["reduced_residual"],
    }, hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--level", type=int, default=1)
    parser.add_argument("--rise-mm", type=int, default=3)
    args = parser.parse_args()
    if args.level < 1:
        parser.error("--level must be positive")

    cli = args.binary.resolve()
    if not cli.is_file():
        parser.error(f"native CLI not found: {cli}")
    study = load_study("s2_iga_study", HERE / "run.py")
    hexstudy = load_study("s2_hex_reference", HERE / "hex_reference.py")
    recipe_path = ROOT / "tests/references/g2_cross_sketches.recipe.json"
    model = model_from_cli(cli, args.rise_mm)
    targets, target_sha = enriched_targets()
    additional_targets, ten_chart_counts, additional_sha = additional_iga_targets()
    matched_target, matched_sha = matched_iga_target()

    # Reproduce the existing ten-chart IGA result at the same material and load.
    iga = study.run(model, 14, args.rise_mm)
    iga_counts = iga_dof_counts(study, model)
    source_nodes, source_cells, *_ = hexstudy.mesh(model, args.level, args.rise_mm)
    hex8_nodes, hex8_cells = source_nodes, source_cells
    hex20_nodes, hex20_cells = hexstudy.elevate_hexes(source_nodes, source_cells)
    if not np.array_equal(hex8_cells, hex20_cells[:, :8]):
        raise RuntimeError("HEX8 and HEX20 do not share the same element corner mesh")
    if not np.array_equal(hex8_nodes, hex20_nodes[: len(hex8_nodes)]):
        raise RuntimeError("HEX20 changed the HEX8 corner geometry")

    hex8 = brief_hex(
        hexstudy.solve_mesh(model, args.level, args.rise_mm, quadratic=False),
        len(hex8_nodes),
    )
    hex20 = brief_hex(
        hexstudy.solve_mesh(model, args.level, args.rise_mm, quadratic=True),
        len(hex20_nodes),
    )
    budget_search, budget_candidates = find_dof_brackets(
        model, hexstudy, targets, args.rise_mm
    )
    exact_options = {"circumferential": 3, "nv_by_col": [10, 1, 3], "height": 1}
    exact_x, exact_c, _, exact_traces, _ = hexstudy.mesh(
        model, 1, args.rise_mm, **exact_options
    )
    exact_topology = topology_counts(exact_c)
    exact_hex20_dofs = free_dofs(exact_x, exact_c, True, hexstudy)
    if exact_hex20_dofs["free_dofs"] != matched_target["free_vector_dofs"]:
        raise RuntimeError(
            "preselected HEX20 budget no longer matches the IGA free DOFs"
        )
    exact_candidate = {
        "options": exact_options,
        "free_dofs": exact_hex20_dofs["free_dofs"],
        "nodes": exact_hex20_dofs["nodes"],
        "total_dofs": exact_hex20_dofs["total_dofs"],
        "elements": len(exact_c),
        "partial_trace_quads": len(exact_traces),
        "topology": exact_topology,
    }
    exact_key = (
        tuple(exact_options["nv_by_col"]),
        exact_options["circumferential"],
        exact_options["height"],
    )
    budget_candidates["configs"][exact_key] = exact_candidate
    budget_search["p322_v6w2_C0"] = {
        "target": matched_target,
        "match_status": "exact",
        "selected": [exact_candidate],
        "selection_note": "configuration selected before displacement evaluation; integer layout 3, [10,1,3], 1 has 11,112 measured free HEX20 vector DOFs",
    }
    candidate_records = budget_candidates["candidate_records"]
    additional_scan = {
        name: {
            **target,
            **nearest_dof_scan(candidate_records, target["free_vector_dofs"]),
        }
        for name, target in additional_targets.items()
    }
    ten_chart_scan = {
        str(target): {
            "source": "parent-requested ten-chart IGA DOF target",
            **nearest_dof_scan(candidate_records, target),
        }
        for target in ten_chart_counts
    }
    additional_scan["p322_v6w2_C0"] = {
        **matched_target,
        "exact_match_found": True,
        "selected_hex20_divisions": exact_options,
        "selected_hex20_free_dofs": exact_hex20_dofs["free_dofs"],
        "dof_difference": 0,
    }
    budget_measurements = {}
    for candidate in budget_candidates["configs"].values():
        options = candidate["options"]
        X, C, _, traces, _ = hexstudy.mesh(model, 1, args.rise_mm, **options)
        quad_x, quad_c = hexstudy.elevate_hexes(X, C)
        topo = topology_counts(C)
        row = hexstudy.solve_mesh(
            model,
            options["height"],
            args.rise_mm,
            quadratic=True,
            mesh_options=options,
        )
        measured = brief_hex(row, len(quad_x))
        measured.update(
            {
                "division_parameters": options,
                "corner_elements": len(C),
                "partial_interface_quads": len(traces),
                "corner_topology": topo,
                "relative_volume_error_to_exact_union": measured["volume_mm3"]
                / float(targets["C0"]["source_volume_mm3"])
                - 1.0
                if "source_volume_mm3" in targets["C0"]
                else None,
            }
        )
        key = tuple(options["nv_by_col"]) + (
            options["circumferential"],
            options["height"],
        )
        budget_measurements[key] = measured
    for continuity, data in budget_search.items():
        target = data["target"]["free_vector_dofs"]
        measured_rows = []
        for selected in data["selected"]:
            options = selected["options"]
            key = tuple(options["nv_by_col"]) + (
                options["circumferential"],
                options["height"],
            )
            result_row = dict(budget_measurements[key])
            result_row["free_dof_difference"] = (
                result_row["free_displacement_dofs"] - target
            )
            result_row["free_dof_difference_percent"] = (
                100.0 * result_row["free_dof_difference"] / target
            )
            reference_ux = data["target"].get("loaded_end_traction_weighted_mean_ux_mm")
            result_row["mean_ux_difference_from_iga_mm"] = (
                None
                if reference_ux is None
                else (
                    result_row["loaded_end_traction_weighted_mean_ux_mm"] - reference_ux
                )
            )
            result_row["mean_ux_relative_difference_from_iga"] = (
                None
                if reference_ux is None
                else (
                    result_row["loaded_end_traction_weighted_mean_ux_mm"] / reference_ux
                    - 1.0
                )
            )
            measured_rows.append(result_row)
        data["measured_hex20_bracket"] = measured_rows
    force = float(iga["force_N"])
    result = {
        "schema": "s2_tilted_hex8_hex20_comparison_v1",
        "source": {
            "recipe": "tests/references/g2_cross_sketches.recipe.json",
            "recipe_sha256": hashlib.sha256(recipe_path.read_bytes()).hexdigest(),
            "native_binary": cli.name,
            "native_binary_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
            "family": "inclined_crossed_tools_v1",
        },
        "conditions": {
            "transverse_rise_mm": args.rise_mm,
            "young_modulus_MPa": 1000.0,
            "poisson_ratio": 0.3,
            "clamp": "x=0, all displacement components",
            "load": "1 MPa traction in global +x on x=100 end",
            "resultant_force_N": force,
            "observable": "traction-weighted mean global ux on loaded end",
            "sampling": "same traction-consistent end integral; force-weighted displacement",
        },
        "iga_ten_chart": {
            **iga_counts,
            "loaded_end_traction_weighted_mean_ux_mm": iga["loaded_end_mean_ux_mm"],
            "quadrature_order": 14,
            "relative_equilibrium_residual": iga["reduced_residual"],
            "continuity_jump_mm": iga["max_solved_jump_mm"],
        },
        "iga_enriched_22_chart_targets": {
            "source": "docs/results/crossed-interface-iga/smooth_bending.json",
            "source_sha256": target_sha,
            "geometry_charts": 22,
            "field_degree": [3, 2, 2],
            "field_spans_vw": [2, 2],
            "transverse_rise_mm": 3,
            "integration_cells": 88,
            "common_observable": "same traction-weighted loaded-end mean global ux",
            "same_material_loading_clamp": True,
            "targets": targets,
            "preselected_exact_budget_case": {
                "source": "docs/results/crossed-interface-iga/matched_iga.json",
                "source_sha256": matched_sha,
                "target": matched_target,
            },
        },
        "same_mesh_comparison": {
            "level": args.level,
            "same_corner_mesh": True,
            "element_corner_connectivity_identical": True,
            "hex8": hex8,
            "hex20": hex20,
            "hex20_minus_hex8_ux_mm": hex20["loaded_end_traction_weighted_mean_ux_mm"]
            - hex8["loaded_end_traction_weighted_mean_ux_mm"],
            "relative_hex20_minus_hex8_ux": hex20[
                "loaded_end_traction_weighted_mean_ux_mm"
            ]
            / hex8["loaded_end_traction_weighted_mean_ux_mm"]
            - 1.0,
            "iga_dof_match": False,
            "dof_match_note": "This baseline mesh is not DOF-matched; separate low-budget searches below bracket both enriched IGA targets.",
        },
        "equal_budget_search": {
            "parameters": {
                "circumferential_divisions": [1, 2],
                "wing_a_profile_divisions": [1, 13],
                "core_profile_divisions": [1, 3],
                "wing_b_profile_divisions": [1, 4],
                "height_divisions": 1,
            },
            "candidate_policy": "Integer chart-following meshes; retain only one connected cell complex, no overused faces, and complete paired partial-interface quadrilaterals. Actual HEX20 free DOFs are counted after the x=0 clamp.",
            "search_counts": {
                k: v
                for k, v in budget_candidates.items()
                if k not in ("configs", "candidate_records")
            },
            "targets_and_brackets": budget_search,
            "additional_conforming_iga_targets": {
                "source": "docs/results/crossed-interface-iga/conforming_bending_study.json",
                "source_sha256": additional_sha,
                "scan": additional_scan,
            },
            "additional_ten_chart_counts": ten_chart_scan,
        },
        "limits": [
            "The hexahedral boundary is chordal; it is not the exact rational spline boundary.",
            "HEX8 and HEX20 share the same corner mesh, but HEX20 adds midside nodes on its edges.",
            "The coarse equal-budget meshes use fewer chart-profile divisions and have chordal boundary-volume error; their response is an equal-DOF diagnostic, not a convergence qualification.",
            "The equal-budget search checks cell connectivity, face incidence and partial-interface pairing, but does not prove exact geometry conformity at every chart seam.",
            "These low-cost tilted cases do not establish three-dimensional convergence.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "hex8_free_dofs": hex8["free_displacement_dofs"],
                "hex20_free_dofs": hex20["free_displacement_dofs"],
                "iga_free_dofs": iga_counts["vector_free_after_clamp"],
                "iga_ux_mm": iga["loaded_end_mean_ux_mm"],
                "hex8_ux_mm": hex8["loaded_end_traction_weighted_mean_ux_mm"],
                "hex20_ux_mm": hex20["loaded_end_traction_weighted_mean_ux_mm"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
