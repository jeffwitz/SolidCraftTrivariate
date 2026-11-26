"""Replay the directional bending controls and representative enriched trace spaces."""

import argparse
import json
from pathlib import Path

import numpy as np

import conforming_bending_study as fitted
from hex_bending_reference import model_from_cli
import hex_reference
import partial_bending_study as partial
import smooth_bending


def read(name):
    return json.loads(Path(__file__).with_name(name).read_bytes())


def verify(work, force, volume, residual):
    assert np.isfinite(work) and work > 0
    np.testing.assert_allclose(force, 40.22222222222222, rtol=2e-8)
    np.testing.assert_allclose(volume, 6799.6047258023, rtol=2e-8)
    assert np.isfinite(residual) and residual < 2e-8


def main(binary):
    model = model_from_cli(binary.resolve(), 3)
    assert len(model["blocks"]) == 10
    base, locations = fitted.fitted_geometry(model, 3)
    assert len(base) == 22
    expected = read("conforming_bending_study.json")["cases"]
    for label in ("p311", "p321", "p312", "p322"):
        old = expected[label]
        row = fitted.solve_case(base, locations, old["degrees_uvw"], 1, 1, 8)
        verify(
            row["strain_energy_N_mm"],
            row["loaded_force_x_N"],
            row["geometry_volume_mm3"],
            row["reduced_residual"],
        )
        assert row["vector_dofs_free"] == old["vector_dofs_free"]
        np.testing.assert_allclose(
            row["loaded_end_mean_ux_mm"], old["loaded_end_mean_ux_mm"], rtol=2e-7
        )
    measured = []
    for old in read("smooth_bending.json")["runs"]:
        nv, nw = old["spans_vw"]
        run = smooth_bending.solve(binary, 3, nv, nw, 8)
        for row, ref in zip(run["fields"], old["fields"], strict=True):
            verify(
                row["energy_N_mm"],
                run["force_N"],
                run["volume_mm3"],
                row["relative_residual"],
            )
            assert (row["continuity"], row["free_vector_dofs"]) == (
                ref["continuity"],
                ref["free_vector_dofs"],
            )
            np.testing.assert_allclose(row["mean_ux_mm"], ref["mean_ux_mm"], rtol=2e-7)
        measured.append(run)
    assert (
        measured[0]["fields"][0]["free_vector_dofs"]
        == measured[1]["fields"][1]["free_vector_dofs"]
    )
    saved = read("partial_bending_study.json")
    for nv, nw in [(1, 1), (2, 2)]:
        old = next(
            r
            for r in saved
            if (r["rise_mm"], r["refine_v"], r["refine_w"]) == (3, nv, nw)
            and r["quadrature_order"] == 8
        )
        row = partial.solve_case(model, 3, nv, nw, 8)
        verify(
            row["strain_energy_Nmm"],
            row["force_N"],
            row["volume_mm3"],
            row["reduced_equilibrium_residual"],
        )
        assert row["free_vector_dofs"] == old["free_vector_dofs"]
        assert row["rigid_modes"] == 6
        assert row["affine_energy_error"] < 1e-8
        assert row["max_partial_face_geometry_mismatch_mm"] < 1e-8
        assert row["max_solved_partial_interface_jump_mm"] < 1e-8
        assert row["trace_operator_heldout_defect"] < 1e-8
        np.testing.assert_allclose(
            row["loaded_end_mean_ux_mm"], old["loaded_end_mean_ux_mm"], rtol=2e-7
        )
    recorded = read("matched_iga.json")
    matched = fitted.solve_case(base, locations, [3, 2, 2], 6, 2, 8)
    mesh = hex_reference.solve_mesh(
        model,
        1,
        3,
        quadratic=True,
        mesh_options={"circumferential": 3, "nv_by_col": (10, 1, 3), "height": 1},
    )
    assert matched["vector_dofs_free"] == mesh["free_dofs"] == 11112
    verify(
        matched["strain_energy_N_mm"],
        matched["loaded_force_x_N"],
        matched["geometry_volume_mm3"],
        matched["reduced_residual"],
    )
    np.testing.assert_allclose(
        matched["loaded_end_mean_ux_mm"], recorded["loaded_end_mean_ux_mm"], rtol=2e-7
    )
    np.testing.assert_allclose(mesh["force_N"], matched["loaded_force_x_N"], rtol=2e-8)
    hex_report = read("hex_bending_reference.json")["equal_budget_search"]
    exact = hex_report["targets_and_brackets"]["p322_v6w2_C0"]
    reference = exact["measured_hex20_bracket"][0]
    assert reference["free_displacement_dofs"] == mesh["free_dofs"]
    np.testing.assert_allclose(
        mesh["ux_mm"], reference["loaded_end_traction_weighted_mean_ux_mm"], rtol=2e-7
    )
    assert mesh["min_gauss_det"] > 0
    print(
        "S2 directional enrichment, strong partial traces and equal-budget C0/C1 checks passed."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    main(parser.parse_args().binary)
