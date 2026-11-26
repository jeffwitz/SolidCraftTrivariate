# Public bore propagation controls

Reference revision: `e5bd5c6`. Geometry changes are confined to the sketch partition and its existing CLI. The deterministic public study in `tests/references/propagated_bores.py` compares ordinary, stacked, diagonal, multi-level and axial bores on a filleted L, including a rotated arm assignment, an unequal-width control and two distinct stacks. All specified tools are physically disjoint and their bore centers are inside material. This is a development study, independent of the externally held private corpus.

For each identical recipe, record failure, fallback or accepted compact path. Give each CLI command 120 s. Read back the serialized model with the CLI audit and require a positive certified determinant bound. Verify full oriented homogeneous traces independently using exact rational arithmetic, require a one-to-one labeling of exterior faces, and reread STEP in FreeCAD/OCCT. Require one valid closed solid with one shell, adaptive volume within 2e-9 relative of the analytic volume, zero disagreements on 160 pseudorandom probes and 16 near-wall probes per bore, and finite sampled condition numbers on the closed 5 × 5 × 5 grid. Count every geometrically valid input, including refusals and timeouts. Invalid accepted geometry fails qualification. Compare counts and quality on common successes. There is no asserted Euler lower bound for these 3D complexes.

## Reproduction

Build the baseline with `cargo build --release --locked --offline` in a checkout of `e5bd5c6`, and build the candidate with the same command in its own checkout. Build the mass-properties helper using the command in the [orthogonal STEP readback](../../tests/references/orthogonal_step_freecad_result.md). Run the following in the candidate checkout with a Python interpreter that can import FreeCAD and Part. Each call writes the input recipes, native models, VTU samples, STEP, BREP and per-case scores to its own directory.

```sh
python3 tests/references/propagated_bores.py "$BASELINE_CLI" \
    build/solidcraft-bridge/occ-volume build/propagation-study/baseline
python3 tests/references/propagated_bores.py target/release/solidcraft-trivariate \
    build/solidcraft-bridge/occ-volume build/propagation-study/candidate
```

The baseline model in each admitted case can also be passed directly to the candidate's `audit` command to check serialized-model compatibility.
