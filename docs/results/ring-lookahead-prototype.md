# Bounded planar ring lookahead prototype

## Pre-registration

- Starting revision: `72e8e19c1bd4b25b41610a11d930f38543e40d40` from the isolated `main` checkout.
- Scope: add an opt-in candidate selector for convex sketch ring partitions; no default planner or recipe-format change.
- Frozen inputs: `seeded_convex_14`, `seeded_convex_04`, and `seeded_convex_26` from `docs/results/sketch-planner-polygon-evaluation/cases.json` (generator seed `20261126`).
- Candidate generation: compare the existing bounded host-separator alternatives and ring dynamic program globally rather than selecting tied-count alternatives by the extruded 3D score. The lookahead variant adds support offsets 1/4 and 3/4 (in addition to 1/3, 1/2, and 2/3), a host-search budget of 128 recursive visits and at most eight host plans. It also tries at most 128 exact shared-junction slides per recipe: for a non-input/non-guide graph vertex lying strictly between collinear neighbors, move it toward either neighbor by 1/16, 1/8, 1/4, or 1/2 of that segment. Slides are rejected before construction if they leave any source segment containing the original vertex.
- Primary objective: minimize material block count. Secondary diagnostic: minimize the maximum 2D Jacobian condition number, sampled on each cell's closed 3×3 parameter grid. Deterministic first-candidate tie break within relative `1e-12`.
- Required gates: `partition_from_graph` exact rational area/volume and complete incidence checks; positive determinant lower bound on each emitted block; no input refusals. No candidate may replace a lower-block-count candidate.
- Baseline: default `partition` at the same starting revision. Report baseline and prototype block count, method, planar score, certificate lower bound, and elapsed time. The planar score is a diagnostic, not a continuous certificate.
- Exit: prototype API remains explicitly opt-in; report all three inputs, including no-change or worse outcomes, and run targeted Rust checks plus `make check` if the prototype is retained.

## Results

Frozen source: `f406011839fdf6dd92f3d7857c53b2b35c90231e`. Release CLI SHA-256: `e2ac9c06f1148c914f53c13ad05e4f71c6900847c758647bcd7e0470cbd55bfc`. Build command: `cargo build --release --locked --offline` with Rust `1.99.1` on Linux x86_64, Intel Core Ultra 5 228V. The isolated tracked worktree was clean at the source commit.

The paired public inputs are the frozen cases `seeded_convex_04`, `seeded_convex_14`, and `seeded_convex_26` from `cases.json` (SHA-256 `f10964bb0417a702a394f66aa96e148b6ead67685a3a9f9f541cea8d6ab98e7b`), plus the specifically authorized uploaded `mixed_sizes.recipe.json` (SHA-256 `9464ebe901b9b3ddb22bba0991760a65e0ad0b0706af9523fae9b57076eec585`). The prior `mixed_sizes.model.json` is SHA-256 `949649c718a79485a1e1147737e088253d46bea07b355db727a7a4556731397c`.

All four inputs built and audited in all three modes. The independent exact section oracle, incidence/Euler checks, area coverage and volume checks passed in every row. Each mode produced the same block count for a given input. Default output was byte-identical to the clean `72e8e19` baseline model on all four inputs, including the uploaded mixed-size reference.

| Case | Blocks | Automatic planar κ | Rerank planar κ | Lookahead planar κ | Automatic → lookahead det. lower bound (mm³) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `seeded_convex_04` | 12 | 123.9416 | 41.1111 | 31.6034 | 180960 → 113100 |
| `seeded_convex_14` | 30 | 229.3883 | 229.3883 | 31.0621 | 5544.609 → 29311.868 |
| `seeded_convex_26` | 33 | 126.5424 | 75.9425 | 75.9425 | 22986.389 → 12131.667 |
| `mixed_sizes` | 14 | 558.9996 | 558.9996 | 45.4834 | 18.743 → 507.962 |

The planar κ values are the independent `planar_spectral_kappa` summary computed by `/tmp/sct-quality-review/scripts/evaluate_polygon_extrusion_quality.py` (SHA-256 `018a75fcf2c379012e25afee2678d858105cfe29b845d3812fb630fe0fddec14`). Lower is better. It is a sampled diagnostic; it does not replace the positive closed-domain determinant certificate. The clean-baseline/default byte comparisons were: `_04` `af61158d…404ce`, `_14` `7f1d7e1a…60d3d`, `_26` `93eaac7a…4f083`, and `mixed_sizes` `949649c7…1397c`, identical between the clean baseline binary and the prototype's automatic mode.

Single warm-cache build/audit wall times (seconds; build and audit measured separately):

| Case | Automatic build/audit | Rerank build/audit | Lookahead build/audit |
| --- | ---: | ---: | ---: |
| `seeded_convex_04` | 0.205 / 0.122 | 0.354 / 0.160 | 0.359 / 0.179 |
| `seeded_convex_14` | 11.845 / 5.994 | 25.269 / 12.893 | 25.710 / 12.618 |
| `seeded_convex_26` | 6.697 / 3.337 | 13.291 / 6.661 | 14.787 / 7.375 |
| `mixed_sizes` | 0.923 / 0.444 | 1.750 / 0.860 | 1.809 / 0.899 |
| **Total** | **19.670 / 9.897** | **40.664 / 20.574** | **42.665 / 21.071** |

The timing table is one run per cell after the release build; it is descriptive, not a controlled performance estimate. The bounded lookahead roughly doubles build-plus-audit time on this four-case targeted sample. It improved planar κ on all four cases; reranking alone improved two and had no change on `_14` or `mixed_sizes`. The expanded search added a further improvement on `_04`, `_14`, and `mixed_sizes`, and matched reranking on `_26`.

Reproduction used the ordinary CLI with the same recipe and model paths for every mode:

```sh
cargo build --release --locked --offline
solidcraft-trivariate build CASE.recipe.json CASE.model.json
solidcraft-trivariate audit CASE.model.json
```

For the ablations, set the recipe's top-level `partition_strategy` to `ring_planar_rerank` or `ring_planar_lookahead`. Omit the field for automatic baseline. The independent exact checks were run by `docs/results/sketch-planner-polygon-evaluation/run.py` (`eval_one`; SHA-256 `fe8771b7d3745fba5ff07d15390771410a11e0371567f267196a160c2b044ed4`). The focused `cargo test --locked --offline --test sketch_generic_separators` passed 5 tests. The frozen isolated source had one strict Clippy `type_complexity` finding in the candidate callback signature; that type-only cleanup was handled separately in the integration branch and was not part of the measured binary.
