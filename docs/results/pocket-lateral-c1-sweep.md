# Lateral C0/C1 sweep on the blind bilinear pocket

## Result scope

This is an independent sparse SciPy solver for a restricted anti-plane Ritz problem on the polygonal `sketch_layer_partition_v1` pocket. It is a consistent subproblem of 3D isotropic linear elasticity for `u=(0,0,w(x,y))`, constant through physical z. It is not unrestricted 3D elasticity and it is not the native Rust inherited-space solve. The published curved-pocket model and its prior mechanics numbers are separate: that model's maps are not bilinear.

The source is the exact 5+4 quadrilateral pocket partition rebuilt by the repository CLI from [the recipe](../../tests/references/blind_bilinear_pocket.recipe.json). The 5 lower charts define the section; 4 upper charts repeat the outer ring. All 9 maps have degree `[1,1,1]`, unit weights, and the CLI audit reports volume `46,800 mm^3` with a positive determinant lower bound of approximately `1,200 mm^3`. The reduction integrates the physical thickness exactly: 20 mm on the four outer section charts and 14 mm on the footprint chart. Thus the pocket remains blind, with a 14 mm floor and a 6 mm upper opening. Four wall tractions account for the removed upper material.

The field is `w=sin(pi*x/60) sin(pi*y/40)` mm, with `mu=1 MPa`; its body force is `-mu*Delta(w)`. Outer stock side faces are clamped, the 4 pocket walls receive manufactured natural traction, and top/floor/bottom have zero anti-plane traction. The tensor-product spline mesh uses n subdivisions in each parametric direction of each real chart, p=3 or 4, and knot multiplicity p for C0 or p-1 for C1. C0 edge coefficients are identified using the 8 actual oriented in-plane interfaces. C1 adds physical normal-gradient equations on every tangential knot span of all 8 interfaces; tangential continuity follows from the common trace. The equations cover the footprint edges, slanted radial seams, and junctions. Bilinear-map Jacobians are included in stiffness, loads, and norms.

The C1 coefficient equations use floating-point rank-revealing reduction and are checked at denser independent points. They are a numerical continuity qualification, not a symbolic certificate. The solver is sparse direct; the table gives actual free scalar DOFs after C1 constraints. No artificial three-component DOF count is reported.

## Error and convergence table

Errors are relative to independent high-order quadrature of the analytic MMS in the physical 3D energy and L2 norms, including the exact thickness weight. Rates are `log2(error at n / error at 2n)` and are diagnostics, not preregistered pass thresholds.

| p | n | Continuity | Free scalar DOFs | Relative energy error | Relative L2 error | Energy rate | L2 rate | q16/q24 |
|---:|---:|:---:|---:|---:|---:|---:|---:|:---:|
| 3 | 1 | C0 | 40 | 1.9868e-2 | 5.5390e-3 | — | — | pass |
| 3 | 2 | C0 | 169 | 2.3992e-3 | 3.2381e-4 | 3.050 | 4.096 | pass |
| 3 | 4 | C0 | 697 | 3.0016e-4 | 2.0419e-5 | 2.999 | 3.987 | pass |
| 3 | 8 | C0 | 2,833 | 3.7478e-5 | 1.2821e-6 | 3.002 | 3.993 | pass |
| 3 | 1 | C1 | 12 | 3.4344e-2 | 1.2998e-2 | — | — | pass |
| 3 | 2 | C1 | 64 | 7.1304e-3 | 1.5862e-3 | 2.268 | 3.035 | pass |
| 3 | 4 | C1 | 280 | 6.7198e-4 | 6.7840e-5 | 3.408 | 4.547 | pass |
| 3 | 8 | C1 | 1,192 | 7.4877e-5 | 3.6996e-6 | 3.166 | 4.197 | **fail** |
| 4 | 1 | C0 | 73 | 2.4635e-3 | 3.4801e-4 | — | — | pass |
| 4 | 2 | C0 | 305 | 1.6952e-4 | 1.5617e-5 | 3.861 | 4.478 | pass |
| 4 | 4 | C0 | 1,249 | 1.0804e-5 | 5.1385e-7 | 3.972 | 4.926 | pass |
| 4 | 8 | C0 | 5,057 | 6.8178e-7 | 1.6064e-8 | 3.986 | 4.999 | **fail** |
| 4 | 1 | C1 | 37 | 4.9745e-3 | 9.7186e-4 | — | — | pass |
| 4 | 2 | C1 | 164 | 4.1713e-4 | 5.2754e-5 | 3.576 | 4.203 | pass |
| 4 | 4 | C1 | 680 | 2.3450e-5 | 1.6126e-6 | 4.153 | 5.032 | pass |
| 4 | 8 | C1 | 2,792 | 1.2382e-6 | 4.3151e-8 | 4.243 | 5.224 | **fail** |

At fixed p and n, C1 uses fewer DOFs but is less accurate than C0 in every row of this study. For p=4,n=8 it reduces the free dimension from 5,057 to 2,792 (45%) while the energy error rises from `6.82e-7` to `1.24e-6` (1.82×). The table exposes error against actual free dimension for each refinement, so readers can compare the C0 and C1 sequences without conflating equal h with equal DOFs. This tradeoff does not establish that strong C1 is best for other loads or unrestricted 3D fields.

## Qualification and failed gates

The original observable threshold is unchanged: relative q-to-q change must be at most `1e-9`. The initial n=1 q12/q16 gate is stored in the machine-readable result. It passes for p=3/C0, p=3/C1 and p=4/C1; p=4/C0 fails at `1.9254e-9`. The full q16/q24 sweep passes 13/16 states. Its three failures are p=3/C1/n=8 (`5.06e-9`), p=4/C0/n=8 (`1.60e-9`), and p=4/C1/n=8 (`1.13e-7`). Their absolute error changes are small, but the preregistered relative criterion governs and the states remain failed.

The p=3/C1/n=8 case had two additional diagnostics, q24/q32 (`2.40e-9`) and q32/q40 (`8.44e-9`). Both failed; no further retries were accepted. The q16/q24 absolute changes in relative energy error and relative L2 error were `3.79e-13` and `4.69e-15`, respectively. Free equilibrium residuals were approximately `3.4e-14`, so the failed quadrature gate is not explained by an inadequate linear solve residual alone. Exact MMS norms changed by at most `1.02e-14`; smallest stiffness eigenvalues were positive (`0.0291` to `0.6114`); free residuals ranged from `9.3e-16` to `1.15e-13`; dense physical C1 jump checks were at most `6.89e-14`.

As a separate post-run load-sign and geometry control, the exact MMS was inserted into the reduced weak form. Body work plus the four pocket-wall tractions matches exact strain norm squared to `4.42e-15` at q24 and `4.02e-15` at q48. The values are 97.6403565038 N·mm for body work, 8.67293255882 N·mm for pocket-wall work, and 106.31328906265 N·mm for exact strain norm squared. This check validates the blind-pocket reduction and wall-traction orientation independently of the discrete solves.

Two earlier implementation failures are preserved separately in [the failure record](pocket-lateral-c1-sweep-initial-failures.json): C1 rows first missed tangential h-knots, and wall tractions were first integrated across spline breaks. Those measurements are excluded. The final implementation samples every tangential span for both constraints and wall loads.

## Reproduction

The CLI rebuilds the bilinear model from the versioned recipe; the result records its SHA-256 as `a1ba1d319dcc9c691534dfa888e5a155a76e481012f6aee79404ab6d2856614a`.

```sh
target/release/solidcraft-trivariate build \
  tests/references/blind_bilinear_pocket.recipe.json \
  build/pocket-lateral-c1/bilinear-pocket.model.json
python3 -m unittest discover -s tests/references -p test_pocket_lateral_c1_sweep.py
python3 tests/references/pocket_lateral_c1_sweep.py \
  --output build/pocket-lateral-c1/replay.json
```

The full command exits nonzero while the three q16/q24 states remain failed; the JSON is written before that exit. The additional q24/q32 and q32/q40 diagnostic is retained in [the retry record](pocket-lateral-c1-sweep-quadrature-retries.json). The run used Python 3.12.3, NumPy 2.4.4, and SciPy 1.17.1. Public source baseline: [fbe8723](https://github.com/jeffwitz/SolidCraftTrivariate/tree/fbe8723). Corresponding local checkout: `1c63c1280e1bef16fcf9a533b6fe70d6e0a8d158`.
