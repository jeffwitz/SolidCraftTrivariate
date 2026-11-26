# P322 analysis on the ten-map S2 crossed union

## Protocol and scope

This study keeps the ten native trivariate maps from the S2 union and their original partial wing-to-core interfaces. It changes only the analysis space: rational Bernstein degree `(3,2,2)` and real, C0 h-refinement in the chart parameters `v` and `w`. The source geometry is elevated homogeneously in `v,w`; subdivision uses de Casteljau restrictions. The mapped Jacobian determinant is checked at every volume quadrature point.

The independent model comes from `tests/references/g2_cross_sketches.recipe.json` and the production `crossed-boolean ... union` CLI for each transverse rise. The tested rises are `+3, 0, -3 mm`. Material is homogeneous isotropic, `E=1000 MPa`, `nu=0.3`. The left outer wing face is fixed in all three displacement components; the right outer wing face receives `1 MPa` axial traction in +x; all other exterior faces are traction-free. The x resultant is `40.2222222222 N`; geometry volume at q14 is `6799.60472580227 mm³`.

The volume basis has 4×3×3=36 scalar coefficients per element. C0 identities are made across real `v,w` h-interfaces and the seven full chart neighbours; the two collapsed polar lines are single-valued. The four partial wing/core faces use enriched trace collocation. At each sampled u-line, the implementation splits at the wing's h-knots and at every preimage of a master `w` knot. It partitions u at crossings of each slave w-knot boundary with master knot planes, then uses 21 u nodes per interval (greater than the cleared-denominator bound `3+3*pw+6=15`) and four w nodes per resulting segment (`pw+1=3` minimum). The partial constraint null space uses SVD with `rcond=1e-12`; the JSON also records the rank under `rcond=1e-10` as a stability check. This is a numerical sampled coupling, not a formal exact mortar certificate. Independent off-grid traces and the physical face maps are tested separately.

## Results

| rise (mm) | v intervals | w intervals | free vector DOFs | trace rank | mean loaded-end ux (mm) | energy (N mm) |
|---:|---:|---:|---:|---:|---:|---:|
| +3 | 1 | 1 | 594 | 66 | 0.2780814182 | 5.59252630 |
| +3 | 2 | 1 | 1152 | 66 | 0.2943430322 | 5.91956543 |
| +3 | 4 | 1 | 2268 | 66 | 0.3012717741 | 6.05891011 |
| +3 | 1 | 2 | 996 | 108 | 0.2784237913 | 5.59941180 |
| +3 | 2 | 2 | 1926 | 108 | 0.2948563732 | 5.92988928 |
| 0 | 1 | 1 | 666 | 42 | 0.0845637513 | 1.70067100 |
| 0 | 2 | 1 | 1224 | 42 | 0.0862122229 | 1.73382359 |
| 0 | 4 | 1 | 2340 | 42 | 0.0872109396 | 1.75390890 |
| −3 | 1 | 1 | 594 | 66 | 0.2780814175 | 5.59252628 |
| −3 | 2 | 1 | 1152 | 66 | 0.2943430319 | 5.91956542 |
| −3 | 4 | 1 | 2268 | 66 | 0.3012717741 | 6.05891012 |

The `+3` ten-map `v2,w2` result (`0.29485637 mm`) is within 0.82% of the independently assembled 22-map P322 `v2,w2` matching-face result (`0.2972940 mm`, `4128` free vector DOFs, `5.978913 N mm`). That comparison uses the same native geometry, material, loading, and clamp. The existing coarser P311 ten-map result is `0.08004945 mm` at either tilt; P322 changes compliance materially, so the old low-degree result is not a reliable stiffness estimate for the tilted specimen. The zero-rise control stays near `0.0872 mm` after v refinement, which separates the tilt response from the straight geometry case.

Numerical checks across the listed cases: six free-body rigid modes; affine uniaxial patch energy within `1.6e-12` relative; partial-face geometry mismatch below `5.7e-14 mm`; independent held-out unit-space trace defect below `2.4e-15`; solved partial-face displacement jump below `1.8e-14 mm`; reduced equilibrium residual below `4.8e-10`. q8/q14 for `+3,v1,w1` change mean ux by `2.19e-9 mm`; for `+3,v4,w1` the change is `7.83e-9 mm`. The reported trace defects are collocation diagnostics and do not convert the sampled coupling into an exact proof.

## Reproduce

From the repository root:

```sh
cargo build --locked --offline
OPENBLAS_NUM_THREADS=1 python3 docs/results/crossed-interface-iga/partial_bending_study.py \
  target/debug/solidcraft-trivariate \
  docs/results/crossed-interface-iga/partial_bending_study.json \
  --plan docs/results/crossed-interface-iga/partial_bending_plan.json
```

`partial_bending_plan.json` freezes the 13 runs, including the q14 checks. The script invokes the native CLI to regenerate every ten-map model; it does not load a manually edited approximation. `partial_bending_study.json` contains DOF and rank accounting, held-out trace probes, interface geometry mismatch, residuals, force, volume, displacement, and energy for each run.
