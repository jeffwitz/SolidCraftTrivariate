# Box minus centered cylinder reference

This control case uses the stock `[0,40] × [0,40] × [0,12]` mm and a
z-directed cylinder centered at `(20,20)` mm, radius `5` mm, with tool interval
`[-2,14]` mm. The Boolean recipe describes the two primitive operands separately. For this through-cut, the native constructor reuses the four-sector sweep. The control recipe describes the same solid using a rectangular sketch and a bore feature.

From the repository root, replay both builds, their native audits, and the
independent exact-trace/interface and numerical-volume checks with:

```sh
cargo build --offline --locked
python3 scripts/check_box_cylinder.py --cli target/debug/solidcraft-trivariate
```

The checker requires both versioned recipes. It rebuilds models into `build/`,
checks the exact rational circle identity on each wall trace, checks common
homogeneous control nets on all interfaces, and compares the two sets of
control maps after ignoring block labels. The Jacobian integral uses a 20-point
Gauss-Legendre rule on each parametric axis and is compared with
`19200 - 300*pi` mm³. Exact rational trace checks do not depend on sampling; the
quadrature value is an independent numerical cross-check, not a proof of
volume.

`classic-independent-audit.json` records the model and recipe SHA-256 values,
the exact checks, the quadrature result, and the route-to-route map comparison.


[Inspect the two grids side by side](https://jeffwitz.github.io/SolidCraftTrivariate/box-cylinder.html). At half height, the four degree-(2,1,1) maps each connect one straight side of the square to a rational quarter-circle. Their seams meet the bore at rational 3–4–5 offsets. The colors in the section plot identify curves of constant u and constant v.

The [Boolean recipe](../../site/box-cylinder/boolean.recipe.json), [classic sketch recipe](../../site/box-cylinder/classic.recipe.json), [geometry audit](../../site/box-cylinder/audit.json), [STEP](../../site/box-cylinder/result.step) and [CAD readback](../../site/box-cylinder/step-audit.json) accompany the viewer. To regenerate the figure, run `python3 scripts/plot_box_cylinder_section.py` with NumPy and Matplotlib. To repeat the CAD readback, build the [OCCT mass helper](../../../tests/references/orthogonal_step_freecad_result.md), set `OCC_VOLUME_BIN` to its executable and run `freecadcmd tests/references/box_cylinder_step.py`.
