# Exact spline parent support: S1

## Scope and admission

`nurbs_sketch_strip_v1` retains two synchronized rational sketch curves and
their source knot intervals. `spline_support::recover_nurbs_strip_parent`
reconstructs their complete tensor-product volume with degrees `(p,1,1)` and
the original source knots. It admits a model only after rebuilding it from its
embedded recipe and comparing every Bézier block with exact knot insertion in
the original homogeneous gauge, up to one positive constant block gauge.

The mandatory cubic recipe is
`tests/references/nurbs_strip.recipe.json`. It has degree `(3,1,1)`, control
shape `5 × 2 × 2`, 20 parent controls, two source spans `[-2,1]` and `[1,4]`,
32 block-control occurrences, and 28 distinct controls after the C0 interface
plane is shared. The interior knot `1` is simple, so the source support carries
C2 continuity in the actual source parameter. The support keeps the actual
intervals; it does not replace them with unit spans.

The representation distinguishes the parent, each source span, its block,
the affine restriction, the homogeneous gauge, and continuity evidence.
Restrictions include axis permutations and reversal flags. `map_parameters`
and `derivative_scale` expose the affine map and its signed metric factor.
Serialization uses decimal-string numerator and denominator pairs for exact
homogeneous values. The stable support ID is tied to a deterministic FNV-1a
recipe fingerprint; this token is for lineage and is not a cryptographic hash.

## Independent checks

Rust tests verify the 20-control parent, exact lineage, serialization,
orientation/permutation, changed-recipe revision, forged-block rejection,
source reversal, and nonuniform spans `[-5,0]` and `[0,4]` with derivative
factors 5 and 4.

The independent Python oracle evaluates the retained source curves directly
with `fractions.Fraction` and Cox-de Boor basis functions. It compares each
parent boundary against the serialized block on `2p+1` distinct exact
parameters per span. The cross-multiplied coordinate identity has degree at
most `2p`; those exact equalities therefore establish the whole-span rational
map identity. It separately checks first derivatives with the signed source
interval length. For the mandatory recipe this yields 56 exact position
comparisons and 16 derivative comparisons. The Python comparison reads the native build and exported parent, checks the stored weights, intervals, gauges and orientation against the recipe, and computes its own reference values.

Run the frozen cubic and nonuniform checks with:

```sh
python3 scripts/check_spline_parent.py --cli target/release/solidcraft-trivariate
python3 scripts/check_spline_parent.py --cli target/release/solidcraft-trivariate --nonuniform-source
python3 scripts/check_spline_parent.py --cli target/release/solidcraft-trivariate --reverse-source
cargo run --locked -- spline-eligibility build/model.json
cargo run --locked -- recover-spline-parent build/model.json build/parent.json
cargo test --locked --offline --lib spline_support
```

## Conversion eligibility

The first `structured_tricubic_m2_v1` profile requires one regular, injective
tensor parent covering the complete material. A valid model outside that
profile remains valid. `eligibility_report` keeps every exact shared-face
incidence and labels its established continuity as C0; it does not infer C1 or
C2 from face coincidence.

The box-minus-through-cylinder fixture has four sector blocks and four retained
joins. Its positive-radius cutter strictly crosses both stock end planes, so
its material has a through tunnel and cannot be the injective image of one
tensor cube. The exact NURBS cross-bore fixture has 26 blocks, 36 retained
joins, and four retained source curves. Those curves are not a complete 3D
parent; the report leaves each join at C0 unless stronger source evidence is
available. Both cases are valid source geometry and return
`unsupported_profile` for conversion.

## Regional source columns

`recover_translation_columns` proves a prismatic column on a block when one
local direction has degree one and every corresponding homogeneous endpoint
control satisfies `H1 - H0 = t W` with equal positive weights. The translation
vector `t` may point in any physical direction. The proof is exact over the
whole parametric block; a failed block is reported locally and does not affect
validity of the source body.

On `box_minus_cylinder.recipe.json`, all four annular sector blocks recover
their exact base source columns and common translation `(0,0,12) mm`. On
`sketch_boolean_blind_pocket.recipe.json`, all nine blocks recover; five use
`(0,0,6) mm` and four use `(0,0,4) mm`. On `nurbs_cross_bores.recipe.json`, 18
of 26 blocks satisfy the selected z-column identity; eight curved blocks are
explicitly retained as non-prismatic. These are regional chart results, not
one parent covering the Boolean body. `scripts/check_spline_eligibility.py`
recomputes the control-column identities with Python `Fraction` and separately
checks the exact through-cylinder preconditions.

The multi-support results are read-only eligibility reports. They do not rename
blocks, change face selections, or replace the accepted source model.

A single audited affine parallelepiped is also recoverable. The compiler checks every homogeneous control, constant positive weight and a strictly positive oriented determinant. For example, the box operand in the tall Boolean suite has one complete regular parent. Recovery and admission to a uniform tricubic starting mesh are reported separately. Nonuniform source spans can retain a valid exact parent. A through bore supplies a topological obstruction to one globally regular injective cube chart. For other families with no successful witness the status remains unsupported by this recovery procedure.
