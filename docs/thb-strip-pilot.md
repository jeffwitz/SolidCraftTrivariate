# Exact interior THB strip pilot

## Scope and acceptance

This pilot starts from the audited cubic NURBS strip parent and keeps its
geometry map unchanged. It constructs a tricubic field space in three recorded
stages: the source parent, an exactly degree-elevated initial space, a uniformly
prepared space, and one strictly interior dyadic THB refinement. It does not
construct a T-mesh or claim support for arbitrary multipatch bodies.

The default prepared-cell mark is half-open `[[2,6],[1,3],[1,3]]` in the
prepared span grid, with tensor storage `(u * n_v + v) * n_w + w` (`w` fastest).
The builder rejects boundary-touching/empty marks, marks with no fully supported
fine function, noncubic source profiles, reversed or uncertified geometry
orientation, and inputs exceeding checked allocation budgets. A separate
alternative interior mark is covered by the Rust tests. Source-mesh admission
uses the same exact `initial_mesh_admission` contract as the canonical parent
eligibility report; this pilot does not redefine that decision.

## Exact construction

Knot insertion and degree elevation are performed with rational arithmetic.
The support and THB truncation are formed from exact tensor extraction rows.
Every local leaf verifies Bernstein partition of unity. Every prepared-to-THB
nested-space identity is checked row by row, and the prepared-to-local
coefficient transfer is checked for full column rank. All four homogeneous
geometry coordinates are also compared exactly on every leaf: direct parent
restriction, prepared tensor extraction, and THB extraction must agree.
The 544 local basis functions also receive an exact rank check when expanded in
the 2,299-function global fine tensor basis.

The source extrusion orientation is certified on every source knot span. Its
signed Jacobian numerator is represented in Bernstein form and recursively
subdivided with exact midpoint de Casteljau arithmetic until every coefficient
is strictly positive. Failure to certify within the depth budget is an explicit
refusal. Runtime sampling uses cached f64 extraction rows per leaf; rational
operations are not performed at quadrature points.

## Prepared-space comparison

Each leaf exposes both `sample` (local THB field) and `sample_prepared` (the
untruncated prepared field) on the identical geometry chart. This permits the
mechanical solver to assemble both operators on the same local partition and
quadrature points. The prepared zero-bottom lift uses the documented w-fast
index order and constrains all `w=0` tensor controls.

The JSON document emitted by the CLI records exact knot vectors, active IDs,
half-open marked cell ranges, polynomial and rational transfers, global field
weights, local extracted rows and weights, and homogeneous geometry controls.
Rational numbers are serialized as decimal-string numerator/denominator pairs.
Regenerate the package and exact transfer metadata from the retained model and
plan with:

```sh
cargo run --release --locked --offline -- refine-thb-strip \
  build/thb-night/model.json tests/references/thb_interior.plan.json \
  build/thb-night/final-package.json
cargo run --release --locked --offline -- export-thb-extraction \
  build/thb-night/final-package.json build/thb-night/final-extraction.json
```

## Reference fixture

For `tests/references/nurbs_strip.recipe.json`, the measured cardinalities are:

| Stage | Scalar DOFs / elements |
|---|---:|
| Source parent | 20 |
| Degree-elevated initial tricubic space | 80 |
| Uniform prepared space (8×4×4 spans) | 539 |
| Uniform global fine comparison space | 2,299 |
| Local THB space (539 coarse + 5 fine functions) | 544 |
| Extracted leaves (112 coarse + 128 refined) | 240 |

The exact prepared-to-THB transfer has rank 539, equal to the full prepared
column rank; the local functions have rank 544 when embedded in the 2,299
function uniform fine basis. This records strict nesting while retaining five
additional local fine modes. The Python `Fraction` oracle independently reconstructs basis
values and gradients from the knot vectors and exact extraction data. Runtime
mechanical residuals and energy comparisons are separate evidence, not implied
by the exact algebraic checks.
