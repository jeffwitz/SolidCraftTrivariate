# Tall-cylinder Boolean operations

Baseline: main `19375aa`. A is the original 40 × 40 × 12 mm stock; B is a vertical cylinder of radius 5 mm centered at (20,20), from z = −12 to z = 24 mm. This experiment covers the two operands, their regularized union, intersection, both directed differences and symmetric difference.

| State | Expected volume (mm³) | Face-connected components |
|---|---:|---:|
| A | 19,200 | 1 |
| B | 900π | 1 |
| Union | 19,200 + 600π | 1 |
| Intersection | 300π | 1 |
| A minus B | 19,200 − 300π | 1 |
| B minus A | 600π | 2 |
| Symmetric difference | 19,200 + 300π | 3 |

Acceptance criteria: reproduce the seven native models, retain the classical four maps of A minus B, verify coverage against the two primitive membership predicates above, within and below the stock, and check volume from the exported rational control nets to relative tolerance 1e−9. Check that each NURBS representation extracts to its source Bézier cells. The two cut disks in B minus A are disconnected. For the symmetric difference the three material interiors are separate; their closures touch along the two circular rims. It must be displayed with that qualification.

The suite extends the box/cylinder constructor, the exact gallery patch maps and the existing box/cylinder page. Tangency, a blind tool, and inclined axes are outside this experiment.

## Measured representations

| State | Bézier cells | NURBS volumes | Shared Bézier faces |
|---|---:|---:|---:|
| A | 1 | 1 | 0 |
| B | 15 | 2 (closed ring) or 3 (half-rings) | 34 |
| Union | 19 | 7 | 42 |
| Intersection | 5 | 3 | 8 |
| A minus B | 4 | 2 | 4 |
| B minus A | 10 | 6 | 16 |
| Symmetric difference | 14 | 8 | 20 |

Each disk section has a regular square core and four surrounding cells. The disk sectors reuse the homogeneous circular traces of the classical bore. In the union, five axial NURBS volumes continue through all three height intervals. The remaining two volumes contain the four original bore sectors. Every source Bézier cell is recovered by extraction from its assigned NURBS volume. These counts describe the published construction.

The cylinder supports a joint angular and axial fusion. The two-volume cover identifies both angular end traces of the annulus. The three-volume cover has two half-rings. [Geometry, seam topology and IGA qualification](closed-ring-iga.md)

## Reproduce the gallery

Build the Rust CLI with `cargo build --release --locked --offline`. For each state, the recipe, native model, native audit and exact NURBS patch map are under `docs/site/box-cylinder/tall/`. For example:

```sh
site=docs/site/box-cylinder/tall
cli=target/release/solidcraft-trivariate
mkdir -p build
"$cli" build "$site/union.recipe.json" build/tall-union.json
"$cli" audit build/tall-union.json
cargo build --release --locked --offline --example cylinder_tensor_covers
python3 scripts/build_box_cylinder_tall_patches.py --check
python3 scripts/check_tall_booleans.py
```

The independent check integrates the oriented cross-section of each serialized straight extrusion from its homogeneous boundary curves and multiplies by its height. It compares occupancy and overlap multiplicity against the Boolean truth table at 3,179 off-boundary probes per state, including points 0.01 mm to either side of the cylinder and both stock caps. The [manifest](../../site/box-cylinder/tall/manifest.json) records the expected volumes, observed cell and patch counts, and face-connected components. The [interactive page](https://jeffwitz.github.io/SolidCraftTrivariate/box-cylinder.html#tall-cylinder) can show a vertical section through the center.

## Construction diagnostic

The first disk candidate combined a transformed primitive sector with the bore's angular control row. It failed the positive-Jacobian gate (`boolean_chart_jacobian_not_certified`) on the cylinder. The box alone passed. The failure was retained as a construction diagnostic; the gate and the volume tolerance were kept fixed.
