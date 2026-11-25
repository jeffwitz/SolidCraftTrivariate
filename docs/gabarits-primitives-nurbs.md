# From sketch primitives to NURBS volumes

The generator chooses construction lines when the primitive is created. It preserves boundary spans, planned junctions, and the homogeneous net for each section patch. The operation then chooses the sweep direction and spans. The counts below describe patches emitted by these library conventions. A single NURBS sketch curve may contain multiple spans.

| Filled primitive | Boundary spans | Geometric junctions | Section patches | Blocks for straight extrusion | Blocks for full turn about Y |
|---|---:|---|---:|---:|---:|
| Rectangle | 4 degree-1 lines | 4 G0 corners | 1 | 1 | 4 |
| Disk | 4 rational degree-2 quarters | 4 G2-or-higher joins | 5 | 5 | 20 |
| Filled ellipse | 4 rational degree-2 quarters | 4 G2-or-higher joins | 5 | 5 | 20 |

The number of boundary NURBS surfaces differs from the number of volume blocks. The extruded rectangle has 6 exposed block faces. The extruded disk and ellipse each have 14: 4 side walls and 5 patches on each of the two caps. For an off-axis revolution, all three templates have 16 exposed faces, one for each boundary-span × quarter-turn pair. Thus the torus has 16 boundary faces for 20 material blocks.

![Section patches extracted from the controls produced by the three templates](figures/gabarits-primitives-nurbs.png)

The figure samples the boundary of each block on the bottom face of the exported extrusion. The labels count blocks in the corresponding extrusion and revolution models. To regenerate it from the CLI outputs:

```sh
cargo build --release --locked --offline
mkdir -p build/templates
python3 tests/references/primitive_sweep_cli.py target/release/solidcraft-trivariate build/templates/models build/templates/checks.json
python3 tests/references/plot_primitive_templates.py build/templates/models build/templates/sections.png
```

For the disk, lines connect the four vertices of a diamond at half-radius to the four cardinal points on the circle. The diamond gives a regular central patch; four side regions connect its sides to the boundary. Thus the center has a nondegenerate surface patch. The ellipse uses an affine image of the same subdivision. For the rectangle, the two side parameters directly define the interior patch.

Each quarter circle comes from the exact rational parameterization `(1−t², 2t)/(1+t²)`, its rotations, and direction changes needed to match weights at joins. These coefficients have an integer homogeneous representation. The constructor knows the junctions before producing volume controls. Tests independently check jets on emitted coefficients, incidences, and determinant sign on every closed cube.

On the unit disk, each side patch has one side of the diamond, two radial segments, and one quarter-circle as its boundary. This boundary is simple and traverses the corresponding sector boundary exactly once. With a strictly positive determinant, the topological degree in the sector is one: every interior point has a preimage. The four sectors and diamond have disjoint interiors. A translation and two positive dilations carry this partition to the ellipse. Positive extrusion length or a strictly positive revolution radius preserves regularity in each angular sector.

For extrusion, a section patch extends across the full positive length. For revolution, each patch is crossed with four quarter-turns. The sketch axis is Y and the radial coordinate must remain strictly positive. A disk offset from the axis therefore gives a solid torus. Crossing or touching the axis needs another template because angular coordinates degenerate on the axis.

The table counts apply to an isolated filled primitive. Removing the primitive from a part also requires maps for the **material around the hole**. For two intersecting removals, lines induced in each sketch depend on the relative profile positions and extrusion directions. The template preserves its original junctions, then the preplanner adds required intersections. [Sketch prepartitioning](prepartition-esquisse-nurbs.md)

The Rust API is `Primitive::plan(operation)`. The plan provides `contour`, `junctions`, `section_patches`, and `block_count()` before any 3D construction. `plan.build()` instantiates the volume. Sketch dimensions are capacity-checked rationals. [Code](../src/primitive_template.rs) · [Tests](../tests/primitive_template.rs)

For the 3 mm extrusion or revolution radius 6 mm from the centroid, all six volumes integrated at Gauss order 12 agree with area × length or Pappus' theorem to at most 5.66 × 10⁻¹⁵ relative. The Bernstein certificate establishes a strictly positive determinant on all 55 blocks across the six volumes. Measured incidences are 0 and 4 interfaces for the extruded and revolved rectangle, and 8 and 52 for the disk and ellipse.

## Replaying the solid torus

The recipe provides the disk center and radius and requests a turn about Y. The planning command reports the expected 20 blocks before construction. [Torus recipe](../tests/references/primitive_circle_revolution.recipe.json) · [Six CLI trials](results/primitive-sweep-cli.json)

```sh
mkdir -p build/templates
cargo run --locked --offline -- plan-sketch tests/references/primitive_circle_revolution.recipe.json
cargo run --locked --offline -- build tests/references/primitive_circle_revolution.recipe.json build/templates/torus.json
cargo run --locked --offline -- audit build/templates/torus.json
cargo run --locked --offline -- export-step build/templates/torus.json build/templates/torus.step
```
