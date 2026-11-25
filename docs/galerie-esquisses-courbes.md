# Curved sketch gallery

The six recipes reuse the blind-pocket plate, two closed cavities, two C-sections, and L-shaped prism from the polygonal series. The Rust constructor preserves the block and interface identities from the sketch plan. It replaces trilinear maps with Bézier volumes of degree (2, 2, 1), using integer homogeneous controls with a common weight. These polynomial maps are special cases of NURBS. Arcs and surfaces appear in model coordinates and in the STEP export.

## Construction rules

The rectangle for each operation defines the center and dimensions of its curved profile. For a through-hole, the constructor moves the four corners toward the center by a factor of 2/3. It connects successive corners with four quadratic arcs. Each arc's middle control is offset outward by half the transverse dimension of the reduced rectangle. Parametric tangents agree at the junctions. The hole's X and Y envelope remains that of the source rectangle.

For a blind pocket or closed cavity, the corresponding factors are 4/5 and 1/4. The four arcs have distinct tangents at their junctions. This construction preserves a regular floor block, including at its corners. On each part, the outer side with the smallest ordinate is also curved. Its middle control is moved outward by one eighth of the profile height in the XY plane. The geometric sagitta at the arc midpoint is half that offset.

Edge controls are shared by adjacent maps. The center control of the biquadratic net is obtained by Coons interpolation. The sections are then extruded in the direction specified by the recipe. For the C, the top face is translated by (3, −2, 7) mm.

| Configuration | Blocks | Interfaces | Exact volume (mm³) |
|---|---:|---:|---:|
| Blind pocket | 9 | 16 | 48,976 |
| Centered cavity | 14 | 28 | 247,952 / 25 |
| Offset cavity | 14 | 28 | 247,952 / 25 |
| C with two through-holes | 9 | 10 | 45,878 / 27 |
| L with through-hole and cavity | 26 | 51 | 49,792 / 9 |
| C with hole shifted by 1 mm | 9 | 10 | 45,878 / 27 |

The block-count lower bound applies to partitions shared across layers, with whole-face interfaces and four prescribed junctions per profile. These four junctions may be smooth for a through-hole. Retaining them as planar vertices is part of the counting assumption. [Combinatorial proof](partition-poche-comparaison.md)

## Replaying a part

From the repository root, this example rebuilds the C with two holes, runs the native audit and independent oracle, then writes STEP.

```sh
cargo build --release --locked --offline
mkdir -p build/curved-gallery
cli=target/release/solidcraft-trivariate
recipe=docs/site/data/c_section_two_holes/c_section_two_holes.recipe.json
model=build/curved-gallery/c.model.json
"$cli" build "$recipe" "$model"
"$cli" audit "$model"
python3 tests/references/curved_sketch_gallery.py --recipe "$recipe" --model "$model" \
  --expected-block-count 9 --output build/curved-gallery/c.audit.json
"$cli" export-step "$model" build/curved-gallery/c.step
```

The Rust audit replays the recipe, compares oriented traces, and certifies positive Jacobians on every closed parametric cube. The Python oracle checks arc controls on the walls, bounds the determinant with rational arithmetic, and integrates volume in two independent ways. One integrates map determinants. The other integrates sketch contours with Green's formula and sums layer volumes. The gallery provides each recipe, model, audit, and STEP for download. [Gallery data](site/data)

All six STEP exports were re-imported with FreeCAD and OCCT 7.8.1. The reader finds one valid solid per file, with two closed shells for parts with internal cavities and one for the other three. Relative volume errors against exact integration are below 10⁻¹⁴. The oracle's `--step build/curved-gallery/c.step` option adds this check when the FreeCAD and Part Python modules are available.
