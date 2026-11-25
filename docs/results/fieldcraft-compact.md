# Local partition of the FieldCraft brackets

The baseline is revision `aa94b1b`. Its five bracket models each have 48 blocks and 74 interfaces. Each arm's perforated exterior face contains 24 panels.

## Qualification criteria

- Rebuild all five solids from the same recipes. Preserve the L outline, R6 fillet, and the positions and diameters of both holes. Keep old models readable and auditable.
- For a rectangular prism with a through bore, connect each outside side directly to a quarter-circle. The target is four through-blocks. For the bracket, try four blocks per arm and two in the bend.
- Verify positive determinant on each closed cube. Check homogeneous-coordinate interfaces, cylindrical walls, and planar footprints.
- Compare each STEP file with its corresponding historical solid using adaptive integration and material classification, including points ±0.02 mm from each wall. Use the geometric thresholds from the preceding replay.
- Count panels on perforated exterior faces directly from the maps. Fewer blocks do not prove minimality among all possible 3D partitions.

## Geometric result

| Measure per bracket | Previous grid | Local partition |
|---|---:|---:|
| Trivariate blocks | 48 | 10 |
| Block interfaces | 74 | 11 |
| Exterior block faces | 140 | 38 |
| Maximum sampled Jacobian condition number, each of five brackets | 2,048 | 13.33 |
| Certified determinant lower bound, canonical timeline import (mm³) | 0.00854492 | 67.7333 |
| Panels on each planar perforated face x = 0 and y = 0 | 24 | 5 |

The condition numbers above are the maximum of σmax(DX)/σmin(DX) over all blocks on the closed tensor grid {0, ¼, ½, ¾, 1}³. Reproduce them with `export-geometry-vtu MODEL 4 OUTPUT.vtu` and take the maximum of the `condition_number` point array. The determinant bounds come from continuous certificates.

The five variants share this topology. Four maps in each arm connect the circle directly to the boundary of the rectangular arm section. The fifth panel on the large face belongs to the bend. The two bend maps are shared by both arms. The eleven interfaces compare homogeneous control nets, with degree elevation when needed. The rectangular single-bore prism has four blocks and four interfaces for each tested extrusion axis X, Y, and Z. Reflected and translated L profiles were also rebuilt in all four plane orientations.

The eight cylindrical bore faces are exported on STEP `CYLINDRICAL_SURFACE` supports. Recognition checks the quarter-circle coordinates and weights, radius, and extrusion direction. The first attempt represented the same surfaces as rational B-splines in STEP. On the Ø10/Ø5 variant, FreeCAD/OCCT classified eight material points as void on a plane tangent to the small hole. Coordinates and file hashes from that first attempt are preserved in the [counterexample](../../tests/counterexamples/fieldcraft_compact_x80_y40_cylindrical_step_supports.json).

FreeCAD/OCCT reread five valid one-shell solids with 38 faces each. The 7,680 points in the two grids, 80 radial probes, and 40 probes on tangent planes agreed with the references. Across the five STEP pairs, the largest absolute adaptive-integration discrepancies were `2.19 × 10⁻¹¹ mm³` in volume and `1.64 × 10⁻⁷ mm²` in area.

## Slender and unequal-width controls

The [slender plate](../../tests/references/slender_plate_host.recipe.json) measures 100 × 12 × 12 mm and has a through-bore of radius 3 mm centered at x=80, y=6. The whole-section candidate fails its spatial regularity check. Four sectors in a 12 × 12 mm host plus two affine boxes give six blocks and six full-face interfaces. Relative to revision `8eb6ea0`, the sampled maximum κ drops from 96 to 6.17 on the same 125-point grid, and the certified determinant lower bound rises from 0.0100 to 13.9629 mm³. Both partitions give material volume 14,400 − 108π mm³.

Adding a third bore of radius 2 mm at x=54 mm on the east arm of the original x64-y64 FieldCraft bracket gives 14 blocks and 16 interfaces, compared with 64 and 99 on the old grid. Sampled κ falls from 1,024 to 8.53, and the certified determinant lower bound from 0.0004605 to 2.10286 mm³. The [three-bore regression](../../tests/regressions/fieldcraft_three_bores.recipe.json) stores the physical dimensions in millimeters. To reproduce the comparison with the archived grid, apply the same edit to the x64-y64 recipe: append a Y-axis bore at [864,64,320] with radius 32 inside its source model, keeping the 1/16 affine transformation.

The [unequal-width bracket](../../tests/references/unequal_fillet_cross_bores.recipe.json) has 12 mm and 20 mm arms, an R4 bend and two transverse R2 bores. It replays on the orthogonal grid with 48 blocks and 74 interfaces. Its material volume is 30,912 − 176π mm³. The plan records why the compact candidate was inapplicable. Independent STEP readback yielded one valid closed solid for the slender plate and the unequal bracket, with adaptive volumes agreeing with their analytic targets. [STEP replay and measurements](../../tests/references/orthogonal_step_freecad_result.md)

## Analysis on the partition

On the original bracket, the unrefined Q2 solve reproduces a 2,000 N resultant over 320 mm² with 171 scalar degrees of freedom, compared with 705 previously. At quadrature order 16, strain energy is 135.445 N·mm with ten blocks and 144.124 N·mm with the archived timeline partition. The latter has κ = 64 on the 125-point grid. The scaled-coordinate archive used in the five-variant quality table has tighter bore guards, κ = 2,048 and energy 143.725 N·mm on the same physical geometry. The 6.02% difference between the compact and timeline energies compares two different analysis spaces. The [energy refinement study](fieldcraft-bracket-pa5-mechanics.md) compares both discretizations under the same load. Local stress convergence requires a separate check.

## Files and replay

[Download the five brackets and control prism](fieldcraft-replay/fieldcraft-compact.zip) · [Open the bracket in the gallery](https://jeffwitz.github.io/SolidCraftTrivariate/?model=fieldcraft_x64_y64)

From a repository clone, these commands rebuild the control prism:

```sh
mkdir -p build
cargo run -- build tests/references/four_sector_plate.recipe.json build/plate.model.json
cargo run -- audit build/plate.model.json
cargo run -- export-step build/plate.model.json build/plate.step
```

For a bracket, use the recipe `docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.recipe.json`. The builder chooses the partition from the sketches. The 48-block models remain auditable using the plan stored in their files.
