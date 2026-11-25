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
| Panels on each planar perforated face x = 0 and y = 0 | 24 | 5 |

The five variants share this topology. Four maps in each arm connect the circle directly to the boundary of the rectangular arm section. The fifth panel on the large face belongs to the bend. The two bend maps are shared by both arms. The eleven interfaces compare homogeneous control nets, with degree elevation when needed. The rectangular single-bore prism has four blocks and four interfaces for each tested extrusion axis X, Y, and Z. Reflected and translated L profiles were also rebuilt in all four plane orientations.

The eight cylindrical bore faces are exported on STEP `CYLINDRICAL_SURFACE` supports. Recognition checks the quarter-circle coordinates and weights, radius, and extrusion direction. The first attempt represented the same surfaces as rational B-splines in STEP. On the Ø10/Ø5 variant, FreeCAD/OCCT classified eight material points as void on a plane tangent to the small hole. Coordinates and file hashes from that first attempt are preserved in the [counterexample](../../tests/counterexamples/fieldcraft_compact_x80_y40_cylindrical_step_supports.json).

FreeCAD/OCCT reread five valid one-shell solids with 38 faces each. The 7,680 points in the two grids, 80 radial probes, and 40 probes on tangent planes agreed with the references. Across the five STEP pairs, the largest absolute adaptive-integration discrepancies were `2.19 × 10⁻¹¹ mm³` in volume and `1.64 × 10⁻⁷ mm²` in area.

## Analysis on the partition

On the original bracket, the unrefined Q2 solve reproduces a 2,000 N resultant over 320 mm² with 171 scalar degrees of freedom, compared with 705 previously. At quadrature order 16, strain energy is 135.445 N·mm with ten blocks and 144.124 N·mm with the old partition. The 6.02% difference compares two different analysis spaces. Geometric qualification leaves convergence of stress and energy at higher refinement open.

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
