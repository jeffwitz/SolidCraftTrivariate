# Exact NURBS strips for layered polygon sketches

Baseline: kernel and published gallery in GitHub revision e0121602. The experiment is restricted to positive certified extruded trilinear cells. The native planner, assembler and sketch interchange implement the two recipe modes. Existing recipes keep their representation. The independent extraction oracle must recover every source cell once and all oriented traces. Reject unknown modes and edited patch controls during replay. Compare patch and cell counts, sampled quality and conforming scalar degrees of freedom over the 46 polygon gallery recipes, including the through octagon, diagonal holes and two-depth pocket. The octagon target is four patches and eight exact spans.

## Build and inspect a patch layout

Set the optional recipe field `nurbs_partition` to `assemble` to keep the current partition, or `optimize` to compare alternative partitions of the same sketch. Both modes emit explicit patches. An omitted field retains the previous output format. The following commands reproduce the octagon comparison from the repository root.

```sh
mkdir -p build
cargo build --release --locked --offline
cli=target/release/solidcraft-trivariate
for mode in assemble optimize; do
  jq --arg mode "$mode" '. + {nurbs_partition: $mode}' \
    docs/site/data/bench-poly-through_octagon_high_aspect/*.recipe.json > "build/octagon-$mode.recipe.json"
  "$cli" plan-sketch "build/octagon-$mode.recipe.json"
  "$cli" build "build/octagon-$mode.recipe.json" "build/octagon-$mode.model.json"
  "$cli" audit "build/octagon-$mode.model.json"
done
```

Each entry in `nurbs_patches` contains its three degrees, three clamped knot vectors, ordered homogeneous control sections, and the ordered source cell IDs. Consecutive sections define one Bézier span. Its four control points at each end are ordered by radial coordinate, then extrusion coordinate. The knot domain along the strip is [0, span count]. Interior knots have multiplicity one at degree one and permit corners. The other two parameters range from zero to one. To transfer a design-face label onto a span, use its source cell ID and its `source_xy_rotations` entry. For patch coordinates (s,r), rotations 0 through 3 give the source coordinates (s,r), (1-r,s), (1-s,1-r), and (r,1-s), respectively. The extrusion coordinate is preserved.

Before serializing a patch, the builder extracts every span and checks all eight ordered controls against the source under one positive homogeneous scale. It verifies that each source cell occurs exactly once. Model audit also replays the planner and checks the knots and controls. Regularity and coverage are inherited from the certified source partition. A corner at an internal knot has one-sided Jacobians. Interfaces across different patches use the source span incidence, so downstream programs can identify subfaces even when whole patch sides do not match one-to-one. The existing STEP and analysis exporters continue to use the extracted cell representation.

Patch count measures the number of volumes for design operations. Span count governs elementwise quadrature. Grouping a fixed partition with C0 knots preserves its conforming Q1 space. A different choice of separators can change the cell count, analysis space and conditioning. A minimum over all possible partitions is outside the search contract.

## Candidate search

The assembler considers straight lines through sketch-cell edges. A strip can follow a straight side while its opposite side bends at internal knots. Matching successive cross-sections must coincide exactly, up to the common scale of their source nets. A candidate must fit the integer exchange format. The assembler examines both in-plane axes and all four positive rotations. It orders strips by descending length and then by source index, accepts disjoint strips, and retains unused cells as singletons. This deterministic search is deliberately bounded; arrangements of strips that turn on both sides remain outside it.

The opt-in planner compares certified plain, guided, rectangular-host, convex-host, radial and ring proposals when applicable. Its ranking is patch count, followed by Bézier cell count and then sampled 3D spectral condition number on a 3×3×3 grid. Explicit guides and forced triangulation go through their existing planner. The independent measurement below uses a denser grid and also measures in-plane distortion. To replay from the repository root:

```sh
python3 scripts/check_sketch_nurbs.py --cli target/release/solidcraft-trivariate
```

## Results on the thickened polygon gallery

| 46 sketches | Assemble the existing partition | Select separators by patch count |
|---|---:|---:|
| NURBS volumes | 789 | 707 |
| Extracted Bézier cells | 1,218 | 1,343 |
| Span interfaces verified by the independent checker | 2,083 | 2,355 |
| Distinct geometric Q1 corner nodes, summed per model | 3,204 | 3,450 |

The selector reduced patch count in 19 cases and tied in 27. It used additional cells in 18 cases. Exact rational integration matched the published reference volume for every case in both modes. Every span reproduced its source control net and declared orientation. The independent checker compared complete oriented interface nets. An initial oracle normalized homogeneous controls before changing their order and interchanged two rotations. Reordering first and normalizing second restored the rotation-invariant comparison.

| Case | NURBS volumes, assemble → optimize | Cells | Scalar Q2 DOFs from analysis export |
|---|---:|---:|---:|
| Through octagon | 4 → 4 | 8 → 8 | 144 → 144 |
| Two diagonally separated diamonds | 4 → 4 | 10 → 10 | 165 → 165 |
| Blind cuts at two depths | 36 → 24 | 36 → 96 | 416 → 1,082 |

The octagon groups the original eight cells as 3 + 3 + 1 + 1. Its certified determinant lower bound exceeds 579.99 mm³ on each unit span. Its exact volume is 50,880 mm³. FreeCAD reread both STEP exports as a single closed valid solid, each measuring 50,879.999999999985 mm³. STEP exports retain the existing boundary subdivision.

On the 5×5×5 diagnostic grid the worst 3D condition number per case increased on 13 parts and decreased on five. The largest ratio was 9.80 on seeded convex case 01. The two-depth pocket went from 6.59 to 15.5. Patch count alone is therefore an unsuitable default cost measure for analysis. Use the two recipe modes to compare design complexity and analysis cost on a given sketch. Q1 coordinate counts in the table are geometric diagnostics; the Q2 column comes from the actual analysis packages.

[Per-case counts, volume and quality](results/sketch-nurbs-comparison.json)
