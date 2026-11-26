# SolidCraftTrivariate

A native Rust volumetric kernel for building and auditing rational trivariate maps, exporting STEP, and preparing their IGA analysis spaces. The geometry workflow starts from a sketch or recipe, derives its partition lines, then builds and audits the volumes in the same crate. Python supplies independent oracles and historical references. [Documentation index](docs/README.md)

## Polygonal sketch input contract

The native `sketch_layer_partition_v1` input is a set of straight-edged polygon contours in the XY plane, with coordinates in integer millimeters. Each loop has 3–64 distinct vertices, a unique identifier, and must be simple: self-crossings, repeated vertices, and contacts between contours are rejected. The outer loop must enclose every cut. Without Boolean operation history, cut loops must also be mutually disjoint and non-nested; with history, nested loops are allowed, but contours still may not touch or cross. A recipe supports at most 12 cuts, an integer height from 1 to 4096 mm, integer cut intervals satisfying `0 ≤ z0 < z1 ≤ height`, and an optional integer top shift of at most 4096 mm per axis. Coordinates and construction-guide endpoints are limited to ±4096 mm; at most 64 construction guides are allowed. Clockwise or counterclockwise input winding is accepted and normalized by the planner.

These are polygonal line segments, not sampled approximations of curves: this input family does not accept arcs or NURBS contours. Automatic partitioning tries certified strategies for the supplied geometry and can still reject a valid contour set when no supported partition is found. Passing one strategy does not imply that every simple polygon is supported. Curved profiles use their separate NURBS recipe families.

The [audit gallery](https://jeffwitz.github.io/SolidCraftTrivariate/) opens four collections: Shapes, Compare, Fillets and chamfers, and Boolean operations. The Shapes explorer contains 68 entries. Browse the 45 polygon benchmark pieces, the overlapping-projections diamond diagnostic, the ten bore propagation controls and the unequal-width two-bore bracket alongside six curved-profile parts and five FieldCraft variants. Each entry includes the native model, recipe, audit and STEP export. [Open the east-wing bore stack](https://jeffwitz.github.io/SolidCraftTrivariate/?model=bench-bracket-east-stack) · [Polygon benchmark](docs/convex-host-rings.md) · [Bracket measurements](docs/results/propagated-bores.md)

The operation studies include the [box-minus-cylinder comparison with four classical sectors](docs/results/box-cylinder/README.md), the [fixed-partition rim models and STEP checks](docs/results/fixed-rims/README.md), the [crossed-tool Boolean constructor and reproduction commands](docs/crossed-tool-boolean.md), and the [IGA and hexahedral comparison at the crossed interface](docs/results/crossed-interface-iga/hexahedral-comparison.md). The [planar partition quality study](docs/results/planar-external-review/README.md) includes its independent input corpus and replay script.

The automatic polygonal planner derives separators, quadrilateral rings and pocket floors directly from the contours. On the 45-piece polygon corpus it selects 1,208 blocks in total, including 910 on the 30 seeded convex sketches. All 45 published recipes and 12 additional sketches passed independent geometry and volume checks. On eight of the 45 cases the sampled 3D condition number increased. Across the current corpus the maximum in-plane condition number is 229.4 on seeded case 14. The maximum 3D value is 892.9 on seeded case 09. [In-plane distortion and extrusion aspect](docs/results/polygon-extrusion-quality/README.md). [Current measurements](docs/convex-host-rings.md) · [Partition selection](docs/architecture.md#planner)

The orthogonal sketch planner builds all ten public bore-propagation controls. The east and north bore stacks each use 17 blocks. The five original FieldCraft variants use ten blocks each. [Bracket measurements](docs/results/propagated-bores.md) · [Downloadable FieldCraft geometry](docs/results/fieldcraft-replay/README.md)

## Quick start: sketch, plan, volumes

Requirements: Rust/Cargo, Python 3.10+, Git, Bash, and Make. From the repository root, the initial Rust workflow does not depend on upstream clones or Python oracles:

```sh
cargo build --release --locked --offline
make demo-native
```

`make check` runs the local Rust and Python tests, demos, and published-geometry replay without requiring any upstream checkout. To also require and validate all six pinned upstream references listed in `references.lock`, run `make bootstrap` once and then `make check-strict`. The standalone `make refs` target checks those clones without installing them.

`plan-sketch` reports contour junctions and cutting levels induced by other operations before blocks are created. In this example, a periodic cubic profile is subtracted by extrusion from a transverse rectangle. The plan is saved as inspectable JSON:

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" plan-sketch tests/references/g2_cross_sketches.recipe.json > build/demo/g2-plan.json
python3 -m json.tool build/demo/g2-plan.json
```

The recorded report contains eight G2 junctions on the periodic profile, G0 corners on the rectangle, four exact crossings at parameters `7/2`, `13/2`, `5/2`, and `15/2`, and transverse levels at 4 and 8 mm. These values are specific to this recipe. Knots inserted to represent a curve and crossings between two sketches retain distinct provenance.

Build the volumetric representation from the same recipe, then audit and export the solid:

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" build tests/references/g2_cross_sketches.recipe.json build/demo/g2-cross.json
"$cli" audit build/demo/g2-cross.json
"$cli" export-step build/demo/g2-cross.json build/demo/g2-cross.step
```

For this recipe, the recorded audit reports 32 blocks, 46 shared interfaces, 100 exterior faces, a determinant lower bound of `5.729862499996598 mm³`, and two through bores. The bound certifies local regularity of the maps on their domains; by itself it does not establish global fidelity for a different recipe.

### Build rings from polygonal sketch contours

For convex stock with convex holes or blind pockets, the automatic planner derives candidate separator lines, matching ring vertices and pocket floors directly from the outlines. The `construction` field can be omitted. Inspect the selected rational sketch graph and its axial layers before exporting the volume.

```sh
recipe=tests/regressions/overlapping_projections_diamonds.recipe.json
target/release/solidcraft-trivariate plan-sketch "$recipe"
target/release/solidcraft-trivariate build "$recipe" build/demo/diamonds.json
target/release/solidcraft-trivariate audit build/demo/diamonds.json
```

Layered polygon recipes also support exact NURBS strip assembly and patch-aware selection of sketch separators. [Recipe modes and benchmark](docs/sketch-nurbs-partition.md)

The two diamond holes overlap in both axis projections. Their profiles are disjoint and can be separated by an oblique chord. [Planner and measured results](docs/convex-host-rings.md)

### Import a SolidCraft timeline

`import-solidcraft-sketch` reads a canonical timeline exported from SolidCraft; it does not launch or execute the user's script. The qualified subset is an orthogonal polygon sketch, its extrusion, a fillet, and two through bores on distinct axes. A versioned fixture lets you replay import, audit, and STEP export:

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" import-solidcraft-sketch tests/references/solidcraft_sketch_bracket.timeline.json build/demo/bracket.json
"$cli" audit build/demo/bracket.json
"$cli" export-step build/demo/bracket.json build/demo/bracket.step
```

The timeline above builds the compact bracket with 10 blocks, 11 shared interfaces, 38 exterior faces and two bores. Its certified determinant lower bound is 67.7333333333291 mm³. Across the five FieldCraft brackets, the worst sampled Jacobian condition number falls from 2,048 on the archived 48-block grid to 13.33 on the compact partition (125 points per block, including its corners). This measures geometric distortion. The [mechanical refinement study](docs/results/fieldcraft-bracket-pa5-mechanics.md) reports energy versus degrees of freedom for each analysis space. Details of the canonical input and atomic rejection are in [the timeline import test](tests/solidcraft_sketch_features.rs). For other SolidCraft script formats, see `import-solidcraft` and its documented subsets in the [campaign results](docs/campagne-20261125.md).

## Construction families and scope

Versioned JSON recipes express operations and dimensions; the audit rebuilds and checks candidate maps before publication. Implemented families include affine pairs, annuli, brackets, axial and moving bores, multi-loop NURBS contours, intersecting NURBS sketches, and rectangle/circle/ellipse primitives made by extrusion or revolution. Each family has its own admission assumptions. Success on a template is not a universal guarantee for every geometry.

The `plan-sketch` preplanner handles junctions, subdivisions, and intersection events required by specific families. Bernstein certificates apply to the domains and methods named in each report. Classifications, grids, and samples remain measured diagnostics; they are not presented as proofs of coverage or injectivity. Rejections and counterexamples are retained under `tests/counterexamples`, `tests/regressions`, and the versioned reports.

The primitive catalog records conventions and replayed results: [NURBS templates](docs/gabarits-primitives-nurbs.md). For imported curves and cross-operations, see [NURBS sketch prepartitioning](docs/prepartition-esquisse-nurbs.md) and the [multi-loop NURBS cycle](docs/cycle-nurbs-multiboucles.md). These documents distinguish measurements for each recipe from properties proved over its domain. The corrected specification remains normative; adaptation notes and campaign logs provide context, recipes, and results without changing its requirements.

The native demo also includes a versioned prismatic pocket:

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" build tests/references/sketch_boolean_blind_pocket.recipe.json build/demo/blind-pocket.json
"$cli" audit build/demo/blind-pocket.json
```

This fixture produces 9 blocks and a certified determinant lower bound of `191.99999999999991 mm³`. The broader reference-test recipe (`sketch_boolean_oracle.py`) uses a different domain and operation history; its values do not replace the demo's.

## Other workflows

### Mechanical analysis

The FieldCraft PA5 driver rebuilds the original two-bore bracket from its sketch recipe, checks its native audit, and runs the prescribed-load elasticity problem. The two pierced arms use four blocks each and the bend uses two. The driver checks the 320 mm² load area, 2,000 N resultant and the volume at higher quadrature orders.

```sh
cargo build --release --locked --offline
python3 scripts/fieldcraft_bracket_pa5.py --max-level 1
```

The output under `build/fieldcraft-bracket-pa5/` includes the rebuilt model, solve results and energy-versus-DOF plot. The [mechanics report](docs/results/fieldcraft-bracket-pa5-mechanics.md) also records earlier energy measurements for comparison.

### Editing bore depth

The polygonal `closed` recipe in the corpus can replay an axial pocket and audit its edit:

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
python3 -c 'import json; print(json.dumps(next(r["recipe"] for r in json.load(open("tests/references/interval_bores_corpus.json")) if r["id"]=="closed")))' > build/demo/cavity.recipe.json
"$cli" build build/demo/cavity.recipe.json build/demo/cavity.json
cat > build/demo/edits.json <<'JSON'
{"schema":"interval_bore_edits_v1","steps":[{"bore":"pilot","z_interval_mm":[0,9],"require_features":["bore/pilot/upper"]}]}
JSON
"$cli" edit-bore-intervals build/demo/cavity.json build/demo/edits.json build/demo/blind.json
"$cli" audit build/demo/blind.json
"$cli" export-step build/demo/blind.json build/demo/blind.step
```

If a required selection disappears, the edit fails before publishing the target. See the [independent edit oracle](tests/references/interval_bore_edits.py).

### Python reference and evaluation

Historical Python references can be replayed independently of the Rust crate:

```sh
python3 -m trivariate capabilities
python3 -m trivariate build tests/references/affine_pair.json build/affine_pair.json
python3 -m trivariate audit build/affine_pair.json
python3 -m trivariate evaluate tests/references/quarter_annulus.json 0.5 0.5 0.5
```

On a fresh clone, `make bootstrap` retrieves the six oracle dependencies pinned in `references.lock`; afterward, `make check-strict` requires their exact locked checkouts and runs the local checks offline. For Rust tests alone: `cargo test --locked --offline`.

## Documents and results


- [Architecture and contracts](docs/architecture.md)
- [L1 evaluator and results](docs/L1.md)
- [Luna development roadmap](docs/LUNA.md)
- [Development instructions](AGENTS.md)
- [Corrected scientific specification](docs/kernel-requirements.md)
- [Reference bracket construction](docs/archive/fieldcraft-design-study.md)
- [Adaptation and review organization](docs/archive/development-workflow-study.md)
- [Campaign log and measured results](docs/campagne-20261125.md)
- [Interactive gallery](https://jeffwitz.github.io/SolidCraftTrivariate/)
