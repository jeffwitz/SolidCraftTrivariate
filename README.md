# SolidCraftTrivariate

A native Rust volumetric kernel for building and auditing rational trivariate maps, exporting STEP, and preparing their IGA analysis spaces. The geometry workflow starts from a sketch or recipe, derives its partition lines, then builds and audits the volumes in the same crate. Python supplies independent oracles and historical references. [Documentation index](docs/README.md)

The gallery contains six curved-profile solids and five FieldCraft bracket variants. Each entry shows a solid section and its sketch arcs, with the recipe, model, audit, and STEP file available for download. [Interactive gallery](https://jeffwitz.github.io/SolidCraftTrivariate/) · [Construction and checks](docs/galerie-esquisses-courbes.md) · [FieldCraft bracket results](docs/results/fieldcraft-compact.md)

The polygonal planner compares local separators with the historical partition on the same recipe. In the recorded corpus, the blind pocket uses 9 instead of 17 volumes and the two-hole C profile uses 9 instead of 35. Both meet the topological lower bound within the conforming partition class studied. Map conditioning also affects analysis. [Algorithm](docs/architecture.md#planner) · [Corpus, proofs, and measurements](docs/partition-poche-comparaison.md) · [Sketch-based construction](docs/partition-esquisses.md)

## Quick start: sketch, plan, volumes

Requirements: Rust/Cargo, Python 3.10+, Git, Bash, and Make. From the repository root, the initial Rust workflow does not depend on upstream clones or Python oracles:

```sh
cargo build --release --locked --offline
make demo-native
```

To install the six references listed in `references.lock` and run the full reference and test suite, run `make bootstrap` once, then `make check`. `make check` checks that those clones are present; it does not install them implicitly.

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

The Rust solver takes an audited model and a versioned job. This example imports the historical bracket template, solves its reference job, and exports the field and stress for ParaView:

```sh
mkdir -p build/demo
cli=target/release/solidcraft-trivariate
"$cli" import-solidcraft upstream/fieldcraft/tests/oracle/filleted-bracket-fillet-stress/design.json build/demo/bracket.json
python3 -c 'import json; print(json.dumps(json.load(open("tests/references/mechanical_jobs/bracket.json"))["job"]))' > build/demo/job.json
"$cli" solve-job build/demo/bracket.json build/demo/job.json build/demo/result.json
"$cli" export-field-vtu build/demo/result.json 4 build/demo/bracket-field.vtu
```

The historical FieldCraft bracket template has 36 blocks and 55 interfaces. The sketch-based builder reproduces the bracket and its four two-bore variants with 10 blocks and 11 interfaces. Each pierced arm uses four blocks and the bend uses two. [Models, STEP, and checks](docs/results/fieldcraft-compact.md). The result retains displacements in mm, reactions, work-energy balance, and per-block stresses at interfaces. The solver handles small-strain linear isotropic elasticity. Comparisons and limits for each experiment are recorded in the [campaign log](docs/campagne-20261125.md).

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

On a fresh clone, `make bootstrap` retrieves the six oracle dependencies pinned in `references.lock`; afterward, `make check` runs offline. For Rust tests alone: `cargo test --locked --offline`.

## Documents and results


- [Architecture and contracts](docs/architecture.md)
- [L1 evaluator and results](docs/L1.md)
- [Luna development roadmap](docs/LUNA.md)
- [Development instructions](AGENTS.md)
- [Corrected scientific specification](CdC_noyau_trivarie_sans_trim_corrige.md)
- [Reference bracket construction](LienSolidecraft.md)
- [Adaptation and review organization](ArtcratAdaptation.md)
- [Campaign log and measured results](docs/campagne-20261125.md)
- [Interactive gallery](https://jeffwitz.github.io/SolidCraftTrivariate/)
