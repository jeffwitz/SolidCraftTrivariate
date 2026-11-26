# Integrating SolidCraft and FieldCraft Around a Trivariate Kernel

## Archival scope

This study examined the `bracket-hex` example: a filleted, two-bore bracket built in SolidCraft, partitioned into **30 cells**, and meshed by FieldCraft. It assessed how that example could seed a trivariate kernel. The central architectural change is to preserve the blocks as native geometry—with spline maps and interfaces—from part construction onward.

The source review found useful construction and meshing components, but also a clear boundary: FieldCraft currently produces meshes from the blocks, while SolidCraft retains B-rep bodies. A native trivariate model must sit between feature history and derived representations.

## 1. The example and what it demonstrates

### The bracket partitioned into 30 blocks

Relevant sources:

- [Construction and partition recipe](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket-blocks/recipe.json).
- [Image of the 30-block partition](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/docs/media/bracket-partition.png).
- [Code preparing the `bracket-hex` calculation](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/engine/src/sample.rs#L556-L580).
- [Multiblock mesher `hex_mesh_cells`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L1-L138).

The part is an L profile **60 × 50 mm**, with **8 mm** legs, extruded by **40 mm**. It has a **radius-6 mm** inside fillet and two **diameter-8 mm** bores, one in each leg. Their axes differ, so the finished part is beyond a simple section extruded in one direction.

The recipe prescribes **29 `solid.split_body` operations**:

| Step | Resulting partition |
|---|---|
| Splits at z=10 and z=30 | Three widthwise slices |
| Splits in the legs and a diagonal split in the bend | Eight cells in each slice |
| Splits around each bore in the middle slice | Each cell containing a bore becomes four cells arranged around it |

The total is

\[
8+(8+3+3)+8=30\ \text{cells}.
\]

**The partition is explicitly recorded in the recipe:** cut planes, their order, and successive body names are supplied. That makes the example useful as a starting case because the block-construction decisions are already available for inspection.

The application opens the example with:

```sh
fieldcraft --sample bracket-hex
```

The documentation then describes a calculation with **5,778 hex20 elements** at a requested size of 3 mm. The **30 partition blocks** must be distinguished from the thousands of elements used to discretize them. [FieldCraft example](https://tangency-labs.github.io/fieldcraft/examples.html?utm_source=chatgpt.com)

### What FieldCraft does with the partition

The organization in `cells.rs` is directly relevant:

1. Recognize each cell as a **logical cube**.
2. Identify shared vertices, edges, and faces.
3. Assign compatible subdivision counts across blocks.
4. Construct each edge once.
5. Construct each face grid once.
6. Fill each interior using **trivariate transfinite interpolation**.
7. Check the resulting mesh.

A shared face therefore has one grid used by both cells. Once topology is established, nodes are shared by identity. The corpus fixture records **30 bodies and 45 interfaces**. [Multiblock algorithm](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L1-L32) · [Fixture](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket-blocks/fixture.json).

The reusable idea is **a common topology followed by common geometry on interfaces**.

### An automatic-decomposition example

The second example is **`crank-arm`**, a crank. After partitioning along its symmetry plane, `mesh.hex.generate` automatically decomposes its section into regions suited to structured meshing.

The implementation is [`hex/decompose.rs`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/decompose.rs#L1-L48). It separates holes, proposes splits around them, and seeks cuts from re-entrant corners and fillets.

Two limits inform the proposed strategy:

- partitions are **internal to the mesher** and are not retained in the geometric model;
- the meshed variant is adjusted to remain sweepable: some fillets and chamfers are removed, and some geometric encounters are changed to avoid degeneracies. [Description of `crank-arm`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/docs/examples.md#L142-L165).

Start with **the controlled bracket partition**, then use the crank to extend automation.

## 2. Is the constraint sufficient?

It is sufficient to define the scope of an initial kernel, but should be stated more precisely:

> The kernel constructs solids from regular trivariate blocks whose interfaces are explicitly shared and parametrically compatible.

“Decomposable into hexahedral-meshable solids” is a good starting point, but four properties remain to be established:

| Property | Contribution to the kernel |
|---|---|
| Each block has cube topology | A three-coordinate parameter domain |
| Each block has a regular, injective spline map | Volumetric geometry usable for analysis |
| Neighboring blocks have compatible traces | Well-defined geometric joins and field assembly |
| Blocks form a partition of the solid | Complete coverage without interior overlap |

The key distinction is **block topology versus parameterization**. Six logical faces provide a construction structure; they do not alone guarantee a good map \(X(u,v,w)\).

One instructive detail is that the 30 bodies have **188 B-rep faces**, whereas 30 logical cubes have 180 sides. FieldCraft can group multiple geometric faces into one composite logical face, for example near a tangent join between a plane and fillet. [Fixture](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket-blocks/fixture.json) · [Composite-face construction](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L9-L15).

Thus a block face may be a spline surface spanning multiple knot intervals. If a logical cell needs further subdivision to obtain regular maps, the kernel must be able to produce that subdivision explicitly.

**The 30-cell plan is a first construction target, with map regularity as an acceptance criterion.** Thirty is a target to verify, not a reason to retain a defective map.

## 3. The decisive change in FieldCraft: retain the maps

### Current state

The output types are explicit:

- `HexMesh` contains nodes, connectivity, boundary faces, and a report;
- `CellsMesh.blocks` contains the subdivision counts in each block direction.

These objects contain no volumetric control net, NURBS weights, or knot vectors. [Sources: `HexMesh`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/mod.rs#L213-L223) · [`CellsMesh`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L52-L72).

FieldCraft does have continuous maps for finite elements and analytic evaluators for some constructions. **What it lacks is retention of design blocks as spline volumes independent of a mesh.** For example, [`transfinite_volume`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/quad/tfi.rs#L204-L240) fills a point array from boundary grids. It is a useful formula reference, but its result remains a grid.

### Proposed objects

| Proposed object | Responsibility |
|---|---|
| `BlockPlan` | Block-construction recipes, incidences, orientations, and provenance |
| `TrivariateComplex` | Constructed spline maps, canonical interfaces, boundary, and validity |
| `AnalysisSpace` | Analysis bases, refinement, DOFs, constraints, and quadrature |

The geometric core would include, for example:

```rust
// Proposed schema: these types do not yet exist in the projects.
struct Trivariate {
    degree: [usize; 3],
    knots: [Vec<f64>; 3],
    homogeneous_control_points: TensorGrid<[f64; 4]>,
}

struct NativeSolid {
    blocks: Vec<Trivariate>,
    interfaces: Vec<SharedInterface>,
    design_faces: Vec<DesignFace>,
    provenance: ConstructionHistory,
    validity: ValidityReport,
}
```

The key point is that **requested analysis subdivisions no longer define geometry**. Construct maps first. They can then be exactly refined, sampled for display, or used to derive an FE mesh.

### What can be reused from the organization

The multiblock mesher already has useful mechanisms for sharing an edge and face, converting face orientations between blocks, and propagating subdivision constraints. For example, [`FaceGrid::map_from`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L939-L971) handles grid-orientation correspondences; the same logic can be adapted to parametric directions.

Spline compatibility must additionally cover **degrees, knots, multiplicities, weights, and trace coefficients**. Equal point counts on both sides are insufficient.

**For the first demonstrator, use whole-face interfaces with identical traces.** This matches FieldCraft's multiblock mode, which already rejects cells meeting without sharing complete faces. Nested interfaces specified in the project requirements can follow. [Source: `hex_mesh_cells` contract](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/cells.rs#L100-L103).

## 4. Constructing the bracket directly in trivariate form

Turn the existing recipe into a **parametric recipe for native blocks**. Today it creates a body, drills it, then splits it. In the new kernel, features should produce or rebuild a block set whose union represents the result.

### A deliberately small initial operator set

| Native operator | Use |
|---|---|
| Planar quadrilateral and extrusion | Straight regions of the legs |
| Quadrilateral with circular boundary and extrusion | Regions around bores |
| Section with circular blend and extrusion | Filleted bend |
| Isoparametric split | Exact subdivision of an existing block |
| Assembly by shared faces | Compose the solid |
| Block-family reconstruction | Update a hole, thickness, or radius |

Extrusion is especially favorable. For a rational section \(\mathbf S(u,v)\), construct

\[
\mathbf X(u,v,w)=\mathbf S(u,v)+wh\,\mathbf d.
\]

If section control points are \(\mathbf P_{ij}\) with weights \(\omega_{ij}\), the two volumetric control layers are

\[
\mathbf P_{ij0}=\mathbf P_{ij},\qquad
\mathbf P_{ij1}=\mathbf P_{ij}+h\mathbf d,
\]

with equal weights on both layers and degree 1 in the extrusion direction. This directly produces a native NURBS volume whose definition is independent of element count.

### Bores become block constructions around voids

For each bore, FieldCraft's partition provides four regions around the hole. Use this as a template:

1. Define the bore arcs exactly.
2. Define the outer boundaries of the region reserved for the bore.
3. Construct four section patches.
4. Extrude them through the leg thickness.
5. Share their interfaces with neighboring blocks.

For this admissible family, “make a hole” becomes **a local volumetric reconstruction with known topology**. The bore wall is a restriction of volume maps and is retained in the kernel.

The exterior faces of these blocks and their neighbors must be parameterized together. Even a geometrically planar interface can have different parameterizations on its two sides, so the canonical shared face must be construction data.

### FieldCraft arcs require a real adaptation

In [`mesh/block.rs`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/block.rs#L177-L188), arcs are evaluated with sine and cosine of an angular parameter. Construct the same geometric arcs using rational quadratic segments instead. Their parameterization differs, so boundary correspondences must be recalculated.

Transfinite formulas remain useful for more general blocks, but spline coefficients must be constructed from compatible boundaries, then weights and Jacobians checked. **Applying TFI to projected points does not automatically produce an exact NURBS map.**

### Start with a qualified fillet family

The bracket's inside fillet comes from a rounded planar section extruded through the part. It is a good first native operation. Start with:

> A circular blend in a section, paired with a specified volumetric partition and explicit validity conditions.

Spatial fillets between arbitrary surfaces should follow later. This sequence delivers useful design operations while keeping volume construction controlled.

## 5. Replacing Solidcore in SolidCraft

### The facade exists, but replacement is not transparent

SolidCraft has a `solidcraft-kernel` crate that groups common operations. However, its `Body` type explicitly wraps a B-rep record and exposes, among other methods:

```rust
pub fn record(&self) -> &Brep
pub fn smlib(&self) -> Result<Arc<solidcraft_smlib::Brep>>
```

Higher layers also use Solidcore types. The change therefore affects the geometric model and its consumers, not just a Cargo dependency. [Definition of `Body`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/kernel/src/body.rs#L165-L183) · [Public accessors](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/kernel/src/body.rs#L799-L808).

### Feature evaluation is the right entry point

SolidCraft already has serialized features, including sketches, extrusions, revolutions, and fillets; parameters and expressions; feature identities; a timeline with references; and evaluation dispatch. Integration belongs around `FeatureKind` and `eval_feature_kind`, where construction intent, profiles, and parameters are still known. [Feature definitions](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/doc/src/document.rs#L887-L934) · [Evaluation](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/doc/src/eval.rs#L14094-L14131).

The restricted contract at that point should be:

> A resolved design operation receives parameters and semantic references, then returns a valid trivariate solid or an explained refusal.

Begin with the operations needed for the bracket. This allows the editor and history to remain in use while the native body gains its own representation.

### Solidcore contains volume types, but not an operational volume kernel

Solidcore does contain a `BSplineVolume` type. Its source explains that the port follows an SMLib build without `USE_VOLUMES`, so volume evaluation is unavailable. `evaluate_simple` explicitly returns `BSplineVolumeUnsupported`. [Implementation of `BSplineVolume`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/surface/b_spline_volume.rs#L31-L38) · [Evaluation](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/surface/b_spline_volume.rs#L307-L324).

There are structures to inspect, but no operational volumetric kernel to activate. Elementary routines do exist for B-spline basis evaluation, curve derivatives, and control-net refinement. They can be qualified and potentially reused behind the new API. [Basis functions](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/nurbs/knots.rs#L444-L512) · [Refinement](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/solidcore/crates/solidcore/src/nurbs/crv_adv/p2.rs#L440-L479).

Replacing Solidcore as the master model remains compatible with reusing selected mathematical components.

### Preserve naming and rendering while adapting their identities

SolidCraft's naming already carries provenance such as feature, face role, subdivision, and copy. [Persistent names](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/doc/src/persistent_name.rs#L1-L10).

Keep three concepts distinct:

- **design face:** “bore wall,” “loaded face,” or “support face”;
- **block face:** parametric boundary of a map;
- **internal interface:** join between two maps.

A bore wall can comprise four block faces while retaining one design identity. Boundary conditions attach to the design identity.

Rendering is easier to connect: SolidCraft already displays meshes with face identifiers. Triangles can be generated from the trivariate complex boundary while much of the display engine is retained. [`SceneMesh`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/viewport/src/scene.rs#L24-L40).

## 6. Using FieldCraft for IGA

FieldCraft is a useful host because geometric dependencies are concentrated in `fieldcraft-geometry`, separately from the mesher. [Geometry dependencies](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/geometry/Cargo.toml#L11-L29) · [Mesher dependencies](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/Cargo.toml#L11-L17).

The adapter would need a native volumetric output. Current output exposes primarily triangulated surfaces, discretized curves, and selected analytic surfaces; it does not carry general NURBS maps.

### IGA also requires solver integration

The current solver is organized around finite-element families. Shape functions are selected by element type, and some structures fix the maximum number of local functions at 27. A tensor-product cubic basis already uses 64 per element. [`ElementType`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/fem/src/element.rs#L57-L74) · [`MAX_NODES`](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/fem/src/shape.rs#L24-L34).

Use **Bézier extraction** with weights, extraction operators, and control-DOF connectivity. This could reuse parts of the mechanics and linear algebra while adapting basis and Jacobian evaluation, local support sizes, interface assembly, boundary conditions, and evaluation at physical result points.

The analysis space could begin with \(C^0\) joins between blocks and spline continuity inside each block. That would already provide an end-to-end native volumetric design and analysis pipeline.

### Existing checks are useful but apply to a different object

FieldCraft checks element Jacobians using Bernstein representations and subdivision. This approach is relevant, but it applies to polynomial hex8/hex20 element maps and uses `f64` coefficients. [Jacobian check](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/quality/jacobian.rs#L1-L15) · [Bernstein coefficients](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/quality/bernstein.rs#L42-L50).

Extend this logic to rational maps using the numerical enclosure specified in the project requirements, and check global partition validity separately. **Sampled hex quality remains useful as a diagnostic; geometric validity belongs to the native maps themselves.**

## 7. Recommended starting strategy

Build **a restricted standalone kernel**, followed by adapters to SolidCraft and FieldCraft. This advances the geometric core without requiring an immediate rewrite of the full CAD application.

### Step 1 — An elementary native construction

First result:

- a rational four-sided section;
- its extrusion;
- evaluation of \(X\) and \(DX\);
- boundary traces;
- exact subdivision;
- validity check;
- save and reload.

Then assemble a few blocks sharing faces.

**Success criterion:** changing display or analysis resolution leaves geometry unchanged.

### Step 2 — The FieldCraft bracket as demonstrator

Reuse its dimensions, bores, fillet, and partition plan, then express blocks directly with the new constructors. Replace body names inherited from successive splits with explicit identities such as leg block, bend block, bore sector, and slice.

The recipe provides an independent volume check. From its dimensions:

\[
\begin{aligned}
V
&=8(60+50-8)\,40\\
&\quad+6^2\left(1-\frac{\pi}{4}\right)40
-2\pi\,4^2\,8\\
&=34080-616\pi\\
&\approx 32144.7789254\ \mathrm{mm}^3.
\end{aligned}
\]

This matches the value recorded in the recipe. It does not alone establish shape equality, but is a useful analytic reference.

**Success criterion:** reproduce the part with native maps, explicit interfaces, and a boundary matching the prescribed primitives.

### Step 3 — Parametric modification

Vary leg thickness, extrusion width, bore diameters and positions, and fillet radius. The kernel must determine whether the block plan remains valid, needs reconstruction, or has left its admissible domain.

The history then has a precise function: carry parameters, dependencies, and entity provenance. The partition becomes an explicit derived construction with its own validity conditions.

### Step 4 — IGA on the same model

Reuse the example's load case: steel, a fixed end of one leg, and a 1 kN force at the other end. Compare displacements, energy, reactions, and stress evaluated at the same physical locations. Treat FieldCraft values as an initial comparison alongside convergence studies.

The README reports 0.1% agreement for the documented tip displacement. The inspected tests use other thresholds and do not all recompute the tetrahedral reference. Treat those values as project-reported results to reproduce, with their precise scope. [Bracket convergence test](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/testkit/tests/bracket_blocks.rs).

### Step 5 — Automate partition proposals

At this stage, `hex/decompose.rs` becomes a basis for proposing block plans: hole separation, bore sectors, cuts from re-entrant corners, fillet handling, and ranking of candidate partitions.

Its decision criterion needs adaptation: current evaluation depends partly on requested mesh size and can retain a free-mesh solution. The new constructor must propose an admissible partition for spline maps and preserve it independently of analysis refinement. [Decomposer principles and evaluation](https://github.com/tangency-labs/fieldcraft/blob/35e834652d8548384c2da8ef1b28d588c9af2ffd/crates/mesh/src/hex/decompose.rs#L19-L48).

## 8. Implications for the proposed architecture

| Existing base | Potential reuse | Project-specific work |
|---|---|---|
| SolidCraft | Sketches, parameters, features, commands, provenance | Native feature evaluation and new body type |
| FieldCraft multiblock | Incidences, orientations, subdivision constraints, examples | Persistent plans and spline-map interfaces |
| Transfinite constructions | Formulas and sampled reference | Rational construction and map checks |
| Solidcore NURBS components | Selected elementary routines after qualification | Volume complex, constructors, and native operations |
| FieldCraft mechanics | Mechanics infrastructure, units, solvers, visualization | IGA bases, extraction, assembly, adapted post-processing |

One practical point concerns the GUI: SolidCraft's `ui-kit`, `look`, and `theme` components have a license distinct from the geometric code. Engine reuse and visual-identity reuse should therefore be reviewed separately when defining a distributable product. [UI license](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/LICENSE-BRAND).

The recommended starting point is the **30-cell bracket partition stored as a persistent geometric object**. It would provide a first kernel able to construct, modify, display, and analyze a perforated, filleted part while remaining trivariate.

The review covered FieldCraft at commit `35e8346`, its SolidCraft v0.42.0 pin, and the later SolidCraft architecture at commit `765f26e`. It inspected source, recipe, data structures, and test assertions; it did not execute the mesher or solver because Rust and Cargo were absent in the reviewed environment.

The initial milestone is therefore precise: construct this bracket directly with rational maps, shared faces, and IGA on the same model. It is rich enough to test the architecture while limiting the operation set to constructions that can be controlled.
