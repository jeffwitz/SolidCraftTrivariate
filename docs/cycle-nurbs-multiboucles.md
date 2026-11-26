# Multi-loop NURBS sketches and crossed extrusions

Historical campaign record. For the current sketch planner, models and reproduction commands, use the [documentation index](README.md).

Starting point: `81cb61d` on `campagne-nurbs-multiboucles`. This cycle follows the completed ten-hour campaign and covers multi-contour sketches, automatic partitioning, pockets and IGA spaces.

## M1 — strip contours, through-hole and pocket

The `nurbs_multiloop_v1` family accepts `bounds_xy_mm:[[xmin,xmax],[ymin,ymax]]`, `height_mm` and `loops`. Each loop has an `id`, `lower` and `upper` curves in `nurbs_sketch_strip_v1` format, and `z_interval_mm:[lo,hi]`. Each pair of curves is joined at both ends by vertical segments. The outer contour is rectangular; interior loops must lie strictly inside it and have strictly separated x-projections. This initial domain includes non-conical side profiles. The planner orders physical events and extends separators through free intervals; cells follow the layers induced by all cut starts and ends. Piece names do not affect the partition.

The Rust constructor is in `src/nurbs_multiloop.rs`; the CLI integrates construction and persistence auditing. An independent geometric oracle evaluates the source B-splines at their original knots and integrates cross-sectional areas to compute removed volumes by axial interval. Acceptance criteria: STEP volume within 1e-8 relative, wall positions within 1e-10 mm, and agreement on inside/outside/boundary classifications. Prospective cases cover the same loop as a through-hole and blind pocket, a second differently shaped loop at another level, and input permutation/renaming. Check material invariance, pocket-floor and wall provenance, and target preservation on refusals. Reject contours touching the stock, overlapping x-projections and degenerate intervals. Report planned decisions, admissions and refusals separately.

The IGA oracle checks face DOF identity, affine-field and gradient reproduction, and equilibrium of constant stress through volume and boundary forces. The relative residual threshold is 1e-8 at order 16. Numerical volume or conditioning estimates are tied to their quadrature points and order. Run `make check` and strict Clippy before integration; record revisions and measured results in the journal.

## Agent roles

Session metadata and the effective model identify three distinct `gpt-6-luna` sessions.

| Role | Session |
|---|---|
| Rust construction | `01a21208-822a-7f02-8a5d-a90895f7a829` |
| Geometry and STEP oracle | `01a21208-b7a7-79d3-93c9-bdc4262b1012` |
| IGA oracle | `01a21208-e10c-78d1-a96e-9423ec00f071` |

## First prospective profile

The initial second-loop profile had three successive poles at x = 70 mm. Its x-derivative vanishes at the interior knot t = 2, so it fails the regularity criterion. The original recipe is retained in `tests/counterexamples/nurbs_multiloop_stationary.json`. The replacement uses abscissae 55, 62, 67, 70, 74 and 80 mm. The independent integral predicts 112940.85164966024 mm³ of material for the replacement mixed recipe; preregistered thresholds are unchanged.

The constructor was integrated at `8025838`. On an index snapshot, `make check` ran 395 Rust and 58 Python tests; strict Clippy passed. The mixed model has 45 blocks, 75 interfaces and 120 exterior faces. Three isolated tamperings (homogeneous control, loop name and face label) each produce `multiloop_stored_partition_mismatch`.

### M1 qualification

The through-hole, pocket, combined case, and permuted/renamed contours produce 12, 27, 45 and 45 blocks, respectively. Each STEP reread by OCC 7.8.1 is a valid closed solid. The mixed-case volume is 112940.85164966017 mm³ versus 112940.85164966024 mm³ from direct integration of the sketch curves. All 1,515 classifications agree; the largest measured deviation at the nine checked points per face is 2.85e-14 mm. All four invalid inputs preserved the prior target byte for byte. The first volume-check attempt evaluated Gaussian points outside the unit cube; both reports preserve the diagnosis and correction. [STEP qualification](results/nurbs-multiloop-delivery-20261125.json) · [Initial diagnostic](results/nurbs-multiloop-delivery-initial-20261125.json) · [Profile reference](results/nurbs-multiloop-geometry.json)

The mixed case has 1,800 scalar Q3 degrees of freedom, 1,200 interface nodal matches and 2,880 affine-field/gradient checks. The relative equilibrium residual at q16 is 2.65e-12. A historical test applied weights twice to homogeneous coordinates; after correction its residual is 1.96e-12. Both M1 oracles used the binary with SHA-256 `72d0b6163bcd2eb5468b8783b9abaf52fc562cd87c5fffe7cc821889c9304040`. [IGA qualification](results/nurbs-multiloop-iga.json) · [Historical replay](results/nurbs-patch-equilibrium-20261125.json) · [Earlier double-weight diagnostic](results/nurbs-patch-equilibrium-20261125-double-weight-diagnostic.json)

## M2 — assemble loops from independent edges

Constructor baseline: `8025838`. A `nurbs_edge_loops_v1` recipe contains the M1 stock and height, an unordered `edges:[{id,curve}]` array and `cuts:[{id,seed_edge,z_interval_mm}]`. Each cut names one contour edge; the constructor must recover the remaining edges, both side profiles and their directions from exact homogeneous endpoints. Initial scope: disjoint components in x, each a simple four-edge loop with two linear vertical closures and two NURBS side curves synchronizable as in M1. Edge and operation names or order must not determine the partition.

Planned files were `src/nurbs_edge_loops.rs`, family integration, targeted tests and recipes/oracles under `tests/references`. Prospective cases: two free loops, a top-open pocket and through-hole, permuted edges, reversed orientations, changed identifiers and changed gauges. Compare independent volume within 1e-8 relative, wall positions within 1e-10 mm, pocket floors and conforming traces; check affine fields and Q3 equilibrium at q16 within 1e-8. Reject open chains, branches, oblique closures and two operations on one component, preserving the prior target on every refusal. Interior edge intersections and arbitrary edge counts are out of scope. At Jeff’s request, M2 is suspended after three targeted Rust tests and construction of a 45-block prospect; independent qualification remains open.

## S1 — rectangular plate with a closed G2 NURBS hole

Starting revision: `58c4409`. The priority is junction verification. Dedicated paths: `src/g2_plate.rs`, `tests/g2_plate.rs`, `tests/references/g2_plate*` and `docs/results/g2-plate*`. The three Luna roles above cover Rust construction, B-rep comparison and IGA/curve-jet oracles.

The prospect is an extruded plate with two through-cuts: a sharp-cornered rectangle and a periodic, rational, cubic, non-conic B-spline with simple knots and positive weights, without an added closing segment. All representations use the same recipe and poles. At every smooth junction, including closure, check position (1e-10 mm), oriented unit tangent (1e-10) and geometric curvature vector (1e-10 mm⁻¹). A G1-but-not-G2 witness must pass the first two checks and fail the third. Also measure junctions between extruded faces: normal and cross-sectional curvature, plus curvature along generators. Plate edges remain sharp, so the assembled analysis space is C0.

First use SolidCore’s actual B-rep backend on the sketch and extrusion, recording revision, command, geometry-access path and outcome. OCCT is an allowed fallback only after a documented SolidCore failure. Compare volume to independent curve integration (1e-8 relative), walls to the extruded curves and classifications in both voids and the material. Check traces, affine reproduction and Q3 equilibrium at q16 (1e-8 relative). Preserve any failure of regularity, closure or comparison with its details.

## S2 — intersecting through-cuts along X and Z

Starting point: `755f6d3`. S2 is the primary demonstration; S1 is a parallel control. The stock is 100 × 50 × 12 mm. The S1 periodic profile cuts through along Z. A rectangular profile in YZ, y = 719/36 to 1081/36 mm and z = 4 to 8 mm, cuts through along X. The tools’ intersection must be positive and smaller than either tool. Their four intersections with constant-y planes lie inside spans at parameters 2.5, 3.5, 6.5 and 7.5. Maps split at these events must cover retained material in all three Z layers.

Dedicated files: `src/g2_cross_plate.rs`, `tests/g2_cross_plate.rs`, and recipes, oracles and reports prefixed `g2_cross` or `g2-cross`. Compare each tool, their intersection and union, and remaining material to independent volumes within 1e-8 relative. On the 3D intersection, distinguish four generators and horizontal arcs at z = 4 and 8. Probe both sides of seams and exposed faces. Retain S1’s smooth-jet and affine-reproduction/equilibrium thresholds; measure normal discontinuity where walls meet. Replay reversed and cyclically shifted profiles, knot insertion and affine knot rescaling; report every admission/refusal and geometric invariance. The exploratory y = 20 and 30 mm case has intersections without known dyadic roots: record an indeterminate result or refusal and verify target preservation.

Construction owns material maps and certification. Geometry first runs the two S2 cuts through SolidCore and records its command, revision and resulting B-rep; OCCT is permitted only after a documented failure. The jet/IGA role compares source curves, represented walls and junctions, then checks the C0 space. S1 results stay separate. The S2 figure is reserved for the primary thread.

The public contract is in `tests/references/g2_cross_sketches.recipe.json`. Both cuts use the same NURBS-sketch format, frame, direction and length. The YZ rectangle is a closed degree-1 loop with homogeneous poles; the XY profile is periodic with variable weights. Before extrusion, `plan-sketch` must expose jet-classified junctions, multiplicities, closure, induced events and provenance IDs. `g2_cross_plate.recipe.json` preserves the originally agreed numerical dimensions for oracle replay. A higher-degree test inserts a straight subdivision; the four corners must remain the only true corners. Native primitives provide partition data through their creation template; imported curves enter through geometric curve checks.

### S2 construction results

The main model has 32 blocks, 46 interfaces and 100 boundary faces. The initial audit gives a common determinant lower bound of 5.7298625 mm³. Order-32 map integration gives 18836.7560490, 15526.8831762 and 18836.7560490 mm³ in the three layers. The four seam generators span z = 4–8 mm; horizontal arcs lie at both ends. The walls meet sharply along these seams. [Map-volume integration](results/g2-cross-map-volume.json) · [Maps and curvature figure](figures/plaque-deux-trous-g2.png)

The 32 retained NURBS wall faces provide 64 exact endpoint-jet checks. Fraction-based comparison with de Boor derivatives recovers position, oriented tangent and curvature at twelve distinct positions. The top and bottom layers each retain twelve smooth junctions. In the middle layer, six junctions remain smooth and four endpoints meet the rectangle. The face generators are straight. [Exact trivariate-wall jet checks](results/g2-cross-wall-continuity.json) · [G1 witnesses and G2 rounded rectangle](results/g2-plate-jets.json)

The first B-rep engine tested was SolidCore 1.0.1. Two rational sweeps followed by kernel booleans give tool volumes 3489.731852993524 and 4022.222222222223 mm³, intersection 712.3493494134794 mm³ and union 6799.604725802266 mm³. The 60 × 50 × 12 stock minus this union has 22 faces, 66 edges, 44 vertices and volume 53200.39527419775 mm³. Trivariate-map quadrature gives 53200.395274197734 mm³. B-rep classification and differential properties are assessed separately. [SolidCore test](results/g2-cross-brep-solidcore-20261125.json)

The native recipe defines the rectangle once as a primitive. Its template generates and persists the NURBS contour, marks all four corners and supplies cutting levels. Imported curves instead pass through span extraction and geometric-jet checks. `plan-sketch` exposes events before `build`. All six primitive pairs (rectangle, circle or ellipse × extrusion or revolution) pass construction, audit, STEP export and analysis-space export. [Primitive-template catalogue](gabarits-primitives-nurbs.md) · [CLI run](results/primitive-sweep-cli.json) · [NURBS sketch prepartitioning](prepartition-esquisse-nurbs.md)

Revision `23d1810` reproduces the delivered model bytes and the STEP DATA section. Reparameterizing knots by t ↦ 7t + 11 preserves block coefficients, topology and design faces exactly. Five tested variants (reversal, cyclic shift, knot insertion, elevated rectangle and straight-side subdivision) pass audit; knot insertion increases the block count to 35. Order-32 volumes match the original. Ordinary y = 20 and 30 mm levels have four real intersections at irrational parameters; the CLI rejects them and preserves the output SHA-256. [Variants, roots and IGA](results/g2-cross-junctions.json) · [Delivery check](results/g2-cross-delivery.json) · [Model](results/g2-cross-plate.model.json) · [Trivariate STEP](results/g2-cross-trivariate.step)

The assembled C0 Q3 space has 1,360 scalar DOFs. Its 736 trace matches and affine-field reproduction at 2,048 points pass. Relative equilibrium residual at order 16 is 2.87 × 10⁻¹⁵. The spatial certificate concerns Jacobian sign on the blocks. Global coverage of the middle layer is supported by layer volumes and independent solid comparison. Generalization to arbitrary contour arrangements remains open.

Numerical map inversion produces 40 classifications consistent with contour parity and both cuts, including the seams and four offsets on each side of every generator and arc segment. None is indeterminate. [Coordinates, classifications and method](results/g2-cross-map-seams.json)

At the native-solid checkpoint, 284 of 296 classifications agree. Both methods place (75, 34, 2) on the wall. In the twelve retained disagreements, SolidCore reports “on boundary” while the oracle reports material; they lie near the rectangle-plane crossings. Requested tolerance is 2 × 10⁻⁷ mm, so pointwise qualification fails on this set. The earlier test, which also included four oracle errors at extrema, remains archived. [296 native classifications and jets](results/g2-cross-brep-probes296.json) · [Initial diagnosis](results/g2-cross-brep-classification-initial-20261125.json) · [SolidCore STEP](results/g2-cross-solidcore.step)

Thirty-two SolidCore B-rep jet measurements bracket the eight periodic knots on both sides at z = 2 and 10 mm. Maximum position difference is 3.05 × 10⁻¹⁴ mm; maximum difference between source curvature and signed normal curvature times the measured normal is 3.77 × 10⁻¹⁴ mm⁻¹. The minimum dot product of the oriented source tangent and normalized surface derivative is 0.9999999999999998. These numerical probes 10⁻⁵ in parameter distance from knots complement rational checks at the delivered map endpoints. Revision `23d1810` passes 415 Rust tests, 58 Python tests, `make check` demonstrations and strict Clippy.

### Classification and STEP reread follow-up

A supplementary sweep of 96 points uses four offsets along the profile normal and the classifier’s requested tolerance (2 × 10⁻⁷ mm). At 10⁻⁵ mm, SolidCore marks all 12 material and 12 void points as boundary. At 10⁻⁴ and 10⁻³ mm it marks 4 of the 12 expected-void points as material; at 10⁻² mm it marks 3 of 12 as material. These observations do not support attributing all disagreements solely to classifier tolerance. [Sweep and coordinates](results/g2-cross-brep-solidcore-20261125.json)

SolidCore built the B-rep and performed its booleans. An independent documentary reread of its STEP with FreeCAD/OCCT finds a valid closed solid with 22 faces and volume 53200.395273882794 mm³. It agrees on 293 of the original 296 classifications. Three wall points at z = 8.00001 mm are reported outside; the minimum distance to reread faces is 10⁻⁵ mm. [Documentary reread](results/g2-cross-solidcore-reread-20261125.json)

The eleven sweep points classified as material by the native kernel but void by the oracle were also checked against the STEP reread in FreeCAD/OCCT. All eleven are void there at the original tolerance. The native-classifier disagreement remains reproducible. [Targeted check on the same STEP](results/g2-cross-solidcore-reread-20261125.json)

SolidCore’s reader loaded both STEP files as one solid, one shell and two regions. It measured 53200.39527419774095 mm³ in each; five probes (material, each individual hole, their intersection and the boundary) agree. For the trivariate STEP, the reader reports reversing 50 faces to reconcile loop direction with `same_sense`. Follow-up on orientation and native classifications is suspended at Jeff’s request. [SolidCore reread](results/g2-cross-trivariate-solidcore-reader-20261125.json)

## S3 — blind pocket on the longitudinal face of S2

Starting point: `441ec30`. The reference input is the delivered S2 model, SHA-256 `19bda241a12bb0cae9e87ad7fa2e02a3e71c45be119115aec9fd8aeb0b48847c`. The sketch attaches to source `rect_x/zmin` on z = 4 mm. The planner extension permits a loop to overlap the stock boundary and the NURBS hole. A rectangular pocket cuts 2 mm deep and ends in material at z = 2 mm. Partitioning uses the maps and existing events. The three Luna roles cover Rust construction/CLI, an independent source oracle, and map/IGA checks.

Before qualification, require strictly positive footprint area in both the existing material and void. Measure the truncated old seam, retained seam, new floor, exposed walls and incremental volume separately. Refusal for an unsupported face, invalid depth or out-of-domain contour must preserve the original model. Thresholds: relative volume and Q3-to-q16 equilibrium 1e-8, position 1e-10 mm, and exact agreement with declared interfaces. A positive Jacobian certificate establishes local regularity only; STEP exports the resulting V-rep.

The prospect was frozen before calculation. The XY sketch frame has origin (0, 0, 4) mm, +X/+Y axes and cut direction −Z. The control rectangle spans x = 0–75 and y = 719/36–25 mm; it opens at x = 0, while its x = 75 edge lies in the NURBS void on this band. The floor remains at z = 2. On the left NURBS wall, the old seam portion from y = 719/36 to 25 must move to z = 2; the portion from y = 25 to 1081/36 remains at z = 4. Footprints in material and void, added volume, old seam and relocated seam are separate checks.

For the principal prospect, the team chose a rectangle strictly inside the stock while retaining the x = 0 rectangle as a control. It spans x = 55–75 and y = 5893/266–25 mm. Its lower y bound is the source NURBS point at t = 13/4, whose x-coordinate is 2513/38 mm. The pocket crosses the old seam portion between parameters 3 and 13/4. Construction adds an x = 55 separator and subdivides the span. Measure the old floor at x < 55 and the remaining arc after the cut. Count the open-side witness separately.

The source-curve oracle predicts 31.4269014385 mm² of material footprint and 25.4903917946 mm² already in the void; expected incremental volume is 62.85380287694 mm³. The old left seam is 10.1198358902 mm long. The two portions remaining at z = 4 measure 5.0599179451 and 2.2099209957 mm; the portion relocated to z = 2 measures 2.8499969494 mm. The laterally open control removes 665.12299798998 mm³. These values come from the source contour before comparison with S3 maps. [Source oracle and recipe](results/g2-face-pocket-geometry-20261125.json)

### S3 delivery

Revision `44cc7dc` audits the embedded S2 model and decodes its maps and design faces. Support is determined from faces carrying floor provenance. On the left wing, exact rational roots of Y − yW determine contour subdivisions; Bernstein coefficients on the right wing certify that the opposite rectangle edge stays in the void over the whole strip. The constructor propagates subdivisions through the maps, cuts to the requested depth and labels new faces. This family accepts a rectangle on the floor’s XY plane, extrusion along −Z and positive integer depth that leaves the floor in the lower layer.

| Result | Inner pocket | Open-side control (x = 0) |
|---|---:|---:|
| Blocks / interfaces / exterior faces | 93 / 191 / 176 | 86 / 177 / 162 |
| Expected removal (mm³) | 62.85380287694 | 665.12299798998 |
| Old left seam retained at z = 4 (mm) | 7.26983894078 | 5.05991794508 |
| Seam relocated to z = 2 (mm) | 2.84999694939 | 5.05991794508 |

For the inner pocket, order-32 map integration gives 53137.54147132078 mm³ remaining. The three new walls cover 49.95488721805 mm² and the floor 31.42690143847 mm². The x = 75 side lies in the void and creates no wall. On both models, surface areas agree within 2 × 10⁻¹³ mm² and arc lengths within 4.5 × 10⁻¹⁵ mm of the source oracle. [Face and arc comparison](results/g2-face-pocket-faces.json) · [Model sections](figures/s3-face-pocket-sections.png)

Changes to x and depth, plus a new cut at t = 25/8 (y = 48689/2070 mm), pass independent volume comparison. Seven refusals, including a modified parent and a straight edge extending beyond the void, preserve the previous model byte for byte. The first candidate recognized the recorded dimensions but failed the open-side witness; its initial diagnosis is retained. [Variants](results/s3-pocket-variants.json) · [Refusals and prints](results/s3-pocket-refusals.json) · [First open-side refusal](results/s3-first-open-side-refusal.json) · [First variant refusals](results/s3-first-variant-refusals.json)

The HEX8 pocket mesh at the z = 4 level has 7,422 nodes and 5,952 cells. Determinants are checked over each full cell using the rounded VTK coordinates. The common NURBS-map regularity lower bound is 2.18330380277 mm³. HEX8 checks and the map certificate concern their respective representations. [HEX8 report](results/hex8-witness.json) · [Public gallery](https://jeffwitz.github.io/SolidCraftTrivariate/) · [Gallery Chrome check](https://github.com/jeffwitz/SolidCraftTrivariate/actions/runs/45606030882)

[S3 recipe](../tests/references/s3_face_pocket.recipe.json) · [Exact model](results/s3-face-pocket.model.json) · [STEP](results/s3-face-pocket.step) · [Open-side control model](results/s3-face-pocket-open-side.model.json) · [Delivery probes](results/s3-pocket-delivery.json)

The main C0 Q3 space has 3,349 scalar DOFs. Its 191 interfaces give 3,056 exact trace matches. Affine field and gradient reproduce at 5,952 points; relative equilibrium residual at q16 is 2.95 × 10⁻¹⁵. All nine targeted classifications on new walls and old seams agree. All 864 probes in the old maps remain in the retained material; these probes do not measure removed volume. Relative order-32 volume errors are 2.74 × 10⁻¹⁶ overall and at most 3.89 × 10⁻¹⁶ per layer. `make check` passes 426 Rust tests, 58 Python tests and demonstrations; strict Clippy passes all targets. [IGA checks and classifications](results/g2-face-pocket-iga.json)

The open-side control has 86 blocks, 3,094 scalar Q3 DOFs, 2,832 trace matches and 5,504 affine value/gradient checks. Its relative q16 equilibrium residual is 2.71 × 10⁻¹⁵. All eight edge probes agree, including the opening and old cross-sectioned edge. The 846 probes in old maps outside the removed envelope remain in material. Order-32 remaining volume is 52535.27227620775 mm³. [Open-side IGA report](results/g2-face-pocket-open-side-iga.json)

The public deployment at `870fcf6` includes the constant-axial-parameter section rendering and final figure. The deployed renderer and figure match the versioned files. Desktop and mobile Chrome checks pass. [Pages deployment](https://github.com/jeffwitz/SolidCraftTrivariate/actions/runs/45606991235) · [Browser check and captures](https://github.com/jeffwitz/SolidCraftTrivariate/actions/runs/45606964966)
