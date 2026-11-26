# Architecture and contracts

## Native Rust kernel

The `solidcraft-trivariate` crate exposes `RationalBezierMap::new` and `evaluate` from `src/lib.rs`. New production functionality goes through this crate. Models and Python tests from the early increments provide references for qualifying ports.

## Planner

For a rectangular prism with a through bore, the orthogonal planner first tries four rational maps. Each connects one full side of the rectangle to one quarter of the circle. For an equal-thickness bracket with a longitudinal fillet and one independent bore in each arm, it tries connecting both arms directly to the two bend maps. A candidate must pass wall, incidence, and Jacobian certificates. For a slender pierced plate, the next candidate encloses the bore in a local square and fills the remaining stock with affine boxes. Groups of parallel bores can have separate hosts in the same arm. For multiple bores it first tries strips separated between bore centers, followed by disjoint square hosts if no strip layout is available. An unequal-width filleted L uses rectangular bend charts. Its pierced wings can also connect to the bend with local hosts. New plans record the reason for falling back, and audit checks it on replay. The serialized model records the strategy so that audit can replay historical partitions too. [FieldCraft case](results/fieldcraft-compact.md)

The polygonal `sketch_layer_partition_v1` path takes a stock outline, cut outlines, axial material/void intervals, and explicit sketch guides. Under automatic strategy, the builder compares separator proposals. It constructs maps for each admitted proposal and verifies the exact material area of every layer, positive Jacobian on the closed domain, and complete-face incidences. It then selects the proposal with the fewest blocks. For equal block counts, it prefers the lower maximum sampled Jacobian condition number on the 3 × 3 × 3 closed grid per block. At relative differences within 1e−12 it retains proposal order. With rectangular stock it also tries strips separating convex four-sided cuts along X or Y; for one convex hexagonal cut it tries two six-spoke rings with an interior cap diagonal. These candidates must pass the same exact coverage, interface and regularity checks. The serialized plan records the method, segments and their provenance, and layer counts. [Construction](../src/sketch_layers.rs) · [Plan reader and replay](../src/sketch_layers_json.rs)

For a strictly convex stock outline around one strictly convex footprint, the builder also casts exact rays from the mean of the inner vertices through every corner of both polygons. Consecutive hits give conforming quadrilaterals. A through-hole leaves the inner polygon empty. A blind pocket needs the same inner boundary on its floor. The builder tries alternate pairings of radial triangles, then subdivides the triangles with shared edge midpoints if necessary. The trial has a combined limit of 48 input corners. A second candidate separates groups of convex cutters inside quadrilateral hosts, using opposite-edge segments and recursive splits. The chord search samples 16 intervals per host side, allows 64 recursive nodes and caps the ring at 32 vertices. Coordinates on the trial lattice use scales 1, 2, 4 and 8. Both candidates compete with the historical plans on the finished solid. [Radial construction](../src/sketch_radial.rs) · [Separator construction](../src/sketch_quad_separators.rs)

For convex polygonal stock with multiple disjoint cuts the automatic planner also separates the holes into convex hosts. A dynamic program walks each host and its hole boundary, allowing extra rational vertices on host sides. Both incident hosts use the same points on a separator. Blind cuts contribute floor cells which reuse the wall vertices. The proposal participates in the same whole-solid selection and historical replay. [Native construction and evaluation](convex-host-rings.md)

An explicit recipe `partition_strategy: "ring_planar_rerank"` or `"ring_planar_lookahead"` opts into experimental planar ring selection. Both keep block count primary and reconstruct every proposal through the exact coverage, incidence, and determinant checks. The lookahead mode adds bounded exact support offsets and shared-junction slides while preserving sketch/guide segments. Neither strategy changes recipes that omit `partition_strategy`; the planar condition score is a sampled diagnostic, not a certificate.

On a bracket assembled from local hosts, two opposing faces can cover the same plane with different extents. The planner recognizes affine rectangular traces from their homogeneous control nets, projects hanging-edge coordinates into the parent parameter intervals, and passes them through already matched faces using their exact axis permutations and reversals. It splits rational control nets in homogeneous coordinates. Different bore axes may share a wing when a slab perpendicular to both directions separates the full tools with strict clearance. The implementation allows eight closure rounds and 1,024 blocks. The plan records the local rational cut parameters, originating faces and per-round child counts. A finished candidate must pass the same whole-face and closed-domain checks as the historical plan. The independent bracket study counts each input, including construction refusals. [Exact trace discovery](../src/mapped_seams.rs) · [Subdivision](../src/mapped_cut_propagation.rs) · [3D evaluation protocol](results/propagated-bores-protocol.md)

A sketch already divided into quadrilaterals can retain its convex faces when every edge matches its neighbor. The pocket plan connects the four corners of the footprint to the stock corners. For some simple contours with an even number of vertices, another method cuts off convex quadrilaterals along visible diagonals. It then attaches each quadrilateral hole to the four corners of its host cell. Separators end at the corresponding vertices. This deterministic method can fail to find a valid decomposition. The orthogonal grid and triangulation followed by quadrangulation remain available for other configurations. [Arrangement](../src/sketch_arrangement.rs) · [Quadrilateral conformity](../src/sketch_quads.rs) · [Local diagonals](../src/sketch_quad_separators.rs)

Pocket identifiers describe the wall and layer relative to the floor. In other guided partitions, a cell is identified by its canonical rational corners and axial interval. These keys protect unchanged cells from global renumbering. Coordinate edits can change geometric keys of incident cells. Matching two plans with different topologies and general selection of the number of cells incident to each corner remain open.

Cross-axis operations and NURBS contours use dedicated paths. `import-solidcraft-sketch` reads the bracket feature history and builds the model from its sketches. `plan-sketch` makes intersection events explicit for supported NURBS profiles. Specialized templates and mechanical models provide shape, cost, and conditioning references. The pocket and C comparison uses complete-face connections and a partition shared across layers. The proposed five-volume case with T-junctions requires a different connection contract for this family. [Comparison and bounds](partition-poche-comparaison.md) · [Demo workflows](../README.md)

### Selection and editing objective

The target contract takes outlines, a 3D builder (extrusion, revolution, etc.), and planned treatments as input. Among candidates that satisfy shape fidelity, coverage, connections, and regularity on the closed domain, it favors the fewest logical volumes. For an edit, it first tries to preserve existing connectivity. At equal volume count, it favors identity stability, then compares spans, control points, and analysis cost. Bézier extraction used to certify or integrate a span does not create a new design-volume identity.

Each template must state its topology class, regularity assumptions, and knot-placement rules. For a revolved profile, the planner checks where parametric corners of the maps land on the contour after an edit. If two exterior faces of one tensor cell meet along an edge that has become smooth on the skin, the determinant vanishes along their shared edge. This rules out that particular candidate. In that case, the plan must record why it was rejected before considering a topology change. Reports distinguish an admitted candidate, the best candidate tested, and a proven minimum. This rule guides planner development.

## Published models

The gallery, standalone replay files and compact FieldCraft bundle carry models rebuilt from their source recipes. `make check-published-geometry` reconstructs them and compares them to the published data, including native STEP contents after normalization of the exchange timestamp. The separately identified SolidCore STEP is an independent CAD reference.

## Historical Python reference model

The versioned recipe is the construction intent. `trivariate.build.build(recipe)` produces an in-memory candidate model. `trivariate.audit.audit(model)` takes only that model and reconstructs expected coefficients from the declared dimensions. The CLI serializes JSON after validation to a temporary file in the destination directory; rename is atomic. Of concurrent valid publications, the last one replaces the previous file.

The closed `affine_pair_v1` family encodes two maps over `[0,1]^3`, in millimeters. Each axis has degree 1 and knots `[0,0,1,1]`. The eight coefficients `(wx,wy,wz,w)` are ordered with k varying fastest, then j, then i. All weights are 1. The `left` and `right` identities come from the two template positions. The `u=1` side of `left` and `u=0` side of `right` reference the same interface oriented by `(v,w)`.

For positive integer widths `a,b`, the maps are `(au,v,w)` and `(a+bu,v,w)`. Exact equality of all coefficients with these maps guarantees positive determinant on the closed domain, blockwise injectivity, equal traces, and a partition of `[0,a+b] × [0,1]^2`. Integer arithmetic in the prototype makes this exact for this family. The oracle returns `valid`, `invalid`, or `indeterminate` for an unknown schema.

## Planned extension

| Layer | Data and responsibility |
| --- | --- |
| BlockPlan | Resolved features, parameters with units, provenance, and expected incidences. |
| TrivariateComplex | Native rational maps, orientations, and boundary. `BlockSide` denotes a logical side, `TopologicalFace` a shared incidence, and `DesignFace` a user-selected identity. |
| Verification | Proofs per block, per interface, and for the complete partition. Record budget, rejection reasons, and bounds. |
| AnalysisSpace | Bases, refinement, trace weights, H¹ assembly, constraints, and quadrature from accepted geometry. |
| Adapters | SolidCraft supplies feature intent. FieldCraft supplies the corpus and finite-element comparisons. Rendering samples maps. |

The first curved connection covers complete faces. Sides have explicit orientation; degrees, knots, multiplicities, weights, and homogeneous coefficients must agree after reorientation and common normalization. Nested connections and partial boundaries require later qualification.

## Local references

- SolidCraft and its `solidcore` subtree provide the body interface, feature evaluation, and coordinate rules. FieldCraft provides `tests/corpus/filleted-bracket-blocks/{recipe,fixture}.json` and `crates/testkit/tests/bracket_blocks.rs`.
- AutoCraft, craft-contrib, Craft Anything, and ArtCraft document reviews, tests, commands, provenance, and shared instructions. Craft Anything and ArtCraft clones use sparse checkouts for the paths under study. The six revisions are listed in `references.lock`.
- Test references use the Python standard library. Upstream Rust applications have their own dependencies and build cycles; they will be compiled when adapter work begins.

## Implemented Rust workflow

The CLI reads versioned JSON recipes. `interchange::build` routes each family to its builder and audit. Audit rechecks homogeneous coefficients, parameter domains, and family-specific incidences. STEP export assembles oriented boundary faces and their components, including cavity shells. [Routing](../src/interchange.rs) · [Boundary](../src/boundary_shell.rs)

Loads and supports use design-face identities. The Q2 space has its own numbering and trace constraints. The dyadic plan belongs to the analysis job. The lifting operator removes dependent equations before solve. Results retain the plan, blocks, nodal loads, displacements, reactions, and work balance. [Mechanical jobs](../src/mechanical_job.rs) · [Constraints](../src/constraint_system.rs)

Job v5 controls quadrature on rational faces. With `derivative: "freeze_primal"`, differentiation reuses each face's primal subdivisions and counts extra panels against the budget. Job v6 accepts a spatially affine body force and a `body_quadrature` rule with minimum and maximum orders, relative tolerance, and point cap. The estimate compares the complete nodal load, volume, resultant, and moment at two orders. Counts expose rejected orders and selection cost. Convergence between rules remains a numerical indicator. [Surface](../src/adaptive_traction.rs) · [Volume](../src/volume_load.rs)

A geometric interval-edit session preserves layers whose geometry and material mask remain compatible. The bracket mechanical session caches local contributions with homogeneous controls, integration parameters, and material properties. Both publish a new state only after admission checks. Reports give each experiment's cache scope, cold checks, and timings. [Geometric editing](../src/interval_bores_session.rs) · [Mechanical editing](../src/bracket_mechanics_session.rs) · [Campaign report](bilan-campagne-20261125.md)
