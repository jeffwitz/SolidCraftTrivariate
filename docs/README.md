# Documentation

Start with the [project README](../README.md) to build the Rust CLI and run a sketch through planning, construction and audit. The [interactive gallery](https://jeffwitz.github.io/SolidCraftTrivariate/) lets you inspect published models and download their inputs.

## Design and implementation

| Question | Where to look |
| --- | --- |
| Which original chart partitions support edge treatments? | [FieldCraft and S3 fixed-partition rim experiments, qualification and STEP files](results/fixed-rims/README.md) |
| How are crossed-tool Boolean results reproduced? | [Constructor and oracle](crossed-tool-boolean.md), [interface mechanics](results/crossed-interface-iga/README.md) |
| What are the requirements and acceptance criteria? | [Kernel specification](kernel-requirements.md) |
| How does a sketch become a volumetric model? | [Architecture](architecture.md), [partitioning from construction lines](partition-esquisses.md), [prepartitioning NURBS sketches](prepartition-esquisse-nurbs.md) |
| Which native profiles and imported curves are supported? | [Primitive templates](gabarits-primitives-nurbs.md), [NURBS strip extrusion](nurbs-sketch.md), [crossed NURBS profiles](nurbs-cross.md), [intersecting cylinders](steinmetz.md) |
| How is a partition selected? | [Current convex-host planner and measured results](convex-host-rings.md), [pocket partition comparison and lower bounds](partition-poche-comparaison.md), [curved sketch corpus](galerie-esquisses-courbes.md), [FieldCraft bracket comparison](results/fieldcraft-compact.md) |
| Which properties have certificates? | [Planar maps](planar-certificate.md), [annular maps](annular-certificate.md), [early interface contract](L1.md) |
| Which edge treatments preserve the block topology? | [Four-block annulus](annulus-four-block-topology.md), [annulus edge states](annulus-original-edge-states.md), [bracket study](results/bracket-native-feasibility.md), [S3 pocket study](results/s3-native-edge-feasibility.md) |

## Experiments and reproduction

Run `make check-published-geometry` to rebuild the gallery downloads, standalone native outputs, FieldCraft bundle and reference recipes with the current kernel. [Reference recipe inventory](results/reference-recipe-replay.md)

The [planar quality report](results/polygon-extrusion-quality/README.md) separates in-plane distortion, near-flat corners, and extrusion aspect on the 45 published polygon models. The [current polygon study](convex-host-rings.md) covers 45 corpus cases and 12 supplemental sketches. To refresh all 68 models, build the release CLI and run `python3 scripts/gallery_catalog.py`. In a Python environment with the FreeCAD modules available, `python3 scripts/check_gallery_geometry.py --read-step` additionally verifies that every native STEP imports as one valid solid with closed shells. The [gallery regeneration guide](site/README.md) describes the inputs and published files.

- The [FieldCraft replay guide](results/fieldcraft-replay/README.md) links the downloadable bundles. The [STEP comparison](results/fieldcraft-kernel-replay.md) and [PA5 energy-versus-DOF study](results/fieldcraft-bracket-pa5-mechanics.md) document different checks on the bracket.
- The [local fillet cost report](results/local-fillet-cost-20261125.md) records the cost of refinement. The [directional quadrature protocol](lot-m12-oracle.md) and [cylinder qualification protocol](lot-p6b-oracle.md) define independent numerical checks.
- The SolidCraft import protocols cover [polygon extrusions](sc5-protocol.md), [circular bores](lot-sc6.md), [oriented sketch planes](lot-sc7.md), [blind bores](lot-sc8.md) and [counterbores](lot-sc9.md). The corresponding kernel protocols cover [interval bores](lot-g21.md) and [counterbores](lot-g23.md).

## Project history and contributor guidance

The [campaign summary](bilan-campagne-20261125.md) gives the outcomes of the initial ten-hour run. The [chronological log](campagne-20261125.md) explains its intermediate results and failures. The [detailed run ledger](campaign-run-ledger.md) retains the individual protocols, implementation revisions, budgets and outcomes. Later work on sketch loops, crossed extrusions and face pockets is recorded in the [multi-loop campaign](cycle-nurbs-multiboucles.md).

The [research positioning](positionnement-campagne-20261125.md), [development workflow study](archive/development-workflow-study.md) and [FieldCraft and SolidCraft design study](archive/fieldcraft-design-study.md) document the project's starting assumptions. For current development conventions, see the [repository instructions](../AGENTS.md) and [Luna guide](LUNA.md).

The [blind sketch-planner evaluation protocol](sketch-planner-blind-evaluation.md) includes the paired comparison against revision `e5bd5c6`. The [bracket propagation study](results/propagated-bores.md) records block counts, map quality and independent STEP readback for ten controls.
