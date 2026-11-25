# Transverse bores defined by NURBS sketch curves

## Lot contract

Start from the campaign branch after `ce372ff`. Preserve the Steinmetz tests. The input contains two sketches in the XY and ZY planes. Each closed profile has two NURBS flanks and two constant-ordinate closing segments. The first sketch is removed along Z, the second along X. Homogeneous control coordinates and knots come from the recipe; rational decimal data may use a common integer scale. Stock bounds are in mm.

The flanks are clamped NURBS of degree 1 through 3, with positive weights, oriented in increasing y. In a band where multiple flanks are active, put their bases on a common knot partition and degree. Their homogeneous parameterization of y must then be identical up to scale. The flanks of one bore share a parameter interval. The two bores may start and end at different levels. Their ranges must use globally compatible parameters for the physical event ordering.

## Criteria fixed before testing

- Test multi-span rational cubic curves with an interior simple knot, varying weights, and different flanks. Compare against an independent evaluation of the input NURBS. Test rejection of desynchronization or crossing flanks and preservation of the target file.
- Check weight signs, strict growth of y, and ordering of all four boundaries in each section over the entire domain. Reference volume is the integral of material lengths in x and z multiplied together. Compare with independent quadrature to 1e-9 relative. Position and incidence tests supplement this check.
- Preserve design-curve degrees and homogeneous coefficients. Reject exact-arithmetic overflow. Distinguish capabilities actually verified in the CLI, STEP export, and analysis export.

## Constructing the crossing with smooth poles

Recipe `nurbs_cross_oval_v1` provides the stock, center, and two NURBS quarter-contours, one in XY and one in ZY. Reflections close both sections. The constructor checks reflection-axis tangents, weights, control ordering, and stock clearance. After knot insertion and exact degree elevation, both meridians share the ordinate and weight at each parameter. If their points are (a(t), e(t)) and (b(t), e(t)), their common 3D edge is exactly (a(t), e(t), b(t)) in local coordinates. Widths a and b may vary independently.

The constructor inserts lines in both sketches at the same knot event. For a single-span curve, it splits at the midpoint parameter. In each octant, it extends both meridians to an interior rectangular guard, then extends that shared surface along X or Z to the stock. The corner block inherits the guard parameterization. The two directions meet on a whole face. Each octant has two blocks per span and one corner block. Topology changes at the poles are handled by the junctions between octants.

## Measurements and reproduction

The flat-ended crossing has 26 blocks and 36 interfaces. Volume by block integration is 50,835.02450412389 mm³. Independent integration of the source curves evaluated with de Boor gives 50,835.02450413309 mm³. On each of the first four spans, seven rational points give exact rank six for the coefficients of a conic. No conic contains these spans. The oracle also checks 702 points in the material.

For the two smooth-pole profiles, the construction has 40 blocks and 68 interfaces. The degree-three IGA package has 1,588 scalar DOFs with C0 coupling. At Gauss order 24, its volume is 63,390.98047047747 mm³, compared with 63,390.98047048173 mm³ from inversion of the source curves and section integration. STEP re-import finds a closed solid with 104 faces. Adaptive OCCT mass integration at tolerance 1e-12 gives 63,390.980470481736 mm³. [STEP checks](results/nurbs-oval-step-20261125.json) · [IGA checks](results/nurbs-oval-space-20261125.json)

The explicit variant `nurbs_cross_bores_bounded_v1` accepts independently parameterized ordinates on monotone flanks with flat ends. It certifies each chord using the control-point envelope, overlays physical levels, and checks stock clearance. The requested tolerance bounds each flank's lateral displacement. A conservative upper bound on 3D boundary displacement is twice that tolerance. At tolerance 0.25 mm, the witness has 82 blocks and a measured maximum lateral error of 0.147 mm.

```sh
mkdir -p build/nurbs-demo
cargo run --locked -- build tests/references/nurbs_cross_oval.recipe.json build/nurbs-demo/oval.json
cargo run --locked -- audit build/nurbs-demo/oval.json
cargo run --locked -- export-step build/nurbs-demo/oval.json build/nurbs-demo/oval.step
cargo run --locked -- export-analysis build/nurbs-demo/oval.json 16 build/nurbs-demo/oval-space.json
```

Monotonicity and synchronization are sufficient conditions for this constructor. A sketch with multiple free loops, oblique axes, or a blind pocket intersecting the crossing requires another decomposition. In general, two rational surfaces can intersect along a curve of positive genus. Exact representation of that edge by a finite chain of regular NURBS cannot then be required. [Katz and Sederberg, *Genus of the intersection curve of two rational surface patches*](https://www.sciencedirect.com/science/article/pii/0167839688900064)

Replay the sweep of sixteen parameter sets with `python3 tests/references/nurbs_oval_sweep.py target/release/solidcraft-trivariate build/nurbs-oval-sweep.json`. Eight two-span recipes have 40 blocks; eight three-span recipes have 56 blocks. On this corpus, the maximum relative volume error at integration order 16 is 4.50e-12. [Sweep results](results/nurbs-oval-sweep-20261125.json)
