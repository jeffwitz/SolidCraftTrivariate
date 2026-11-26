# Exact axial knot removal and independent field evaluation

Two extrusion knots at the box end planes can be removed from the geometry of the tall cylinder. The same planes remain available as source subfaces in the union. The Rust assembler rebuilds the tensor control net from the source Bézier maps and a layout that supplies degrees, cell placements and orientations. It elevates degrees up to three, matches homogeneous traces up to positive scale, and supports exact Bézier extraction on the assembled C⁰ grid. The tall-box Boolean release pipeline runs this assembler on its proposed layouts and rechecks the published nets. Layout selection remains a separate step.

| Geometry | Design volumes | Geometry spans, before → after | Stored homogeneous controls, before → after | Retained interpatch subfaces |
|---|---:|---:|---:|---:|
| Cylinder, core and closed ring | 2 | 15 → 5 | 108 → 54 | 12 |
| Cylinder, core and two open strips | 3 | 15 → 5 | 116 → 58 | 18 |
| Box–cylinder union | 7 | 19 → 9 | 172 → 106 | 30 |

The two-volume cover also retains three self-seam records. The checker reinserts the original axial knots using rational arithmetic and compares all four components of every original control. It also checks rational positions, patch ownership, exterior faces, interfaces, source spans and field knots against the original cover. The current reduction applies to a clamped degree-one extrusion direction identified by the source parameter map. The compact adapter accepts signed 64-bit integer knots and homogeneous coefficients and requires positive weights. [Geometry measurements](geometry.json)

## Cylinder torsion

The independent elasticity calculation evaluates positions, Jacobians and load geometry from the compact two-volume cover. It keeps the original field knots, rational weights, 15 integration cells and trace constraints, including the closed ring. With order-16 quadrature, both runs have 204 vector degrees of freedom, including 153 free. Both give a strain energy of 0.01767145867653 N·mm.

| Quantity | Relative change |
|---|---:|
| Stiffness, Frobenius norm | 2.58 × 10⁻¹⁶ |
| Load, Euclidean norm | 0 |
| Solved displacement, Euclidean norm | 5.37 × 10⁻¹⁴ |
| Strain energy | 3.89 × 10⁻¹⁵ |

The comparison tests the geometry change at fixed discretization. No reduction in solve cost or field dimension was measured. The union has been checked geometrically. [Torsion measurements](torsion.json) · [Separate compositional Boolean experiment](../composed-axial/README.md)

## Reproduce

Run `make check-axial-knot-reduction` for the exact geometry replay or `make check-cylinder-iga` for the geometry and mechanical replays (NumPy and SciPy required). Generated layouts, assembled covers and compact geometry are written under `build/axial-knot-reduction/`.

The public CLI operations are `assemble-patches MODEL LAYOUT OUTPUT` and `reduce-axial-knots PATCHES OUTPUT`. The compact document stores original subfaces and source spans alongside separate geometric knots and original analysis knots. The reference field additionally uses the original cover's weights and global trace numbering.
