# Separation chords and ring quality

The external `mixed_sizes` sketch contains a hexagonal through-hole and two triangular through-holes. The reference release reproduces the model supplied by the reviewer byte for byte. Its worst planar Jacobian condition number occurs at the triangle corner (112, 75) mm. An internal separator junction and the adjacent hole edge form an angle of 179.79° there.

![The same sketch and its partition before and after a shared junction moves](mixed-sizes.svg)

| Strategy | Volume blocks | Maximum sampled planar κ |
| --- | ---: | ---: |
| Production default | 14 | 559.00 |
| Rerank the original ring candidates | 14 | 559.00 |
| Explore additional chords and slide an unpinned junction | 14 | 45.48 |

The winning partition preserves the three input contours and the material volume of 269,820 mm³. The certified lower bound on the volume-map Jacobian determinant rises from 18.74 to 507.96 mm³. The minimum sampled corner sine rises from 0.00371 to 0.06556. The planar condition numbers are numerical measurements on the material maps; the positivity bound is a separate continuous certificate.

The experimental selector first minimizes the number of material blocks. At equal count it compares ring maps by their worst planar condition number on a closed 3 × 3 grid, then tries bounded rational displacements of constructed junctions along a straight edge. Every face incident on the junction receives the same new point. Source vertices and guide endpoints are fixed. An independent evaluation samples the exported volume maps on a 5 × 5 × 5 grid.

The corpus includes the 45 original sketches, 12 published supplemental sketches, and the 11 sketches from the independent reviewer. Only `mixed_sizes` from the final group was used to develop the correction. The two suites have different independent checks and are reported separately.

See the [decision protocol](protocol.md) and [production baseline audit](../polygon-extrusion-quality/README.md).
