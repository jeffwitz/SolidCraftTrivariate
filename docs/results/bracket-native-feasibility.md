# Fixed-topology native edge treatment: bracket witness

## Scope and starting data

- Starting revision: `d25b6b6e17f1c95b8371a9ae714298fd87861677`.
- Files changed by this investigation: this report only.
- Live source model: `docs/site/data/models/filleted-bracket.model.json`.
- Native source family: 36 blocks, 55 complete-face interfaces, 106 exterior block sides; the fixed-family bracket audit checks these counts.
- Independent sharp source and edge oracle: `tests/references/calculix_bracket.py` and `tests/references/edge_treatments_brep.py`.
- Operation: the independent all-edge fillet uses radius 0.25 mm. This note tests whether the same block and side incidences can represent that operation. It does not replan or modify the source model.
- Admission thresholds: retain the 36/55/106 block-side inventory, keep every block Jacobian strictly positive on its closed cube, and match the independent all-23-edge treatment. The local smoothness condition is equality of tangent planes along the treated edge.
- Inputs rejected by this scope: changing a block identity, face incidence, or exterior-side identity; a zero-Jacobian block; or treating an untested edge as covered.
- Exit criterion: a native construction and audit for all 23 selected edges, or one exact fixed-topology obstruction proving that the all-edge result is inadmissible under the unchanged inventory. The latter is established below.

## Identified same-block physical crease

Block `arm_x_30_50_z0_10` has degree `[1,1,1]` and the exact Cartesian extent
`x=[30,50]`, `y=[0,8]`, `z=[0,10]` mm. Its homogeneous control weights are all one.
The two block sides

- `axis=1, upper=false` (`v=0`, plane `y=0`), and
- `axis=2, upper=false` (`w=0`, plane `z=0`)

are both exterior sides. Their common block edge is exactly
`(x,0,0)`, `30 <= x <= 50`, of length 20 mm. An exact whole-face comparison of this JSON against every other block face found no matching face for either side, confirming that both are exterior rather than internal interfaces.

This subedge is part of independent sharp-source crease #2: the oracle records the complete physical edge `(x,0,0)`, `0 <= x <= 60`, length 60 mm. At its midpoint `(30,0,0)`, the all-edge fillet STEP has shell distance `0.10355339059327377` mm. This equals the distance from the sharp corner to a radius-0.25 mm quarter-circle centered at `(30,0.25,0.25)`:

`0.25 * (sqrt(2) - 1) = 0.10355339059327376... mm`.

The oracle output records one valid, closed solid after STEP re-import, 44 faces, and volume `32134.92052955058 mm^3`. Source edge #2 is classified as a crease; the witness confirms that this crease is displaced by the fillet at the same location as the native block subedge.

## Fixed-incidence obstruction

Let `F(u,v,w)` be the regular trivariate map of `arm_x_30_50_z0_10`, with axes aligned to x, y, z. Along its common edge `v=w=0`, the two exterior boundary faces have tangent planes

`T_v = span(F_u, F_w)` for `v=0`, and `T_w = span(F_u, F_v)` for `w=0`.

If the edge is made smooth while these two sides remain the incident physical boundary faces, their tangent planes must agree. Since the edge tangent `F_u` is nonzero, equality implies `F_v` lies in `span(F_u,F_w)`, hence
`det(F_u,F_v,F_w)=0` on that edge. This contradicts the positive-Jacobian requirement for a regular block map on the closed parameter cube.

Therefore this particular physical crease cannot become a smooth fillet edge while retaining this block's two adjacent exterior sides as the same incident patches and retaining a regular trivariate map. Keeping these incidences permits an alternative experiment with a singular boundary parameterization, as on the original annulus. A separate cylindrical boundary patch with altered incidence is another construction. Neither alternative is a certified complete model of this bracket. Knot insertion alone cannot change the fact that these two cube sides meet at one edge.

This is a local, exact obstruction for the requested fixed topology. It does not rule out a new conforming partition, nor does it establish that every one of the other selected edges has the same obstruction.

## Coverage and disposition

The independent OCCT source has 23 selected physical creases out of 27 enumerated edges (2 are periodic seams and 2 are smooth). Its all-edge radius-0.25 mm fillet is a valid closed STEP solid. Against that denominator, this investigation establishes:

- 1/23 selected edges has a same-block, fixed-incidence obstruction, witnessed by source edge #2 and `arm_x_30_50_z0_10`;
- 22/23 selected edges have no individual native feasibility result in this scoped pass;
- 0/23 selected edges has a qualified native fillet construction in the existing block/incidence topology.

Because one required edge is obstructed under the fixed topology, the all-23-edge native operation is not admissible with the existing blocks and incidences. The independent OCCT result is evidence for the target geometry only; it is not a native trivariate construction or proof of native coverage.

## Reproduction and limits

The block and sides are read directly from `docs/site/data/models/filleted-bracket.model.json`. The fixed source inventory is enforced by `src/bracket_audit.rs::audit_bracket`. The independent edge classification, all-at-once fillet, STEP re-import, and midpoint distance witness are in `tests/references/edge_treatments_brep.py` and `docs/results/edge-treatments/reference.json` under `source_models.bracket.operations.fillet.midpoint_witnesses`.

The midpoint witness is a sampled distance, not a domain-wide deviation certificate. The differential obstruction above is exact and does not depend on that sample or on the fillet volume. No viewer geometry is generated by this report.
