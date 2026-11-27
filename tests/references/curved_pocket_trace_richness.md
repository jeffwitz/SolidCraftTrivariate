# Curved pocket physical trace ranks

Reproduce from the repository root with:

```sh
python3 scripts/check_curved_pocket_traces.py \
  --model docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json \
  --output tests/references/curved_pocket_trace_richness.json
```

The source is the nine-map curved pocket. Its five lower-layer XY maps form
eight complete oriented interior section seams; the four wall maps are checked
to preserve the same XY controls with height. Cubic and quartic tensor fields
use endpoint-Hermite spans with C1 internal knots. Seam values and tangential
derivatives are identified in the common physical edge parameter. Physical
transverse-gradient equations use exact rational chart Jacobians and retain
every unisolvence endpoint, including repeated vertex rows.

The JSON reports exact rational row ranks for isolated two-chart seams and the
whole five-chart layout. It also reports observable projection ranks using
`rank([A; P]) - rank(A)`: `P` is the sequence of value Bernstein coefficients,
the first-side physical-gradient numerator samples with the known Jacobian
determinant omitted, or both together. The five-chart projections are measured
after value/tangent gluing, after the local seam, after the other seven seams,
and after all seams. Prime-field ranks are retained as an independent check on
the exact elimination.

## Recorded ranks

The five-patch scalar space starts from five independent chart spaces and
identifies the value and tangent trace on each seam. The table reports its
dimension before transverse coupling, the exact physical transverse rank, and
the remaining section-space dimension after all eight seams.

| Degree `p` | Spans per axis `n` | After value/tangent gluing | Transverse rank | After all physical C1 seams |
|---:|---:|---:|---:|---:|
| 3 | 1 | 52 | 36 | 16 |
| 3 | 2 | 136 | 68 | 68 |
| 3 | 4 | 424 | 140 | 284 |
| 3 | 8 | 1,480 | 284 | 1,196 |
| 4 | 1 | 89 | 44 | 45 |
| 4 | 2 | 260 | 84 | 176 |
| 4 | 4 | 872 | 172 | 700 |
| 4 | 8 | 3,176 | 348 | 2,828 |

The four outer/radial seams have isolated two-patch transverse ranks
`5, 9, 17, 33` for cubic fields and `6, 11, 21, 41` for quartic fields as
`n=1,2,4,8`. The four footprint seams give `6, 10, 20, 40` and
`7, 12, 24, 48`, respectively. These are exact ranks after matching the
value/tangent trace; all vertex equations are included.

For each edge, the observable projection ranks are also tabulated in the JSON.
Adding only the other seven seam constraints leaves the ranks of edge value,
transverse numerator, and their joint projection unchanged from the
value/tangent-conforming baseline. Adding the edge's own C1 equations reduces
the transverse and joint ranks; imposing all eight seams gives the same
per-edge projection dimensions as imposing that edge alone. This is a rank
statement, not an equality of the projected subspaces.

This is a section-space trace audit. It does not include the axial field
factor, clamp, elasticity solve, geometry coverage proof, or convergence claim.
