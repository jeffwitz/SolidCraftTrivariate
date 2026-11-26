# Bracket fillet geometry and rank audit

## Disposition

The x64 FieldCraft fillet example emits a fixed-partition candidate with 10
logical charts, 30 span elements, and 53 declared interfaces. The candidate
sets `singular_boundary_qualified=true` and reports determinant lower bound
zero. It is not a strictly positive-Jacobian chart on the closed cubes: the
independent exact Bernstein audit finds zero-Jacobian contact rails on the
first and last 0.25 mm portions of both 8 mm bores. The toroidal surface
passes the exact round invariant on all 16 angular patches. These checks
qualify the intended geometry and the location of its rank loss; they do not
make the volume maps regular on their closed domains.

## Exact surface and knots

Both source bores have radius 4 mm and axial bounds `[0,8] mm`. The fillet
radius is one quarter-millimetre tick, `r=0.25 mm`. The fixed parameter breaks
are `[0,1,31,32]/32`, corresponding to physical stations `[0,0.25,7.75,8]
mm. Each mouth uses four quarter-circle angular charts, so the four physical
mouth edges are represented by 16 fillet patches.

For each emitted end patch, the oracle clears the rational denominator and
checks this torus identity as a complete bivariate polynomial:

```text
(X² + Y² + Z² + ((R+r)²-r²) W²)² - 4(R+r)² W² (X² + Y²) = 0,
```

Here `Z` is the axial offset from the meridional-circle center, while `X,Y`
are the transverse offsets from the bore axis. This verifies the complete
0.25 mm constant-radius meridional surface, rather than a few sampled points.
At five rational angular parameters, it also verifies that the meridional
endpoints are exactly `(rho,z)=(R+r,z_end)` at the mouth and `(R,z_end±r)` at
the cylinder tangent, with no transverse shift off the original radial ray.

The independent Bernstein determinant check finds the only zero strata on the
16 fillet patches. In map coordinates `(u,v,w)=(angle, radial, meridional)`,
the first mouth has `det(J)=0` on `v=1,w=0` for every `u`; the other mouth has
`det(J)=0` on `v=1,w=1`. Determinant Bernstein coefficients are nonnegative,
and no other zero strata occur. The eight cylindrical bore spans and six
split elbow spans have no zero strata. The certificate code is the independent
integer implementation in `tests/references/s3_stock_rim.py::certify`.

## Independent volume reference

For the four fillets together, the analytical removal is

```text
4π(2 R r² (1 − π/4) + r³ (5/3 − π/2))
= (101/48)π − (17/32)π²
= 1.3672072038498007 mm³.
```

This is a design reference, not a B-rep volume measurement. The serialized
trivariate retains the exact zero-Jacobian contact rails.

## Reproduction

From the repository root:

```sh
cargo run --offline --example bracket_edge_edit -- fillet
cargo run --offline -- edit-fieldcraft-bore-rims \
  docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json \
  fillet build/fieldcraft-x64-fillet.json
python3 tests/references/bracket_fillet_rank_oracle.py
```

The example writes `build/bracket-edge-edit-fillet.json`; the generic CLI writes
`build/fieldcraft-x64-fillet.json`; the oracle writes
`build/bracket-fillet-independent.json`. It confirms all 30 generic CLI span
maps match the direct example's homogeneous maps exactly after one positive
common scale per map, with identical fixed parameter intervals and parent
topology. The source model SHA-256 is
`82fb0d3f7a0977a47c705886d13632900dc62429be572231fc4f152a38688e09`. The
candidate SHA-256 audited here is
`f6dbe0daa607bc1763158be27f77314e180845dbbc9621d70948683b72e3f1e3`; the
generic CLI file SHA-256 is
`4da023195bf582327c7256e5883ccd322cfc794bc26202ddc81961343199972e`.
