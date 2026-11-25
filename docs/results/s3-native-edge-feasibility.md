# S3 native edge-treatment feasibility at 0.25 mm

## Scope and conclusion

This note concerns the delivered S3 pocket model and its existing 93-block
trivariate partition. The target is a complete treatment of its physical crease
edges with a 0.25 mm fillet or equal-setback chamfer, while retaining that
partition. No remeshing, added blocks, or alternate partition is included.

The exact fixed-parameter treatment is obstructed on a named top-lip edge of
the NURBS through-hole. Its constant-distance contact rail has a genus-3
algebraic graph, so it cannot be represented exactly by a rational Bezier rail
that preserves the source-span correspondence. The witness block has no unused
face side to receive a blend: all six sides are already assigned to two
exterior faces and four interfaces. This does not prove that every possible
global rational repartition is impossible. A certified rational approximation
is feasible in principle, but no such approximation has yet been constructed
or certified in this fixed partition.

The independent OCCT reference does demonstrate a useful separate result: it
applies 0.25 mm operations to all 50 detected S3 creases in one operation and
produces valid closed STEP solids. It is not a native trivariate result and
does not certify the requested radius/setback or chart Jacobians.

## Frozen geometry and named edge

- Start revision: `d25b6b6e17f1c95b8371a9ae714298fd87861677`.
- Parent: the delivered S2 model embedded in
  `tests/references/s3_face_pocket.recipe.json`; source-model SHA-256 is recorded
  in `docs/cycle-nurbs-multiboucles.md`.
- S3 model: `docs/results/s3-face-pocket.model.json`, 93 blocks and 191
  interfaces.
- Operation under study: physical setback/radius `1/4 mm`.
- Named edge witness: the top lip on source span `src0`, in block
  `s3_i_u_zbase_g2_s0_p0_z2_src0`, along local source-span coordinate
  `u in [0,1]`, at `v=1, w=1` (world `z=12 mm`). Its incident block sides are
  `v1 = nurbs_z_hole/wall` and `w1 = stock/zmax`. Both are exterior design
  faces of this same block, so this is a physical crease, not a smooth chart
  seam.

The other four sides are all occupied by existing interfaces:

| Side | Existing neighbor | Consequence |
|---|---|---|
| `w0` | `s3_i_u_zbase_g2_s0_p0_z1_src0`, side `w1` | lower axial layer |
| `u0` | `s3_i_u_zbase_g2_s11_p0_z2_src7`, side `u1` | periodic sector closure |
| `u1` | `s3_i_u_zbase_g2_s1_p0_z2_src1`, side `u0` | next source span |
| `v0` | `s3_o_u_zbase_g2_s0_p0_z2_src0`, side `v1` | radial companion block |

There is no unassigned block side to repurpose as the blend face. Reassigning
one of these sides would remove an existing interface and require a matching
change to its neighbor; that is not a local face reassignment within the
existing partition.

## Exact obstruction on the top-lip edge

The independent symbolic reproducer is
`tests/references/edge_offset_obstruction.py`; its exact output is
`docs/results/edge-offset-obstruction.json`. On the first active source span,
the homogeneous curve is

```text
X(t) = 322 t^3 - 516 t^2 - 30 t + 824
Y(t) = 134 t^3 - 216 t^2 + 30 t + 324
W(t) =   4 t^3 -   6 t^2        + 10,   0 <= t <= 1.
```

Let `D=(X'W-XW', Y'W-YW')`. The exact squared-speed numerator is

```text
D.D = 288 P(t), with the speed denominator W(t)^4,
P(t) = 73 t^8 + 120 t^7 + 156 t^6 - 1776 t^5 + 766 t^4
       - 712 t^3 + 2796 t^2 + 625.
```

Exact polynomial arithmetic gives `gcd(P,P')=1` and `degree(P)=8`. A Sturm
chain has four sign variations at both 0 and 1, so `P` has no real zero on the
source interval; `P(0)=625` and `P(1)=2048`. Thus the real edge is regular over
the full span, while the square root does not simplify algebraically.

For a planar constant-distance rail, the unit normal contains
`D/sqrt(D.D)`. The rail at 0.25 mm therefore contains `sqrt(P(t))` in its
coordinate functions (the nonzero constant factor `sqrt(288)` does not alter
the function-field genus). The graph that retains the source correspondence is
birational to

```text
s^2 = P(t).
```

Because `P` is square-free of degree eight, the smooth projective completion
of this hyperelliptic curve has genus `(8-2)/2 = 3`. A rational reparameterization
of this graph would make its function field rational, contradicting genus 3.
Consequently neither knot insertion nor a finite rational control-net edit can
make this exact, synchronized constant-setback rail rational. The result
applies to this named source span and to a rail required to preserve its source
correspondence. The physical rail itself also has genus three. At an offset point (U,V),
its possible foot parameters x obey two polynomial equations. Set
A=X(x)−U W(x), B=Y(x)−V W(x). Those equations are
A²+B²−W²/16=0 and A Dx(x)+B Dy(x)=0. Their degrees in x are at most
6 and 7. Specialize (U,V) at the actual offset of t=0 on either side.
The degrees remain 6 and 7 and an exact gcd over the resulting quadratic
number field is the monic polynomial x. Hence a nonzero subresultant
minor gives generic gcd degree at most one. The actual foot gives a
common root, so the generic monic gcd over the field of rail coordinates
is linear. Its coefficient recovers t rationally from (U,V). One nonzero
component of the offset vector then recovers the signed square root of
the speed numerator rationally. Thus projection from the genus-three
normal graph to the physical rail is birational. A finite collection of
rational NURBS edge spans cannot trace even an open arc of this rail
exactly, including under a change of parameter. This rules out an exact
constant-setback rail *as an explicit rational patch edge*. An implicit
intersection or a bounded approximation is a different representation.

The same obstruction affects both candidate operations at this edge: on the
planar `stock/zmax` face, both a circular fillet and an equal-setback chamfer
need a contact rail at a fixed physical distance from the curved lip. The
second rail and the blend surface do not remove the non-rational first rail.

## Can a rational approximation be certified?

Yes in principle; it is a different geometry claim from exact radius/setback.
On this compact, regular span, one can subdivide `[0,1]`, approximate
`1/sqrt(P(t))` by positive rational-coefficient polynomials on each subspan,
and form rational Bezier approximations of the contact rail. A certificate
could bound the residual of `P(t) q(t)^2 - 1` with outward-rounded interval or
Bernstein bounds, then turn that into a world-space rail error bound using the
exact `D` and positive `W`. The rational patches would still need exact shared
end traces, a certified positive Jacobian for every modified volume chart, and
vertex/junction closure across all incident edges.

No polynomial degree, approximation tolerance, or topology budget was
pre-registered for that approximate design, and none was changed after seeing
the obstruction. No approximation was fitted here. In particular, the present
fixed partition offers no spare side in the witness block; a certified blend
must be represented by modifying existing face maps and coordinating all
four neighbor interfaces, or else it would require the out-of-scope partition
change. Approximation existence alone is not a qualified S3 treatment.

## Independent full-solid reference (not native qualification)

`tests/references/edge_treatments_brep.py` selected the physical creases of the
independently reconstructed S3 boundary solid. Its inventory has 51 edges: 50
classified creases and one periodic seam. The seam was excluded. At 0.25 mm,
OCCT's all-at-once fillet and chamfer each produce a one-solid, valid, closed
candidate with 101 faces; both exported STEP files re-import as valid, closed
solids.

| Result | Fillet | Chamfer |
|---|---:|---:|
| Re-imported STEP volume (mm^3) | 53130.25366450874 | 53120.45283565439 |
| Source-edge midpoint displacement range (mm) | 0.06789617–0.10355365 | 0.15415830–0.17677756 |
| Midpoints checked | 50 | 50 |

These are shell-displacement witnesses, not measurements that every point on
every new rail is exactly 0.25 mm from its incident source face. The report
also records construction/re-import volume differences of `1.10e-4 mm^3`
(fillet) and `3.72e-5 mm^3` (chamfer). Validity, closure, STEP round-trip, and
volume do not certify a native trivariate partition, exact radius, exact
setback, or positive Jacobians. Full data and STEP hashes are in
`docs/results/edge-treatments/reference.json`.

## Reproducibility and exit status

Relevant commands from the frozen experiment were:

```sh
rtk python3 tests/references/edge_offset_obstruction.py
rtk python3 tests/references/edge_treatments_brep.py S3
```

The first derives the obstruction from the versioned rational source curve. The
second requires the campaign's configured FreeCAD/OCCT environment and writes
its outputs under `build/edge-treatments/`. Existing artifacts are retained;
this report does not claim that those commands were rerun during its drafting.

**Status:** exact constant-distance native rail on the named edge: refused by
an exact genus certificate; same-partition local face reassignment: no spare
side; certified approximate rail and chart-level all-edge fillet/chamfer:
open. The independent OCCT all-crease solid results are positive references,
not substitutes for those open native criteria.
