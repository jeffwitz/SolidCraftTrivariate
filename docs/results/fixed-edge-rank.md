# Fixed-partition edge regularity: bracket and S3

## Question and conclusion

This audit identifies block-parameter edges where two incident sides are both exposed parts of the physical boundary. If those two faces are required to meet smoothly along their complete common edge, a regular trivariate map on the unchanged block must have zero Jacobian on that edge. The conclusion also holds when several such edge requirements meet at a vertex: vertex blending cannot remove the edge-local rank condition. A blend represented by a new face and new incidences is a different topology and is outside this fixed-incidence result.

The independent scan found 36 such **block-parametric edge incidences** in the original 10-block bracket and 85 in the delivered 93-block S3 pocket. These are not deduplicated counts of physical B-rep crease curves; a geometric crease can be represented by multiple chart incidences. Exact rational derivative evaluations at five parameters per incidence classified every sampled pair as physically sharp; none was tangent. Representative pairs are listed below. Their existing certified positive Jacobian bounds show that the current maps are regular, so changing the two incident boundary faces to have equal tangent planes would contradict that regularity requirement.

This is a local obstruction to smoothing the original two faces while retaining their shared block edge. It does not prove that no redesigned topology or added blend patches can produce a globally smooth solid. Such a redesign must replace the original face incidences and be certified separately.

## Frozen inputs

| Part | Model | SHA-256 | Blocks / interfaces / exposed sides | Existing certified Jacobian lower bound |
|---|---|---|---:|---:|
| Original compact bracket | `docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json` | `82fb0d3f7a0977a47c705886d13632900dc62429be572231fc4f152a38688e09` | 10 / 11 / 38 | 67.73333333332907 mm³ |
| S3 face pocket | `docs/results/s3-face-pocket.model.json` | `f5bc0939201da16ef718d2b12ff0d46fae30b0d524b15480d69acaacc0d96140` | 93 / 191 / 176 | 2.1833038027716114 mm³ |

The bounds are read from the corresponding existing audits (`fieldcraft_x64_y64.audit.json` and `s3-pocket-delivery.json`); they are not recomputed by this script. They are positive whole-domain certificates, unlike the local edge samples reported here.

## Explicit local witnesses

The bracket block `wing_east_east` has the edge `u` with `v1,w0`. Its two exposed faces are `bore/east/wall` and `outer/profile`; exact rational derivatives at `u = 0, 1/4, 1/2, 3/4, 1` give a 90° angle between the two tangent planes at each point. This is a physical bore-to-stock crease, rather than an internal chart seam. A second bracket example is `elbow_east_fillet`, edge `u` with `v0,w0`, where `outer/profile` meets `start`; the five sampled tangent-plane angles are also 90°.

For S3, block `s3_i_u_zbase_g2_s0_p0_z2_src0`, edge `u` with `v1,w1`, joins exposed faces `nurbs_z_hole/wall` and `stock/zmax`. This is the top lip of the through-hole against the top stock face. Its exact rational tangent planes are perpendicular at all five sample parameters. The sampled absolute Jacobian determinant magnitudes on this edge range from 211.2 to 384 mm³ in the map's unit-cube coordinates.

The scan found no locally tangent pair among edges with two exposed sides in either model: 36/36 bracket pairs and 85/85 S3 pairs have a nonzero exact cross product of their face normals at at least one tested rational parameter. In these two particular witnesses the angle remains 90° at all five tested parameters. This is sampled geometry evidence for sharpness, not a proof that an edge's tangent angle has a specific value at every parameter.

## Rank argument

Let `F(r,s,t)` be a regular map of a block. Consider the cube edge `r=r0, s=s0`, with `t` varying. The two incident boundary faces have tangent planes

```text
P_r = span(F_s, F_t)  (r = r0)
P_s = span(F_r, F_t)  (s = s0).
```

The edge tangent is `F_t`, and regularity requires it and the two transverse directions to be linearly independent. If `P_r = P_s`, then `F_s` belongs to `span(F_r,F_t)`. Therefore `det(F_r,F_s,F_t)=0`. This contradicts a strictly positive Jacobian lower bound on the closed block. The same argument applies after permuting axes and at either endpoint of the free parameter. It is local to each edge, so satisfying additional tangent constraints at a shared vertex cannot make the determinant nonzero along the edge.

An inserted fillet/blend patch can make the *new* surface chain smooth by replacing the old edge with two new rails and an intervening face. That changes the incident face set and requires new patch maps, interfaces, and positive-Jacobian certificates. It is not evidence that the two original exposed faces can be made tangent while preserving their existing common block edge.

## Reproduction

From the repository root:

```sh
python3 tests/references/fixed_edge_rank.py --output build/fixed-edge-rank.json
```

The script reads the two frozen homogeneous Bezier models and enumerates each block's twelve parametric edges. It keeps an edge only when both incident coordinate faces occur in `topology.exterior`. It evaluates the rational-map derivatives exactly with Python `Fraction` arithmetic at five rational edge parameters, compares the resulting face normals by an exact cross product, and prints floating angle/Jacobian magnitudes as readable diagnostics. It does not use the Rust constructor, a CAD kernel, or sampled geometry to establish the rank implication.

The scan output is regenerated under `build/`; it is intentionally not a separate source of truth. Input SHA values and the findings needed to review the conclusion are recorded above.

## Limits

- Five exact parameter samples can witness non-tangency at those points but do not certify an angle bound over an entire edge.
- The rank implication itself is exact, conditional only on a regular edge tangent and unchanged two-face incidence.
- No alternate partition, blend-patch construction, manufacturing radius, or B-rep validity claim is assessed here.
- The audit does not infer smoothness from shared vertices or from unrelated internal block interfaces.
