# Composed native axial Boolean

## Protocol

- Starting revision: `b16bb585bb9f84413b9b89dc98db1175d3fa39d9`.
- Scope: `examples/composed_axial.rs`, `tests/composed_axial.rs`, and this directory.
- Preconditions: both stages use `reference_rectangle_v1`, axis Z, ordered integer millimeter intervals in the axial Boolean domain; the first result must be nonempty and connected so its occupied set is representable by one `AxialOperand` interval.
- Independent oracle: classify sample points in the analytic product `[-10,10] × [-10,10] × Z`, and compute volume as `400 × (length(union(A,B)) - length(union(A,B) ∩ C))`.
- Exit thresholds: exact interval membership and homogeneous control-net equality after restrictions; absolute analytic-volume error below `1e-10 mm³`.
- Rejected input: a disconnected first-stage union. The demonstration refuses instead of replacing it by its bounding interval, which would fill the gap.

## Construction and lineage

The first stage calls `build_axial_boolean(Union, A, B)` and audits the resulting model. The current public builder accepts `AxialOperand`, not an arbitrary prior model, so the second-stage operand is encoded by deriving its single occupied interval from the audited first result; the selected cells must be contiguous. The second stage then calls and audits `build_axial_boolean(Difference, derived_first_result, C)`.

To verify actual map reuse as well as the builder's interval encoding, the example also constructs the final carried map set by restricting the first-stage blocks directly to the surviving cells. Uncut maps are copied unchanged; a cutter that slices an inherited cell creates an exact affine restriction in z while retaining the source block ID, source interval, result interval, and original A/B ownership. Every carried control net is compared exactly with the corresponding restriction of the second-stage cold build. Final face lineage resolves the second-stage composite operand through those carried maps; faces created by C retain C lineage.

This is a demonstration for products of one fixed section and intervals on Z. Its intervals reduce the Boolean to a one-dimensional occupancy partition. It does not perform generic 3D intersection, and it does not establish geometric intersection for arbitrary sections, axes, or oblique tools.

## Reproduce

From the repository root:

```sh
cargo test --test composed_axial
cargo run --example composed_axial
```

`recipe.json` records four admitted operand triples and one required refusal. The demo's default case is `A=[0,8]`, `B=[4,12]`, `C=[5,9]` mm. Its first result has cells `[0,4]` from A, `[4,8]` from A and B, and `[8,12]` from B. The C cut splits the inherited `[4,8]` cell at 5 mm and removes `[5,9]` overall. The final material is `[0,5] ∪ [9,12]`, volume `3200 mm³`; the three carried maps retain their source identities and match the cold builder restrictions.

The four positive cases have analytic volumes `3200`, `4000`, `3200`, and `4400 mm³`. Tests integrate volume from the carried control-net bounding boxes and query only those carried maps for point coverage; expected material points have coverage exactly one and expected void points coverage zero. Samples lie strictly inside the rectangle and at quarter-mm z coordinates, avoiding cut planes. Face lineage maps source face tags and original A/B owners from first-stage `FaceProvenance`; new C faces remain tagged C. The disconnected case `[0,2] ∪ [5,8]` is rejected before stage two.

## Axial tensile check

Each carried result interval is also used as a two-node, constant-area, linear 1D axial bar element. Elements sharing a z endpoint share one node; every disconnected material component is solved separately with its lower endpoint clamped and a `1000 N` tensile load at its upper endpoint. The cross-section is the `20 × 20 mm` reference rectangle (`A=400 mm²`) and `E=210000 N/mm²`. For every component in all four variants, nodal displacement, tip displacement, and assembled strain energy match `u_tip = F L/(E A)` and `U = F² L/(2 E A)` within `1e-11` relative error.

This is a 1D axial specialization on the carried interval partition. It validates axial composition and load transfer for the constant-area rectangle; it is not a 3D elasticity result and makes no claim about lateral or shear response.

## Scope limits

The current `AxialOperand` API stores one section family and one interval. This example composes a connected first result. The direct map restriction preserves the first model's charts and provenance, but the second-stage canonical builder remains the independent cold comparison. Disconnected unions need a multi-interval operand representation before they can be passed as a single second-stage input.
