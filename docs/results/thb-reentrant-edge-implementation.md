# Reentrant-edge THB implementation record

The frozen physical-error campaign remains exactly as recorded in
[`thb-reentrant-edge-protocol.md`](thb-reentrant-edge-protocol.md). This note
records the Rust basis implementation and the structural width-6 experiment;
it does not alter the frozen boundary lifting or error-comparison protocol.

## Native construction

`src/thb_reentrant_edge.rs` builds a tensor cubic, C2 interior-knot THB basis
on the three positively oriented affine charts of the L section. The prepared
level has 8×8 spans per chart. The width-4 hierarchy refines nested corner
boxes at `[0,1/2]^2`, `[0,1/4]^2`, and `[0,1/8]^2`. Its charts share complete
polynomial traces on the two declared interfaces. The exact field coefficient
maps preserve the global constant and include the lower-level truncation terms;
they are not just an overlay of independently refined cells.

For the width-4 schedule, the active scalar DOFs, integration leaves, and
unsewn source-level function counts are:

| Maximum level | Active scalar DOFs | Leaves | Source-level functions before seam identification |
|---:|---:|---:|---:|
| 0 | 341 | 192 | 363 |
| 1 | 477 | 336 | 315, 192 |
| 2 | 613 | 480 | 315, 144, 192 |
| 3 | 749 | 624 | 315, 144, 144, 192 |

The full uniform level-3 reference has 13,333 scalar DOFs after matching the
same chart traces. Each source-level hierarchy's exact polynomial transfer
has full rank. The tests check prepared-to-level-3 reconstruction and each
successive map L0→L1, L1→L2, and L2→L3 by exact rational Bernstein
coefficients. The two seam trace bases are compared as coefficient vectors in
the common finest knot basis, and their sampled global traces are replayed on
three interior points per seam.

## Exact box-cover and level admissibility

The level-3 width-4 integration cover is replayed on the 64×64 finest grid per
chart with no duplicate or uncovered bin. The exact local extraction shows at
most two consecutive active source levels on every leaf. Its histogram,
reported as `(minimum source level, maximum source level, consecutive level
count, number of contributing levels) → leaf count`, is:

| Level range | Span | Contributing levels | Leaves |
|---|---:|---:|---:|
| 0–0 | 1 | 1 | 144 |
| 0–1 | 2 | 2 | 117 |
| 1–1 | 1 | 1 | 27 |
| 1–2 | 2 | 2 | 117 |
| 2–2 | 1 | 1 | 27 |
| 2–3 | 2 | 2 | 117 |
| 3–3 | 1 | 1 | 75 |

## Wider box admission check

As a separate structural experiment, `build_hierarchy_with_widths(3,
&[6,6,6])` accepts nested boxes `[0,3/4]^2`, `[0,3/8]^2`, and
`[0,3/16]^2`. It gives 1,277 active scalar DOFs and 1,164 leaves against the
same 13,333-DOF uniform level-3 space. Exact partition, direct transfer rank,
cell reconstruction, chart traces, and nonoverlapping coverage all pass. Of
the 1,164 leaves, 597 have one contributing source level and 567 have two;
none exceed the same two-level bound. This is structural box coverage only;
the frozen physical-error campaign did not include a width-6 solve.

The builder refuses levels above three, width vectors of the wrong length,
widths outside 4–8 spans, and any seam or exact-transfer inconsistency. Tests
include `usize::MAX` level inputs to ensure rejection occurs before allocation.

## Verification

The native integration test is `tests/thb_reentrant_edge.rs`; the public Rust
entry points are `build_hierarchy`, `build_hierarchy_with_widths`,
`build_uniform`, `audit_exact_partition`, `audit_exact_transfer`,
`audit_nested_transfers`, `audit_leaf_cover`, and `integration_leaves`.
The extraction API returns Bernstein rows as `bernstein[4*i+j][local_dof]`.
An independent Cox–de Boor/Fraction oracle and the Ritz error curve are
reported separately from this structural record. The level-1 smoke runs with
the normal test suite; the expensive level-3 transfer, width-6, and uniform
extraction qualifications are ignored by default and must be included in the
release check with:

```sh
cargo test --release --locked --offline --test thb_reentrant_edge -- --include-ignored
```
