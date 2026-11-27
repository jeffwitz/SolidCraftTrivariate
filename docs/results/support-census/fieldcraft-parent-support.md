# FieldCraft parent-support and native H1 replay

## Exact support question

This companion examines every full-face interface on the five compact
FieldCraft bracket variants. Each variant has ten blocks and eleven declared
joins (55 interface observations total). The source models are the versioned
`affine_image_v1` gallery representations; the underlying generating profile
family is `orthogonal_sketch_features_v1` with partition method
`four_sector_local_sweeps`.

For each join, the script first elevates the two homogeneous trivariate
Bernstein nets to a common componentwise degree, extracts the oriented face
nets, and solves for a positive constant projective gauge. It then compares
the entire tangential Bernstein control net of each inward normal derivative
from order one through the common normal degree. A candidate requires one
positive local interval-length ratio whose q-th power gives the exact ratio
for every coefficient of the q-th homogeneous jet. This is an exact
necessary condition for inheriting the complete smooth homogeneous jet across
a parent split; it is stricter than position-only face equality.

Result: **0/55** joins pass the complete-jet test. The report stores the
orientation, common degree, gauge, each tested order, and a rational
coefficient witness for every mismatch. The model records block IDs but no
normalized source intervals for these faces, so no interval length or source
span can be claimed as proven. The fresh-release `spline-eligibility` oracle
reports `no_complete_parent_witness` for each variant, with evidence `no
complete source restriction proof in this recovery family`; its
`known_generating_sources` list is empty. All joins remain C0 position matches
with higher-order continuity not certified. Together these results qualify
**zero recoverable complete-jet parent subdivisions**; they do not assert that
arbitrary C0 reassembly into a spline with deliberately repeated knots is
impossible. The separate `multiple_generating_supports` observation concerns
the crossed-bore family, not FieldCraft.

The audit fixture itself includes a positive two-block degree-one parent
split with local interval lengths in ratio 1:2. Its complete homogeneous jet
and inferred first-to-second interval ratio `1/2` are recovered exactly,
while a degree-two control with incompatible second jet is rejected.
This guards against a vacuous all-negative implementation.

## Native C0 H1 space and mechanics

The reusable `nurbs_space::export_audited_multiblock` now accepts any model
that passes the native geometry audit and whole-face incidence decode. The
existing CLI route is `export-nurbs-space MODEL OUTPUT`; the affine FieldCraft
wrapper is now an accepted source family. The field elevates every block to
the maximum geometry degree, here `[2,2,2]`, and identifies the nine trace
coefficients on each oriented interface. For x64/y64 this gives 171 scalar
DOFs over ten blocks and 99/99 exact shared trace-node identifications. Its
affine patch test evaluates four points per block against a full nonsymmetric
affine displacement map. The independent native audit certifies the complete
geometry, not just these points.

The existing Rust `solve-job` is already a native multiblock H1 mechanics
path for the same degree-two C0 field: `analysis_space::build_analysis_space`
elevates to 27 local Q2 nodes per block, validates the shared-face geometry,
and unions the oriented face nodes; `mechanical_job::solve_job` then assembles
the rational Q2 stiffness and solves it by PCG. The new integration test
compares the full local-to-global equivalence partition from that solver to
the degree-harmonized `nurbs_space` export, not merely the DOF count. Their
partitions match exactly. This is why no second mechanics solver was added.

On the published x64/y64 model, the existing q16 steel job reproduces 171
scalar DOFs, integrated volume `32144.77892538882 mm³`, strain energy
`135.44475139018235 N·mm`, relative residual `7.54e-11`, and six rigid
constraints. Relative force and moment imbalances are `2.85e-11` and
`1.05e-11`; work/energy defect is `2.77e-14`. This is a native Q2 mechanics
replay on the same C0 function space, not an independent mechanical
discretization or a proof of stress convergence.

The separate axial inherited-extrusion compiler remains unavailable for the
FieldCraft source: its five regions have certified translations along local
z, y, and x, so one global ascending physical-z plane sequence does not cover
all cuts. The generic audited NURBS H1 export and existing mechanics path do
not require that axial extrusion assumption.

## Reproduce

From the repository root, with the offline Rust toolchain:

```sh
python3 scripts/audit_fieldcraft_parent_support.py \
  --cli target/release/solidcraft-trivariate \
  --output build/support-census/fieldcraft-parent-support.json
python3 -m unittest tests/test_fieldcraft_parent_support.py
cargo test --locked --offline --test fieldcraft_nurbs_h1
target/release/solidcraft-trivariate export-nurbs-space \
  docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json \
  build/support-census/fieldcraft-x64-y64-h1-c0.json
```

The detailed per-interface evidence is regenerable at
`build/support-census/fieldcraft-parent-support.json`; the mechanical result
is regenerable with the q16 job and source feature tags recorded in the Rust
integration test. The report records the executable SHA-256 used for each replay. The signed normal convention covers equal as well as opposite endpoint flags. Each gauge check includes the entire face; higher-order jets that both vanish impose no interval-ratio restriction.
