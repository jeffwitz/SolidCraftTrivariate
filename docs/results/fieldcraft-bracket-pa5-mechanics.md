# FieldCraft `bracket-hex` job on PA5

## Provenance and loading

This replays the original FieldCraft `filleted-bracket-fillet-stress` static load case against the reconstructed 48-block PA5 model in `build/solidcraft_sketch_bracket.json`. Source revision is pinned in `references.lock` (`tangency-labs/fieldcraft` `35e834652d8548384c2da8ef1b28d588c9af2ffd`). Source files: `upstream/fieldcraft/tests/oracle/filleted-bracket-fillet-stress/{recipe.json,design.json,README.md}`.

The SolidCraft sketch profile is ordered `(0,0) -> (60,0) -> (60,8) -> (8,8) -> (8,50) -> (0,50)`, then extruded along +z through 40 mm. Thus historical face `F2:side:l1` is the y=0 foot face, and `F2:side:l5` is the y=50 end face of the upright leg. This run selects these by exact model planes and `outer/sketch_1` provenance: fix all three displacement components on y=0, and apply +x traction on y=50.

The original recipe applies a 2000 N +x force to F2:side:l5. The current native job API accepts a uniform traction vector, so the source resultant is represented by 6.25 MPa +x on the exact 320 mm2 face (x=0..8, z=0..40), giving +2000 N. Steel is E=210000 MPa, nu=0.3. Quadrature order 8; no adaptive refinement. Job input: `tests/references/mechanical_jobs/fieldcraft_bracket_hex_pa5.json`.

## Run and result

Command: `rtk cargo run --quiet -- solve-job build/solidcraft_sketch_bracket.json tests/references/mechanical_jobs/fieldcraft_bracket_hex_pa5.json build/fieldcraft_bracket_hex_pa5.result.json`. The result has 48 geometric blocks, 705 scalar Q2 DOFs (2115 interleaved displacement coefficients), 1749 reduced unknowns, six independent rigid constraints, 298 iterations, and relative residual 9.79e-11.

Selected support patches: 24; loaded patches: 4. The loaded area is 320.00000000000017 mm2; applied resultant is [1999.9999999999995, 0, 0] N. Applied moment about the origin is [0, 40000.000000000015, -99999.99999999996] N mm. Support reaction is [-2000.00000002177, -2.98e-9, 9.92e-11] N; relative force balance 1.10e-11 and relative moment balance 5.81e-12. Strain energy is 144.1284781 N mm; work-energy relative error 3.25e-14.

An independent SciPy assembly of the serialized Q2 blocks and boundary load confirms 320 mm2, +2000 N and the same moment; its load RHS differs from Rust by at most 1.42e-13 N. The independent sparse solve differs from the Rust displacement coefficients by 5.33e-13 relative. This is a qualification of the transferred job and equilibrium on PA5; it does not reproduce FieldCraft’s tet10 stress peak, mesh, or solver. The result JSON is regenerable under `build/`.
