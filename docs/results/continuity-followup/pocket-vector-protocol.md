# Native vector elasticity on the bilinear pocket — prospective supplement

Freeze before the first solve, based on fbe8723. Use the real thickness
profile and source block incidence of the uncurved scaffold of
docs/site/data/plate_blind_pocket/plate_blind_pocket.recipe.json (9 trilinear
blocks, 5 section footprints). Confirm the volume is 46800 mm³.
Use the manufactured displacement (0,0,((x+10)/70)^5 z/20) mm with the
existing isotropic material, bottom clamp, body force and all exact natural
tractions including the pocket walls and floor. The z dependence lies exactly
in the chosen thickness space. Solve for all three displacement components
in the global height-linear bottom-clamped trial space. This fixes the axial
field, separately from full thickness refinement.

Keep C1 inside every tensor section patch in both branches (internal knot
multiplicity p-1). Join the original patches with either C0 or with full
physical C1. At fixed p compare n=1,2,4,8 subdivisions in both section
directions for p=3 and p=4. Construct exact continuity constraints by polynomial
unisolvence, retain all vertex equations, solve the signed rational kernel.
Preserve refusal or failure for any required row. Use q8; independently repeat
the last two n levels at q10 if discrepancies in rates, geometry or work occur.
Require positive geometry, volume defect <1e-10 relative, relative equilibrium
residual <1e-9, virtual-work defect <1e-7. Report both energy and L2 rates,
actual free vector DOFs, and any loss of conditioning. Declare a rate consistent
with order p in energy only if the last measured doubling is within 0.5 of p;
two fine levels alone do not establish an asymptotic theorem. Independent
controls: compare n=1 bilinear cubic C0 with the existing native elevated
space at axial degree one, and check physical C1 against the independent
section experiment. Expanded field coefficients are mechanical unknowns,
while the nine geometric blocks and all material boundaries are fixed.
