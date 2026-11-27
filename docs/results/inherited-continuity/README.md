# Inherited continuity on constructed solids

The campaign starts from the published tree at 3f9d2d915c918121503668691b489a39781c6dd5, the reference tree named by the [original T-spline specification](../../recipe-driven-tsplines.md). GitHub main pointed to that revision at startup. The working integration clone has the same tracked tree. The preexisting working directory contains separate historical experiments. The [frozen protocol](protocol.json) qualifies inherited continuity on the source solids before choosing an adaptive basis.

The bounded exact audit for neighboring Bézier restrictions is documented in [pairwise parent recovery](pair-parents.md). Its pairwise algebraic result is separate from CAD feature provenance and whole-solid certification.

An interface is counted once per model state. A matching homogeneous trace proves C0 compatibility up to a positive constant gauge. Recovering the two neighboring restrictions of one regular extrusion permits additional field smoothness across that interval. Physical C1 between independent supports is a different operation: its field space must satisfy the transverse gradient equations. Neither operation changes material membership or removes semantic load faces.

### Which route to use

- Use the canonical parent basis to restore source continuity on the two-span strip. A simple interior knot in a cubic curve gives C2 in its actual parameter.
- On a layered part, use the inherited axial space first. For a curved multiblock solid, the optional [matched-chart plan](../../../tests/references/matched_chart_jets.plan.json) also tests each lateral face and reports the exact interval ratio for admitted C1 joins. [Results on the three curved parts](#automatic-lateral-selection-on-the-curved-parts)
- Use `refine-thb-strip` for a local increase in displacement degrees of freedom inside a recovered strip. Use `solve-job` on the saved package to select original design faces and solve the prescribed load. [Local refinement and native elasticity](#canonical-parent-field-and-elasticity)

### Arithmetic and quadrature controls

The first all-lateral cubic-section cylinder run exceeded checked i128 denominator storage on the uniformly spaced exact probes. The low-degree controls completed. The revised exact evaluation visits distinct low-denominator rationals first; the polynomial unisolvence bound and required zero residuals stay fixed. The initial cylinder q12/q16 comparison also found a relative strain-energy change of about 1.0e-6 on the original section degrees. Higher quadrature is required for the final comparison at the protocol's 1e-8 threshold.

## Measured effect on the field

The 68 polygonal source recipes have 3,136 complete inter-cell faces. In physical height, 617 joins join two restrictions of the same translation of an entire homogeneous section. With a piecewise quadratic axial field, matching its first derivative on these 617 faces reduces the sum of scalar dimensions from 7,829 to 6,751. The other 2,519 faces are lateral. Area, material membership, load selections and the geometry integration cells are unchanged. [Census by corpus](../support-census/README.md)

The native layer test accepts a section translated with constant horizontal shear per unit height, including oblique polygonal extrusions. It checks the translation on every homogeneous control and the oriented trace on the full face. The field compiler uses true interval lengths in the derivative equations. For the source strip it can also restore a spline in its generating-curve parameter. The ordinary cylindrical sector junction has independent charts on its two sides. A common position or a smooth geometric surface does not establish equal transverse field derivatives.

### Curved cylinder and blind pocket

We prescribe the physical displacement `u=(0,0,((x−xc)/width)² (z−zmin)/height)`. The body force is −div σ and all free boundaries carry σn, including the pocket floor. The base displacement vanishes. We integrate the unrestricted virtual work in the volume as a separate check on body and boundary loads. The table uses relative **energy-norm error** and vector unknowns remaining after essential conditions. All modes use the same source geometry and a quadratic axial field. “Cubic section” elevates both in-plane degrees to three. One rail includes all its restrictions over the height.

| Field | Cylinder free DOFs | Cylinder error | Pocket free DOFs | Pocket error |
|---|---:|---:|---:|---:|
| Source section, C0 | 306 | 26.1604% | 294 | 0.12809% |
| Inherited axial C1 | 204 | 26.1927% | 222 | 0.12813% |
| Plus one physical C1 rail | 144 | 30.7001% | 186 | 0.17473% |
| Plus all physical C1 rails | 36 | 95.3694% | 48 | 50.5679% |
| Cubic section, C0 | 936 | 5.8809% | 600 | 0.12609% |
| Cubic section, inherited axial C1 | 624 | 5.8837% | 456 | 0.12658% |
| Cubic section, plus one rail | 552 | 6.1868% | 411 | 0.13861% |
| Cubic section, plus every rail | 144 | 55.0193% | 180 | 21.3347% |

Each row uses a single in-plane element per source patch. On the cylinder, two spaces with 144 free coefficients give energy errors of 30.7% and 55.0% under different constraints. On the published curved pocket, all nine geometry maps have biquadratic sections that differ from their bilinear corner interpolants. The 21.33% result measures the fully constrained coarse field on those curved charts. A single level cannot diagnose asymptotic locking. The [h-refinement study](../continuity-followup/) uses the separately rebuilt bilinear scaffold to test degree and continuity on a geometry covered by the bilinear AS-G1 construction.

Selective compatibility already has a measurable benefit on the curved pieces: the S3, multiloop and crossed-plate solves use 36%, 45% and 30% fewer free displacement coefficients, respectively, while the largest absolute change in relative energy error is 0.006 percentage points. The load, geometry and degree stay fixed in each pair. [Per-part dimensions and errors](#automatic-lateral-selection-on-the-curved-parts)

For a physical C1 request the Rust compiler uses exact unisolvent gradient equations over the complete lateral face, solves for a signed rational lift, checks the equations and independently checks affine reproduction. The independent Python oracle instead clears the weight denominators and tests the determinant of the three tangents of the physical field graph coefficient by coefficient. On the cubic cylinder, all 1,248 column/layer identities across 24 lateral faces vanish. On the cubic pocket all 732 identities across 12 faces vanish. The unconstrained low-degree negative controls exhibit nonzero jump polynomials on both models. These experiments establish physical first-derivative continuity of the **field**. They do not construct a global smooth coordinate system for the material.

At the prescribed tolerance the strain energy and volume-averaged displacement stabilize between quadrature orders 16 and 24, depending on the case. The build reports retain all orders, solve residuals, stress norms and virtual-work defects. The difficult all-cubic cylinder requires coefficients larger than a two-prime reconstruction can guarantee. Four-prime rational reconstruction, followed by exact residual verification, resolves it with 60 scalar coefficients before clamping. Rank is certified by the modular lower bound together with a verified complementary rational null space. A floating-point pivot heuristic can choose coordinates; it never certifies an interface. [Machine-readable measurements](mechanical-evidence.json)

### Which joins admit an affine transition?

For each lateral face, write the transverse derivative of one chart as α times the transverse derivative of the other plus β times its tangent. The affine test checks the complete polynomial identities after clearing denominators and requires α to retain a nonzero sign on the interval. The normalization fixes the coefficient of the derivative being expressed to one. In the stored face orientation it succeeds on 601/2,519 joins. Testing both chart orders gives 626/2,519 (24.85%), including 265 constant transitions. This sufficient first-order test leaves 1,893 outside its domain. It supplies coefficients for a local transition, not a proof that changing entire parameterizations around every vertex is possible.

Two further mechanics jobs keep only the strong lateral constraints selected by this test and use cubic physical loading `u_z=(1+x/width)³ (z−zmin)/height`. With quadratic in-plane and axial spaces, the seeded oblique pentagonal pocket has 342 free DOFs and 0.1879% energy error after restoring its axial field. Adding the 12 selected lateral-face constraints gives 180 free DOFs and 0.3031% error. The triangular blind pocket changes from 366 DOFs and 0.1894% to 303 DOFs and 0.1945% across five selected joins. These two post-census probes include an independent exact polynomial jet check. They are demonstrations of a selective field space, with no corpus-wide optimality claim.

## Decision after the interface census

1. Restore a proven common support before imposing additional field constraints. Keep boundary names and source restrictions available even when an analysis knot disappears.
2. Keep separate-support C0 as the default for elasticity. Select a physical C1 join explicitly where the model needs one, and compare the resulting error and free dimension. For planar transitions prioritize reconstruction or design of compatible chart coordinates at sketch time. Measure vertex-cycle compatibility and approximation power before treating a collection of individual transitions as one reparameterized solid.
3. For adaptive field enrichment on a recovered tensor chart, use the hierarchical pilot to test a genuinely local change at fixed geometry. The retained 2,519 lateral joins and S2's partial footprints need interface-specific treatment whichever interior basis is chosen. A single-parent T-mesh cannot smooth those joins. A THB pilot also leaves them intact. Future unstructured T-splines on extraordinary charts or partial faces would require an additional admission contract.

The hierarchical construction uses the truncation definition of Giannelli, Jüttler and Speleers. The analyses of adaptive hierarchical B-splines and analysis-suitable T-splines motivate separate alternatives. Neither paper certifies the present three-dimensional elasticity implementation. Any future T-mesh profile must use one stated basis and boundary convention throughout. [THB basis definition](https://doi.org/10.1016/j.cagd.2012.03.025) · [Hierarchical IGAFEM](https://arxiv.org/abs/1701.07764) · [T-spline IGAFEM](https://arxiv.org/abs/1910.01311)

## Reproduction

```sh
make check-inherited-corpus
```

The native `export-parent-field` reports the polynomial and physical coefficient maps of the recovered NURBS parent. On each span the relation is C = diag(w_span)⁻¹ E diag(w_parent). `solve-parent-field` uses these canonical parent basis functions. `export-inherited-space` builds the signed field extraction for a requested list of joins. `audit-inherited-space` replays its model and plan. `solve-inherited-space` verifies that replay before integration. `solve-inherited` provides the separate torsion and axial manufactured load cases. The original source model remains the owner of geometry and face provenance.

```sh
cli=target/release/solidcraft-trivariate
mkdir -p build/parent
"$cli" build tests/references/nurbs_strip.recipe.json build/parent/model.json
"$cli" recover-spline-parent build/parent/model.json build/parent/support.json
"$cli" export-parent-field build/parent/model.json build/parent/field.json
python3 scripts/check_parent_field.py
```

## Canonical parent field and elasticity

The native strip has 32 block Bernstein control occurrences, 28 globally identified C0 coefficients and 20 canonical NURBS coefficients. The simple cubic knot restores C2 in the source curve parameter. At q8 the independently assembled C0 stiffness and load pull back to the canonical parent with relative defects below 1e−14. The 60-dimensional free-body stiffness has exactly six numerical zero eigenvalues at the preregistered threshold. Its first positive eigenvalue is 0.005576 times the largest. All six analytical rigid modes have normalized residuals below 6e−16; the clamped matrix admits Cholesky factorization. [Spectrum and affine energy](parent-spectrum.json)

The local test uses the [frozen marking and forcing protocol](local-refinement-protocol.md). The external material boundary and the two source cells retain their original identifiers. A saved package embeds the source and refinement request; replay reconstructs its basis, integration leaves and boundary lineage before assembly.

```sh
cli=target/release/solidcraft-trivariate
mkdir -p build/parent
"$cli" build tests/references/nurbs_strip.recipe.json build/parent/model.json
"$cli" refine-thb-strip build/parent/model.json tests/references/thb_interior.plan.json build/parent/local.json
"$cli" check-thb-galerkin build/parent/local.json 4
make check-local-thb
```

The two source spans now have **20 geometry controls**, 80 coefficients after field degree elevation to tricubic, and 539 scalar functions on the prepared tensor mesh. Marking 16 of its 128 cells introduces five local scalar functions. Its interior fine box is 8 × 4 × 4 spans wide. With four spans per cubic basis support, (8−3)(4−3)(4−3) = 5 complete fine supports fit inside. The resulting 544-function basis integrates on 240 leaves, with 112 unchanged coarse leaves and 128 fine leaves. Uniform refinement throughout would give 2,299 scalar functions. All 160 exterior leaf sides retain a source design face. At q4 the common-partition Galerkin relative defects are 9.54e−17 for stiffness and 2.66e−15 for load. After applying the feature clamp they are 9.76e−17 and 2.69e−15 respectively. The transfer commutes exactly with the essential-condition lift. The prepared-to-local map has exact rank 539, and all 544 local functions are linearly independent in the uniform fine space. An independent Cox–de Boor evaluator also checks shape values, derivatives and geometry on 6,480 samples of the 240 leaves. [Independent basis evaluation](thb-function-values.json)

| Frozen Gaussian load | Free vector DOFs, prepared → local | Energy-norm error, prepared → local | Converged quadrature |
|---|---:|---:|---:|
| Center inside refinement | 1,386 → 1,401 | 25.6753% → 21.2992% | 10 |
| Center outside refinement | 1,386 → 1,401 | 28.5439% → 28.5430% | 12 |

The first load gains 17.04% in relative energy-norm error with 15 additional free displacement coefficients. The same mark gains just 0.00339% for the off-center load. The maximum converged free residual is below 1e−11 and the independent load/virtual-work defect is below 5.5e−12. Both fields evaluate on exactly the same geometry and quadrature points. The [quadrature histories](thb-evidence.json) include error, physical observables, work and timing.

The ordinary mechanical entry point accepts a saved local field and a prescribed full-feature clamp, surface tractions and constant body load. The result stores all displacement coefficients in the replayed field numbering.

```sh
"$cli" solve-job build/parent/local.json tests/references/thb_end_traction.job.json build/parent/solution.json
```

For the 1 MPa axial traction the native solve returns 9.19306904 N on the selected end. The opposite support reaction balances it to less than 1e−12 N per component. Across all 1,632 displacement coefficients the maximum difference from the exactly prolonged affine field is 7.16e−12 mm. An unknown design-face name is rejected with the prior result intact. [Feature-load audit](traction-audit.json)

## Acceptance against the revised scope

| Requirement | Evidence and scope |
|---|---|
| A01, A02 — parent and continuity | Whole-net reconstruction with the original knot intervals, sign and positive gauge; nonuniform and reversed source tests. Parent-coordinate C2 on the strip, physical gradient checks on selected multipatch joins. |
| A03, A05 — basis evaluation and unity | Independent Cox–de Boor and Bernstein evaluation on the canonical parent; independent hierarchical active-set, direct basis, derivatives and mass replay. |
| A04 — admissibility | The selected hierarchical profile admits one interior dyadic level and exact parent lineage. T-mesh anchors, boundary extension, and a general three-dimensional T-junction closure algorithm are deferred. |
| A06, A07 — preservation and transfer | Separate four-component homogeneous data and weight-scaled physical extraction. Inclusion checked on every common integration leaf. The local stiffness/load test adds an independent numerical transfer check. |
| A08, A09 — local support and integration | A box of 4 × 2 × 2 prepared cells is split in all three coordinates. The two exterior thickness layers remain coarse. Integration respects source and field span boundaries. |
| A10, A11 — mechanics | Native dynamic assembly, six rigid modes on the strip, clamped positivity, Galerkin consistency and native body/surface loads. The frozen C2 cube convergence criterion on levels 2–4 fails. See the cube results below. |
| A12, A13 — selection and admission | Compact source and plan replay, complete design-face coverage, separate source validity and field eligibility, rejection of damaged lineage and invalid marks before publication. |
| A14 — regressions | Existing source models, Q2 mechanics, exact axial knot reduction and crossed-interface mechanics retain their independent replay gates. |

The 68-model census measures continuity within extrusion columns. The projected lateral constraints implement physical C1 of the selected displacement space with unchanged geometry. A patchwise reparameterization, a general multi-level error estimator, an unstructured T-spline and a general mortar solver remain outside this qualified path.

## Automatic lateral selection on the curved parts

The optional native `matched_chart_jets_v1` plan checks all four homogeneous components of the transverse derivative on every lateral face. After orienting the two nets, the ratio of normal interval lengths must be one positive constant across the entire face. The field compiler then matches the corresponding derivatives of the weighted displacement numerator. A rational coefficient oracle independently recovers the same faces and ratios. A second oracle checks the physical gradient of every affected basis function as a polynomial identity. These witnesses establish local chart compatibility. Historical common-parent provenance requires separate evidence.

| Part | Whole faces | Inherited axial / compatible lateral / retained C0 | Free vector DOFs, C0 → axial → selected | Relative energy error, C0 → selected |
|---|---:|---:|---:|---:|
| Crossed G2 plate | 46 | 16 / 4 / 26 | 1,236 → 924 → 864 | 5.4335% → 5.4367% |
| S3 face pocket | 191 | 57 / 16 / 118 | 2,679 → 1,866 → 1,710 | 3.6405% → 3.6417% |
| NURBS multiloop plate | 75 | 21 / 26 / 28 | 1,056 → 792 → 576 | 1.3936% → 1.3995% |

The three complete-face physical checks cover 48, 231 and 416 column/layer identities respectively, all exactly zero. Each mechanics comparison preserves the geometry and polynomial degrees. The off-axis bores and cut faces keep their original boundary tags. These results are for the prescribed manufactured displacement; other loads can excite different modes. [Curved-part numerical histories](curved-chart-evidence.json)

```sh
target/release/solidcraft-trivariate export-inherited-space docs/results/s3-face-pocket.model.json tests/references/matched_chart_jets.plan.json build/s3-matched-field.json
python3 scripts/check_curved_supports.py
```

The ordinary `solve-job` route also uses this saved S3 field. A unit axial stress acts on the stock top, the bore's two axial faces and the pocket floor, with the stock bottom clamped and Poisson ratio zero. Against the independently integrated source volume, the strain energy is 2,656.87707357 N·mm (relative defect 1.78e−13). The maximum error across 2,232 extracted local displacement controls is 1.49e−10 mm. The axial load is 4,709.18901225 N, with a relative global force imbalance of 6.69e−13. [Solution and reaction audit](s3-stress-audit.json)

## Smooth-cube convergence decision

On the affine unit cube, the smooth axial sine solution gives maximum-continuity cubic rates of 2.3884 in energy norm and 2.9415 in displacement L² between two and four cells per direction. These miss the preregistered lower limits of 2.4 and 3.2. At four cells the C²/C⁰ energy-error ratio is 2.096, above the allowed 1.5. Both spaces pass the quadrature stability check. The separate four-to-eight-cell replay gives rates of 3.1225 and 4.2599 for C² and 2.9852 and 4.0316 for C⁰. The original S5 convergence gate remains failed. [Frozen protocol and interpretation](../support-census/cube-convergence-results.md) · [Native results](cube-convergence.json)

## Which pairs recover a whole common map?

On the 68 polygon models, 881 of 3,136 face-adjacent pairs recover a single tensor Bézier parent. Splitting each parent reproduces both complete input nets exactly. The 881 include all 617 axial joins and 264 lateral joins (131 along the first chart's u-axis and 133 along v). The counts describe individual interfaces; two successful pairs may share a cell. A maximal simultaneous aggregation requires a separate compatibility check.

Of the nine axial Boolean states, 20 of 132 face pairs recover such a parent. The source NURBS strip illustrates a different case. Its two cubic spans share a true NURBS parent with a simple interior knot and C² continuity. They fail the stricter test for two restrictions of one Bézier polynomial. On S3, 62 pairs have a whole-net parent, while 73 interfaces support the measured C¹ field. Each test therefore answers a different design or analysis question. [Per-model polygon and Boolean evidence](../support-census/pair-parent-expanded.json) · [Curved models and FieldCraft](pair-parents.md)
