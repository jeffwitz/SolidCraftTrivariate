# Codex specification for recipe derived trivariate T-splines

**Project:** [SolidCraftTrivariate](https://github.com/jeffwitz/SolidCraftTrivariate)
**Specification date:** 27 November 2026
**Reviewed reference:** `3f9d2d915c918121503668691b489a39781c6dd5`
**Proposed repository destination:** `docs/recipe-driven-tsplines.md`
**Assignment:** implement the mandatory S0–S6 delivery below, using the existing native kernel and constructive recipes.

## 1. Objective and architectural decision

Extend the existing constructive trivariate kernel so that recipes preserve their generating spline supports, recover continuity inherited from those supports, and compile an admitted support into an analysis-suitable trivariate T-spline representation with exact local refinement.

The first delivery must demonstrate the complete chain:

> Existing recipe → exact parent support → inherited spline continuity → native trivariate T-spline basis → local three-dimensional refinement → exact geometry and field transfer → native elasticity calculation.

The project already constructs related Bézier maps from known curves, intervals, orientations and geometric operations. Use that constructive information. A general inverse reconstruction from an unordered collection of blocks is unnecessary for the first admitted families.

The geometric model remains volumetric throughout construction. Each admitted material cell represents its entire parameter domain. Preserve the existing guarantees concerning fidelity, regularity, injectivity, coverage and semantic boundaries. The T-spline layer extends the representation and analysis spaces; it does not replace these guarantees.

The long-term objective is one logical volumetric model carrying adaptive spline spaces, with artificial C0 interfaces removed wherever the construction supports the requested continuity. This delivery must already remove such interfaces on a qualified existing family. Arbitrary smooth gluing between independent supports is a separately qualified extension.

### 1.1 Meaning of completion

S0–S6 are implementation work, not a request for another proposal. Complete the native functionality, CLI, independent checks, mechanical demonstration and documentation. A prerequisite milestone alone does not complete the assignment. If a genuine mathematical or access blocker prevents a mandatory case, preserve the reproducible evidence and report precisely which acceptance criterion remains open.

In particular, completion requires:

- an actual T-mesh, anchors, local knot vectors and basis functions;
- exact Bézier extraction and refinement operators;
- a rational curved example with nonconstant weights;
- refinement localized in all three coordinates;
- use of that basis by the native mechanical assembly;
- explicit accounting of any remaining interfaces and their continuity.

## 2. Starting point and required reuse

Read `AGENTS.md`, `docs/architecture.md`, `docs/LUNA.md` and the relevant contracts in `docs/kernel-requirements.md` before editing. Inspect the actual working revision: the reference above is a reviewed baseline, not an instruction to reset a newer checkout. Preserve unrelated user changes.

Repository instructions require Rust production functionality, independent oracles, English code and documentation, versioned formats, audited atomic adoption and `make check` at integration. Follow these requirements. Add the new spline profile explicitly to the architecture and requirements rather than silently changing the meaning of existing H1/C0 contracts. [R1–R2]

### 2.1 Existing capabilities

| Existing implementation | Role in this delivery |
| --- | --- |
| `src/nurbs_curve.rs` | Exact curve extraction and source knot intervals. |
| `src/nurbs_extrusion.rs`, `src/nurbs_sketch_json.rs` | Paired source curves, exact strip volumes, source intervals and original recipe data. |
| `src/nurbs_cross.rs` | Synchronization, restriction and degree elevation of generating curves. |
| `src/primitive_template.rs`, `src/g2_cross_plate.rs`, `src/sketch_partition_plan.rs` | Constructive plans, orientations, event provenance and generating-curve information. |
| `src/sketch_junction.rs` | Exact rational curve jets and G-continuity classification. |
| `src/nurbs_tensor_assembly.rs` | Existing native C0 tensor assembly, homogeneous gauge alignment, degree elevation and exact extraction. |
| `src/nurbs_knot_reduction.rs` | Existing exact degree-one axial knot removal, with subfaces, lineage and original analysis knots retained. |
| `src/incidence.rs`, `src/mapped_seams.rs`, `src/mapped_cut_propagation.rs` | Interface orientation, parametric restrictions and subdivision conventions. |
| `src/spatial_certificate.rs` and family audits | Geometric certification and admissibility. |
| `src/adaptive_partition.rs`, `src/adaptive_transfer.rs` | Useful provenance, budget and transfer conventions; the existing Q2 refinement is a separate space. |
| `src/constraint_system.rs`, `src/elasticity.rs`, `src/mechanical_job.rs` | Boundary constraints, mechanics and result conventions. |
| `docs/results/crossed-interface-iga/smooth_bending.py` | Independent experimental reference for weighted C1 coefficient relations inside charts. |

Reuse these capabilities instead of duplicating them. Keep the native tensor assembler's C0 contract and the axial reduction API working. Its supplied tensor layout is not a general parent parameterization: `assemble_cover` currently expects unit cell intervals. Preserve actual source intervals in the new support contract. [R3–R5]

The current `AnalysisBlock` and `elasticity::BasisSample` specialize the native analysis path to 27 Q2 functions. The new path must introduce a variable-size extracted-element interface; merely exporting T-spline geometry while retaining independent Q2 fields would leave the assignment incomplete. [R6]

### 2.2 Existing progress that must remain visible

The axial reduction reports already distinguish geometric compaction from field dimension. The C1 bending experiment already distinguishes continuity inside charts from C0 chart-to-chart joins. Preserve these distinctions in the new results. Do not present either capability as newly invented or as an existing general T-spline implementation. [R4–R5]

## 3. First admitted T-spline profile

Implement the named profile `structured_tricubic_m2_v1`:

- dimension three and degree `(3,3,3)`;
- one regular, injective tensor parent covering the entire admitted material domain;
- an initially uniform tensor parameter mesh, allowing a declared affine rescaling of each coordinate;
- simple interior knots and open boundary treatment;
- globally structured analysis-suitable T-splines;
- cyclic directional bisection with the documented admissibility closure, using `m=2`;
- a positive rational denominator preserved as a function by refinement.

Use Morgenstern's globally structured three-dimensional construction for the mesh and independence contract. Use the boundary/active-function conventions and refinement Algorithm 2.1 of Gantner–Praetorius. Their Section 3.3.2 establishes the required nesting, citing Morgenstern's 2017 thesis, Theorem 5.4.12. Their Section 4.1 describes the fixed-positive-denominator rational setting. Record the precise source versions and algorithm-to-code mapping. [R9–R10]

The 2016 article alone is not the nesting reference. Its independence result and the later nesting statement have distinct roles. The elliptic convergence analysis in the reference is also not a validation of this project's elasticity implementation.

The following are separate extension profiles: nonuniform initial meshes beyond the admitted affine normalization, repeated interior knots, periodic supports, extraordinary topology, and smooth joins between distinct parent charts. Preserve their existing geometry and analysis paths. Report conversion eligibility without marking an already valid source solid invalid.

Geometry of degree below three may be elevated exactly for this profile. Admission additionally requires that all four homogeneous geometry components belong to the initial spline space. Degree elevation alone cannot turn a general C0 patch assembly into a smooth admissible parent.

A two-dimensional T-spline extruded through the full height is a useful intermediate case. The final local-refinement demonstration must have different refinement patterns at different positions along the third coordinate.

## 4. Constructive support and lineage contract

### 4.1 Separate identities

Distinguish:

1. the logical material body;
2. the generating support and its revision;
3. a restriction of that support used by a source block;
4. a design face or semantic subface;
5. a T-mesh cell and a basis-function anchor;
6. an integration or certification cell.

Extracting a support into more Bézier cells must not create new design-volume identities. Refinement must not silently rename a selected loading face.

### 4.2 Minimum records

Introduce equivalent typed records, adjusting names to the existing codebase when appropriate:

| Record | Required information |
| --- | --- |
| `TensorSupport` | Stable ID, source revision, degrees, actual knot vectors, homogeneous controls, parameter domain, generating recipe and certificate references. |
| `ParametricRestriction` | Parent ID, three exact intervals, axis permutation, reversals and the explicit affine coordinate map. |
| `SourceSpan` | Source block ID, restriction, homogeneous gauge, generating operation and event provenance. |
| `ContinuityEvidence` | Interface, coordinate frame, direction, claimed order, method, assumptions and verification status. |
| `InterfaceRole` | Same-parent subdivision, join between supports, intentional material interface, boundary or periodic identification. |
| `SplineVolume` | Support, basis profile, active basis records, global homogeneous coefficients, extraction data and semantic-face references. |

Existing recipes retain their schemas. Store new support metadata in a versioned associated representation or in a backward-compatible versioned extension. Legacy inputs without this metadata must still build through their existing path. A known recipe adapter may reconstruct the metadata exactly; it must then verify it.

Use canonical exact rational serialization. If a new encoding is needed, represent numerator and positive denominator as decimal integer strings, not binary floating-point numbers. Respect current checked arithmetic and storage limits.

### 4.3 Required identity

For every admitted source cell:

\[
F_e = F_{\mathrm{parent}}\circ\rho_e.
\]

Prove this by exact restriction, orientation and coefficient comparison after a documented positive constant homogeneous gauge. This is a sufficient initial admission test. Do not claim that failure of this particular test proves that two rational maps represent different physical solids.

Choose one canonical homogeneous gauge for the parent at admission and keep it fixed throughout refinement. Record each source block's positive constant gauge relative to that parent. A local comparison gauge must not redefine the parent's denominator or the global field weights.

Keep the real parent interval lengths. Differentiation through `rho_e` must include their scale factors and orientation. Replacing every source interval by an unrelated unit interval destroys information required for derivative compatibility.

A construction support may extend beyond the final material in future Boolean adapters. That fact must never activate the absent material. The first T-spline profile requires a complete parent material domain; partial parent occupation belongs to a later explicit profile, without introducing cut-cell integration implicitly.

## 5. Continuity and exactness requirements

### 5.1 Recover inherited continuity first

Prioritize reusing known parent knots and controls. Reconstructing the original parent from its recipe is preferable to solving a general inverse fitting problem. General knot removal is useful only where parent data are unavailable and the exact relation can be established.

Classify interfaces before acting:

- same-parent subdivisions can inherit the parent's continuity;
- joins between coordinated constructions need the appropriate transverse relations;
- joins with insufficient evidence retain their existing behavior or return an explicit conversion refusal;
- physical interfaces keep the continuity required by the mechanical model.

Curve G2 classification is not a certificate of C2 parametric volume continuity. A smooth source direction does not certify a transverse join between different regions. The continuous target is qualified per interface and per direction.

### 5.2 Physical field continuity

For a claimed physical C1 scalar field across a chart interface, the required relations are equality of traces and equality of physical gradients:

\[
\widehat u_1=\widehat u_2,
\qquad
DF_1^{-T}\nabla\widehat u_1
=DF_2^{-T}\nabla\widehat u_2,
\]

evaluated at matched physical points. For a common regular parent, use its coordinate transitions to establish this relation. Test the relations for the basis space, not only for one displacement solution.

For the first tricubic profile, preserve the C2 continuity associated with simple interior knots. The minimum advertised smoothness for the admitted physical interior is C1. A claim of higher physical continuity must include the corresponding regularity of the geometry map and inverse.

### 5.3 Geometric certificates

Spline basis admissibility, positive Jacobian, injectivity and material coverage are separate properties. Exact restriction of an already certified parent may inherit its geometric evidence, with the chain-rule factor from `D rho_e` recorded. If controls or parameterization change, provide the corresponding new evidence.

Arithmetic overflow or exhausted certification budgets must produce a qualified unavailable/indeterminate result. They are not proofs of negative Jacobian or geometric incompatibility. Sampling remains a diagnostic, not a continuous-domain certificate.

## 6. Basis, extraction and rational transfer conventions

Fix the following matrix orientation in code and tests, or document one consistent equivalent.

### 6.1 Polynomial extraction

Let `N_A` be the nonrational T-spline blending functions. On an integration cell `e`, let `B_e,a` be the tensor Bernstein basis and define:

\[
N_A|_e=\sum_a E^e_{aA}B_{e,a}.
\]

`E_e` has Bernstein rows and active-global-function columns. If `H_A=(w_A P_A,w_A)` are homogeneous geometric controls, then:

\[
H^e=E_e H,\qquad w^e=E_e w.
\]

Publish degrees, active global IDs, the sparse extraction operator, local parameter restriction and any source gauge. Verify the four homogeneous components against the existing recipe maps exactly on a common extraction partition.

### 6.2 Rational field extraction

The rational functions are:

\[
R_A=\frac{w_A N_A}{W},\qquad W=\sum_Aw_A N_A>0.
\]

For physical field coefficients `q_A`, extraction into local rational Bernstein coefficients uses:

\[
C^e=\operatorname{diag}(w^e)^{-1}E_e\operatorname{diag}(w),
\qquad q^e=C^e q.
\]

Require `w_A > 0` and every extracted `w_e,a > 0` before applying these diagonal inverses; a positive denominator function alone does not establish positivity of every coefficient. Keep the polynomial extraction `E_e` and weighted extraction `C_e` distinct. For a scalar local operator, use `C_e^T K_e C_e` with the matching geometric map. For three-component elasticity with component-interleaved numbering, use `C_e,v = C_e tensor I_3` and `C_e,v^T K_e C_e,v`. An equivalent direct evaluation of active global functions is acceptable.

The implementation must expose the complete basis, including functions needed to represent nonzero geometry and boundary traces. Do not accidentally use only the homogeneous-Dirichlet subspace from a PDE reference as the geometric basis.

### 6.3 Exact refinement

For a refinement, return a sparse matrix `P` with new rows and old columns satisfying:

\[
N^{\mathrm{old}}_j=\sum_iP_{ij}N^{\mathrm{new}}_i,
\qquad H_{\mathrm{new}}=P H_{\mathrm{old}},
\qquad w_{\mathrm{new}}=P w_{\mathrm{old}}.
\]

The rational field prolongation is:

\[
T=\operatorname{diag}(w_{\mathrm{new}})^{-1}
P\operatorname{diag}(w_{\mathrm{old}}),
\qquad q_{\mathrm{new}}=Tq_{\mathrm{old}}.
\]

Preserve the denominator as the same function, even though its coefficients change. Require the positive-weight contract of the admitted profile. Avoid arbitrary weight renormalizations that alter the represented function.

Check every new coefficient `w_new,i > 0` exactly before forming `T`. Positivity of the preserved denominator function alone is insufficient to justify an inverse diagonal of its coefficients.

Prove these identities using exact extraction on a common Bézier partition. Numerical SVD and point samples can diagnose defects; they must not define the production basis or certify nesting. Use constructive knot operations or the local refinement relations of the selected family. A small local system of exact Bézier coefficient identities is also acceptable for computing a transfer column: bound it by the affected supports, solve it with checked rational arithmetic, and verify every coefficient identity independently. A global sampled fit is not an acceptable substitute.

Geometry compaction and changing the analysis space remain distinct operations. Removing a field knot generally reduces its space and does not preserve an arbitrary old field. Exact refinement must preserve every field in the old space. Approximate coarsening is outside S0–S6.

## 7. T-mesh and refinement implementation

Represent exact index coordinates, parameter coordinates, levels, active cells, anchors and local knot vectors. Use the selected cyclic bisection schedule and its degree-dependent neighborhood closure. Do not replace that schedule with an unrestricted octree rule under the same profile name.

Implement the boundary extension and active-anchor conventions consistently. Test face, edge and corner cases independently. Audit the final structure in addition to recording the refinement history.

Each refinement result must record:

- explicitly marked cells;
- cells additionally refined by admissibility closure;
- old and new basis counts;
- changed anchors and local knot vectors;
- affected basis-support region;
- exact homogeneous and rational-field transfer operators;
- numerical and arithmetic budgets consumed;
- unchanged semantic references and certificate dependencies.

A support can change beyond the cells directly split. Define the affected region from both closure and changed basis supports. Test stability outside that region instead of promising that every unmarked cell is untouched.

The T-mesh is not necessarily the integration partition. Build integration cells resolving all polynomial breaks of active basis functions, geometry and boundary traces. Keep temporary common partitions used for proofs separate from the production basis. A global tensor overlay is acceptable as a small independent oracle; it must not become the hidden production refinement algorithm.

Use explicit cell-count, level, support-count and arithmetic budgets. A refused operation preserves the previous accepted model and output file. Report the limiting resource and the requested refinement; never silently replace local refinement by global uniform refinement.

## 8. Extracted-element interface and mechanics

Introduce an element contract with variable-length active function IDs, values and gradients. Geometry evaluation, field evaluation and quadrature are separate responsibilities.

At a common integration point, the interface supplies:

- physical position, Jacobian and determinant from the admitted geometry;
- active scalar basis values and physical gradients;
- mappings to global scalar and vector coefficients;
- boundary/subface provenance when applicable.

Keep the current Q2 implementation available through an adapter or shared assembly core. Extend the smallest necessary assembly layer; do not replace unrelated solvers, sensitivity code or UI adapters in this task.

For the new path, provide native Rust stiffness and load assembly, essential constraints, solve and result export. Reuse established units: mm, MPa and N. Python may provide an independent reference solve, not the only working new solver.

T-spline coefficients are not generally nodal values. Apply boundary conditions through the qualified trace representation and existing semantic selections. The mandatory first mechanical cases use complete semantic faces and exactly representable prescribed traces. Preserve current partial-face cases as regressions; unsupported new partial-trace configurations must be identified explicitly.

Evaluate geometry and field at the same physical integration point. When their parameter domains differ, compose the recorded restrictions. Do not substitute a compact geometric knot vector for the field knot vector without an explicit space change.

## 9. Implementation lots and exit criteria

### S0 — Baseline and executable contract

Record the actual starting revision, worktree state, relevant files, source versions, initial corpus, thresholds, budgets and independent oracle design. Inspect newer work before adding modules.

Add this delivery profile to the architecture and a concise machine-readable corpus manifest. Record which assertions are exact, which follow from the selected admissibility theorem, and which are numerical diagnostics.

Run or reuse a documented baseline of the required integration checks. Preserve any unrelated pre-existing failure with its reproduction. Avoid repeated full campaigns after every small edit.

**Exit:** scope, required recipes, thresholds and provenance are fixed before implementation measurements.

### S1 — Parent supports and inherited continuity

Implement a common support/lineage representation and the adapter for `nurbs_sketch_strip_v1`. Reuse the source NURBS data retained in `nurbs_sketch_json.rs`.

The mandatory first case is `tests/references/nurbs_strip.recipe.json`. Its curves have degree three and knot vector:

```text
[-2, -2, -2, -2, 1, 4, 4, 4, 4]
```

They share five positive weights and compatible generating coordinates. Reconstruct the exact parent volume of degree `(3,1,1)` with a `5 × 2 × 2` control grid. The interior source knot is simple and carries C2 continuity. [R7]

Expected structural accounting before any tricubic elevation:

| Representation | Homogeneous controls or occurrences |
| --- | ---: |
| Two separate Bézier blocks | 32 occurrences |
| C0 tensor assembly sharing the interface plane | 28 controls |
| Original parent NURBS support | 20 controls |

These are different representation counts, not a solver-performance result. Verify the counts and the exact extracted maps. Preserve semantic faces and source-span identities.

Also test actual nonuniform source intervals, using an admitted synchronization case such as `nurbs_cross_bores.recipe.json`. This tests lineage and derivative scaling; it does not automatically admit that entire multi-support solid to the first T-spline profile.

**Exit:** exact parent recovery, C2 source evidence, reversed/permuted restrictions, nonuniform metric handling, round trip and forged-lineage rejection.

### S2 — Variable-size extracted spaces

Implement the extracted-element and rational-evaluation contracts. Support both the legacy Q2 adapter and the recovered parent tensor space.

Build an independent evaluator using the original curve/tensor definition, rather than calling the new extraction code. Check values, derivatives, rational weights and physical gradients.

**Exit:** the Q2 adapter reproduces legacy assembly at identical quadrature and numbering; the parent space reproduces coordinates, affine physical displacement and rigid motions. Geometry-only representation changes preserve a fixed field discretization.

### S3 — Native T-spline basis and extraction

Implement `structured_tricubic_m2_v1`: mesh, boundary treatment, anchors, local knot vectors, basis evaluation, audit and extraction.

Convert the recovered strip support to the initial tricubic tensor case by exact degree elevation in the linear directions. Its homogeneous geometry must remain exactly representable. With the original two spans and open knot vectors, the complete initial tricubic tensor basis has `5 × 4 × 4 = 80` scalar functions before boundary constraints. Keep the compact 20-control parent available as provenance and a geometry reference.

For this fixture, the affine normalization `u_index = (u_source + 2) / 3` maps the two source spans to `[0,1]` and `[1,2]`. The initial index domain is therefore `[0,2] × [0,1] × [0,1]`. Record this map and its derivative factor explicitly. The 80-function count concerns the initial complete tensor-limit basis, not a promise about the count after local refinement.

**Exit:** tensor-limit equivalence, exact homogeneous reconstruction, positive denominator, full boundary basis and structural admissibility checks. S3 alone is not a completed local T-spline delivery.

### S4 — Local 3D refinement and exact transfer

Implement admissibility closure, basis updates, exact extraction updates and the operators `P` and `T`.

Use a reproducible refinement request on the curved strip parent after sufficient exact uniform pre-refinement to leave a remote unaffected region. Freeze that preparation before measuring locality. Apply successive local marks through the cyclic directions so the final mesh cannot be expressed as a refined section extruded uniformly over the full height.

Test an interior request, a boundary-face request, an edge request and a corner request. Include adjacent marks whose closure interacts, duplicate marks and budget exhaustion.

**Exit:** genuine 3D T-junction structure, exact geometry and denominator preservation, transfer of arbitrary old fields, composition of transfers, stable remote data and atomic refusal. Compare local and uniform refinement costs at comparable resolution without imposing an invented universal gain.

### S5 — Native mechanics and controlled comparisons

Run the new basis through native three-dimensional elasticity. Use three complementary tests:

1. An affine physical field on the curved rational parent, with compatible exact boundary data and loads.
2. A non-affine manufactured elasticity solution on a regular affine cube, with analytically derived body force and boundary data, to assess convergence independently of geometric approximation.
3. A curved recipe-based mechanical case using the same material, loading and boundary definitions for the C0, restored tensor-spline and locally refined T-spline spaces.

A localized manufactured solution or a frozen geometric marking request may drive the first adaptivity study. A new general residual estimator is not required for this delivery.

For the smooth manufactured cube case, preregister at least three uniform refinement levels, the error norms, quadrature policy and expected asymptotic behavior of the tricubic space. Fix the finite-level convergence acceptance bands in S0. Passing S5 requires those checks to pass; a table of measurements without a convergence decision is insufficient. Do not impose the same rate criterion on arbitrary local refinement sequences or singular physical solutions.

Compare energy error, displacement error and a stated stress norm where a reliable reference exists. Report degrees of freedom, integration work, matrix nonzeros, runtime and memory separately. Higher continuity alone is not a guarantee of lower error at fixed spans.

Use existing HEX20 comparisons when the physical case matches. Keep the current crossed-interface studies reproducible. Do not attribute old reference numbers from a different geometry or field space to the new demonstration. A superiority claim over FEM requires an actual non-IGA comparison at the claimed accuracy and cost.

**Exit:** new native solve, mechanical identities, measured convergence, rational curved case and reproducible cost/error tables.

### S6 — CLI, persistence and release evidence

Expose the common implementation through a versioned CLI and library API. Proposed command names are:

```text
compile-spline RECIPE PLAN OUTPUT
audit-spline MODEL
refine-spline MODEL REQUEST OUTPUT
export-spline-analysis MODEL JOB OUTPUT
solve-spline MODEL JOB OUTPUT
```

These are proposed additions, not commands assumed to exist at the reviewed baseline. Follow current CLI conventions when selecting final names. Store the profile, source identity, admission evidence, refinement history and basis/extraction definitions needed for an independent replay.

Provide a strict smooth-conversion mode for the admitted parent case. For a complex input outside that mode, return an eligibility report naming unsupported interfaces while preserving the accepted source model. Exercise this reporting on two existing multi-support recipes, including a box–cylinder Boolean case and a synchronized NURBS construction, frozen in S0. Identify inherited parent restrictions, every retained C0 join, and the exact reason a full smooth conversion is outside the admitted profile. This is a required integration check, not a claim that these complete solids have been converted to the first single-parent T-spline profile.

Use the existing audited publication and atomic replacement pattern. Test that invalid requests leave an existing destination byte-for-byte unchanged. [R8]

Add `check-spline-parent`, `check-tspline` and `check-tspline-mechanics` targets, or equivalent clearly documented targets. Put the bounded mandatory regression set into `make check`; retain larger performance sweeps as explicit commands.

**Exit:** the complete recipe-to-solve path is callable and replayable; mandatory integration checks pass and reports identify the actual tested revision.

## 10. Acceptance matrix

| ID | Required evidence | Failure example that must be detected |
| --- | --- | --- |
| A01 | Exact reconstruction of the strip parent and source blocks | Wrong source interval or altered interior control. |
| A02 | Continuity in parent coordinates and corresponding physical field regularity | Ignoring unequal interval lengths or equating G2 with parametric C2. |
| A03 | Independent direct-basis versus extracted evaluation | Transposed extraction or swapped parameter axes. |
| A04 | Admissible T-mesh, valid anchors and local knot vectors | Unsupported imported connectivity or broken boundary extension. |
| A05 | Full basis and rational partition of unity | Omitting boundary functions or inconsistent weights. |
| A06 | Exact four-component geometry and denominator preservation | Matching sampled positions while changing the rational map. |
| A07 | Exact nesting and weighted arbitrary-field transfer | Applying the homogeneous operator directly to physical coefficients. |
| A08 | Local refinement genuinely varying in the third coordinate | Refinement of a 2D section repeated through the full height. |
| A09 | Correct integration partition | Integrating across an unresolved polynomial break. |
| A10 | Native elasticity, affine reproduction and six rigid modes | Missing modes, extra mechanisms or incorrect physical gradients. |
| A11 | Galerkin transfer identities and numerical convergence | Similar end displacement hiding an incorrect space transfer. |
| A12 | Stable semantic selections and audited atomic updates | Lost loading face or partially overwritten rejected output. |
| A13 | Exactness/validity/eligibility statuses kept separate | Treating an arithmetic limit as invalid source geometry. |
| A14 | Legacy tensor, knot-reduction and mechanics regressions | New T-spline path altering accepted old recipes. |

### 10.1 Exact tests

Use exact arithmetic for knot coordinates, restrictions, extraction coefficients, homogeneous comparisons and refinement identities in the admitted rational scope. These tests require equality, not a floating-point tolerance.

Verify the basis-function inclusion column by column on a common extraction partition. This proves transfer of every old field; random vectors are additional diagnostics. Verify composition by comparing operators or their exact functional action after bringing outputs to a common basis.

Check continuity relations on every affected interface and at intersections of interface constraints. Use rational derivative identities with actual interval factors. A finite list of evaluation points alone is insufficient.

### 10.2 Numerical thresholds

The following are initial thresholds for bounded, normalized regression cases. Freeze the case scales and error normalizations in the S0 manifest. If a threshold is unsupported by the selected arithmetic or conditioning, establish the reason before the qualification run and version the protocol; do not relax it after a failed run just to obtain acceptance.

| Diagnostic | Initial criterion |
| --- | --- |
| Direct versus extracted positions | Maximum error divided by a declared characteristic length ≤ `1e-12`. |
| Basis partition of unity | Maximum absolute defect ≤ `1e-12`; exact identity where provided algebraically. |
| Values and scaled first derivatives | Maximum normalized defect ≤ `1e-11`. |
| Same-space stiffness and load comparison | Relative Frobenius/Euclidean defect ≤ `1e-10`, at identical quadrature and compatible numbering. |
| Free-system equilibrium | Relative residual ≤ `1e-9`, with denominator and zero-load convention recorded. |
| Affine patch field | Normalized displacement/gradient defect ≤ `1e-9` on the selected well-conditioned cases. |
| Quadrature replay | Change of the reported energy and displacement observable ≤ `1e-8` before attributing changes to the space. |

Report actual residuals, not only pass flags. Scale derivative errors by the stated parameter and physical lengths. Use absolute tolerances for zero reference quantities.

For free-body mechanics on one connected regular test domain, verify the six known rigid-mode vectors and inspect the small-system spectrum with a preregistered scaled threshold. After sufficient essential constraints, check the expected positive definiteness. Do not use a floating-point rank test as the proof of T-spline basis independence.

### 10.3 Galerkin consistency

With vector prolongation `T_v = T tensor I_3`, check:

\[
K_{\mathrm{old}}=T_v^T K_{\mathrm{new}}T_v,
\qquad f_{\mathrm{old}}=T_v^T f_{\mathrm{new}}.
\]

For this algebraic check, integrate both spaces on the same common cells with the same quadrature nodes. Check unconstrained matrices first, then compatible boundary liftings. This avoids mistaking different quadrature errors for a transfer defect. Separately assess quadrature convergence and approximation error.

### 10.4 Negative cases

Include missing or falsified parent lineage; wrong axis orientation; unequal intervals incorrectly treated as equal; nonpositive denominator; arithmetic overflow; budget exhaustion; unsupported repeated knots or periodic topology; and a valid C0 join whose transverse continuity is not established.

A refused smooth conversion must preserve the original accepted geometry. Distinguish demonstrated incompatibility, unsupported profile and indeterminate computation. A required positive case that is refused counts as a delivery coverage failure.

## 11. Implementation organization

Suggested module boundaries, to adapt after inspecting current symbols:

| Concern | Suggested modules |
| --- | --- |
| Parent supports and provenance | `spline_support.rs`, `spline_lineage.rs` |
| Inherited continuity | `spline_continuity.rs` and recipe adapters |
| Generic element interface | `extracted_space.rs`, `spline_evaluation.rs` |
| T-mesh and admissibility | `tspline_mesh.rs`, `tspline_audit.rs` |
| Basis and Bézier extraction | `tspline_basis.rs`, `tspline_extraction.rs` |
| Refinement and transfer | `tspline_refinement.rs` |
| Versioned exchange and CLI | `spline_json.rs`, additions to `main.rs` and crate exports |

The module names are guidance, not a requirement to create empty abstractions. Reuse checked rational operations and sparse structures already available. Avoid coupling the geometry library to the gallery or UI.

Document the source, version, license and numerical contract before copying upstream code. Independent implementation of the cited mathematics is acceptable. Do not change the project license or add an unqualified geometry backend as part of this task.

Use independent Python oracles for source-curve evaluation, exact coefficient identities and small mechanical comparisons. They must not obtain their expected result by calling the production function under test. A high-resolution global overlay can be an oracle on small cases, with its cost clearly separated from production timings.

## 12. Deliverables and execution protocol

Deliver:

1. Native Rust library capabilities and versioned CLI commands for S0–S6.
2. The source-support and T-spline data contracts, including matrix orientations and rational conventions.
3. A frozen corpus manifest containing positive cases, negative cases, thresholds and budgets.
4. Independent oracles and bounded regression tests.
5. An English report under `docs/results/recipe-driven-tsplines/` with exactness, continuity, refinement and mechanics results.
6. Reproduction commands, environment, source recipes and tested revisions.
7. Clear remaining-profile coverage: multiple parents, periodic joins, repeated knots, Boolean changes and any remaining C0 interfaces.

Keep regenerable outputs under `build/`. Track only the intended recipes, concise reports and regression assets. If published geometry changes, rebuild it from its source recipe and run the published-geometry checks.

Run the existing integration gates appropriate to affected paths, including:

```sh
make check
make check-cylinder-iga
make check-crossed-bending
```

`make check` already includes the published-geometry path at the reviewed baseline. Avoid repeating it needlessly. Also run the new parent, T-spline and mechanical targets. Record dependencies or unavailable external references honestly; do not report unexecuted checks as passed.

For numerical studies, report source/compact geometry controls, logical supports, T-mesh cells, integration cells, scalar functions, vector coefficients, free DOFs, nonzeros, quadrature points, assembly time, solve time and peak memory separately. Include closure amplification and coverage failures. Separate warm-cache and cold-build measurements when comparing updates.

### 12.1 Suggested final Codex report

- Implemented scope and reference commit.
- Reused capabilities and actual new modules.
- Required cases and acceptance criteria, with exact/numerical status.
- Continuity removed, continuity retained and the reason for each retained join.
- Geometry preservation and rational field-transfer results.
- Mechanical error/cost comparisons with their physical and numerical controls.
- Commands run, actual results and remaining reproducible blockers.

## 13. Extensions after the mandatory delivery

Extend one construction family at a time once S0–S6 pass:

- additional parent adapters for primitive sweeps and synchronized NURBS constructions;
- admission of nonuniform initial meshes under an explicit verified algorithm;
- periodic annular supports with qualified seam continuity;
- joins between different supports with explicit transition and transverse-derivative contracts;
- integration of the admitted spline compilation into Boolean/edit sessions, preserving existing lineage and semantic selections;
- intentional material interfaces, independent refinement on multiple supports and qualified trace coupling;
- error-driven refinement, coarsening and design sensitivities with their own transfer and differentiability contracts.

For Boolean operations, preserve existing parent relations where they remain valid and identify newly created joins. Do not claim automatic global smoothness after an arbitrary Boolean sequence. The useful next result is a growing set of constructive recipes for which that smoothness and refinement are actually qualified.

## 14. References and evidence

All repository references below are pinned to the reviewed commit. Codex must inspect later changes before implementation.

- **R1. Repository workflow.** [AGENTS.md](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/AGENTS.md).
- **R2. Existing contracts.** [Architecture](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/docs/architecture.md), [kernel requirements](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/docs/kernel-requirements.md).
- **R3. Native tensor assembly.** [nurbs_tensor_assembly.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/nurbs_tensor_assembly.rs).
- **R4. Exact axial reduction.** [nurbs_knot_reduction.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/nurbs_knot_reduction.rs), [reported geometry and torsion checks](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/docs/results/axial-knot-reduction/README.md).
- **R5. Existing continuity experiment.** [smooth_bending.py](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/docs/results/crossed-interface-iga/smooth_bending.py), [bending-enrichment report](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/docs/results/crossed-interface-iga/bending-enrichment.md).
- **R6. Current native analysis interface.** [analysis_space.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/analysis_space.rs), [elasticity.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/elasticity.rs), [nurbs_space.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/nurbs_space.rs).
- **R7. Mandatory parent-recovery case.** [nurbs_strip.recipe.json](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/tests/references/nurbs_strip.recipe.json), [nurbs_extrusion.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/nurbs_extrusion.rs), [nurbs_sketch_json.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/nurbs_sketch_json.rs).
- **R8. Publication and integration commands.** [main.rs](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/src/main.rs), [Makefile](https://github.com/jeffwitz/SolidCraftTrivariate/blob/3f9d2d915c918121503668691b489a39781c6dd5/Makefile).
- **R9. Morgenstern, P. (2016).** Globally Structured Three-Dimensional Analysis-Suitable T-Splines: Definition, Linear Independence and m-Graded Local Refinement. SIAM Journal on Numerical Analysis 54(4), 2163–2186. [DOI 10.1137/15M102229X](https://doi.org/10.1137/15M102229X), [open author text](https://arxiv.org/abs/1505.05392). Basis: structured three-dimensional admissibility, independence and refinement. Record the exact version used.
- **R10. Gantner, G.; Praetorius, D. (2020).** Adaptive IGAFEM with optimal convergence rates: T-splines. Computer Aided Geometric Design 81, 101906. [DOI 10.1016/j.cagd.2020.101906](https://doi.org/10.1016/j.cagd.2020.101906), [published open text](https://pure.uva.nl/ws/files/60294158/1_s2.0_S0167839620300935_main.pdf). Implementation references: Sections 2.1–2.3, Algorithm 2.1, nesting in Section 3.3.2, rationalization in Section 4.1. Its nesting discussion cites Morgenstern (2017), Theorem 5.4.12; the thesis itself is not treated here as independently inspected evidence.
- **R11. Beuchler, S.; Hiniborch, R.; Morgenstern, P. (2024).** Deal.t: an implementation of multivariate analysis suitable T-splines within the deal.ii framework. Engineering with Computers 40, 3901–3928. [DOI and open article](https://doi.org/10.1007/s00366-024-02002-1). Implementation reference for data organization, local functions and extraction; copying software requires a separate license review.

## 15. Startup instruction for Codex

Read this specification and the repository instructions. Inspect the current revision and identify any capability added after the reviewed baseline. Implement S0–S6 in dependency order, beginning with exact parent recovery for `nurbs_strip.recipe.json`. Reuse the existing native tensor assembly, axial reduction, geometry audits and semantic boundaries. Deliver a real rational trivariate T-spline basis with admissible local 3D refinement, exact homogeneous and physical-field transfer, and a native elasticity demonstration. Validate each milestone with independent evidence, then run the required integration gates. Report the actual implemented scope and results against the acceptance matrix.
