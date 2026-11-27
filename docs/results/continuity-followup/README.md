# Continuity, approximation and local refinement

This follow-up revisits three observations from the [inherited-space campaign](../inherited-continuity/README.md) and tests a wider local THB hierarchy. Geometry, scalar coefficients, free displacement coefficients and integration cells are reported separately.

## What the earlier measurements establish

On S3, the NURBS multiloop plate and the crossed G2 plate, selecting compatible axial and lateral chart derivatives reduced free vector DOFs by 36%, 45% and 30%. The relative energy-norm errors changed by at most 0.006 percentage points on the prescribed load. The selection preserves the original maps and design-face labels.

The earlier 21.33% pocket error came from the nine *curved* section maps of the published model, at one in-plane element per patch. The 48-DOF low-degree row and the coarse rational cylinder cannot establish a convergence rate. The original curved-pocket result and its source model remain available.

The normalized affine transition and AS-G¹ admit different data. The former succeeds in either direction on 626 of 2,519 lateral polygonal interfaces. The latter supplies exact linear gluing factors on all 2,519. The 617 axial interfaces have their own classification. These are local compatibility certificates. Approximation on a given solid still depends on the space and its refinement. [Census and exact witnesses](../transverse-as-g1.md)

## Matched dimensions on the smooth affine cube

| Free vector DOFs | Uniform cubic C⁰ / C² spans per axis | Relative energy error C⁰ / C² | C⁰ error divided by C² |
|---:|---:|---:|---:|
| 6,084 | 4 / 10 | 0.00555250 / 0.00073965 | 7.51 |
| 11,520 | 5 / 13 | 0.00288067 / 0.00032368 | 8.90 |

Both pairs pass the prospective 1.5 threshold on the new smooth Gaussian load. The error decreases with increasing dimension in both spaces. Check the quadrature-point counts before using these figures to reason about time or cost. On the original sine field the Q3 C² result at eight spans was already measured, and the original equal-cell S5 gate remains failed. [Protocol](cube-protocol.md) · [Paired native results](cube-matched.json) · [Original sine data and interpretation](../support-census/cube-convergence-results.md)

## Do all lateral C¹ conditions lock the bilinear pocket?

The bilinear control recipe generates nine trilinear solids. They define five distinct section patches, with four slanted radial seams and four footprint edges. The native vector experiment compares two fields on the same geometry, at p=3 and p=4 with n=1, 2, 4, 8 intervals per chart direction. Both have C¹ continuity inside charts and the same bottom-clamped linear height dependence. Only the interchart constraints vary. All three physical displacement components are free. The quintic manufactured field is exactly linear in height; all external tractions, including those on the pocket walls and floor, participate in the assembly.

| Section degree | Interchart continuity | Free vector DOFs at n=8 | Relative energy-norm error at n=8 | Energy / L² order, n=4→8 |
|---:|---|---:|---:|---:|
| 3 | C⁰ | 4,440 | 5.34389e−5 | 2.916 / 3.832 |
| 3 | C¹ | 3,780 | 5.64429e−5 | 2.979 / 3.917 |
| 4 | C⁰ | 9,528 | 3.23534e−7 | 4.009 / 5.004 |
| 4 | C¹ | 8,676 | 3.39158e−7 | 4.052 / 5.070 |

The C¹ energy-error multipliers relative to C⁰ fall from 1.41 to 1.06 for cubic and from 1.24 to 1.05 for quartic. The last observed orders are consistent with p in energy and p+1 in L². This case shows a coarse approximation penalty with recovered order under refinement; it provides no evidence of persistent locking. These are two measured doublings on the fine end of one sequence, with fixed height dependence. The curved gallery model requires its own h-study.

The first cubic n=8 C¹ solve using ordinary B-spline interpolation coordinates reached the PCG iteration limit. A complete-pivot coordinate proposal exceeded the rational reconstruction budget, and a subsequent rational coordinate exchange exceeded its denominator budget. The final implementation expresses the same section space in endpoint Hermite values and derivatives (with one additional middle Bernstein coefficient per quartic interval). Reversing a shared trace reverses its tangential derivative sign. The compact coordinate system allows the native exact signed lift to be used directly. The coarse cubic C⁰ dimension, energy, L² error, volume and strain energy agree with the preexisting elevated-element assembly. All sixteen primary states and eight q10 repeats pass. [Predeclared geometry and thresholds](pocket-vector-protocol.md) · [All native observations](pocket-vector.json) · [Independent scalar experiment and instructions](../pocket-lateral-c1-sweep.md)

The independent sine-loaded pocket uses a scalar displacement constant
through the thickness and a side clamp. It also varies internal knot
multiplicity between the two branches. Its C⁰ and C¹ DOF counts consequently
have different meanings from the native vector counts. At the common
q16/q24 replay, 13 of 16 scalar states meet the 1e−9 relative-observable
gate. Three n=8 rows miss it: cubic C¹ and quartic C⁰/C¹. Their error
measurements remain visible with their failed status. The initial
q12/q16 coarse gate also failed for quartic C⁰. These auxiliary runs
do not enter the vector convergence rates above. [Scalar quadrature retries](../pocket-lateral-c1-sweep-quadrature-retries.json)

## Scalar singularity at a reentrant edge

The section is the union of three unit squares with a 270° corner. The benchmark prescribes the harmonic field with gradient magnitude proportional to r^(−1/3) and projects its outer trace. Local cubic THB refines nested boxes at that corner. The uniform control refines the entire square grid. Both calculations use the same physical problem and the same finest span width in each row.

| Finest spans/side | THB total / free scalar DOFs | Uniform total / free scalar DOFs | Relative H¹ seminorm error, THB / uniform |
|---:|---:|---:|---:|
| 8 | 341 / 261 | 341 / 261 | 2.25750% / 2.25750% |
| 16 | 477 / 389 | 1,045 / 901 | 1.42207% / 1.42204% |
| 32 | 613 / 517 | 3,605 / 3,333 | 0.895823% / 0.895798% |
| 64 | 749 / 645 | 13,333 / 12,805 | 0.564334% / 0.564310% |

The 749-coefficient local space improves on the 3,605-coefficient uniform result of the preceding level. At the finest common resolution its energy error differs from uniform by 0.0043% relative. A supplementary six-span-wide box gives 0.564312% at 1,277 total coefficients (1,161 free). Increasing the number of free coefficients by 80% from the four-span case changes its fine-level error by less than 0.005% relative. Local refinement leaves the far boundary coarse, visible in the relative boundary L² error (1.01e−7 locally and 2.45e−11 uniformly).

Rust assembles the Laplace bilinear form from the THB cell extractions, projects the boundary in its actual trace space and reports all five quadrature orders. On the singular cells the error integration splits the square into two triangles with a cubic radial substitution. The independent implementation performs exact univariate Boehm knot insertions with Fraction arithmetic, truncates the tensor basis on the finest grid, checks trace coefficients and unity, compares the basis and physical gradients at 96 chart probes, then solves an independently assembled sparse tensor Galerkin system in SciPy. The reported H¹ errors use integrated squared gradient differences. This experiment qualifies a scalar section field, a prescribed mark schedule and one edge singularity. [Frozen protocol](../thb-reentrant-edge-protocol.md) · [Exact transfer and admissibility checks](../thb-reentrant-edge-implementation.md) · [Native histories and independent controls](reentrant-ritz.json)

## Reproduction

Run `make check-continuity-followup` from the repository root in a release-capable Rust environment with NumPy and SciPy. The extended suites write their full-precision results under `build/`. The native vector solves, the scalar section controls and the exact geometry classifier each have a different scope and can be replayed independently.
