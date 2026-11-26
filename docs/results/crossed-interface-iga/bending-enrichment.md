# Crossed-tool union under eccentric loading

On the tilted ten-chart union, enriching the longitudinal and thickness degrees from (3,1,1) to (3,2,2) and using four longitudinal spans raises the mean loaded-end displacement from 0.080049 to 0.301272 mm. The fine HEX20 comparator gives 0.307772 mm. A directional control on conforming geometry identifies longitudinal under-resolution in bending as the dominant observed effect.

## Protocol recorded before the enriched solves

The primary observable is external work divided by the applied x force, in mm. Keep the original ten rational maps, material (E = 1,000 MPa, ν = 0.3), clamp at x = 0 and 1 MPa right-end traction. Compare the original ten-map strong partial trace space with the conforming 22-map subdivision of the same material. Keep the two geometries and their quadrature mappings independent of field degree. Report new field knots or subintervals, trace treatment, free vector degrees of freedom, equilibrium and checks away from constraint points.

The primary case has a 3 mm rise. Zero rise is the loading control. Repeat the original degree (3,1,1) and increase v and w to (3,2,2). The conforming control also tests an increase in v or w alone at fixed mesh. Refine in v and w, including the traces. Use the common HEX20 sequence as a numerical reference and report its mesh-to-mesh change. A comparison at equal degrees of freedom requires an actual count match; give both counts for all other comparisons.

The diagnostic is the change in compliance under these interventions. If the unchanged ten-map topology recovers most of the displacement and the two enriched models approach the reference, the original approximation space accounts for much of the observed stiffness. An improvement under p/h enrichment alone cannot uniquely identify shear locking. A residual difference at comparable resolution is a reason to examine trace coupling and geometry approximation. Loaded displacement and global equilibrium alone cannot qualify local stresses.

## Prior observations

The tilted ten-map result with eight v intervals was 0.166888 mm. The coarse 22-map control was 0.083168 mm. HEX20 gave 0.306767 and 0.307772 mm on the two published meshes. [Original geometry and loading](README.md) · [HEX20 study](hexahedral-comparison.md)

A local straight-beam bending mode explains the directional test. In longitudinal and thickness coordinates x and z, the displacement components uₓ = −κxz and u_z = κx²/2 give γₓz = ∂uₓ/∂z + ∂u_z/∂x = 0. Linear longitudinal interpolation can reproduce the first component but cannot reproduce the second on a single span. Raising the longitudinal degree can therefore release spurious shear energy. The tilted assembly also contains a curved core and genuinely three-dimensional stresses, so this example is an interpretation of the directional tests rather than an analytical solution for S2.

## Which direction causes the stiffness?

On the conforming model, elevating the longitudinal degree alone recovers 81% of the coarse displacement shortfall relative to the fine HEX20 result. Elevating thickness alone recovers 0.24%. With both degrees elevated, a fourfold split in v brings the displacement much closer to HEX20 than a fourfold split in w, at similar free-DOF counts. [Directional runs and quadrature replay](conforming_bending_study.json)

| Geometry charts | Field degree (u,v,w) | Spans (v,w) | Free vector DOFs | Mean loaded displacement (mm) |
|---:|---|---|---:|---:|
| 22 | (3,1,1) | (1,1) | 474 | 0.083168 |
| 22 | (3,1,2) | (1,1) | 819 | 0.083711 |
| 22 | (3,2,1) | (1,1) | 762 | 0.265044 |
| 22 | (3,2,2) | (1,1) | 1,302 | 0.281110 |
| 22 | (3,2,2) | (1,4) | 4,542 | 0.281286 |
| 22 | (3,2,2) | (4,1) | 4,200 | 0.303480 |
| 22 | (3,2,2) | (4,2) | 7,620 | 0.303751 |

The q8 → q12 repeat on the (2,2) mesh changes mean displacement by 5.4 × 10⁻⁹ relative. The original exact stock and tool definitions, the resultant 40.222222 N and the loaded-end observable are shared across cases.

## Continuity at a fixed field dimension

The 22-chart control also permits an independent change of field continuity. Across new internal v and w knots the geometry denominator has a continuous first derivative. For the quadratic field, equality of the two derivatives of its homogeneous numerator eliminates the shared control plane by a weight-adjusted average of its neighbors. A sparse prolongation enforces all such relations, including their intersections. Chart-to-chart joins remain C⁰. Each pair uses identical rational geometry, clamps, loads, integration cells and quadrature. The residual of the full set of coefficient constraints on the prolongation is below 10⁻¹². [Solver](smooth_bending.py) · [Measurements](smooth_bending.json)

| v × w spans per chart | Internal field continuity | Free vector DOFs | Integration cells | Mean loaded displacement (mm) | Difference from fine HEX20 |
|---|---|---:|---:|---:|---:|
| 2 × 2 | C⁰ | 4,128 | 88 | 0.297294 | −3.40% |
| 2 × 2 | C¹ | 2,520 | 88 | 0.289699 | −5.87% |
| 3 × 3 | C⁰ | 8,514 | 198 | 0.301743 | −1.96% |
| 3 × 3 | C¹ | 4,128 | 198 | 0.298178 | −3.12% |

At exactly 4,128 free vector DOFs, simple internal knots allow three spans per direction. The displacement difference from fine HEX20 decreases from 3.40% to 3.12%, with 198 integration cells versus 88. At fixed spans, the C¹ space is a proper subspace of C⁰ and gives a stiffer response. These are observable-specific accuracy comparisons. Assembly and solve speed have not been benchmarked.

## The original ten charts with four partial interfaces

The original geometry retains ten charts throughout. The field is elevated to (3,2,2) and its traces are subdivided along with its interior. On each overlapping face, the wing and central basis functions are evaluated at the same physical point. A numerical null space of the overdetermined trace equations gives the admissible displacement coefficients. Additional evaluations test the trace operator and the solved displacement between constraint points. This is a numerical continuity check for the enriched spaces. [Solver](partial_bending_study.py) · [Measurements](partial_bending_study.json)

| Field | v × w spans per chart | Free vector DOFs | Mean loaded displacement (mm) |
|---|---|---:|---:|
| (3,1,1) | 1 × 1 | 222 | 0.080049 |
| (3,2,2) | 1 × 1 | 594 | 0.278081 |
| (3,2,2) | 1 × 2 | 996 | 0.278424 |
| (3,2,2) | 2 × 1 | 1,152 | 0.294343 |
| (3,2,2) | 2 × 2 | 1,926 | 0.294856 |
| (3,2,2) | 4 × 1 | 2,268 | 0.301272 |

On the tilted model, doubling the thickness resolution alone has little effect. At four longitudinal spans the original topology reaches 97.89% of the fine HEX20 displacement. The corresponding conforming 22-chart model reaches 98.61% at 4,200 free vector DOFs. The two spaces use different numbers of unknowns and different interface traces. This comparison isolates a large field-enrichment effect but cannot uniquely assign the remaining difference to strong coupling, discretization of the bulk, or chordal geometry in the HEX20 reference.

On the 2 × 1 and 2 × 2 meshes, sampling 19 or 29 points per longitudinal trace segment and varying the relative singular-value cutoff from 10⁻¹⁰ to 10⁻¹⁴ gives the same scalar constraint ranks (66 and 108). The stability check concerns the ranks for these two meshes.


At zero rise, the same ten-chart (3,2,2) field with four longitudinal spans gives 0.087211 mm against the fine HEX20 value of 0.089171 mm. The original (3,1,1) value was 0.076785 mm. Reversing the rise from +3 to −3 mm preserves the enriched axial response to within 10⁻⁸ mm on all three longitudinal meshes.

Increasing quadrature from q8 to q14 on the (4,1) tilted ten-chart space changes displacement by 7.9 × 10⁻⁹ mm. Reactions and displacement continuity qualify the global solution here; assessing stresses near the partial faces requires a separate local convergence study.

## Equal free dimensions against HEX20

The comparison uses the +3 mm tilted geometry, 40.222222 N load, identical clamp and traction-weighted mean displacement. For this additional test, the mesh pair was chosen by free dimension before solving for its displacements. The 22-chart rational model uses the same material domain as the ten-chart union. The HEX20 mesh follows its charts with straight boundary chords. [IGA result](matched_iga.json) · [HEX meshes and results](hex_bending_reference.json)

| Field | Free vector DOFs | Integration cells or elements | Volume (mm³) | Mean loaded displacement (mm) | Difference from fine HEX20 |
|---|---:|---:|---:|---:|---:|
| IGA (3,2,2), C⁰, 6 × 2 | 11,112 | 264 | 6,799.605 | 0.305584 | −0.71% |
| HEX20, chart-following | 11,112 | 522 | 6,779.590 | 0.306089 | −0.55% |
| Fine HEX20, n=2 | 125,004 | 8,184 | 6,794.582 | 0.307772 | reference |

The two equal-budget displacements differ by 0.000504 mm (0.165% of the matched HEX20 result). The matched HEX20 represents 0.294% less volume. The fine HEX20 displacement differs by 0.33% from the published n=1 HEX20 value. These comparisons establish agreement for global compliance in this case. They do not rank the two methods for stress accuracy or runtime.

## Reproduction

From the repository root with NumPy and SciPy installed, run `make check-crossed-bending`. It rebuilds the native geometry and re-solves selected directional, partial-interface and equal-dimension continuity cases. The commands below replay the larger studies. The hexahedral comparison also uses SymPy for the original ten-chart trace calculation.

```sh
cargo build --release --locked --offline
mkdir -p build
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
study=docs/results/crossed-interface-iga
binary=target/release/solidcraft-trivariate
python3 "$study/conforming_bending_study.py" "$binary" build/s2-conforming.json
python3 "$study/partial_bending_study.py" "$binary" build/s2-partial.json --plan "$study/partial_bending_plan.json"
python3 "$study/smooth_bending.py" "$binary" build/s2-smooth.json
python3 "$study/hex_bending_reference.py" "$binary" build/s2-hex-budget.json
```
