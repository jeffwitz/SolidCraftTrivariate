# A boundary-singular IGA chart against regular P2 finite elements

Run `OPENBLAS_NUM_THREADS=1 python3 run.py` with NumPy and SciPy to regenerate results.json. This standalone reference calculation uses antiplane linear elasticity with unit shear modulus, per unit extrusion length. Both solvers use the same triangle with vertices (−1,0), (1,0), (0,1), homogeneous prescribed displacement on x+y=1, and manufactured body force and tractions on the other faces. The exact out-of-plane displacement is (1−x−y) exp(0.6x+0.4y).

The FEM mesh uses regular affine triangles with quadratic Lagrange elements. The IGA chart maps the unit square by F(s,t)=(s−t,st), with unit weights and quadratic tensor-product B-spline analysis. It is injective even on the boundary. At (0,0) its two parameter tangents are nonzero and collinear. Its determinant s+t vanishes there, giving a singular longitudinal line after extrusion. The worst stiffness density for bounded parametric derivatives grows as 1/(s+t), which has finite integral 2 log(2) over the square. The corner cell uses a Duffy transform on two triangles. The solve uses Gauss order 10 and error evaluation uses order 20. Reassembly at order 18 checks quadrature sensitivity.

| Intervals | FEM free DOFs | FEM relative energy norm error | IGA free DOFs | IGA relative energy norm error |
|---:|---:|---:|---:|---:|
| 2 | 10 | 3.499% | 12 | 2.528% |
| 4 | 36 | 0.949% | 30 | 0.626% |
| 8 | 136 | 0.244% | 90 | 0.156% |
| 16 | 528 | 0.0615% | 306 | 0.0390% |

Both methods reproduce an affine displacement to relative energy norm error below 2e−14. The manufactured solution converges at about second order in the energy norm. At 16 intervals, the spectral condition numbers of the reduced unpreconditioned stiffness matrices are about 2450 for FEM and 13432 for IGA. These metrics compare different bases and dimensions. For IGA at eight intervals, the shear at x=0, y=1e−8 is −0.66738 versus the exact −0.60000. Global energy convergence coexists with poor local stress very near the rail. This experiment establishes integrability for the specified map and basis. Each fillet chart requires its own integrability argument and numerical qualification under the intended load.
