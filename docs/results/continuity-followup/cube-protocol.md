# Matched-dimension cube study — prospective protocol

Base: public main fbe8723. Keep the original sine study and its failed
equal-cell gate unchanged. Freeze the following before solving this load.

Use the existing three-dimensional small-strain elasticity assembly on the
affine unit cube, lambda = 1.5 MPa and mu = 1 MPa. Clamp the whole bottom,
apply the analytic body force and exact natural tractions for displacement
(0, 0, z exp(-0.5 sum_i ((x_i-c_i)/l_i)^2)) mm, with
c=(0.31,0.57,0.68) mm and l=(0.30,0.41,0.36) mm. This profile varies in all
three directions. Use the existing local_gaussian_mms_v1 implementation.
Check analytic stress divergence by independent finite differences.

Compare uniform cubic C0 subdivisions n=4,5 with C2 subdivisions n=10,13.
Each pair has exactly the same one-dimensional basis count b=3 n_C0+1
=n_C2+3, the same b^3 scalar coefficients and 3 b² (b-1) free vector
coefficients. Keep integration cell count, CPU time and accuracy separate.
Use q=6 per direction; repeat the first pair at q=8. Solve and integrate
in 3D. Preserve all failures. Reject invalid b, quadrature or levels
before reserving the grid.

Prospective positive result: C2 relative energy-norm error smaller than C0
by at least 1.5 times at **both** dimensions; both energy errors decrease
with increasing dimension. Report L2 errors without a comparative threshold.
Require integrated volume within 2e-11 of one, residual <1e-9, relative
virtual-work defect <1e-8, q6/q8 relative change <1e-4 in each error.
Check counts, positivity of the sampled geometry Jacobian, and affine
residuals <2e-11. These are accuracy-per-coefficient checks for this load;
they supply no universal convergence rate, efficiency or time-to-solution claim.

Reproduce with the release cube integration test named matched_dimension_gaussian.
