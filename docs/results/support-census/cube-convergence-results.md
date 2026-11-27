# Uniform cubic cube convergence: measured result

## Reproduction

Protocol: [`cube-convergence-protocol.md`](cube-convergence-protocol.md).
Environment: Rust/Cargo 1.99.1, offline release build. Reproduce the recorded
replay with:

```sh
cargo test --release --locked --offline --test cube_convergence cube_convergence_replay -- --ignored --nocapture
```

The complete replay, including every observation, quadrature defect, analytic
diagnostic and covariance, is regenerated at
`build/cube-convergence/replay.json`. The original strict gate remains a
separate ignored test and deliberately fails on the preregistered C2 rate and
matched-error-ratio criteria; its measurements are not replaced by the
supplement below.

## Primary q8 observations

Errors are relative energy norm and displacement L2. The same one-cell
solution serves as both spaces because both have the same tensor-cubic basis on a single span.

| Space | n | Scalar DOFs | Energy error | L2 error |
|---|---:|---:|---:|---:|
| C2 maximum | 1 | 64 | 1.15255397e-1 | 4.66439651e-2 |
| C2 maximum | 2 | 125 | 1.64674732e-2 | 3.41120343e-3 |
| C2 maximum | 4 | 343 | 3.14514602e-3 | 4.44045812e-4 |
| C0 Bernstein | 1 | 64 | 1.15255397e-1 | 4.66439651e-2 |
| C0 Bernstein | 2 | 343 | 1.16181513e-2 | 2.17706260e-3 |
| C0 Bernstein | 4 | 2197 | 1.50075674e-3 | 1.33169299e-4 |

| Space | Energy order, n2 to n4 | L2 order, n2 to n4 | Rate verdict |
|---|---:|---:|---|
| C2 maximum | 2.38842030 | 2.94150036 | **fail** |
| C0 Bernstein | 2.95261847 | 4.03104946 | pass |

The C2 finest-level energy error is 2.09570674 times the C0 error, above the
preregistered 1.5 limit. The strict gate therefore fails even though both
errors decrease monotonically, the C2 space uses fewer DOFs, and all geometry,
affine patch, residual, virtual-work and volume checks pass. Across all six
rows the volume error is at most `1.97e-13 mm^3`, relative residual at most
`9.81e-12`, and virtual-work defect at most `2.93e-13`.

At n=4 the q8/q12 relative defects are `2.44e-10` and `2.83e-14` for C2 L2
and energy, and `4.37e-10` and `3.32e-14` for C0 L2 and energy. Volume agrees
within `1.59e-13 mm^3`.

## Frozen n4 to n8 supplement

The supplement uses the same q8 problem and the unchanged preregistered rate
bands. Both spaces pass these rate bands:

| Space | n4 DOFs | n8 DOFs | Energy order | L2 order | Verdict |
|---|---:|---:|---:|---:|---|
| C2 maximum | 343 | 1331 | 3.12254774 | 4.25990659 | pass |
| C0 Bernstein | 2197 | 15625 | 2.98520318 | 4.03158704 | pass |

This asymptotic result does not change the failed primary gate. It suggests
the C2 rate failure at n2 to n4 was pre-asymptotic for this sine field, but it
does not meet the separate n4 C2/C0 error-ratio requirement.

## Analytic diagnostics and scope

For the exact axial sine field, the integrated body-force resultant is
`7*pi = 21.9911485751 N`, the top traction resultant is
`-3.5*pi = -10.9955742876 N`, the bottom support reaction is the same
negative value, and exact strain energy is `7*pi^2/8 = 8.63590385095 N mm`.
These are analytic equilibrium/energy diagnostics, not fitted solver targets.
The sample covariance of `log(DOF)` and `log(error)` for `(n=1,2,4)` is:

| Space | Energy covariance | Displacement L2 covariance |
|---|---:|---:|
| C2 maximum | -1.50328110 | -1.93715874 |
| C0 Bernstein | -3.83388204 | -5.17500185 |

This measurement qualifies the affine unit cube, the specified isotropic
material, sine field, loading, and uniform cubic spaces only. It does not
qualify curved geometry, stress accuracy, or arbitrary boundary conditions.
