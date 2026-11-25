# M12 directional surface quadrature oracle

## Frozen scope

This oracle qualifies the optional `surface_quadrature.split = "directional"`
mode on the v5 mechanical-job CLI. The production splitter estimates error in
the face coordinates C, U, V, and UV, then bisects the dominant axis when its
estimate is more than twice the other axis estimate; otherwise it bisects both.
The estimates remain diagnostics and are not mathematical error bounds.

Inputs are generated from the frozen G20 circular-bore corpus, plus two
single-axis shear variants `[32,0]` and `[0,32]` of the extreme-shear recipe.
The cases cover a straight bore, an oblique bore, x-only/y-only/mixed shear,
the bore wall, a selected end plane, and a rectangular exterior face. Loads
include constant world traction, pressure, and a mixed off-diagonal stress
tensor. The bore chart tests arc/extrusion directions; the end plane tests
arc/radial directions. The rectangular face supplies an affine chart control.
Existing M11 materials, supports, and mechanics are retained. The only
quadrature change is the explicit directional split mode.

## Numerical protocol

- Every case is run twice from a fresh CLI process and a fresh output path:
  first with `max_panels = 4096`, then with `max_panels = 65536`.
- Both attempts use `relative_tolerance = 1e-9`, `max_depth = 8`, and the
  preregistered load integration `order = 4`. A refusal at 4096 remains a
  recorded refusal; the 65536 run is a separate attempt, not a replacement.
- Each output starts with a sentinel before invocation. Every failed attempt
  must preserve its sentinel byte-for-byte. Diagnostics and process wall time
  are archived, including budget refusals.
- The independent load oracle evaluates serialized homogeneous Bernstein
  control nets and rational surface Jacobians in Python/NumPy. Orders 96 and
  128 are compared before order 128 is used as the reference. It reconstructs
  the full nodal vector, area, resultant, and moment on all selected physical
  faces. Surface work is checked against the reported traction work using the
  solved displacement.
- Acceptance requires the q96/q128 diagnostic and every CLI-vs-reference
  relative discrepancy in area, full nodal vector, resultant, moment, and work
  to be at most `1e-8`. For a near-zero resultant or moment, the scale is the
  integrated absolute physical force or moment, respectively.
- A case is qualified only if the job is solved, all reference comparisons
  pass, the solver's true residual is at most `1e-10`, and its force and moment
  balance outputs are finite. A process exit code alone is not qualification.
- For overall corpus qualification, the 65536 attempt for every case must
  qualify. A refusal at 4096 is accepted as a verified budget outcome only
  when its diagnostic is explicit and the sentinel remains unchanged; it is
  still counted and shown in the report.

## Reporting

`tests/references/directional_surface_load.py` records each selected face,
load law, budget, panel count, deepest split, per-face `[du,dv]` accepted grid,
maximum estimate, reference errors, refusal diagnostic, and output-preservation result. Both local chart axes must be represented by at least one accepted split in the corpus. The result report
is `docs/results/directional-surface-load-20261125.json`. The corpus includes
both attempted budgets in its denominator; a budget refusal is not hidden by a
later successful run. The reported convergence comparison is empirical and
does not extend the certificate's mathematical domain.
