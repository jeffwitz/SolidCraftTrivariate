# Exact planar Bézier regularity certificate

## Contract

`certify_planar(degrees, controls)` accepts two nonnegative integer degrees and
a rectangular tensor grid of homogeneous triples `(wx, wy, weight)`. The
control list is flattened with the second parameter index varying fastest.
Each scalar must be exactly an `int`, finite `float`, or `Fraction`; floats are
interpreted as their exact binary rational values. Boolean values and implicit
coercions are rejected. Every control weight must be strictly positive.

The result is one of:

- `valid`: an exact positive lower bound for the physical determinant on the
  entire closed unit square;
- `invalid`: malformed input, a nonpositive weight, a zero parametric degree,
  or a determinant that is nonpositive at a parameter-square corner;
- `indeterminate`: the four corner values are positive but the global
  Bernstein coefficient lower bound is not strictly positive.

The certificate establishes local regularity and positive orientation only.
It says nothing about global injectivity, self-overlap, boundary matching, or
coverage of a target domain.

## Exact polynomial construction

Let `H=(A,B,W)` be the tensor-product homogeneous polynomial. The physical map
is `(A/W,B/W)`. Its oriented Jacobian determinant is

\[
  \det D(A/W,B/W)
  = \frac{\det(H_u,H_v,H)}{W^3}.
\]

All input values are converted to `Fraction` before arithmetic. Derivative
Bernstein coefficients use `d*(c[i+1]-c[i])`. Products use the exact
Bernstein product identity in each parameter:

\[
  (fg)_k=\sum_{i+j=k} f_i g_j
  \frac{\binom{m}{i}\binom{n}{j}}{\binom{m+n}{k}}.
\]

For input bidegree `(du,dv)`, the numerator is represented exactly at
bidegree `(3du-1,3dv-1)` when both degrees are positive. No evaluation grid,
floating tolerance, or sampled determinant is used. The Bernstein convex-hull
property gives `det(H_u,H_v,H) >= min(coefficients)` on the closed square.
Likewise `W <= max(weight controls)`; positivity of all control weights gives
`W>0`. Therefore a positive numerator coefficient minimum divided by the cube
of the maximum weight is a conservative lower bound for the physical
determinant.

The four corner values equal the corresponding corner Bernstein coefficients.
Any nonpositive corner is an exact counterexample to strict positivity and is
reported `invalid`. If the corners are positive but another coefficient is
nonpositive, this coefficient envelope is inconclusive and the result is
`indeterminate`. A future exact subdivision may tighten the envelope; it must
cover the full closed square and preserve rational arithmetic. No threshold is
adjusted to force acceptance.

## E1 usage

For each of the four pierced-plate sectors, take the `k=0` controls from its
degree-(2,1,1) block and retain `(wx,wy,weight)` in `(i,j)` order. Certify all
four patches independently. A positive result on each patch proves local
regularity for each planar sector. It does not certify the full extruded
blocks, their interfaces, or the plate's global topology.
