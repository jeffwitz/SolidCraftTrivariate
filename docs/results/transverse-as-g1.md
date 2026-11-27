# Exact AS-G1 census for bilinear extruded sections

## Preregistered contract

- Starting revision: `1c63c1280e1bef16fcf9a533b6fe70d6e0a8d158`.
- Scope: `src/transverse_transitions.rs`, `tests/transverse_transitions.rs`,
  `scripts/census_transverse_transitions.py`, and this report plus the native
  coverage inventory.
- Inputs: lateral interfaces of degree-`[1,1,1]` extruded blocks with one
  shared positive homogeneous weight per block. Other degrees or projective
  weights are explicitly unsupported. Axial faces are outside the section
  test.
- Exact admission conditions: linear alpha factors; quadratic or lower total
  beta; relatively prime alpha factors; positive alpha product on the closed
  edge; nonzero section Jacobians throughout the closed edge; and an exact
  decomposition of total beta into two linear side factors. Checked rational
  arithmetic uses `i128`; overflow is a refusal, not an approximate result.
- Independent oracle: separately implemented Python `Fraction` polynomial
  arithmetic, compared record-for-record against the native CLI.
- Negative cases: unsupported family or weights, non-orientable data, singular
  edge jets, failed linear beta split, and rational overflow. Existing
  constant/affine transition outcomes remain a separate diagnostic.
- Exit criteria: include the 45 historical polygon models and 23 supplemental
  or held-out models; achieve exact Python/native report parity; pass positive
  and negative Rust tests. No tracked historical census output is regenerated.

## Constructive factors and proof

For the inward section derivatives `n₁(t)`, `n₂(t)` and common edge tangent
`t(t)`, the implementation forms the exact determinant identity

```text
det(n₂,t) n₁ + det(t,n₁) n₂ + det(n₁,n₂) t = 0.
```

It reduces the three polynomial coefficients by their common gcd and solves
exactly for linear side factors in

```text
beta = alpha_first * beta_second + alpha_second * beta_first.
```

The report contains the two alpha polynomials, both linear side beta
polynomials, total beta, common factor, coprimality and orientation results,
closed-edge regularity, and whether the split identity was verified. Polynomial
coefficients are ascending powers of the oriented edge parameter `t ∈ [0,1]`
and use reduced rational strings. The native implementation verifies the
reconstructed beta split; the Python oracle independently reconstructs it and
checks the determinant identity.

This is a per-interface algebraic certificate. It does not establish vertex
compatibility across a whole layout, a global reparameterization, or the
approximation properties of a particular spline space. The pre-existing
constant/affine transition fields retain their old meaning.

## Corpus replay

The support census has 68 models: 45 historical polygon models plus 23
supplemental or held-out models. Python/native reports matched exactly on all
68 models and 3,136 interfaces. The 2,519 lateral interfaces all receive an
AS-G1 construction; 617 axial interfaces are counted separately as
nonlateral. All 45 historical models remain included and have positive lateral
proofs. No refused or unsupported record was removed from the denominator.

Executed from the repository root:

```sh
cargo test --locked --offline --test transverse_transitions
cargo build --locked --offline
python3 scripts/census_transverse_transitions.py --cli target/debug/solidcraft-trivariate
```

The full `check-inherited-corpus` target still runs this census before its
downstream join-mechanics check. Saved and published measurements were not
edited.

## Literature

The condition follows the published AS-G1 definition: linear alpha and side
beta data, relatively prime alpha factors, positive alpha product, and the
geometric continuity identity. The source identifies piecewise bilinear
multi-patch parameterizations as one subclass of AS-G1 geometries:

M. Kapl, G. Sangalli and T. Takacs, “Isogeometric analysis with C¹
functions on planar, unstructured quadrilateral meshes,” *SMAI Journal of
Computational Mathematics* S5 (2019), 67–86, DOI
[10.5802/smai-jcm.52](https://doi.org/10.5802/smai-jcm.52), Definition 3.1 and
equations (3.3)–(3.5).
