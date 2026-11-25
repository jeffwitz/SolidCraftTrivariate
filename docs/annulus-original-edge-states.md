# Fillets and chamfers on the sharp annulus partition

The prototype reuses the four `full_annulus_v1` volumes, their identities, four interface faces, and the dyadic angular coefficients of the sharp part. The displayed variants use the same triquadratic basis, knot vector in each direction, and weights. Each volume has 63 poles (3 × 3 × 7). Additional knots are internal to the volumes. All three states retain four sectors, with three axial integration spans per sector.

| State | Gauss-16 volume (mm³) | Regularity |
|---|---:|---|
| [Sharp](https://github.com/jeffwitz/SolidCraftTrivariate/blob/38825ceaf8d70141181e48ee3d4cc06387d87deb/docs/site/data/models/annulus-four-sharp.json) | 301.592894744620 | Positive determinant on closed spans |
| [Four chamfers](https://github.com/jeffwitz/SolidCraftTrivariate/blob/38825ceaf8d70141181e48ee3d4cc06387d87deb/docs/site/data/models/annulus-four-chamfer.json) | 299.236700254428 | Positive determinant on closed spans |
| [Four fillets](https://github.com/jeffwitz/SolidCraftTrivariate/blob/38825ceaf8d70141181e48ee3d4cc06387d87deb/docs/site/data/models/annulus-four-fillet-singular.json) | 300.581607414644 | Positive in material, zero on four contact circles. IGA not qualified |

## Fillet map

In the lower span, s is the radial parameter (0 outside, 1 inside), t is the local axial parameter, and R=1/4 mm. Set a=2s−1, D=1+t², d=R(1−t)²/D, b=1−d, E=t(1−t)/D, and k=R/2. The nominal radius and height are

    ρ = 3 − a b
    z = 2R t²/D + k(1−a²)E.

On both walls a=±1, the k term vanishes. The functions (1−t²)/D and 2t/D trace an exact quarter-circle in the nominal meridian profile. Revolution uses the original dyadic angular approximation. The upper span is symmetric to the lower span about z=4 mm. In the middle span, the profile is rectangular.

The k term regularizes **the interior of the end faces**. Up to a positive constant, the meridional factor of the oriented Jacobian is b z_t + 2 a² k b′ E. Here b>0, b′≥0, E≥0, and

    D² z_t = 4R t + k(1−a²)(1−2t−t²).

For t>0 this last term is strictly positive because 1−2t−t² ≥ −2t and 4R>2k. At t=0 it is positive whenever |a|<1. Thus the only zeros are the two circles at z=0 and their symmetric counterparts at z=8, at nominal radii 2.25 and 3.75 mm. The angular factor remains strictly positive as in the original part. To verify coverage, at fixed radius solve a=(3−ρ)/b(t). The permitted t values form an interval because b is increasing, and z is strictly increasing there in the interior by the Jacobian formula. The boundary describes the rounded profile, and the four angular sectors each cover their quadrant once.

A cube edge between two tangent physical faces forces this rank loss in these four sectors. An IGA field test and qualification of this singular-map family remain to be done. The earlier report [on the strict Jacobian condition](annulus-four-block-topology.md) explains why the regular variant uses a different partition.

## Reproducible checks

The independent check operates directly on gallery models. On every span and angular row, nine exact-fraction comparisons determine componentwise identity of the homogeneous biquadratic polynomials (numerators and denominator) with the analytic formula. The four interface faces agree coefficient by coefficient. Gauss quadrature orders 12 and 16 agree with the analytic sharp volume or treated STEP reference to within 10⁻¹⁰ mm³.

```sh
cargo test --offline --example annulus_original_edges
cargo run --offline --example annulus_original_edges -- fillet > build/annulus-original-fillet.json
python3 tests/references/annulus_original_edge_states.py
```

[Rust constructor](../examples/annulus_original_edges.rs) · [Independent oracle](../tests/references/annulus_original_edge_states.py)
