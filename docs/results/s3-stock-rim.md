# Fixed-partition treatment of the two S3 stock rims

Starting revision: 34b9720. Target the eight straight horizontal outer edges of the delivered 100 × 50 × 12 mm S3 stock, using 0.25 mm setback or nominal radius. The four straight sections of each perimeter meet with mitres. Retain all original block IDs, complete-face neighbor relations and internal feature boundaries. Reject other stock conventions or inconsistent charts. Use exact rational tensor spans grouped under their parent blocks. Test incidence against the sharp parent, quadrature volume against the rectangular section formula, and determinant signs over whole closed cells with integer Bernstein coefficients. Classify boundary rank loss explicitly. The Rust entry point is exposed through the native CLI.

## Results

All three states have the original 93 logical charts, 191 logical interfaces
with their original orientations, and 176 logical boundary sides. Each state
uses 145 tensor Bézier spans. The inserted axial stations at z = 0.25 and
11.75 mm affect the cap layers only. In each stock-wall row, a joint
x/y contraction gives the outer rectangle the same inset on all four sides.
The opposite guide row retains its xy coordinates. All material holes,
pocket walls and the floor retain their original physical geometry. At
rectangular corners the two modified walls meet in a mitre. The vertical
stock edges and the edges of the internal cuts remain untreated.

| State | Spans with positive determinant on the closed cube | Other spans | Removed material, mm³ |
|---|---:|---|---:|
| Sharp | 145 | 0 | 0 |
| Chamfer | 145 | 0 | 18.708333333333 |
| Fillet | 119 | 26 positive in the interior, singular along their cap contact edge | 8.035585080112 |

The oracle expands the homogeneous determinant in the tensor Bernstein
basis with integer coefficients. For each open face, edge, interior and
vertex, it checks that the corresponding restriction has a positive
coefficient or vanishes identically. All coefficients are nonnegative.
The wall sections satisfy the exact 0.25 mm setback line or circle equation,
including the arc endpoints. For height-dependent setback d the cross-section
area removed is 300d − 4d². Integration over both ends agrees with direct
quadrature of the emitted charts to within 1e−7 mm³ in removed volume at
orders 12 and 16. This checks the geometric fillet; closed-domain strict
Jacobian positivity for this fillet fails at the cap contact rails.

## Replay

From the repository root:

```sh
cargo build --release --locked --offline
mkdir -p build/s3-rim
for state in sharp chamfer fillet; do
  target/release/solidcraft-trivariate edit-s3-stock-rim docs/results/s3-face-pocket.model.json "$state" "build/s3-rim/$state.json"
done
python3 tests/references/s3_stock_rim.py docs/results/s3-face-pocket.model.json build/s3-rim
```

The same CLI accepts a lowercase `.step` output extension. Re-reading all
three native STEP files with FreeCAD/OCCT produced one valid closed solid
with one shell and 228 STEP faces per file. The default OCCT volume integral
returned 53,137.542055, 53,118.833696 and 53,129.506445 mm³ respectively
(sharp, chamfer, fillet). The native quadrature values are 53,137.541471,
53,118.833138 and 53,129.505886 mm³. STEP face count includes the axial
knot spans. The 93 logical trivariate maps are preserved in the JSON export.

A separate shell probe checked material and void at ±0.001 mm from each
predicted wall contact at three depths, on both ends and near two mitred
corners. All 72 probes per state agreed. Reproduce with
`SCT_STEP_DIR="$PWD/build/s3-rim" freecadcmd tests/references/s3_stock_rim_step.py`
after exporting the three STEP files. The two treated STEP files and CAD
observations are in [fixed-rims](fixed-rims/).

The Rust CLI also performs a pre-publication sign check. For a non-contact
span it uses the closed-cube spatial certificate. For a fillet wall span
it verifies the exact cap plane and transverse tangency using integer
homogeneous coordinates, then bounds every Bernstein coefficient away
from the identified edge strictly above zero. The closed-domain determinant
lower bound reported for the full fillet is zero. Sharp and chamfer give
a positive common lower bound of 0.0170570 mm³ in local span coordinates.
