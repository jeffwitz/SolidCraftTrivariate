# Pairwise Bézier parent recovery

## Protocol

The audit asks whether two adjacent, oriented tensor Bézier blocks are exact restrictions of one tensor Bézier parent. It processes declared whole-face interfaces in model order. Inputs must pass the model audit, have positive homogeneous weights, and use degrees supported by the bounded native arithmetic. The two traces must agree after exact Bernstein degree elevation up to one positive constant homogeneous gauge.

For each candidate pair, the first normal jet determines one positive interval-width ratio shared by the complete face. The implementation reconstructs the unsplit parent controls with exact rational arithmetic, then splits that parent with de Casteljau and compares every homogeneous coefficient of both blocks. Opposite nonzero first jets with zero sum are a mathematical refusal; zero jets or arithmetic-budget failures are indeterminate. The independent Python oracle uses `Fraction` and repeats orientation, elevation, gauge, full normal-jet checks, and both complete de Casteljau restrictions.

This establishes pairwise mathematical parent support only. It does not identify a historical CAD feature, certify the whole solid, or prove global regularity, coverage, or disjoint interiors. Only complete declared faces are considered; partial-face matching is outside this audit.

## Corpus results

`make check-inherited-corpus` runs the native report and the independent oracle for each row below. The runner creates its scratch directory before building the generated multiloop model. Reports and oracle details are regenerable under `build/pair-parent-corpus/`.

| Source | Interfaces | Exact parents recovered | Recovered on lateral axes `u`/`v` | Refused | Indeterminate |
|---|---:|---:|---:|---:|---:|
| S3 face pocket | 191 | 62 | 5 | 129 | 0 |
| NURBS multiloop | 75 | 29 | 8 | 46 | 0 |
| G2 crossed plate | 46 | 16 | 0 | 30 | 0 |
| Five FieldCraft models | 55 | 0 | 0 | 55 | 0 |
| **Total** | **367** | **107** | **13** | **260** | **0** |

The independent oracle replayed all 107 recovered parents exactly. Every native full-normal-jet candidate reconstructed both complete nets; none failed the exact parent check. All 55 retained FieldCraft interfaces were negative controls and refused recovery. The compact frozen result, including source model and report SHA-256 values, is [pair-parents.json](pair-parents.json). The runner also writes a regenerable copy at `build/pair-parent-corpus/summary.json`.

## Reproduction

```sh
cargo build --release --locked --offline
python3 scripts/check_pair_parent_corpus.py --cli target/release/solidcraft-trivariate
```

The focused Rust tests are `cargo test --locked --offline --test pair_parent`. They include unequal split intervals, positive gauge changes, swapped and reversed orientations, degree elevation, a false higher-order control, malformed faces, nonpositive weights, and the opposing-slope refusal case.
