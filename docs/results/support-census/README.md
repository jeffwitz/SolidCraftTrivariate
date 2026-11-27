# Exact support and interface census

## Scope and identity rules

The polygon corpus contains 45 historical cases, 12 supplements, and 11 held-out recipes. Recipes are deduplicated by SHA-256 over canonical JSON; every original corpus membership remains attached to the record. In this snapshot there are 68 source entries and 68 unique recipe hashes, with no duplicate groups.

For each built model, the census lists every declared topological interface with its two block sides and stored orientation. It extracts the homogeneous face control nets, elevates each side exactly with rational arithmetic, applies the orientation, and tests equality up to one positive constant gauge. This proves the shared geometric face trace represented by those nets. It does not prove that the adjacent maps came from one parent spline support, nor does it imply G1/G2 geometry continuity or continuity of an analysis field.

The independent S2 crossed-union study is recorded separately because its model format does not declare the partial faces as ordinary whole-face topology links. For each of three rise states it has seven complete chart neighbors with exact homogeneous face traces, four partial wing-to-core footprints, and two collapsed polar-line constraints. The last item is a line constraint, not a face, and is excluded from interface counts. The four partial links retain their documented preimage-aware trace treatment and its limitations; they are not relabeled as exact whole-face supports.

## Measurements

| Corpus | Recipes / unique hashes | Native interfaces | Exact homogeneous traces | pz=2 C0 scalar DOFs, summed by model | pz=2 inherited C1 scalar DOFs, summed by model | Inherited joins |
|---|---:|---:|---:|---:|---:|---:|
| Historical polygon cases | 45 / 45 | 2,074 | 2,074 | 5,092 | 4,412 | 385 |
| Supplements | 12 / 12 | 729 | 729 | 1,750 | 1,430 | 183 |
| Held-out recipes | 11 / 11 | 333 | 333 | 987 | 909 | 49 |
| **Total** | **68 / 68** | **3,136** | **3,136** | **7,829** | **6,751** | **617** |

All 68 recipes built and passed native audit. Both native field packages compiled for all 68 models with the current release CLI. The C0 reference plan uses `axial_degree=2`, `axial_continuity=0`, `restore_strip=false`, and no transverse-interface requests. The inherited candidate uses the same plan with `axial_continuity=1`. The degree-of-freedom totals are sums over separate models, not a single combined analysis space. Per recipe, the JSON records both DOF counts, inherited-join and retained-C0-join counts, and both package paths. No native compilation refusal was omitted from the denominator.

The separately measured named states are:

| State group | Model states | Declared interfaces | Exact full-face traces |
|---|---:|---:|---:|
| Tall-cylinder Boolean states | 7 | 124 | 124 |
| Box–cylinder routes | 2 | 8 | 8 |
| Compact FieldCraft variants | 5 | 55 | 55 |
| S2 crossed union at rises −3, 0, +3 mm | 3 | 33 total: 21 full + 12 partial | 21 full; 12 partial retained as partial |

The fourteen native model states and three S2 states are reported as state observations. The table deliberately does not merge repeated templates across variants into a claimed count of unique geometric interfaces.

## Reproduce

From this checkout, rebuild the release binary and run the census:

```sh
cargo build --release --locked --offline
python3 scripts/census_support_interfaces.py \
  --cli target/release/solidcraft-trivariate \
  --work build/support-census \
  --output build/support-census/census.json
python3 -m unittest tests/test_support_census.py
```

The detailed output is regenerable under `build/support-census/`: canonical recipe copies, one native model and audit per unique recipe, C0/C1 inherited-space packages, and `census.json`. The JSON records corpus manifests and hashes, model hashes, CLI hash and source revision, all 3,136 interface records, all 68 native C0/C1 observations, plus the named-state interface records. It also retains all refusals and mismatches if a future replay produces them.

The parent/support classification is intentionally narrow: an exact face trace is a representation fact. Same-parent continuity evidence, full support recovery, and partial-interface coupling remain separate claims with their own cited evidence.

The separate [FieldCraft parent-support and native H1 replay](fieldcraft-parent-support.md) audits all 55 interfaces of its five compact variants, includes a positive parent-split control, and compares the degree-harmonized C0 field with the existing mechanics space.

The [compact role index](role-index.json) records every polygon interface by ordinal in its fingerprinted source model. It also groups all exterior sides by their source feature. The native field evidence and both directional first-jet tests are cross-referenced during `make check-inherited-corpus`.
