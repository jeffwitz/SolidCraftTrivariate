//! Sparse exact elimination of homogeneous jet equations on scalar field coefficients.
use crate::constraint_system::Ratio;
use std::collections::BTreeMap;

pub(crate) type Row = BTreeMap<usize, Ratio>;
pub(crate) type Result<T> = std::result::Result<T, &'static str>;

pub(crate) fn add(row: &mut Row, i: usize, v: Ratio) -> Result<()> {
    let value = row
        .get(&i)
        .copied()
        .unwrap_or(Ratio::zero())
        .checked_add(v)?;
    if value == Ratio::zero() {
        row.remove(&i);
    } else {
        row.insert(i, value);
    }
    if row.len() > 4096 {
        return Err("indeterminate: jet row budget");
    }
    Ok(())
}
fn axpy(row: &mut Row, k: Ratio, other: &Row) -> Result<()> {
    for (&i, &v) in other {
        add(row, i, k.checked_mul(v)?)?;
    }
    Ok(())
}

/// All row indices here are pivot priorities (lower first). Retain the original
/// equations for a separate exact residual check after backward substitution.
pub(crate) fn nullspace(n: usize, equations: &[Row]) -> Result<(Vec<Row>, Vec<usize>)> {
    match direct_nullspace(n, equations) {
        Err(reason) if reason.contains("arithmetic budget") || reason.contains("overflow") => {
            modular_nullspace(n, equations)
        }
        result => result,
    }
}
fn direct_nullspace(n: usize, equations: &[Row]) -> Result<(Vec<Row>, Vec<usize>)> {
    let mut pivots: BTreeMap<usize, Row> = BTreeMap::new();
    for original in equations {
        let mut row = original.clone();
        while let Some((&i, &c)) = row.first_key_value() {
            if i >= n {
                return Err("constraint index outside field");
            }
            if let Some(other) = pivots.get(&i) {
                axpy(&mut row, Ratio::zero().checked_sub(c)?, other)
                    .map_err(|_| "indeterminate: jet forward elimination arithmetic budget")?;
            } else {
                let inv = Ratio::new(c.denominator(), c.numerator())?;
                for v in row.values_mut() {
                    *v = v
                        .checked_mul(inv)
                        .map_err(|_| "indeterminate: jet pivot normalization arithmetic budget")?;
                }
                pivots.insert(i, row);
                break;
            }
        }
    }
    let free: Vec<_> = (0..n).filter(|i| !pivots.contains_key(i)).collect();
    let mut lift = vec![Row::new(); n];
    for (j, &i) in free.iter().enumerate() {
        lift[i].insert(j, Ratio::one());
    }
    for (i, row) in pivots.iter().rev() {
        let mut terms = Row::new();
        for (&j, &v) in row.range(i + 1..) {
            axpy(&mut terms, Ratio::zero().checked_sub(v)?, &lift[j])?;
        }
        lift[*i] = terms;
    }
    for equation in equations {
        let mut residual = Row::new();
        for (&i, &v) in equation {
            axpy(&mut residual, v, &lift[i])?;
        }
        if !residual.is_empty() {
            return Err("nonzero exact jet residual");
        }
    }
    Ok((lift, free))
}

// A modular rank is a lower bound on the characteristic-zero rank when the
// input denominators are units. If reconstruction gives a verified rational
// kernel of complementary dimension with an identity free submatrix, both
// bounds are sharp. Modular agreement alone is never used to admit a field.
const P1: u64 = (1u64 << 61) - 1;
const P2: u64 = 2_305_843_009_213_693_921;
fn product<const P: u64>(a: u64, b: u64) -> u64 {
    (u128::from(a) * u128::from(b) % u128::from(P)) as u64
}
fn inverse<const P: u64>(mut a: u64) -> u64 {
    let mut power = P - 2;
    let mut result = 1;
    while power > 0 {
        if power & 1 == 1 {
            result = product::<P>(result, a);
        }
        a = product::<P>(a, a);
        power >>= 1;
    }
    result
}
fn residue<const P: u64>(x: Ratio) -> Result<u64> {
    let den = x.denominator().rem_euclid(i128::from(P)) as u64;
    if den == 0 {
        return Err("indeterminate: modular denominator");
    }
    Ok(product::<P>(
        x.numerator().rem_euclid(i128::from(P)) as u64,
        inverse::<P>(den),
    ))
}
type ModularRow = BTreeMap<usize, u64>;
fn mod_sub<const P: u64>(row: &mut ModularRow, scale: u64, other: &ModularRow) -> Result<()> {
    for (&j, &v) in other {
        let old = row.get(&j).copied().unwrap_or(0);
        let value = (old + P - product::<P>(scale, v)) % P;
        if value == 0 {
            row.remove(&j);
        } else {
            row.insert(j, value);
        }
    }
    if row.len() > 4096 {
        return Err("indeterminate: modular jet row budget");
    }
    Ok(())
}
fn modular_basis<const P: u64>(
    n: usize,
    equations: &[Row],
) -> Result<(Vec<ModularRow>, Vec<usize>)> {
    let mut pivots: BTreeMap<usize, ModularRow> = BTreeMap::new();
    for e in equations {
        let mut r = ModularRow::new();
        for (&i, &c) in e {
            if i >= n {
                return Err("invalid jet index");
            }
            let v = residue::<P>(c)?;
            if v != 0 {
                r.insert(i, v);
            }
        }
        while let Some((&i, &v)) = r.first_key_value() {
            if let Some(p) = pivots.get(&i) {
                mod_sub::<P>(&mut r, v, p)?;
            } else {
                let inv = inverse::<P>(v);
                for c in r.values_mut() {
                    *c = product::<P>(*c, inv);
                }
                pivots.insert(i, r);
                break;
            }
        }
    }
    let free: Vec<_> = (0..n).filter(|i| !pivots.contains_key(i)).collect();
    let mut reduced = vec![ModularRow::new(); n];
    for (j, &i) in free.iter().enumerate() {
        reduced[i].insert(j, 1);
    }
    for (&i, row) in pivots.iter().rev() {
        let mut expansion = ModularRow::new();
        for (&j, &v) in row.range(i + 1..) {
            mod_sub::<P>(&mut expansion, v, &reduced[j])?;
        }
        reduced[i] = expansion;
    }
    Ok((reduced, free))
}

fn reconstruct(a: u64, b: u64) -> Result<Ratio> {
    // CRT and Euclidean rational reconstruction. Both primes and their
    // product fit signed i128; every Euclidean coefficient is bounded by p.
    let modulus = i128::from(P1) * i128::from(P2);
    let delta = (b + P2 - a % P2) % P2;
    let x =
        i128::from(a) + i128::from(P1) * i128::from(product::<P2>(delta, inverse::<P2>(P1 % P2)));
    let bound = (modulus / 2).isqrt();
    let (mut m, mut r, mut t, mut u) = (modulus, x, 0i128, 1i128);
    while r > bound {
        let q = m / r;
        (m, r) = (r, m - q * r);
        (t, u) = (u, t - q * u);
    }
    if u == 0 || u.abs() > bound {
        return Err("indeterminate: two-prime rational reconstruction bound");
    }
    let result = Ratio::new(r, u)?;
    if residue::<P1>(result)? != a || residue::<P2>(result)? != b {
        return Err("indeterminate: rational reconstruction residues");
    }
    Ok(result)
}
fn modular_nullspace(n: usize, equations: &[Row]) -> Result<(Vec<Row>, Vec<usize>)> {
    let (a, free) = modular_basis::<P1>(n, equations)?;
    let (b, other) = modular_basis::<P2>(n, equations)?;
    if free != other {
        return Err("indeterminate: incompatible modular pivot profiles");
    }
    let (c, fc) = modular_basis::<2_305_843_009_213_693_907>(n, equations)?;
    let (d, fd) = modular_basis::<2_305_843_009_213_693_723>(n, equations)?;
    if fc != free || fd != free {
        return Err("indeterminate: incompatible four-prime pivot profiles");
    }
    let mut cache = BTreeMap::new();
    let mut lift = Vec::with_capacity(n);
    for (row, (a, b)) in a.iter().zip(&b).enumerate() {
        let mut values = Row::new();
        // Columns vanishing at both primes are checked by the exact kernel
        // identity below. Both modular pivot profiles give the rank lower bound.
        let keys: std::collections::BTreeSet<_> = a
            .keys()
            .chain(b.keys())
            .chain(c[row].keys())
            .chain(d[row].keys())
            .copied()
            .collect();
        for i in keys {
            let x = *a.get(&i).unwrap_or(&0);
            let y = *b.get(&i).unwrap_or(&0);
            let residues = [
                x,
                y,
                *c[row].get(&i).unwrap_or(&0),
                *d[row].get(&i).unwrap_or(&0),
            ];
            let value = if let Some(&v) = cache.get(&residues) {
                v
            } else {
                let small = reconstruct(x, y).ok().filter(|&r| {
                    residue::<2_305_843_009_213_693_907>(r) == Ok(residues[2])
                        && residue::<2_305_843_009_213_693_723>(r) == Ok(residues[3])
                });
                let r = match small {
                    Some(v) => v,
                    None => crate::spline_reconstruction::recover(residues)?,
                };
                cache.insert(residues, r);
                r
            };
            values.insert(i, value);
        }
        lift.push(values);
    }
    for equation in equations {
        let mut residual = Row::new();
        for (&i, &v) in equation {
            axpy(&mut residual, v, &lift[i])?;
        }
        if !residual.is_empty() {
            return Err("indeterminate: reconstructed jet kernel has nonzero exact residual");
        }
    }
    Ok((lift, free))
}

/// Choose an interpolatory coordinate system for coupled transverse jets. Floating
/// arithmetic proposes a column permutation only; rank, identities and acceptance
/// still come entirely from exact elimination of the original equations.
pub(crate) fn conditioned_kernel(n: usize, equations: &[Row]) -> Result<(Vec<Row>, Vec<usize>)> {
    let mut pivots: Vec<(usize, Vec<(usize, f64)>)> = Vec::new();
    let mut chosen = vec![false; n];
    for row in equations {
        let mut r = vec![0.; n];
        let mut norm = 0_f64;
        for (&i, &v) in row {
            r[i] = v.numerator() as f64 / v.denominator() as f64;
            norm = norm.max(r[i].abs());
        }
        if norm == 0. {
            continue;
        }
        for (i, other) in &pivots {
            let v = r[*i];
            if v != 0. {
                for &(j, x) in other {
                    r[j] -= v * x;
                }
                r[*i] = 0.;
            }
        }
        let Some((i, &v)) = r
            .iter()
            .enumerate()
            .filter(|(i, _)| !chosen[*i])
            .max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))
        else {
            continue;
        };
        if !v.is_finite() || v.abs() < 1e-10 * norm {
            continue;
        }
        chosen[i] = true;
        pivots.push((
            i,
            r.iter()
                .enumerate()
                .filter_map(|(j, &x)| if x != 0. { Some((j, x / v)) } else { None })
                .collect(),
        ));
    }
    let order: Vec<_> = pivots
        .iter()
        .map(|(i, _)| *i)
        .chain((0..n).filter(|&i| !chosen[i]))
        .collect();
    let mut ids = vec![0; n];
    for (j, &i) in order.iter().enumerate() {
        ids[i] = j;
    }
    let transformed = equations
        .iter()
        .map(|r| r.iter().map(|(&i, &v)| (ids[i], v)).collect())
        .collect::<Vec<_>>();
    let (lift, free) = match nullspace(n, &transformed) {
        Ok(result) => result,
        Err(_) => {
            if let Ok(basis) = nullspace(n, equations) {
                return Ok(basis);
            }
            // An interpolatory basis can have large denominators even when
            // the same space has a small basis in other free coordinates.
            // Alternate reproducible pivot orders, retaining exact checks.
            for seed in 1..=8 {
                let mut perm: Vec<_> = (0..n).collect();
                let mut state = seed as u64;
                for i in (1..n).rev() {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    perm.swap(i, ((state >> 32) as usize) % (i + 1));
                }
                let mut index = vec![0; n];
                for (j, &i) in perm.iter().enumerate() {
                    index[i] = j;
                }
                let matrix = equations
                    .iter()
                    .map(|r| r.iter().map(|(&i, &v)| (index[i], v)).collect())
                    .collect::<Vec<_>>();
                if let Ok((basis, roots)) = nullspace(n, &matrix) {
                    return Ok((
                        index.iter().map(|&i| basis[i].clone()).collect(),
                        roots.into_iter().map(|j| perm[j]).collect(),
                    ));
                }
            }
            return Err("indeterminate: exact jet kernel reconstruction budget");
        }
    };
    let free = free.into_iter().map(|j| order[j]).collect();
    Ok((ids.into_iter().map(|i| lift[i].clone()).collect(), free))
}

/// Eliminate inherited axial conditions first. Each transverse row is then
/// restricted exactly to those roots. Prefer the short section equations
/// over their redundant axial interpolants to retain a well-scaled basis.
pub(crate) fn joint_kernel(
    n: usize,
    axial: &[Row],
    transverse: &[Row],
) -> Result<(Vec<Row>, Vec<usize>)> {
    let (first, roots) = nullspace(n, axial)?;
    if transverse.is_empty() {
        return Ok((first, roots));
    }
    if roots.len() > 12000 {
        return Err("indeterminate: transverse elimination budget");
    }
    let restrict = |row: &Row, lift: &[Row]| -> Result<Row> {
        let mut output = Row::new();
        for (&j, &v) in row {
            axpy(&mut output, v, &lift[j])?;
        }
        Ok(output)
    };
    let mut projected = transverse
        .iter()
        .map(|r| restrict(r, &first))
        .collect::<Result<Vec<_>>>()?;
    projected.sort_by_key(|r| {
        (
            r.len(),
            r.values().map(|v| v.denominator()).max().unwrap_or(1),
        )
    });
    let (second, independent) = conditioned_kernel(roots.len(), &projected)?;
    let lift = first
        .iter()
        .map(|r| restrict(r, &second))
        .collect::<Result<Vec<_>>>()?;
    for r in axial.iter().chain(transverse) {
        if !restrict(r, &lift)?.is_empty() {
            return Err("nonzero composed jet residual");
        }
    }
    Ok((lift, independent.into_iter().map(|i| roots[i]).collect()))
}
