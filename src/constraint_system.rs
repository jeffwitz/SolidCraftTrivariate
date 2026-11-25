//! Exact global reduction of rational scalar trace equations.

use crate::face_constraints::Row;
use std::collections::{BTreeMap, VecDeque};

/// A canonical rational number with checked `i128` arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Ratio {
    numerator: i128,
    denominator: i128,
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn from_magnitudes(
    negative: bool,
    numerator: u128,
    denominator: u128,
) -> Result<Ratio, &'static str> {
    if denominator == 0 {
        return Err("zero rational denominator");
    }
    let divisor = gcd_u128(numerator, denominator);
    let n = numerator / divisor;
    let d = denominator / divisor;
    if d > i128::MAX as u128 {
        return Err("rational denominator exceeds i128");
    }
    let signed_n = if negative {
        if n == (1u128 << 127) {
            i128::MIN
        } else if n <= i128::MAX as u128 {
            -(n as i128)
        } else {
            return Err("rational numerator exceeds i128");
        }
    } else {
        i128::try_from(n).map_err(|_| "rational numerator exceeds i128")?
    };
    Ok(Ratio {
        numerator: signed_n,
        denominator: d as i128,
    })
}

impl Ratio {
    /// Construct a reduced rational with positive denominator.
    pub fn new(numerator: i128, denominator: i128) -> Result<Self, &'static str> {
        if denominator == 0 {
            return Err("zero rational denominator");
        }
        let negative = (numerator < 0) ^ (denominator < 0);
        from_magnitudes(
            negative,
            numerator.unsigned_abs(),
            denominator.unsigned_abs(),
        )
    }

    /// Return the reduced numerator.
    pub fn numerator(self) -> i128 {
        self.numerator
    }

    /// Return the positive reduced denominator.
    pub fn denominator(self) -> i128 {
        self.denominator
    }

    /// Return rational zero.
    pub const fn zero() -> Self {
        Self {
            numerator: 0,
            denominator: 1,
        }
    }

    /// Return rational one.
    pub const fn one() -> Self {
        Self {
            numerator: 1,
            denominator: 1,
        }
    }

    fn is_zero(self) -> bool {
        self.numerator == 0
    }

    fn checked_neg(self) -> Result<Self, &'static str> {
        Self::new(
            self.numerator
                .checked_neg()
                .ok_or("rational negation overflow")?,
            self.denominator,
        )
    }

    pub(crate) fn checked_add(self, other: Self) -> Result<Self, &'static str> {
        let common = gcd_u128(self.denominator as u128, other.denominator as u128) as i128;
        let left_factor = other.denominator / common;
        let right_factor = self.denominator / common;
        let left = self
            .numerator
            .checked_mul(left_factor)
            .ok_or("rational addition overflow")?;
        let right = other
            .numerator
            .checked_mul(right_factor)
            .ok_or("rational addition overflow")?;
        let numerator = left
            .checked_add(right)
            .ok_or("rational addition overflow")?;
        let reduce_by = gcd_u128(numerator.unsigned_abs(), common as u128) as i128;
        let denominator = right_factor
            .checked_mul(other.denominator / reduce_by)
            .ok_or("rational denominator overflow")?;
        Self::new(numerator / reduce_by, denominator)
    }

    fn checked_sub(self, other: Self) -> Result<Self, &'static str> {
        self.checked_add(other.checked_neg()?)
    }

    pub(crate) fn checked_mul(self, other: Self) -> Result<Self, &'static str> {
        let left_cancel = gcd_u128(self.numerator.unsigned_abs(), other.denominator as u128);
        let right_cancel = gcd_u128(other.numerator.unsigned_abs(), self.denominator as u128);
        let a = self.numerator / left_cancel as i128;
        let d = other.denominator / left_cancel as i128;
        let c = other.numerator / right_cancel as i128;
        let b = self.denominator / right_cancel as i128;
        let numerator = a.checked_mul(c).ok_or("rational multiplication overflow")?;
        let denominator = b.checked_mul(d).ok_or("rational denominator overflow")?;
        Self::new(numerator, denominator)
    }

    fn checked_div(self, other: Self) -> Result<Self, &'static str> {
        if other.numerator == 0 {
            return Err("rational division by zero");
        }
        let numerator_cancel = gcd_u128(
            self.numerator.unsigned_abs(),
            other.numerator.unsigned_abs(),
        );
        let denominator_cancel = gcd_u128(self.denominator as u128, other.denominator as u128);
        let a = self.numerator.unsigned_abs() / numerator_cancel;
        let c = other.numerator.unsigned_abs() / numerator_cancel;
        let b = self.denominator as u128 / denominator_cancel;
        let d = other.denominator as u128 / denominator_cancel;
        let numerator = a.checked_mul(d).ok_or("rational division overflow")?;
        let denominator = b.checked_mul(c).ok_or("rational denominator overflow")?;
        let negative = (self.numerator < 0) ^ (other.numerator < 0);
        from_magnitudes(negative, numerator, denominator)
    }
}

/// One complete equation row mapping a scalar DOF to reduced free roots.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LiftLine {
    /// Original scalar DOF represented by this line.
    pub dof: usize,
    /// `(reduced free-root index, exact coefficient)`, sorted by reduced index.
    pub terms: Vec<(usize, Ratio)>,
    /// Exact affine offset introduced by prescribed values.
    pub offset: Ratio,
}

/// Exact affine lift from reduced free roots to all original scalar DOFs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Lift {
    /// Original free root DOF indices in ascending order.
    pub free_roots: Vec<usize>,
    /// One full representation row per original DOF, ordered by DOF index.
    pub lines: Vec<LiftLine>,
}

type Expansion = BTreeMap<usize, Ratio>;

fn insert_sum(map: &mut Expansion, index: usize, value: Ratio) -> Result<(), &'static str> {
    if value.is_zero() {
        return Ok(());
    }
    let sum = map
        .get(&index)
        .copied()
        .unwrap_or_else(Ratio::zero)
        .checked_add(value)?;
    if sum.is_zero() {
        map.remove(&index);
    } else {
        map.insert(index, sum);
    }
    Ok(())
}

fn add_scaled(
    target: &mut Expansion,
    source: &Expansion,
    scale: Ratio,
) -> Result<(), &'static str> {
    for (&index, &coefficient) in source {
        insert_sum(target, index, coefficient.checked_mul(scale)?)?;
    }
    Ok(())
}

fn expand_dof(
    dof: usize,
    equations: &[Vec<Row>],
    memo: &[Option<Expansion>],
) -> Result<Option<Expansion>, &'static str> {
    let mut canonical: Option<Expansion> = None;
    for row in &equations[dof] {
        if row.denominator <= 0 || row.terms.is_empty() {
            return Err("invalid scalar constraint row");
        }
        let mut sum = 0i128;
        let mut self_coefficient = Ratio::zero();
        let mut rhs = Expansion::new();
        for &(master, numerator) in &row.terms {
            if master >= equations.len() || numerator <= 0 {
                return Err("invalid scalar constraint term");
            }
            sum = sum
                .checked_add(numerator)
                .ok_or("constraint row denominator overflow")?;
            let coefficient = Ratio::new(numerator, row.denominator)?;
            if master == dof {
                self_coefficient = self_coefficient.checked_add(coefficient)?;
            } else {
                let master_expansion = memo[master].as_ref().ok_or("unresolved master DOF")?;
                add_scaled(&mut rhs, master_expansion, coefficient)?;
            }
        }
        if sum != row.denominator {
            return Err("constraint row weights do not sum to denominator");
        }
        let remaining = Ratio::one().checked_sub(self_coefficient)?;
        if remaining.is_zero() {
            if rhs.is_empty() {
                // Exact identity row. It contributes no dependency.
                continue;
            }
            return Err("constraint row imposes a relation between free roots");
        }
        for coefficient in rhs.values_mut() {
            *coefficient = coefficient.checked_div(remaining)?;
        }
        rhs.retain(|_, coefficient| !coefficient.is_zero());
        if let Some(previous) = &canonical {
            if previous != &rhs {
                return Err("conflicting duplicate scalar constraints");
            }
        } else {
            canonical = Some(rhs);
        }
    }

    Ok(canonical)
}

fn substitute_prescribed(
    expansion: &Expansion,
    root_values: &BTreeMap<usize, Ratio>,
) -> Result<(Expansion, Ratio), &'static str> {
    let mut free = Expansion::new();
    let mut offset = Ratio::zero();
    for (&root, &coefficient) in expansion {
        if let Some(&value) = root_values.get(&root) {
            offset = offset.checked_add(coefficient.checked_mul(value)?)?;
        } else {
            insert_sum(&mut free, root, coefficient)?;
        }
    }
    Ok((free, offset))
}

/// Resolve exact rational trace equations and admissible essential values.
///
/// Every equation is `dependent = sum(weight_i/denominator * master_i)`.
/// Duplicate equations are compared only after expanding all
/// masters to a canonical free-root expression. Cycles and constraints that
/// would require elimination among free roots are rejected.
pub fn resolve(
    n: usize,
    equations: &[(usize, Row)],
    prescribed: &[(usize, Ratio)],
) -> Result<Lift, &'static str> {
    let mut by_dof = vec![Vec::<Row>::new(); n];
    for (dependent, row) in equations {
        if *dependent >= n {
            return Err("dependent DOF index is out of range");
        }
        by_dof[*dependent].push(row.clone());
    }

    // Process masters before dependents without consuming a call frame for
    // each level of a long chain. Repeated terms are counted symmetrically.
    let mut dependents = vec![Vec::new(); n];
    let mut pending = vec![0usize; n];
    for (dof, rows) in by_dof.iter().enumerate() {
        for row in rows {
            for &(master, _) in &row.terms {
                if master >= n {
                    return Err("invalid scalar constraint term");
                }
                if master != dof {
                    dependents[master].push(dof);
                    pending[dof] += 1;
                }
            }
        }
    }
    let mut ready: VecDeque<_> = (0..n).filter(|&dof| pending[dof] == 0).collect();
    let mut memo = vec![None; n];
    let mut roots = vec![false; n];
    let mut processed = 0;
    while let Some(dof) = ready.pop_front() {
        let expansion = expand_dof(dof, &by_dof, &memo)?;
        roots[dof] = expansion.is_none();
        memo[dof] = Some(expansion.unwrap_or_else(|| BTreeMap::from([(dof, Ratio::one())])));
        processed += 1;
        for &dependent in &dependents[dof] {
            pending[dependent] -= 1;
            if pending[dependent] == 0 {
                ready.push_back(dependent);
            }
        }
    }
    if processed != n {
        return Err("cyclic scalar constraint graph");
    }

    let mut requested = BTreeMap::<usize, Ratio>::new();
    for &(dof, value) in prescribed {
        if dof >= n {
            return Err("prescribed DOF index is out of range");
        }
        if let Some(previous) = requested.insert(dof, value) {
            if previous != value {
                return Err("conflicting prescribed values");
            }
        }
    }

    let mut root_values = BTreeMap::<usize, Ratio>::new();
    let mut dependent_values = Vec::new();
    for (dof, value) in requested {
        if roots[dof] {
            root_values.insert(dof, value);
        } else {
            dependent_values.push((dof, value));
        }
    }
    for (dof, value) in dependent_values {
        let (free, offset) =
            substitute_prescribed(memo[dof].as_ref().expect("all DOFs expanded"), &root_values)?;
        if !free.is_empty() {
            return Err("prescribed dependent DOF requires elimination or basis adaptation");
        }
        if offset != value {
            return Err("prescribed dependent DOF conflicts with trace equations");
        }
    }

    let free_roots: Vec<usize> = (0..n)
        .filter(|dof| roots[*dof] && !root_values.contains_key(dof))
        .collect();
    let reduced_index: BTreeMap<usize, usize> = free_roots
        .iter()
        .enumerate()
        .map(|(reduced, &root)| (root, reduced))
        .collect();

    let mut lines = Vec::with_capacity(n);
    for (dof, entry) in memo.iter().enumerate() {
        let (expanded_free, offset) =
            substitute_prescribed(entry.as_ref().expect("all DOFs expanded"), &root_values)?;
        let mut terms = Vec::with_capacity(expanded_free.len());
        for (root, coefficient) in expanded_free {
            let reduced = *reduced_index
                .get(&root)
                .ok_or("unresolved nonfree root in lift")?;
            terms.push((reduced, coefficient));
        }
        if let Some(&value) = root_values.get(&dof) {
            if !terms.is_empty() || offset != value {
                // A root expansion is its own unit coefficient and is
                // substituted above, so this is a consistency guard.
                return Err("prescribed root lift is inconsistent");
            }
        }
        lines.push(LiftLine { dof, terms, offset });
    }
    Ok(Lift { free_roots, lines })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ratio(n: i128, d: i128) -> Ratio {
        Ratio::new(n, d).unwrap()
    }

    fn row(terms: &[(usize, i128)], denominator: i128) -> Row {
        Row {
            terms: terms.to_vec(),
            denominator,
        }
    }

    #[test]
    fn chains_are_expanded_and_equation_order_is_irrelevant() {
        let equations = vec![
            (1, row(&[(0, 1), (2, 1)], 2)),
            (3, row(&[(1, 1), (2, 1)], 2)),
        ];
        let lift = resolve(4, &equations, &[]).unwrap();
        let reverse: Vec<_> = equations.iter().cloned().rev().collect();
        assert_eq!(lift, resolve(4, &reverse, &[]).unwrap());
        assert_eq!(lift.free_roots, vec![0, 2]);
        assert_eq!(
            lift.lines[1].terms,
            vec![(0, ratio(1, 2)), (1, ratio(1, 2))]
        );
        assert_eq!(
            lift.lines[3].terms,
            vec![(0, ratio(1, 4)), (1, ratio(3, 4))]
        );
    }

    #[test]
    fn duplicate_valence_constraints_compare_after_exact_expansion() {
        let equations = [
            (2, row(&[(0, 1), (1, 1)], 2)),
            (2, row(&[(0, 1), (0, 1), (1, 2)], 4)),
        ];
        let lift = resolve(3, &equations, &[]).unwrap();
        assert_eq!(lift.free_roots, vec![0, 1]);
        assert_eq!(
            lift.lines[2].terms,
            vec![(0, ratio(1, 2)), (1, ratio(1, 2))]
        );
        assert!(resolve(3, &[(2, row(&[(0, 1)], 1)), (2, row(&[(1, 1)], 1))], &[]).is_err());
    }

    #[test]
    fn identities_are_removed_and_cycles_are_rejected() {
        let identity = resolve(2, &[(1, row(&[(1, 2)], 2))], &[]).unwrap();
        assert_eq!(identity.free_roots, vec![0, 1]);
        assert_eq!(identity.lines[1].terms, vec![(1, Ratio::one())]);
        assert!(resolve(2, &[(0, row(&[(1, 1)], 1)), (1, row(&[(0, 1)], 1))], &[],).is_err());
    }

    #[test]
    fn root_prescriptions_create_offsets_and_dependent_values_are_checked() {
        let equations = [(2, row(&[(0, 1), (1, 1)], 2))];
        let lift = resolve(3, &equations, &[(0, ratio(2, 1))]).unwrap();
        assert_eq!(lift.free_roots, vec![1]);
        assert_eq!(lift.lines[0].offset, ratio(2, 1));
        assert_eq!(lift.lines[2].offset, ratio(1, 1));
        assert_eq!(lift.lines[2].terms, vec![(0, ratio(1, 2))]);

        let prescribed = [(0, ratio(2, 1)), (1, ratio(4, 1)), (2, ratio(3, 1))];
        let fixed = resolve(3, &equations, &prescribed).unwrap();
        assert!(fixed.free_roots.is_empty());
        assert_eq!(fixed.lines[2].offset, ratio(3, 1));
        assert!(resolve(
            3,
            &equations,
            &[(0, ratio(2, 1)), (1, ratio(4, 1)), (2, ratio(2, 1))]
        )
        .is_err());
        assert!(resolve(3, &equations, &[(2, ratio(3, 1))]).is_err());
    }

    #[test]
    fn rational_arithmetic_is_canonical_and_overflow_checked() {
        let reduced = Ratio::new(6, -8).unwrap();
        assert_eq!(reduced, ratio(-3, 4));
        assert_eq!((reduced.numerator(), reduced.denominator()), (-3, 4));
        assert_eq!(Ratio::new(i128::MIN, i128::MIN).unwrap(), Ratio::one());
        assert!(Ratio::new(1, 0).is_err());
        assert!(ratio(i128::MAX, 1).checked_add(Ratio::one()).is_err());
    }
}
