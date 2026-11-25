//! Outward-rounded Bernstein certificate for local regularity of complete
//! rational trivariate Bézier maps on the closed unit cube.
//!
//! This certifies only the sign of the Jacobian determinant. It does not
//! establish injectivity, boundary compatibility, or partition of a solid.

use crate::planar::{binomial_ratio, Interval};

/// Explicit search limits. The leaf budget counts terminal regions in the
/// binary subdivision tree; `max_depth` bounds the number of cuts on a path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CertificateBudget {
    pub max_leaves: usize,
    pub max_depth: u8,
}

impl Default for CertificateBudget {
    fn default() -> Self {
        Self {
            max_leaves: 4096,
            max_depth: 20,
        }
    }
}

/// Final status for the local determinant sign.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CertificateStatus {
    /// Every leaf has a rigorously positive determinant lower bound.
    Valid,
    /// A dyadic corner is rigorously nonpositive, or the input is invalid.
    Invalid,
    /// The requested conclusion was not proved within the arithmetic/budget.
    Indeterminate,
}

/// One exact dyadic interval `[numerator/2^depth,(numerator+1)/2^depth]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DyadicInterval {
    pub numerator: u32,
    pub depth: u8,
}

/// Parameter box for one Bernstein leaf, in `(u,v,w)` order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DyadicBox {
    pub axes: [DyadicInterval; 3],
}

/// A closed-box corner where a nonpositive determinant is proved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonPositiveWitness {
    pub domain: DyadicBox,
    /// `false` selects the lower endpoint; `true` the upper endpoint.
    pub upper_corner: [bool; 3],
    /// Conservative upper bound for `det(DX)` in mm³ at the witness.
    pub determinant_upper_bound_mm3: f64,
}

/// Auditable status, bound, witness localization, and actual search use.
#[derive(Debug, Clone, PartialEq)]
pub struct SpatialCertificate {
    pub status: CertificateStatus,
    /// Conservative positive lower bound, present only for `Valid`.
    pub determinant_lower_bound_mm3: Option<f64>,
    /// Localized nonpositive point, present when a witness is established.
    pub witness: Option<NonPositiveWitness>,
    pub leaves_used: usize,
    pub max_depth_reached: u8,
    pub reason: Option<&'static str>,
}

#[derive(Clone)]
struct TensorBernstein {
    degrees: [usize; 3],
    values: Vec<Interval>,
}

impl TensorBernstein {
    fn flat([i, j, k]: [usize; 3], d: [usize; 3]) -> usize {
        (i * (d[1] + 1) + j) * (d[2] + 1) + k
    }

    fn get(&self, p: [usize; 3]) -> Interval {
        self.values[Self::flat(p, self.degrees)]
    }

    fn derivative(&self, axis: usize) -> Self {
        let mut degrees = self.degrees;
        let degree = degrees[axis];
        degrees[axis] -= 1;
        let mut values = Vec::with_capacity((degrees[0] + 1) * (degrees[1] + 1) * (degrees[2] + 1));
        for i in 0..=degrees[0] {
            for j in 0..=degrees[1] {
                for k in 0..=degrees[2] {
                    let p = [i, j, k];
                    let mut next = p;
                    next[axis] += 1;
                    values.push(
                        self.get(next)
                            .sub(self.get(p))
                            .mul(Interval::exact(degree as f64)),
                    );
                }
            }
        }
        Self { degrees, values }
    }

    fn add(&self, rhs: &Self) -> Self {
        debug_assert_eq!(self.degrees, rhs.degrees);
        Self {
            degrees: self.degrees,
            values: self
                .values
                .iter()
                .zip(&rhs.values)
                .map(|(a, b)| a.add(*b))
                .collect(),
        }
    }

    fn negate(&self) -> Self {
        Self {
            degrees: self.degrees,
            values: self
                .values
                .iter()
                .map(|v| v.mul(Interval::exact(-1.0)))
                .collect(),
        }
    }

    fn product(&self, rhs: &Self) -> Self {
        let degree = std::array::from_fn(|a| self.degrees[a] + rhs.degrees[a]);
        let count = (degree[0] + 1) * (degree[1] + 1) * (degree[2] + 1);
        let mut values = vec![Interval::exact(0.0); count];
        for i in 0..=self.degrees[0] {
            for j in 0..=self.degrees[1] {
                for k in 0..=self.degrees[2] {
                    let a = [i, j, k];
                    for m in 0..=rhs.degrees[0] {
                        for n in 0..=rhs.degrees[1] {
                            for q in 0..=rhs.degrees[2] {
                                let b = [m, n, q];
                                let factor = binomial_ratio(self.degrees[0], i, rhs.degrees[0], m)
                                    .mul(binomial_ratio(self.degrees[1], j, rhs.degrees[1], n))
                                    .mul(binomial_ratio(self.degrees[2], k, rhs.degrees[2], q));
                                let term = self.get(a).mul(rhs.get(b)).mul(factor);
                                let out = [i + m, j + n, k + q];
                                let index = Self::flat(out, degree);
                                values[index] = values[index].add(term);
                            }
                        }
                    }
                }
            }
        }
        Self {
            degrees: degree,
            values,
        }
    }

    fn bisect(&self, axis: usize) -> [Self; 2] {
        let mut children = [
            Self {
                degrees: self.degrees,
                values: vec![Interval::exact(0.0); self.values.len()],
            },
            Self {
                degrees: self.degrees,
                values: vec![Interval::exact(0.0); self.values.len()],
            },
        ];
        let n = self.degrees[axis];
        let other: Vec<usize> = (0..3).filter(|a| *a != axis).collect();
        for a in 0..=self.degrees[other[0]] {
            for b in 0..=self.degrees[other[1]] {
                let mut line = Vec::with_capacity(n + 1);
                for t in 0..=n {
                    let mut p = [0; 3];
                    p[axis] = t;
                    p[other[0]] = a;
                    p[other[1]] = b;
                    line.push(self.get(p));
                }
                for depth in 0..=n {
                    let mut lp = [0; 3];
                    lp[axis] = depth;
                    lp[other[0]] = a;
                    lp[other[1]] = b;
                    children[0].values[Self::flat(lp, self.degrees)] = line[0];
                    let mut rp = [0; 3];
                    rp[axis] = n - depth;
                    rp[other[0]] = a;
                    rp[other[1]] = b;
                    children[1].values[Self::flat(rp, self.degrees)] = line[n - depth];
                    for t in 0..(n - depth) {
                        line[t] = line[t].add(line[t + 1]).mul(Interval::exact(0.5));
                    }
                }
            }
        }
        children
    }

    fn variation_axis(&self) -> usize {
        let mut scores = [0.0_f64; 3];
        for axis in 0..3 {
            if self.degrees[axis] == 0 {
                continue;
            }
            let other: Vec<usize> = (0..3).filter(|a| *a != axis).collect();
            for a in 0..=self.degrees[other[0]] {
                for b in 0..=self.degrees[other[1]] {
                    for t in 0..self.degrees[axis] {
                        let mut p = [0; 3];
                        p[axis] = t;
                        p[other[0]] = a;
                        p[other[1]] = b;
                        let mut q = p;
                        q[axis] += 1;
                        let x = self.get(p);
                        let y = self.get(q);
                        let xm = x.lo * 0.5 + x.hi * 0.5;
                        let ym = y.lo * 0.5 + y.hi * 0.5;
                        scores[axis] = scores[axis].max((ym - xm).abs());
                    }
                }
            }
        }
        (0..3)
            .filter(|a| self.degrees[*a] > 0)
            .max_by(|a, b| scores[*a].total_cmp(&scores[*b]))
            .unwrap_or(0)
    }

    fn corner_nonpositive(&self) -> Option<([bool; 3], Interval)> {
        for mask in 0..8 {
            let p = std::array::from_fn(|axis| {
                if mask & (1 << axis) == 0 {
                    0
                } else {
                    self.degrees[axis]
                }
            });
            let value = self.get(p);
            if value.hi <= 0.0 {
                return Some((std::array::from_fn(|axis| mask & (1 << axis) != 0), value));
            }
        }
        None
    }
}

#[derive(Clone, Copy)]
struct SearchNode {
    domain: DyadicBox,
    depth: u8,
}

/// Certify the Jacobian determinant of one complete homogeneous map.
///
/// Controls are ordered with the last parameter varying fastest and stored as
/// `[X,Y,Z,W]` integer coefficients. Every input is enclosed after exact
/// binary scaling by `2^-63`; Bernstein products and subdivision use outward
/// intervals. The determinant lower bound uses the physical unit represented
/// by the control ratios (mm³ for native bracket blocks). The routine makes no
/// global injectivity claim.
pub fn certify_spatial(
    degrees: [usize; 3],
    controls: &[[i64; 4]],
    budget: CertificateBudget,
) -> SpatialCertificate {
    let invalid = |reason| SpatialCertificate {
        status: CertificateStatus::Invalid,
        determinant_lower_bound_mm3: None,
        witness: None,
        leaves_used: 0,
        max_depth_reached: 0,
        reason: Some(reason),
    };
    let indeterminate = |reason, leaves, depth| SpatialCertificate {
        status: CertificateStatus::Indeterminate,
        determinant_lower_bound_mm3: None,
        witness: None,
        leaves_used: leaves,
        max_depth_reached: depth,
        reason: Some(reason),
    };
    let supported = matches!(degrees, [1, 1, 1] | [2, 1, 1] | [2, 2, 1]);
    let count = (degrees[0] + 1)
        .checked_mul(degrees[1] + 1)
        .and_then(|x| x.checked_mul(degrees[2] + 1));
    if !supported || count != Some(controls.len()) {
        return invalid("unsupported degree family or control count");
    }
    if controls.iter().any(|c| c[3] <= 0) {
        return invalid("homogeneous weights must be strictly positive i64");
    }
    if budget.max_leaves == 0 || budget.max_depth > 31 {
        return indeterminate("invalid search budget", 0, 0);
    }

    let components: [TensorBernstein; 4] = std::array::from_fn(|c| TensorBernstein {
        degrees,
        values: controls
            .iter()
            .map(|point| Interval::from_i64_scaled(point[c]))
            .collect(),
    });
    let derivatives: [[TensorBernstein; 3]; 4] = std::array::from_fn(|component| {
        std::array::from_fn(|axis| components[component].derivative(axis))
    });
    // Matrix rows are homogeneous components (X,Y,Z,W), columns are
    // (H_u,H_v,H_w,H), whose determinant equals W^4 det(DX).
    let matrix: [[TensorBernstein; 4]; 4] = std::array::from_fn(|row| {
        [
            derivatives[row][0].clone(),
            derivatives[row][1].clone(),
            derivatives[row][2].clone(),
            components[row].clone(),
        ]
    });
    let numerator = determinant(&matrix);
    let w = &components[3];
    let w_lower = w.values.iter().map(|x| x.lo).fold(f64::INFINITY, f64::min);
    let w_upper = w
        .values
        .iter()
        .map(|x| x.hi)
        .fold(f64::NEG_INFINITY, f64::max);
    let w_range = Interval {
        lo: w_lower,
        hi: w_upper,
    };
    if !w_range.finite() || w_lower <= 0.0 {
        return indeterminate("positive-weight interval enclosure failed", 0, 0);
    }
    let w4 = w_range.mul(w_range).mul(w_range).mul(w_range);
    if !w4.finite() || w4.lo <= 0.0 {
        return indeterminate("weight fourth-power enclosure failed", 0, 0);
    }

    let start = SearchNode {
        domain: DyadicBox {
            axes: [DyadicInterval {
                numerator: 0,
                depth: 0,
            }; 3],
        },
        depth: 0,
    };
    let mut stack = vec![(numerator, start)];
    let mut leaves = 0usize;
    let mut leaves_created = 1usize;
    let mut max_depth = 0u8;
    let mut lower_numerator = f64::INFINITY;
    while let Some((poly, node)) = stack.pop() {
        max_depth = max_depth.max(node.depth);
        if let Some((corner, value)) = poly.corner_nonpositive() {
            let quotient = Interval::exact(value.hi).div(w4);
            if !quotient.finite() {
                return indeterminate("witness quotient is nonfinite", leaves_created, max_depth);
            }
            let determinant_upper_bound_mm3 = if value.hi == 0.0 { 0.0 } else { quotient.hi };
            return SpatialCertificate {
                status: CertificateStatus::Invalid,
                determinant_lower_bound_mm3: None,
                witness: Some(NonPositiveWitness {
                    domain: node.domain,
                    upper_corner: corner,
                    determinant_upper_bound_mm3,
                }),
                leaves_used: leaves_created,
                max_depth_reached: max_depth,
                reason: Some("a dyadic corner has a nonpositive determinant upper bound"),
            };
        }
        let min_lower = poly
            .values
            .iter()
            .map(|x| x.lo)
            .fold(f64::INFINITY, f64::min);
        if min_lower > 0.0 {
            lower_numerator = lower_numerator.min(min_lower);
            leaves += 1;
            continue;
        }
        if node.depth >= budget.max_depth {
            return indeterminate(
                "depth budget exhausted before positivity was proved",
                leaves_created,
                max_depth,
            );
        }
        if leaves_created >= budget.max_leaves {
            return indeterminate(
                "leaf budget exhausted before positivity was proved",
                leaves_created,
                max_depth,
            );
        }
        let axis = poly.variation_axis();
        if node.domain.axes[axis].depth >= 31 {
            return indeterminate(
                "dyadic domain depth exceeds supported locator",
                leaves_created,
                max_depth,
            );
        }
        let [left, right] = poly.bisect(axis);
        let old = node.domain.axes[axis];
        let mut left_domain = node.domain;
        let mut right_domain = node.domain;
        let Some(left_numerator) = old.numerator.checked_mul(2) else {
            return indeterminate("dyadic locator overflow", leaves_created, max_depth);
        };
        let Some(right_numerator) = left_numerator.checked_add(1) else {
            return indeterminate("dyadic locator overflow", leaves_created, max_depth);
        };
        left_domain.axes[axis] = DyadicInterval {
            numerator: left_numerator,
            depth: old.depth + 1,
        };
        right_domain.axes[axis] = DyadicInterval {
            numerator: right_numerator,
            depth: old.depth + 1,
        };
        let child_depth = node.depth + 1;
        stack.push((
            right,
            SearchNode {
                domain: right_domain,
                depth: child_depth,
            },
        ));
        stack.push((
            left,
            SearchNode {
                domain: left_domain,
                depth: child_depth,
            },
        ));
        leaves_created += 1;
    }
    if leaves == 0 || !lower_numerator.is_finite() {
        return indeterminate(
            "no finite positive Bernstein lower bound",
            leaves_created,
            max_depth,
        );
    }
    let lower = Interval::exact(lower_numerator).div(w4);
    if !lower.finite() || lower.lo <= 0.0 {
        return indeterminate(
            "determinant quotient is not strictly positive",
            leaves_created,
            max_depth,
        );
    }
    SpatialCertificate {
        status: CertificateStatus::Valid,
        determinant_lower_bound_mm3: Some(lower.lo),
        witness: None,
        leaves_used: leaves_created,
        max_depth_reached: max_depth,
        reason: None,
    }
}

fn determinant(matrix: &[[TensorBernstein; 4]; 4]) -> TensorBernstein {
    let mut permutations = Vec::with_capacity(24);
    permute([0, 1, 2, 3], 0, &mut permutations);
    let base_degrees = matrix[0][3].degrees;
    let target: [usize; 3] = std::array::from_fn(|a| 4 * base_degrees[a] - 1);
    let mut result = TensorBernstein {
        degrees: target,
        values: vec![Interval::exact(0.0); (target[0] + 1) * (target[1] + 1) * (target[2] + 1)],
    };
    for permutation in permutations {
        let mut term = matrix[permutation[0]][0].clone();
        for col in 1..4 {
            term = term.product(&matrix[permutation[col]][col]);
        }
        let inversions = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
            .filter(|(i, j)| permutation[*i] > permutation[*j])
            .count();
        result = result.add(&if inversions % 2 == 0 {
            term
        } else {
            term.negate()
        });
    }
    result
}

fn permute(values: [usize; 4], depth: usize, out: &mut Vec<[usize; 4]>) {
    if depth == 4 {
        out.push(values);
        return;
    }
    for i in depth..4 {
        let mut next = values;
        next.swap(depth, i);
        permute(next, depth + 1, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn affine_controls(scale: i64) -> Vec<[i64; 4]> {
        (0..8)
            .map(|i| {
                [
                    ((i as i64 >> 2) & 1) * scale,
                    2 * ((i as i64 >> 1) & 1) * scale,
                    3 * (i as i64 & 1) * scale,
                    scale,
                ]
            })
            .collect()
    }

    fn affine_controls_degree(degrees: [usize; 3]) -> Vec<[i64; 4]> {
        let mut controls = Vec::new();
        for i in 0..=degrees[0] {
            for j in 0..=degrees[1] {
                for k in 0..=degrees[2] {
                    controls.push([
                        2 * (i as i64) / degrees[0] as i64,
                        2 * (j as i64) / degrees[1] as i64,
                        2 * (k as i64) / degrees[2] as i64,
                        2,
                    ]);
                }
            }
        }
        controls
    }

    fn near_singular_map(off_diagonal: i64) -> Vec<[i64; 4]> {
        let scale = 16_384_i64;
        let diagonal = [-5_120_i64, 11_264_i64];
        (0..8)
            .map(|i| {
                let u = ((i >> 2) & 1) as usize;
                let v = ((i >> 1) & 1) as i64;
                let w = (i & 1) as i64;
                [
                    (u as i64) * scale,
                    diagonal[u] * v + 16 * w,
                    off_diagonal * v + diagonal[u] * w,
                    scale,
                ]
            })
            .collect()
    }

    fn contains_dyadic(interval: DyadicInterval, point_num: u32, point_depth: u8) -> bool {
        let common = interval.depth.max(point_depth);
        let lo = (interval.numerator as u64) << (common - interval.depth);
        let hi = ((interval.numerator as u64) + 1) << (common - interval.depth);
        let point = (point_num as u64) << (common - point_depth);
        point >= lo && point <= hi
    }

    #[test]
    fn affine_and_projective_scale() {
        let a = certify_spatial([1, 1, 1], &affine_controls(1), CertificateBudget::default());
        let b = certify_spatial([1, 1, 1], &affine_controls(7), CertificateBudget::default());
        assert_eq!(a.status, CertificateStatus::Valid);
        assert_eq!(b.status, CertificateStatus::Valid);
        let alo = a.determinant_lower_bound_mm3.unwrap();
        let blo = b.determinant_lower_bound_mm3.unwrap();
        assert!(alo > 0.0 && alo <= 6.0);
        assert!(blo > 0.0 && blo <= 6.0);
        assert!((alo - blo).abs() < 1e-10);
        assert_eq!(a.leaves_used, 1);
    }

    #[test]
    fn detects_negative_well_missed_by_quarter_grid_and_localizes_it() {
        let controls = near_singular_map(16);
        let epsilon = 1.0 / 1024.0;
        for numerator in 0..=4 {
            let u = numerator as f64 / 4.0;
            let d = u - 5.0 / 16.0;
            assert!(d * d - epsilon * epsilon > 0.0);
        }
        let cert = certify_spatial(
            [1, 1, 1],
            &controls,
            CertificateBudget {
                max_leaves: 1024,
                max_depth: 16,
            },
        );
        assert_eq!(cert.status, CertificateStatus::Invalid, "{cert:?}");
        let witness = cert.witness.unwrap();
        assert!(witness.determinant_upper_bound_mm3 <= 0.0);
        assert!(
            contains_dyadic(witness.domain.axes[0], 5, 4),
            "{:?}",
            witness.domain.axes[0]
        );
        assert_eq!(
            witness.domain.axes[0],
            DyadicInterval {
                numerator: 4,
                depth: 4
            }
        );
        assert_eq!(
            witness.domain.axes[1],
            DyadicInterval {
                numerator: 0,
                depth: 0
            }
        );
        assert_eq!(
            witness.domain.axes[2],
            DyadicInterval {
                numerator: 0,
                depth: 0
            }
        );
        assert_eq!(witness.upper_corner, [true, false, false]);
        assert!(cert.leaves_used <= 1024 && cert.max_depth_reached <= 16);
    }

    #[test]
    fn certifies_strictly_positive_narrow_minimum_with_bounded_subdivision() {
        let cert = certify_spatial(
            [1, 1, 1],
            &near_singular_map(-16),
            CertificateBudget {
                max_leaves: 4096,
                max_depth: 20,
            },
        );
        assert_eq!(cert.status, CertificateStatus::Valid, "{cert:?}");
        let lower = cert.determinant_lower_bound_mm3.unwrap();
        let exact_min = (1.0_f64 / 1024.0).powi(2);
        assert!(lower > 0.0 && lower <= exact_min);
        assert_eq!(cert.leaves_used, 5);
        assert_eq!(cert.max_depth_reached, 4);
    }

    #[test]
    fn singularity_and_exhaustion_are_indeterminate() {
        let singular = near_singular_map(0);
        let result = certify_spatial(
            [1, 1, 1],
            &singular,
            CertificateBudget {
                max_leaves: 64,
                max_depth: 8,
            },
        );
        assert_eq!(
            result.status,
            CertificateStatus::Indeterminate,
            "{result:?}"
        );
        let limited = certify_spatial(
            [1, 1, 1],
            &near_singular_map(-16),
            CertificateBudget {
                max_leaves: 2,
                max_depth: 1,
            },
        );
        assert_eq!(limited.status, CertificateStatus::Indeterminate);
        assert!(limited.leaves_used <= 2);
        assert!(limited.determinant_lower_bound_mm3.is_none());
    }

    #[test]
    fn rejects_invalid_degrees_counts_and_weights() {
        assert_eq!(
            certify_spatial([1, 2, 1], &vec![[0; 4]; 12], CertificateBudget::default()).status,
            CertificateStatus::Invalid
        );
        let mut controls = affine_controls(1);
        controls.pop();
        assert_eq!(
            certify_spatial([1, 1, 1], &controls, CertificateBudget::default()).status,
            CertificateStatus::Invalid
        );
        let mut controls = affine_controls(1);
        controls[0][3] = 0;
        assert_eq!(
            certify_spatial([1, 1, 1], &controls, CertificateBudget::default()).status,
            CertificateStatus::Invalid
        );
    }

    #[test]
    fn supports_all_three_degree_families() {
        for degrees in [[1, 1, 1], [2, 1, 1], [2, 2, 1]] {
            let controls = affine_controls_degree(degrees);
            let cert = certify_spatial(degrees, &controls, CertificateBudget::default());
            assert_eq!(
                cert.status,
                CertificateStatus::Valid,
                "{degrees:?}: {cert:?}"
            );
            let lower = cert.determinant_lower_bound_mm3.unwrap();
            assert!(lower > 0.0 && lower <= 1.0);
        }
    }

    #[test]
    fn certifies_shear_depending_on_the_third_parameter() {
        let controls = (0..8)
            .map(|i| {
                let u = ((i >> 2) & 1) as i64;
                let v = ((i >> 1) & 1) as i64;
                let w = (i & 1) as i64;
                [4 * u + 2 * w, 4 * v + w, 4 * w, 4]
            })
            .collect::<Vec<_>>();
        let cert = certify_spatial([1, 1, 1], &controls, CertificateBudget::default());
        assert_eq!(cert.status, CertificateStatus::Valid, "{cert:?}");
        let lower = cert.determinant_lower_bound_mm3.unwrap();
        assert!(lower > 0.0 && lower <= 1.0);
        assert_eq!(cert.leaves_used, 1);
    }

    #[test]
    fn certifies_supported_e6_extruded_bracket_blocks() {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        assert_eq!(blocks.len(), 36);
        for block in blocks {
            let cert =
                certify_spatial(block.degrees, &block.controls, CertificateBudget::default());
            assert_eq!(
                cert.status,
                CertificateStatus::Valid,
                "{}: {cert:?}",
                block.id
            );
            assert!(cert.determinant_lower_bound_mm3.unwrap() > 0.0);
        }
    }
}
