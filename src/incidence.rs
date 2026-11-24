//! Exact whole-face incidence for degree-one and degree-two integer Bézier nets.

use super::bracket::BracketBlock;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BlockSide {
    pub block: String,
    pub axis: usize,
    pub upper: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Orientation {
    pub swap: bool,
    pub flip_first: bool,
    pub flip_second: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interface {
    pub first: BlockSide,
    pub second: BlockSide,
    /// Transpose the second trace first, then reverse the requested output axes.
    pub orientation: Orientation,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Incidence {
    pub interfaces: Vec<Interface>,
    pub exterior: Vec<BlockSide>,
}

#[derive(Clone, PartialEq, Eq)]
struct Trace {
    degree: [usize; 2],
    points: Vec<[i128; 4]>,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

impl Trace {
    fn point(&self, i: usize, j: usize) -> [i128; 4] {
        self.points[i * (self.degree[1] + 1) + j]
    }

    fn corner_key(&self) -> Result<[[i128; 4]; 4], &'static str> {
        let [p, q] = self.degree;
        let mut corners = [
            self.point(0, 0),
            self.point(0, q),
            self.point(p, 0),
            self.point(p, q),
        ];
        for point in &mut corners {
            let scale = point.iter().copied().fold(0, gcd);
            *point = point.map(|v| v / scale);
        }
        corners.sort();
        if corners.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err("collapsed face");
        }
        Ok(corners)
    }

    fn oriented(&self, o: Orientation) -> Self {
        let [p, q] = self.degree;
        let degree = if o.swap { [q, p] } else { [p, q] };
        let mut points = Vec::new();
        for i in 0..=degree[0] {
            for j in 0..=degree[1] {
                let a = if o.flip_first { degree[0] - i } else { i };
                let b = if o.flip_second { degree[1] - j } else { j };
                points.push(if o.swap {
                    self.point(b, a)
                } else {
                    self.point(a, b)
                });
            }
        }
        Self { degree, points }
    }

    fn biquadratic(mut self) -> Self {
        for axis in 0..2 {
            if self.degree[axis] == 2 {
                continue;
            }
            let mut degree = self.degree;
            degree[axis] = 2;
            let mut points = Vec::new();
            for i in 0..=degree[0] {
                for j in 0..=degree[1] {
                    let (index, a, b) = if axis == 0 {
                        (i, self.point(0, j), self.point(1, j))
                    } else {
                        (j, self.point(i, 0), self.point(i, 1))
                    };
                    points.push(std::array::from_fn(|k| match index {
                        0 => 2 * a[k],
                        1 => a[k] + b[k],
                        _ => 2 * b[k],
                    }));
                }
            }
            self = Self { degree, points };
        }
        // Two elevations multiply the i64 input bound by at most four.
        let scale = self.points.iter().flatten().copied().fold(0, gcd);
        for p in &mut self.points {
            *p = p.map(|v| v / scale);
        }
        self
    }
}

fn extract(block: &BracketBlock, axis: usize, upper: bool) -> Trace {
    let [p, q, r] = block.degrees;
    let mut points = Vec::new();
    let fixed = if upper { block.degrees[axis] } else { 0 };
    for i in 0..=p {
        for j in 0..=q {
            for k in 0..=r {
                if [i, j, k][axis] == fixed {
                    points.push(block.controls[(i * (q + 1) + j) * (r + 1) + k].map(i128::from));
                }
            }
        }
    }
    let degree = match axis {
        0 => [q, r],
        1 => [p, r],
        _ => [p, q],
    };
    Trace { degree, points }
}

/// Reconstruct incidences within a whole-face template. Proving disjoint
/// interiors and ruling out partial contacts require a geometric family audit.
pub fn enumerate_faces(blocks: &[BracketBlock]) -> Result<Incidence, &'static str> {
    let mut ids = BTreeSet::new();
    let mut buckets = BTreeMap::<_, Vec<(BlockSide, Trace)>>::new();
    for block in blocks {
        if block.id.is_empty() || !ids.insert(&block.id) {
            return Err("invalid block identity");
        }
        if block.degrees.iter().any(|d| !(1..=2).contains(d))
            || block.controls.len() != block.degrees.iter().map(|d| d + 1).product::<usize>()
            || block.controls.iter().any(|p| p[3] <= 0)
        {
            return Err("invalid block net");
        }
        for axis in 0..3 {
            for upper in [false, true] {
                let trace = extract(block, axis, upper);
                let key = trace.corner_key()?;
                let side = BlockSide {
                    block: block.id.clone(),
                    axis,
                    upper,
                };
                buckets.entry(key).or_default().push((side, trace));
            }
        }
    }
    let mut result = Incidence {
        interfaces: Vec::new(),
        exterior: Vec::new(),
    };
    for mut faces in buckets.into_values() {
        faces.sort_by(|a, b| a.0.cmp(&b.0));
        match faces.as_slice() {
            [(side, _)] => result.exterior.push(side.clone()),
            [(a, ta), (b, tb)] => {
                let first = ta.clone().biquadratic();
                let second = tb.clone().biquadratic();
                let mut matches = Vec::new();
                for swap in [false, true] {
                    for flip_first in [false, true] {
                        for flip_second in [false, true] {
                            let orientation = Orientation {
                                swap,
                                flip_first,
                                flip_second,
                            };
                            if first == second.oriented(orientation) {
                                matches.push(orientation);
                            }
                        }
                    }
                }
                if matches.len() != 1 {
                    return Err("nonunique or incompatible trace");
                }
                result.interfaces.push(Interface {
                    first: a.clone(),
                    second: b.clone(),
                    orientation: matches[0],
                });
            }
            _ => return Err("nonmanifold face"),
        }
    }
    result.interfaces.sort_by(|a, b| a.first.cmp(&b.first));
    result.exterior.sort();
    Ok(result)
}
