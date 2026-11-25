//! Construct planar section patches and oriented boundary spans directly from
//! sketch primitives. Cell counts and joins are determined before sweeping.
use crate::{
    bracket::BracketBlock, constraint_system::Ratio as R, nurbs_cross::NurbsCurve2,
    nurbs_curve::normalize, sketch_junction::Continuity,
};

type Error = &'static str;
type H = [R; 3];

#[derive(Debug, Clone, Copy)]
pub enum Primitive {
    Rectangle { min: [R; 2], max: [R; 2] },
    Circle { center: [R; 2], radius: R },
    Ellipse { center: [R; 2], radii: [R; 2] },
}

#[derive(Debug, Clone, Copy)]
pub enum Operation {
    /// Translate the filled XY sketch along positive Z.
    Extrude { length: i64 },
    /// Full turn about sketch Y, with the whole section strictly in x > 0.
    Revolve,
}

#[derive(Debug, Clone)]
pub struct SectionPatch {
    pub degrees: [usize; 2],
    pub controls: Vec<[i64; 3]>,
}

pub struct PrimitivePlan {
    pub contour: NurbsCurve2,
    /// Geometric join at the start of each emitted boundary span.
    pub junctions: Vec<Continuity>,
    pub section_patches: Vec<SectionPatch>,
    operation: Operation,
}

fn int(n: i64) -> R {
    R::new(n.into(), 1).expect("integer")
}

fn contour(spans: &[Vec<H>]) -> Result<NurbsCurve2, Error> {
    let degree = spans[0].len() - 1;
    let mut controls = vec![spans[0][0]];
    let mut knots = vec![0; degree + 1];
    for (k, s) in spans.iter().enumerate() {
        controls.extend_from_slice(&s[1..]);
        knots.extend(std::iter::repeat_n(k as i64 + 1, degree));
    }
    knots.push(spans.len() as i64);
    Ok(NurbsCurve2 {
        degree,
        knots,
        controls: normalize(&controls)?,
    })
}

fn rectangle(min: [R; 2], max: [R; 2], operation: Operation) -> Result<PrimitivePlan, Error> {
    for d in 0..2 {
        if max[d].checked_sub(min[d])?.numerator() <= 0 {
            return Err("rectangle dimensions must be positive");
        }
    }
    let point = |x, y| [x, y, R::one()];
    let [a, b, c, d] = [
        point(min[0], min[1]),
        point(max[0], min[1]),
        point(max[0], max[1]),
        point(min[0], max[1]),
    ];
    Ok(PrimitivePlan {
        contour: contour(&[vec![a, b], vec![b, c], vec![c, d], vec![d, a]])?,
        junctions: vec![Continuity::Corner; 4],
        section_patches: vec![SectionPatch {
            degrees: [1, 1],
            controls: normalize(&[a, d, b, c])?,
        }],
        operation,
    })
}

fn rotate([mut x, mut y, w]: [i64; 3], quarter: usize) -> [i64; 3] {
    for _ in 0..quarter {
        (x, y) = (-y, x);
    }
    [x, y, w]
}

// Both parameter directions run counterclockwise. Alternating the gauges
// makes the endpoint weights match around a closed four-quarter ring.
fn arc(q: usize) -> [[i64; 3]; 3] {
    let base = if q.is_multiple_of(2) {
        [[1, 0, 1], [1, 1, 1], [0, 2, 2]]
    } else {
        [[2, 0, 2], [1, 1, 1], [0, 1, 1]]
    };
    base.map(|p| rotate(p, q))
}

fn ellipse(center: [R; 2], radii: [R; 2], operation: Operation) -> Result<PrimitivePlan, Error> {
    if radii.iter().any(|r| r.numerator() <= 0) {
        return Err("ellipse radii must be positive");
    }
    let transform = |[x, y, w]: [i64; 3]| -> Result<H, Error> {
        let w = int(w);
        Ok([
            center[0]
                .checked_mul(w)?
                .checked_add(radii[0].checked_mul(int(x))?)?,
            center[1]
                .checked_mul(w)?
                .checked_add(radii[1].checked_mul(int(y))?)?,
            w,
        ])
    };
    let mut core = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            core.push(transform([2 - i - j, i - j, 4])?);
        }
    }
    let mut section_patches = vec![SectionPatch {
        degrees: [2, 2],
        controls: normalize(&core)?,
    }];
    let mut boundary = Vec::new();
    for q in 0..4 {
        let mut net = Vec::new();
        for p in [[2, 0, 4], [1, 1, 4], [0, 2, 4]] {
            net.push(transform(rotate(p, q))?);
        }
        // The inner edge has weight four. Give the unit circle the same
        // global gauge so the radial seams agree from both neighboring cells.
        let outer = arc(q)
            .map(|p| transform(p.map(|c| 4 * c)))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        net.extend_from_slice(&outer);
        section_patches.push(SectionPatch {
            degrees: [1, 2],
            controls: normalize(&net)?,
        });
        boundary.push(outer);
    }
    Ok(PrimitivePlan {
        contour: contour(&boundary)?,
        junctions: vec![Continuity::Curvature; 4],
        section_patches,
        operation,
    })
}

impl Primitive {
    pub fn plan(self, operation: Operation) -> Result<PrimitivePlan, Error> {
        match operation {
            Operation::Extrude { length } if length <= 0 => {
                return Err("extrusion length must be positive")
            }
            Operation::Revolve => {
                let min_r = match self {
                    Self::Rectangle { min, .. } => min[0],
                    Self::Circle { center, radius } => center[0].checked_sub(radius)?,
                    Self::Ellipse { center, radii } => center[0].checked_sub(radii[0])?,
                };
                if min_r.numerator() <= 0 {
                    return Err("full-turn template requires positive radial clearance");
                }
            }
            _ => {}
        }
        match self {
            Self::Rectangle { min, max } => rectangle(min, max, operation),
            Self::Circle { center, radius } => ellipse(center, [radius; 2], operation),
            Self::Ellipse { center, radii } => ellipse(center, radii, operation),
        }
    }
}

impl PrimitivePlan {
    pub fn block_count(&self) -> usize {
        self.section_patches.len()
            * match self.operation {
                Operation::Extrude { .. } => 1,
                Operation::Revolve => 4,
            }
    }

    /// Instantiate the already partitioned section. No fitting or search for
    /// knots, seams, or cylinder intersections is needed for an isolated sweep.
    pub fn build(&self) -> Result<Vec<BracketBlock>, Error> {
        let mut blocks = Vec::with_capacity(self.block_count());
        for (section, patch) in self.section_patches.iter().enumerate() {
            let [u, v] = patch.degrees;
            match self.operation {
                Operation::Extrude { length } => {
                    let mut net = Vec::new();
                    for &[x, y, w] in &patch.controls {
                        for z in [0, length] {
                            net.push([int(x), int(y), int(z).checked_mul(int(w))?, int(w)]);
                        }
                    }
                    blocks.push(BracketBlock {
                        id: format!("section_{section}"),
                        degrees: [u, v, 1],
                        controls: normalize(&net)?,
                    });
                }
                Operation::Revolve => {
                    for q in 0..4 {
                        let mut net = Vec::new();
                        for &[x, y, w] in &patch.controls {
                            for [c, s, k] in arc(q) {
                                net.push([
                                    int(x).checked_mul(int(c))?,
                                    int(y).checked_mul(int(k))?,
                                    int(x).checked_mul(int(s))?,
                                    int(w).checked_mul(int(k))?,
                                ]);
                            }
                        }
                        blocks.push(BracketBlock {
                            id: format!("section_{section}_turn_{q}"),
                            degrees: [u, v, 2],
                            controls: normalize(&net)?,
                        });
                    }
                }
            }
        }
        Ok(blocks)
    }
}
