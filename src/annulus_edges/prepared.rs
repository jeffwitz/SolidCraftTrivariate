//! A fixed multi-span NURBS atlas for the gallery annulus dimensions.
//! cargo run --release --example annulus_prepared_volumes -- fillet 2 > model.json
//! Degrees are (meridian, angle, transverse); the last index varies fastest.
#[cfg(test)]
use crate::spatial_certificate::{certify_spatial, CertificateStatus};

type Point = [i64; 2]; // (radius,height) in 1/256 mm
type H = [i64; 3]; // projective meridian, same gauge on every band
const MID: [Point; 4] = [[768, 0], [1024, 1024], [768, 2048], [512, 1024]];
const CORNER: [Point; 4] = [[1024, 0], [1024, 2048], [512, 2048], [512, 0]];
const CENTER: Point = [768, 1024];

use super::EdgeProfile as State;

fn mid(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| (a[i] + b[i]) / 2)
}
fn weighted(p: Point, w: i64) -> H {
    [p[0] * w, p[1] * w, 256 * w]
}
fn towards(c: Point, p: Point) -> Point {
    std::array::from_fn(|i| c[i] + 64 * (p[i] - c[i]).signum())
}
fn segment(a: Point, b: Point, weights: [i64; 3]) -> [H; 3] {
    [
        weighted(a, weights[0]),
        weighted(mid(a, b), weights[1]),
        weighted(b, weights[2]),
    ]
}
fn outer(side: usize, state: State) -> Vec<H> {
    let (start, end, c) = (MID[side], MID[(side + 1) % 4], CORNER[side]);
    let (a, b) = (towards(c, start), towards(c, end));
    let (left, right) = match state {
        State::Sharp | State::Chamfer => {
            let joint = match state {
                State::Sharp => c,
                _ => mid(a, b),
            };
            (segment(a, joint, [4, 4, 5]), segment(joint, b, [5, 6, 8]))
        }
        State::Fillet => {
            // Exact half restrictions of (a,1),(c,1),(b,2).
            let m =
                std::array::from_fn(|i| weighted(a, 1)[i] + weighted(c, 2)[i] + weighted(b, 2)[i]);
            let p = std::array::from_fn(|i| weighted(c, 2)[i] + weighted(b, 4)[i]);
            (
                [weighted(a, 4), weighted(mid(a, c), 4), m],
                [m, p, weighted(b, 8)],
            )
        }
    };
    let mut row = segment(start, a, [4; 3]).to_vec();
    row.extend(&left[1..]);
    row.extend(&right[1..]);
    row.extend(&segment(b, end, [8, 6, 4])[1..]);
    row
}
struct Chart {
    id: String,
    degree: [usize; 3],
    shape: [usize; 3],
    net: Vec<[i64; 4]>,
}

fn charts(state: State, sectors: usize) -> Vec<Chart> {
    assert!(sectors == 2 || sectors == 4);
    let inner = MID.map(|p| mid(p, CENTER));
    let mut result = Vec::new();
    for sweep in 0..sectors {
        let angular_spans = 4 / sectors;
        let mut angles = Vec::new();
        for sub in 0..angular_spans {
            let quarter = sweep * angular_spans + sub;
            // Alternate the homogeneous parametrization with its time reversal.
            // This aligns both Euclidean angular derivatives at each quadrant seam.
            let arc = if quarter & 1 == 0 {
                [[1, 0, 1], [1, 1, 1], [0, 2, 2]]
            } else {
                [[2, 0, 2], [1, 1, 1], [0, 1, 1]]
            };
            for (i, [mut x, mut y, w]) in arc.into_iter().enumerate() {
                if sub > 0 && i == 0 {
                    continue;
                }
                for _ in 0..quarter {
                    (x, y) = (-y, x);
                }
                angles.push([x, y, w]);
            }
        }
        for patch in 0..5 {
            let is_core = patch == 4;
            let rim = if is_core { vec![] } else { outer(patch, state) };
            let depth = if is_core { 9 } else { 2 };
            let mut net = Vec::new();
            for s in 0..9_i64 {
                for &[x, y, w] in &angles {
                    for v in 0..depth {
                        let [r, z, h] = if is_core {
                            // Parallelogram; its four affine sides match the band inner rows.
                            weighted(
                                std::array::from_fn(|i| {
                                    inner[0][i]
                                        + (s * (inner[1][i] - inner[0][i])
                                            + v as i64 * (inner[3][i] - inner[0][i]))
                                            / 8
                                }),
                                4,
                            )
                        } else if v == 0 {
                            rim[s as usize]
                        } else {
                            weighted(
                                std::array::from_fn(|i| {
                                    inner[patch][i]
                                        + s * (inner[(patch + 1) % 4][i] - inner[patch][i]) / 8
                                }),
                                4,
                            )
                        };
                        net.push([r * x, r * y, z * w, h * w]);
                    }
                }
            }
            result.push(Chart {
                id: format!("turn-{sweep}-section-{patch}"),
                degree: [2, 2, if is_core { 2 } else { 1 }],
                shape: [9, angles.len(), depth],
                net,
            });
        }
    }
    result
}

fn knots(n: usize, p: usize) -> Vec<usize> {
    let spans = (n - 1) / p;
    let mut k = vec![0; p + 1];
    for i in 1..spans {
        k.extend(std::iter::repeat_n(i, p));
    }
    k.extend(std::iter::repeat_n(spans, p + 1));
    k
}

#[cfg(test)]
impl Chart {
    fn cells(&self) -> Vec<crate::bracket::BracketBlock> {
        let mut cells = Vec::new();
        let [p, q, r] = self.degree;
        let [nx, ny, nz] = self.shape;
        for i in (0..nx - 1).step_by(p) {
            for j in (0..ny - 1).step_by(q) {
                for k in (0..nz - 1).step_by(r) {
                    let mut controls = Vec::new();
                    for u in i..=i + p {
                        for v in j..=j + q {
                            for w in k..=k + r {
                                controls.push(self.net[(u * ny + v) * nz + w]);
                            }
                        }
                    }
                    cells.push(crate::bracket::BracketBlock {
                        id: format!("{}-{i}-{j}-{k}", self.id),
                        degrees: self.degree,
                        controls,
                    });
                }
            }
        }
        cells
    }
}

use super::{numbers, object, points, string, Json};
pub(super) fn document(state: State, sectors: usize) -> Json {
    let volumes = charts(state, sectors).into_iter().map(|c| {
        let axes = (0..3).map(|a| numbers(knots(c.shape[a], c.degree[a])));
        object([
            ("id", Json::String(c.id)),
            ("degrees", numbers(c.degree)),
            ("shape", numbers(c.shape)),
            ("knots", Json::Array(axes.collect())),
            ("xyzw", points(c.net)),
        ])
    });
    object([
        ("unit", string("mm")),
        ("fastest_index", string("transverse")),
        ("charts", Json::Array(volumes.collect())),
    ])
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::incidence::enumerate_faces;
    #[test]
    fn check_all_three_states_on_the_same_ten_charts() {
        let template = charts(State::Sharp, 2);
        for state in [State::Sharp, State::Chamfer, State::Fillet] {
            let model = charts(state, 2);
            assert_eq!(model.len(), 10);
            let mut lower = f64::INFINITY;
            let mut cells = Vec::new();
            for (a, b) in model.iter().zip(&template) {
                assert_eq!((&a.id, a.degree, a.shape), (&b.id, b.degree, b.shape));
                assert!(a.net.iter().zip(&b.net).all(|(x, y)| x[3] == y[3]));
                cells.extend(a.cells());
            }
            assert_eq!(cells.len(), 128);
            for cell in &cells {
                let cert = certify_spatial(cell.degrees, &cell.controls, Default::default());
                assert_eq!(
                    cert.status,
                    CertificateStatus::Valid,
                    "{} {:?}",
                    cell.id,
                    cert
                );
                lower = lower.min(cert.determinant_lower_bound_mm3.unwrap());
            }
            let inc = enumerate_faces(&cells).unwrap();
            // Four meridian bands each expose four spans on each angular quadrant.
            assert_eq!(inc.exterior.len(), 64);
            assert_eq!(inc.interfaces.len(), 352);
            println!("10 charts, 128 spans, minimum determinant bound {lower:.12}");
        }
    }
}
