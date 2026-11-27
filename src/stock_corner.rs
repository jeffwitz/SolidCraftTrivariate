//! Geometric eligibility and verification for a local treatment of a stock corner.
//!
//! The module is independent of any CAD host. It recognizes a straight, full-span edge at a
//! convex corner between perpendicular planar supports, derives a right-handed local
//! frame, detects nearby trims on the two supporting planes, and bounds the removed volume.

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BodyMetrics {
    pub volume: f64,
    pub closed: bool,
    pub shell_count: usize,
    pub bounds: ([f64; 3], [f64; 3]),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    pub length: f64,
    pub start: [f64; 3],
    pub end: [f64; 3],
    /// Optional midpoint witness, used for support-plane trim checks.
    pub mid: Option<[f64; 3]>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Face {
    /// Unit outward normal.
    pub normal: [f64; 3],
    pub point: [f64; 3],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Treatment {
    Radius(f64),
    EqualDistanceChamfer(f64),
}

impl Treatment {
    fn amount(self) -> f64 {
        match self {
            Self::Radius(x) | Self::EqualDistanceChamfer(x) => x,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Input {
    pub body: BodyMetrics,
    pub selected: Edge,
    pub adjacent_faces: [Face; 2],
    /// Other edges on either adjacent face. The selected edge itself may be included.
    pub support_trims: Vec<Edge>,
}

/// Follow tangent neighbors across split B-rep rims by default. The native host resolves
/// adjacent edges, preserving a continuous treatment around a hole or other split rim.
pub const DEFAULT_TANGENT_CHAIN: bool = true;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Profile {
    /// Remove the square corner outside a quarter-circle of this radius.
    QuarterCircle { radius: f64 },
    /// Remove the equal-leg triangular corner, producing a 45 degree face.
    EqualDistanceChamfer { distance: f64 },
}

/// One directed boundary segment in the local sketch plane. The ordered segments form a
/// closed wire; an arc's `start` must equal the previous segment's end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SketchSegment {
    Line {
        start: [f64; 2],
        end: [f64; 2],
    },
    Arc {
        center: [f64; 2],
        start: [f64; 2],
        end: [f64; 2],
        clockwise: bool,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeEnvelope {
    /// The nominal removed volume for the ideal swept profile.
    pub nominal: f64,
    /// Permitted interval, including the established relative numerical allowance.
    pub minimum: f64,
    pub maximum: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct StockCornerPlan {
    pub origin: [f64; 3],
    pub x_dir: [f64; 3],
    pub y_dir: [f64; 3],
    pub edge_length: f64,
    pub tool_length: f64,
    pub profile: Profile,
    /// Complete, ordered, closed profile in local sketch coordinates, including tool padding.
    pub profile_segments: Vec<SketchSegment>,
    pub nearby_support_trim: bool,
    pub volume: VolumeEnvelope,
    pub original_bounds: ([f64; 3], [f64; 3]),
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}
fn finite_point(p: [f64; 3]) -> bool {
    p.into_iter().all(f64::is_finite)
}
fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let den = dot(ab, ab);
    let t = if den > 0.0 {
        (dot(ap, ab) / den).clamp(0.0, 1.0)
    } else {
        0.0
    };
    distance(p, [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]])
}

fn full_box_perimeter(t: Edge, bounds: ([f64; 3], [f64; 3]), tol: f64) -> bool {
    let chord = distance(t.start, t.end);
    if (t.length - chord).abs() > tol || chord <= tol {
        return false;
    }
    let v = [
        t.end[0] - t.start[0],
        t.end[1] - t.start[1],
        t.end[2] - t.start[2],
    ];
    let Some(axis) = (0..3).find(|&i| v[i].abs() / chord > 1.0 - 1e-9) else {
        return false;
    };
    if (0..3).filter(|&i| i != axis).any(|i| v[i].abs() > tol) {
        return false;
    }
    let (lo, hi) = bounds;
    (t.start[axis].min(t.end[axis]) - lo[axis]).abs() <= tol
        && (t.start[axis].max(t.end[axis]) - hi[axis]).abs() <= tol
        && (0..3).filter(|&i| i != axis).all(|i| {
            (t.start[i] - t.end[i]).abs() <= tol
                && ((t.start[i] - lo[i]).abs() <= tol || (t.start[i] - hi[i]).abs() <= tol)
        })
}

fn nearby_trim(input: &Input, tol: f64, amount: f64) -> bool {
    let e = input.selected;
    input.support_trims.iter().any(|t| {
        if !t.length.is_finite()
            || t.length < 0.0
            || !finite_point(t.start)
            || !finite_point(t.end)
            || full_box_perimeter(*t, input.body.bounds, tol)
        {
            return false;
        }
        let samples = [Some(t.start), Some(t.end), t.mid]
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        input.adjacent_faces.iter().any(|f| {
            samples.iter().all(|p| {
                dot(
                    [p[0] - f.point[0], p[1] - f.point[1], p[2] - f.point[2]],
                    f.normal,
                )
                .abs()
                    <= tol
            }) && samples
                .iter()
                .any(|p| segment_distance(*p, e.start, e.end) - t.length <= amount + tol)
        })
    })
}

/// Analyze eligibility and return the local sweep frame and quantitative checks.
pub fn analyze(input: &Input, treatment: Treatment) -> Result<StockCornerPlan, &'static str> {
    let body = input.body;
    let amount = treatment.amount();
    if !amount.is_finite() || amount <= 0.0 {
        return Err("treatment amount must be positive and finite");
    }
    if !body.closed || body.shell_count != 1 || !body.volume.is_finite() || body.volume <= 0.0 {
        return Err("body must be a closed measured solid");
    }
    let (lo, hi) = body.bounds;
    if !finite_point(lo) || !finite_point(hi) {
        return Err("body bounds must be finite");
    }
    let size = std::array::from_fn::<_, 3, _>(|i| hi[i] - lo[i]);
    if size.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err("body bounds must have positive extents");
    }
    let scale = size.into_iter().fold(1.0_f64, f64::max);
    let tol = scale * 2.0e-8;
    let e = input.selected;
    if !e.length.is_finite() || e.length < 0.0 || !finite_point(e.start) || !finite_point(e.end) {
        return Err("selected edge geometry is invalid");
    }
    let chord_vec = [
        e.end[0] - e.start[0],
        e.end[1] - e.start[1],
        e.end[2] - e.start[2],
    ];
    let chord = distance(e.start, e.end);
    if !chord.is_finite() || chord <= tol {
        return Err("selected edge must be nonzero and straight");
    }
    let raw_dir = chord_vec.map(|x| x / chord);
    let axis = (0..3)
        .find(|&i| raw_dir[i].abs() > 1.0 - 1e-9)
        .filter(|&i| (0..3).filter(|&j| j != i).all(|j| raw_dir[j].abs() <= 1e-9))
        .ok_or("selected edge must be axis-aligned")?;
    // Canonical positive direction makes the frame independent of B-rep edge orientation.
    let mut edir = [0.0; 3];
    edir[axis] = 1.0;
    if (e.length - chord).abs() > tol || (chord - size[axis]).abs() > tol {
        return Err("selected edge must span the full stock length");
    }
    let start = if e.start[axis] < e.end[axis] {
        e.start
    } else {
        e.end
    };
    let end = if e.start[axis] < e.end[axis] {
        e.end
    } else {
        e.start
    };
    if (start[axis] - lo[axis]).abs() > tol || (end[axis] - hi[axis]).abs() > tol {
        return Err("selected edge must span complete stock bounds");
    }
    let trans = (0..3).filter(|&i| i != axis).collect::<Vec<_>>();
    for &i in &trans {
        if (start[i] - end[i]).abs() > tol
            || ((start[i] - lo[i]).abs() > tol && (start[i] - hi[i]).abs() > tol)
        {
            return Err("selected edge must be at an outer stock corner");
        }
    }
    if amount >= size[trans[0]].min(size[trans[1]]) - tol {
        return Err("treatment must fit within both stock dimensions");
    }
    let [f0, f1] = input.adjacent_faces;
    for f in [f0, f1] {
        let n = f.normal;
        if !finite_point(n) || !finite_point(f.point) {
            return Err("adjacent face data must be finite");
        }
        let norm = dot(n, n).sqrt();
        if (norm - 1.0).abs() > 1.0e-8 {
            return Err("adjacent face normals must be unit vectors");
        }
        let ax = (0..3)
            .find(|&i| n[i].abs() > 1. - 1e-9)
            .filter(|&i| (0..3).filter(|&j| j != i).all(|j| n[j].abs() <= 1e-9))
            .ok_or("adjacent faces must be axis-aligned")?;
        if ax == axis {
            return Err("adjacent faces must be transverse to the selected edge");
        }
    }
    let n0 = f0.normal;
    let n1 = f1.normal;
    if dot(n0, n1).abs() > 1.0e-9 {
        return Err("adjacent faces must be perpendicular exterior planes");
    }
    let ax0 = (0..3).find(|&i| n0[i].abs() > 1. - 1e-9).unwrap();
    let ax1 = (0..3).find(|&i| n1[i].abs() > 1. - 1e-9).unwrap();
    if ax0 == ax1 {
        return Err("adjacent faces must be distinct transverse supports");
    }
    for (f, ax) in [(f0, ax0), (f1, ax1)] {
        let bound = if f.normal[ax] > 0. { hi[ax] } else { lo[ax] };
        if (f.point[ax] - bound).abs() > tol || (start[ax] - bound).abs() > tol {
            return Err("selected edge is not on the exterior stock planes");
        }
    }
    let (x_dir, opposite) = if dot(cross(edir, n0), [-n1[0], -n1[1], -n1[2]]) > 1.0 - 1.0e-9 {
        (n0, n1)
    } else if dot(cross(edir, n1), [-n0[0], -n0[1], -n0[2]]) > 1.0 - 1.0e-9 {
        (n1, n0)
    } else {
        return Err("adjacent face orientation does not prove a convex outer corner");
    };
    let y_dir = cross(edir, x_dir);
    if dot(y_dir, [-opposite[0], -opposite[1], -opposite[2]]) <= 1.0 - 1.0e-9 {
        return Err("local frame orientation is not convex");
    }
    let inset = 1.0_f64.min(size[axis] / 4.0);
    let origin = [
        start[0] - edir[0] * inset,
        start[1] - edir[1] * inset,
        start[2] - edir[2] * inset,
    ];
    let profile = match treatment {
        Treatment::Radius(r) => Profile::QuarterCircle { radius: r },
        Treatment::EqualDistanceChamfer(d) => Profile::EqualDistanceChamfer { distance: d },
    };
    // Swept ideal section: square minus a quarter disk for a radius, or an equal-leg
    // right triangle for the 45-degree chamfer.
    let section = match treatment {
        Treatment::Radius(r) => r * r * (1.0 - std::f64::consts::FRAC_PI_4),
        Treatment::EqualDistanceChamfer(d) => 0.5 * d * d,
    };
    let nominal = section * chord;
    let allowance = body.volume * 2.0e-8;
    let nearby_support_trim = nearby_trim(input, tol, amount);
    let pad = inset;
    let profile_segments = match treatment {
        Treatment::Radius(r) => vec![
            SketchSegment::Line {
                start: [-r, 0.0],
                end: [-r, -pad],
            },
            SketchSegment::Line {
                start: [-r, -pad],
                end: [pad, -pad],
            },
            SketchSegment::Line {
                start: [pad, -pad],
                end: [pad, r],
            },
            SketchSegment::Line {
                start: [pad, r],
                end: [0.0, r],
            },
            SketchSegment::Arc {
                center: [-r, r],
                start: [0.0, r],
                end: [-r, 0.0],
                clockwise: true,
            },
        ],
        Treatment::EqualDistanceChamfer(d) => vec![
            SketchSegment::Line {
                start: [-d, 0.0],
                end: [-d, -pad],
            },
            SketchSegment::Line {
                start: [-d, -pad],
                end: [pad, -pad],
            },
            SketchSegment::Line {
                start: [pad, -pad],
                end: [pad, d],
            },
            SketchSegment::Line {
                start: [pad, d],
                end: [0.0, d],
            },
            SketchSegment::Line {
                start: [0.0, d],
                end: [-d, 0.0],
            },
        ],
    };
    Ok(StockCornerPlan {
        origin,
        x_dir,
        y_dir,
        edge_length: chord,
        tool_length: chord + 2.0 * inset,
        profile,
        profile_segments,
        nearby_support_trim,
        volume: VolumeEnvelope {
            nominal,
            minimum: if nearby_support_trim {
                0.0
            } else {
                (nominal - allowance).max(0.0)
            },
            maximum: nominal + allowance,
        },
        original_bounds: body.bounds,
    })
}

/// Check one-solid closure, unchanged stock bounds, and a volume delta inside the plan envelope.
pub fn validate_result(
    before: BodyMetrics,
    after: BodyMetrics,
    plan: &StockCornerPlan,
) -> Result<(), &'static str> {
    if !after.closed || !after.volume.is_finite() || after.volume <= 0.0 {
        return Err("result must remain a closed solid");
    }
    if after.shell_count != 1 {
        return Err("result must have exactly one boundary shell");
    }
    if !finite_point(after.bounds.0) || !finite_point(after.bounds.1) {
        return Err("result bounds must be finite");
    }
    let (lo, hi) = plan.original_bounds;
    let scale = (0..3).map(|i| hi[i] - lo[i]).fold(1.0_f64, f64::max);
    let tol = scale * 2.0e-8;
    if (0..3)
        .any(|i| (after.bounds.0[i] - lo[i]).abs() > tol || (after.bounds.1[i] - hi[i]).abs() > tol)
    {
        return Err("result bounds changed");
    }
    let removed = before.volume - after.volume;
    if !removed.is_finite()
        || removed <= 0.0
        || removed < plan.volume.minimum
        || removed > plan.volume.maximum
    {
        return Err("removed volume is outside the ideal profile envelope");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture(trim: bool) -> Input {
        let selected = Edge {
            length: 10.0,
            start: [0., 0., 0.],
            end: [0., 0., 10.],
            mid: Some([0., 0., 5.]),
        };
        Input {
            body: BodyMetrics {
                volume: 1000.,
                closed: true,
                shell_count: 1,
                bounds: ([0., 0., 0.], [10., 10., 10.]),
            },
            selected,
            adjacent_faces: [
                Face {
                    normal: [-1., 0., 0.],
                    point: [0., 0., 0.],
                },
                Face {
                    normal: [0., -1., 0.],
                    point: [0., 0., 0.],
                },
            ],
            support_trims: if trim {
                vec![Edge {
                    length: 3.,
                    start: [0., 2., 4.],
                    end: [0., 5., 4.],
                    mid: Some([0., 3.5, 4.]),
                }]
            } else {
                vec![]
            },
        }
    }
    #[test]
    fn radius_and_chamfer_share_frame_but_have_distinct_profiles_and_volumes() {
        let i = fixture(false);
        let r = analyze(&i, Treatment::Radius(2.)).unwrap();
        let c = analyze(&i, Treatment::EqualDistanceChamfer(2.)).unwrap();
        assert_eq!(r.edge_length, 10.);
        assert!((r.volume.nominal - 40. * (1. - std::f64::consts::FRAC_PI_4)).abs() < 1e-14);
        assert_eq!(c.volume.nominal, 20.);
        assert_eq!(r.origin, [0., 0., -1.]);
        assert_eq!(r.x_dir, [0., -1., 0.]);
        assert_eq!(r.y_dir, [1., 0., 0.]);
        assert_eq!(r.profile_segments.len(), 5);
        assert_eq!(c.profile_segments.len(), 5);
        for segments in [&r.profile_segments[..], &c.profile_segments[..]] {
            for pair in segments.windows(2) {
                let end = match pair[0] {
                    SketchSegment::Line { end, .. } | SketchSegment::Arc { end, .. } => end,
                };
                let start = match pair[1] {
                    SketchSegment::Line { start, .. } | SketchSegment::Arc { start, .. } => start,
                };
                assert_eq!(end, start);
            }
            let first = match segments[0] {
                SketchSegment::Line { start, .. } | SketchSegment::Arc { start, .. } => start,
            };
            let last = match segments[segments.len() - 1] {
                SketchSegment::Line { end, .. } | SketchSegment::Arc { end, .. } => end,
            };
            assert_eq!(first, last);
        }
        assert!(matches!(r.profile, Profile::QuarterCircle { radius: 2. }));
        assert!(matches!(
            c.profile,
            Profile::EqualDistanceChamfer { distance: 2. }
        ));
    }
    #[test]
    fn reversed_edge_parameterization_keeps_canonical_origin_and_frame() {
        let mut i = fixture(false);
        i.selected.start = [0., 0., 10.];
        i.selected.end = [0., 0., 0.];
        let p = analyze(&i, Treatment::Radius(1.)).unwrap();
        assert_eq!(p.origin, [0., 0., -1.]);
        assert_eq!(p.x_dir, [0., -1., 0.]);
        assert_eq!(p.y_dir, [1., 0., 0.]);
        assert!(p.origin[2] < 0. && p.origin[2] + p.tool_length > 10.);
    }
    #[test]
    fn detects_trim_on_support_plane_and_rejects_off_support_edge() {
        assert!(
            analyze(&fixture(true), Treatment::Radius(1.))
                .unwrap()
                .nearby_support_trim
        );
        let mut i = fixture(true);
        i.support_trims[0].start[0] = 0.1;
        i.support_trims[0].end[0] = 0.1;
        i.support_trims[0].mid = Some([0.1, 3.5, 4.]);
        assert!(
            !analyze(&i, Treatment::Radius(1.))
                .unwrap()
                .nearby_support_trim
        );
    }
    #[test]
    fn rejects_partial_nonconvex_and_oversized_cases() {
        let mut i = fixture(false);
        i.selected.end[2] = 8.;
        i.selected.length = 8.;
        assert!(analyze(&i, Treatment::Radius(1.)).is_err());
        let mut i = fixture(false);
        i.adjacent_faces[1].normal = [0., 1., 0.];
        i.adjacent_faces[1].point = [0., 10., 0.];
        assert!(analyze(&i, Treatment::Radius(1.)).is_err());
        assert!(analyze(&fixture(false), Treatment::EqualDistanceChamfer(10.)).is_err());
    }
    #[test]
    fn rejects_oblique_edge_and_excludes_stock_perimeter_trim() {
        let mut i = fixture(false);
        i.selected = Edge {
            length: 10.,
            start: [0., 0., 0.],
            end: [6., 8., 0.],
            mid: Some([3., 4., 0.]),
        };
        i.adjacent_faces = [
            Face {
                normal: [0.8, -0.6, 0.],
                point: [0., 0., 0.],
            },
            Face {
                normal: [0., 0., 1.],
                point: [0., 0., 0.],
            },
        ];
        assert!(analyze(&i, Treatment::Radius(0.25)).is_err());
        i = fixture(false);
        i.support_trims.push(Edge {
            length: 10.,
            start: [0., 0., 0.],
            end: [0., 10., 0.],
            mid: Some([0., 5., 0.]),
        });
        assert!(
            !analyze(&i, Treatment::Radius(0.25))
                .unwrap()
                .nearby_support_trim
        );
    }
    #[test]
    fn result_checks_bounds_closure_and_volume() {
        let trimmed = fixture(true);
        let plan = analyze(&trimmed, Treatment::EqualDistanceChamfer(2.)).unwrap();
        assert!(validate_result(trimmed.body, trimmed.body, &plan).is_err());
        let input = fixture(false);
        let plan = analyze(&input, Treatment::EqualDistanceChamfer(2.)).unwrap();
        let after = BodyMetrics {
            volume: 980.,
            closed: true,
            shell_count: 1,
            bounds: input.body.bounds,
        };
        assert!(validate_result(input.body, after, &plan).is_ok());
        let mut bad = after;
        bad.closed = false;
        assert!(validate_result(input.body, bad, &plan).is_err());
        bad = after;
        bad.shell_count = 2;
        assert!(validate_result(input.body, bad, &plan).is_err());
        bad = after;
        bad.bounds.1[0] = 9.99;
        assert!(validate_result(input.body, bad, &plan).is_err());
        bad = after;
        bad.volume = 970.;
        assert!(validate_result(input.body, bad, &plan).is_err());
    }

    #[test]
    fn tangent_chain_policy_defaults_on() {
        assert!(DEFAULT_TANGENT_CHAIN);
    }
}
