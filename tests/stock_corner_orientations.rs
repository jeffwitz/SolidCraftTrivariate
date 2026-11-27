//! A consumer of the public kernel API, independent of the CAD host.
use solidcraft_trivariate::stock_corner::{
    analyze, BodyMetrics, Edge, Face, Input, SketchSegment, Treatment,
};

#[test]
fn setbacks_lie_on_the_two_stock_walls_for_all_twelve_edges() {
    let lo = [-5.0, 3.0, 11.0];
    let hi = [7.0, 19.0, 31.0];
    let d = 2.0;
    for axis in 0..3 {
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        for su in [-1.0, 1.0] {
            for sv in [-1.0, 1.0] {
                let mut base = lo;
                base[u] = if su > 0.0 { hi[u] } else { lo[u] };
                base[v] = if sv > 0.0 { hi[v] } else { lo[v] };
                let mut end = base;
                end[axis] = hi[axis];
                let face = |j: usize, sign: f64| {
                    let mut normal = [0.0; 3];
                    normal[j] = sign;
                    Face {
                        normal,
                        point: base,
                    }
                };
                let mut input = Input {
                    selected: Edge {
                        start: end,
                        end: base,
                        length: hi[axis] - lo[axis],
                        mid: None,
                    },
                    body: BodyMetrics {
                        bounds: (lo, hi),
                        closed: true,
                        shell_count: 1,
                        volume: 3840.0,
                    },
                    adjacent_faces: [face(u, su), face(v, sv)],
                    support_trims: vec![],
                };
                for treatment in [Treatment::Radius(d), Treatment::EqualDistanceChamfer(d)] {
                    let plan = analyze(&input, treatment).unwrap();
                    let sweep_sign = plan.x_dir[u] * plan.y_dir[v] - plan.x_dir[v] * plan.y_dir[u];
                    assert_eq!(sweep_sign, 1.0);
                    assert!(plan.origin[axis] < lo[axis]);
                    assert!(plan.origin[axis] + sweep_sign * plan.tool_length > hi[axis]);
                    let mut contacts = Vec::new();
                    for segment in &plan.profile_segments {
                        let local = match segment {
                            SketchSegment::Line { start, .. }
                            | SketchSegment::Arc { start, .. } => start,
                        };
                        let depth = |j: usize, sign: f64| {
                            -sign
                                * (plan.origin[j] - base[j]
                                    + local[0] * plan.x_dir[j]
                                    + local[1] * plan.y_dir[j])
                        };
                        let pair = (depth(u, su), depth(v, sv));
                        if pair.0 >= 0.0 && pair.1 >= 0.0 {
                            contacts.push(pair);
                        }
                    }
                    assert_eq!(contacts.len(), 2);
                    assert!(contacts.contains(&(0.0, d)) && contacts.contains(&(d, 0.0)));
                    input.adjacent_faces.swap(0, 1);
                    assert_eq!(analyze(&input, treatment).unwrap(), plan);
                }
            }
        }
    }
}
