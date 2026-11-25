//! Small three-cell Q2 mesh and exact trace-relation reconstruction utilities.

use solidcraft_trivariate::analysis_space::{AnalysisBlock, AnalysisSpace};
use solidcraft_trivariate::face_constraints::{nested_constraints, FaceNet, Row, Span};
use solidcraft_trivariate::incidence::Orientation;
use std::collections::BTreeMap;

type Rational = (i128, i128);
type FaceBounds = [(Rational, Rational); 2];
type Overlap = (FaceBounds, [Span; 2]);

#[derive(Clone)]
pub struct Relation {
    pub master_id: String,
    pub slave_id: String,
    pub master_side: (usize, bool),
    pub slave_side: (usize, bool),
    pub spans: [Span; 2],
}

#[derive(Clone)]
pub struct TMesh {
    pub space: AnalysisSpace,
    pub relations: Vec<Relation>,
    pub equations: Vec<(usize, Row)>,
}

#[derive(Clone, Copy)]
struct Face {
    block: usize,
    axis: usize,
    upper: bool,
    plane: (i128, i128),
    bounds: [((i128, i128), (i128, i128)); 2],
    degrees: [usize; 2],
    local_dofs: [usize; 9],
    controls: [[i64; 4]; 9],
    area: f64,
}

fn ratio(num: i128, den: i128) -> (i128, i128) {
    let mut a = num.unsigned_abs();
    let mut b = den.unsigned_abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    let d = a.max(1) as i128;
    let sign = if den < 0 { -1 } else { 1 };
    (num / d * sign, den / d * sign)
}

fn face(block_index: usize, block: &AnalysisBlock, axis: usize, upper: bool) -> Face {
    let free = match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let fixed = if upper { 2 } else { 0 };
    let mut controls = [[0_i64; 4]; 9];
    let mut local_dofs = [0; 9];
    for a in 0..3 {
        for b in 0..3 {
            let mut ijk = [0; 3];
            ijk[axis] = fixed;
            ijk[free[0]] = a;
            ijk[free[1]] = b;
            let local = (ijk[0] * 3 + ijk[1]) * 3 + ijk[2];
            let out = a * 3 + b;
            controls[out] =
                block.controls[local].map(|x| i64::try_from(x).expect("small fixture net"));
            local_dofs[out] = block.dofs[local];
        }
    }
    let plane_axis = axis;
    let plane = ratio(
        i128::from(controls[0][plane_axis]),
        i128::from(controls[0][3]),
    );
    let bounds = std::array::from_fn(|tangent| {
        let coordinate = free[tangent];
        let values = controls
            .iter()
            .map(|p| ratio(i128::from(p[coordinate]), i128::from(p[3])));
        let (mut low, mut high) = (None, None);
        for value in values {
            low = Some(low.map_or(value, |x: (i128, i128)| {
                if value.0 * x.1 < x.0 * value.1 {
                    value
                } else {
                    x
                }
            }));
            high = Some(high.map_or(value, |x: (i128, i128)| {
                if value.0 * x.1 > x.0 * value.1 {
                    value
                } else {
                    x
                }
            }));
        }
        (low.unwrap(), high.unwrap())
    });
    let length = |lo: (i128, i128), hi: (i128, i128)| {
        (hi.0 as f64 / hi.1 as f64) - (lo.0 as f64 / lo.1 as f64)
    };
    Face {
        block: block_index,
        axis,
        upper,
        plane,
        bounds,
        degrees: [2, 2],
        local_dofs,
        controls,
        area: length(bounds[0].0, bounds[0].1) * length(bounds[1].0, bounds[1].1),
    }
}

fn compare(a: (i128, i128), b: (i128, i128)) -> std::cmp::Ordering {
    (a.0 * b.1).cmp(&(b.0 * a.1))
}

fn overlap(a: &Face, b: &Face) -> Option<Overlap> {
    if a.axis != b.axis || a.upper == b.upper || a.plane != b.plane {
        return None;
    }
    let mut bounds = [((0, 1), (0, 1)); 2];
    let mut spans = [Span::Full; 2];
    for axis in 0..2 {
        let lo = if compare(a.bounds[axis].0, b.bounds[axis].0).is_gt() {
            a.bounds[axis].0
        } else {
            b.bounds[axis].0
        };
        let hi = if compare(a.bounds[axis].1, b.bounds[axis].1).is_lt() {
            a.bounds[axis].1
        } else {
            b.bounds[axis].1
        };
        if !compare(lo, hi).is_lt() {
            return None;
        }
        bounds[axis] = (lo, hi);
        if lo == a.bounds[axis].0 && hi == a.bounds[axis].1 {
            spans[axis] = Span::Full;
        } else if lo == a.bounds[axis].0 && compare(hi, midpoint(a.bounds[axis])).is_eq() {
            spans[axis] = Span::Lower;
        } else if compare(lo, midpoint(a.bounds[axis])).is_eq() && hi == a.bounds[axis].1 {
            spans[axis] = Span::Upper;
        } else {
            return None;
        }
    }
    Some((bounds, spans))
}

fn midpoint(bounds: ((i128, i128), (i128, i128))) -> (i128, i128) {
    let n = bounds.0 .0 * bounds.1 .1 + bounds.1 .0 * bounds.0 .1;
    let d = 2 * bounds.0 .1 * bounds.1 .1;
    ratio(n, d)
}

fn face_net(face: &Face) -> FaceNet {
    FaceNet {
        degrees: face.degrees,
        controls: face.controls.to_vec(),
    }
}

pub fn rebuild(space: AnalysisSpace) -> Result<TMesh, &'static str> {
    let all_faces: Vec<Face> = space
        .blocks
        .iter()
        .enumerate()
        .flat_map(|(i, block)| {
            (0..3).flat_map(move |axis| [false, true].map(move |upper| face(i, block, axis, upper)))
        })
        .collect();
    let mut relations = Vec::new();
    let mut equations = Vec::new();
    for i in 0..all_faces.len() {
        for j in i + 1..all_faces.len() {
            let a = all_faces[i];
            let b = all_faces[j];
            if a.block == b.block {
                continue;
            }
            let Some((_intersection, spans_a)) = overlap(&a, &b) else {
                continue;
            };
            // Select the containing coarse net as master. Equal faces use a
            // stable block ID order so the scalar constraint graph is acyclic.
            let (master, slave, spans) = if a.area > b.area
                || (a.area == b.area && space.blocks[a.block].id < space.blocks[b.block].id)
            {
                (a, b, spans_a)
            } else {
                let Some((_intersection, spans_b)) = overlap(&b, &a) else {
                    return Err("face overlap is not nested");
                };
                (b, a, spans_b)
            };
            let rows = nested_constraints(
                &face_net(&master),
                &face_net(&slave),
                spans,
                Orientation {
                    swap: false,
                    flip_first: false,
                    flip_second: false,
                },
            )?;
            for (row_index, row) in rows.into_iter().enumerate() {
                equations.push((
                    slave.local_dofs[row_index],
                    Row {
                        terms: row
                            .terms
                            .into_iter()
                            .map(|(m, weight)| (master.local_dofs[m], weight))
                            .collect(),
                        denominator: row.denominator,
                    },
                ));
            }
            relations.push(Relation {
                master_id: space.blocks[master.block].id.clone(),
                slave_id: space.blocks[slave.block].id.clone(),
                master_side: (master.axis, master.upper),
                slave_side: (slave.axis, slave.upper),
                spans,
            });
        }
    }
    relations.sort_by(|a, b| {
        (&a.master_id, &a.slave_id, a.master_side).cmp(&(&b.master_id, &b.slave_id, b.master_side))
    });
    equations.sort_by_key(|(dof, _)| *dof);
    Ok(TMesh {
        space,
        relations,
        equations,
    })
}

pub fn base_three_cell() -> Result<TMesh, &'static str> {
    let boxes = [
        ("master", [0_i64, 2], [0_i64, 4], [0_i64, 2]),
        ("neighbor_low", [2, 4], [0, 2], [0, 2]),
        ("neighbor_high", [2, 4], [2, 4], [0, 2]),
    ];
    let blocks: Vec<_> = boxes
        .iter()
        .enumerate()
        .map(|(bi, (id, x, y, z))| {
            let mut controls = Vec::with_capacity(27);
            let mut dofs = Vec::with_capacity(27);
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        let coord = [
                            x[0] + (x[1] - x[0]) * i as i64 / 2,
                            y[0] + (y[1] - y[0]) * j as i64 / 2,
                            z[0] + (z[1] - z[0]) * k as i64 / 2,
                        ];
                        controls.push([
                            i128::from(coord[0]),
                            i128::from(coord[1]),
                            i128::from(coord[2]),
                            2,
                        ]);
                        dofs.push(bi * 27 + (i * 3 + j) * 3 + k);
                    }
                }
            }
            AnalysisBlock {
                id: (*id).to_owned(),
                controls,
                dofs,
            }
        })
        .collect();
    rebuild(AnalysisSpace {
        ndof: 81,
        blocks,
        boundary_dofs: BTreeMap::new(),
    })
}

#[allow(dead_code)]
pub fn replace_refined(mesh: &TMesh, refine_ids: &[&str]) -> Result<TMesh, &'static str> {
    let mut blocks = Vec::new();
    let mut offset = 0;
    for block in &mesh.space.blocks {
        if refine_ids.contains(&block.id.as_str()) {
            let mut parent_block = block.clone();
            parent_block.dofs = (0..27).collect();
            let local = AnalysisSpace {
                ndof: 27,
                blocks: vec![parent_block],
                boundary_dofs: BTreeMap::new(),
            };
            let refined = solidcraft_trivariate::h_refine::refine_uniform(&local)?;
            for mut child in refined.blocks {
                for dof in &mut child.dofs {
                    *dof += offset;
                }
                blocks.push(child);
            }
            offset += refined.ndof;
        } else {
            let mut copy = block.clone();
            for (local, dof) in copy.dofs.iter_mut().enumerate() {
                *dof = offset + local;
            }
            blocks.push(copy);
            offset += 27;
        }
    }
    rebuild(AnalysisSpace {
        ndof: offset,
        blocks,
        boundary_dofs: BTreeMap::new(),
    })
}
