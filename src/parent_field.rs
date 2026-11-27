//! Canonical rational NURBS displacement field on a recovered tensor support.
//! Cox recursion is restricted to a single knot interval in Bernstein form.
//! E maps polynomial coefficients, C maps physical rational coefficients.
use crate::{
    bracket_json::decode_blocks,
    constraint_system::{Lift, LiftLine, Ratio as Q},
    extracted_space::{AxisExtraction, ExtractedElement},
    incidence::enumerate_faces,
    json::{object, Json},
    section_manufactured::{solve_sampled, IntegrationDomain},
    spline_constraints::{add, Row},
    spline_support::recover_complete_parent,
};
use std::collections::{BTreeMap, BTreeSet};
type Result<T> = std::result::Result<T, &'static str>;

/// Rows index Bernstein coefficients on the requested subinterval. Columns
/// index functions in the original knot vector. A subinterval may be shorter
/// than a knot span. Crossing a nonempty knot interval is refused.
pub(crate) fn axis_rows(knots: &[Q], degree: usize, interval: [Q; 2]) -> Result<Vec<Row>> {
    if !(1..=3).contains(&degree)
        || knots.len() < 2 * (degree + 1)
        || knots
            .windows(2)
            .any(|pair| pair[0].cmp_value(&pair[1]).is_gt())
        || !interval[0].cmp_value(&interval[1]).is_lt()
    {
        return Err("invalid knot interval for polynomial extraction");
    }
    let n = knots.len() - degree - 1;
    let span = (degree..n)
        .find(|&i| {
            !knots[i].cmp_value(&interval[0]).is_gt()
                && !interval[1].cmp_value(&knots[i + 1]).is_gt()
        })
        .ok_or("extraction interval crosses a knot or leaves the active domain")?;
    let mut previous = BTreeMap::from([(span, vec![Q::one()])]);
    for d in 1..=degree {
        let mut next = BTreeMap::new();
        for i in span - d..=span {
            let mut coeffs = vec![Q::zero(); d + 1];
            for (j, left) in [(i, true), (i + 1, false)] {
                let Some(source) = previous.get(&j) else {
                    continue;
                };
                let width = knots[j + d].checked_sub(knots[j])?;
                if width == Q::zero() {
                    continue;
                }
                let affine = interval.map(|t| {
                    if left {
                        t.checked_sub(knots[j])
                    } else {
                        knots[j + d].checked_sub(t)
                    }
                });
                let a = affine[0]?.checked_div(width)?;
                let b = affine[1]?.checked_div(width)?;
                for (k, &value) in source.iter().enumerate() {
                    let low = value
                        .checked_mul(a)?
                        .checked_mul(Q::new((d - k) as i128, d as i128)?)?;
                    let high = value
                        .checked_mul(b)?
                        .checked_mul(Q::new((k + 1) as i128, d as i128)?)?;
                    coeffs[k] = coeffs[k].checked_add(low)?;
                    coeffs[k + 1] = coeffs[k + 1].checked_add(high)?;
                }
            }
            next.insert(i, coeffs);
        }
        previous = next;
    }
    let mut rows = vec![Row::new(); degree + 1];
    for (i, values) in previous {
        for (k, v) in values.into_iter().enumerate() {
            add(&mut rows[k], i, v)?;
        }
    }
    for row in &rows {
        let sum = row.values().try_fold(Q::zero(), |s, &v| s.checked_add(v))?;
        if sum != Q::one() || row.values().any(|v| v.numerator() < 0) {
            return Err("invalid polynomial partition of unity");
        }
    }
    Ok(rows)
}
fn numeric_axis(rows: &[Row]) -> (Vec<usize>, AxisExtraction) {
    let active: Vec<usize> = rows
        .iter()
        .flat_map(|r| r.keys().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut columns = Vec::new();
    let mut coefficients = Vec::new();
    let mut row_offsets = vec![0];
    for r in rows {
        for (&i, &v) in r {
            columns.push(active.binary_search(&i).unwrap());
            coefficients.push(v);
        }
        row_offsets.push(columns.len());
    }
    let axis = AxisExtraction {
        degree: rows.len() - 1,
        basis_count: active.len(),
        columns,
        coefficients,
        row_offsets,
    };
    (active, axis)
}
fn row_json(row: &Row) -> Json {
    Json::Array(
        row.iter()
            .map(|(i, c)| {
                Json::Array(vec![
                    Json::Number(i.to_string()),
                    Json::Array(vec![
                        Json::String(c.numerator().to_string()),
                        Json::String(c.denominator().to_string()),
                    ]),
                ])
            })
            .collect(),
    )
}

pub struct ParentField {
    pub elements: Vec<ExtractedElement>,
    pub ndof: usize,
    pub extraction: Json,
    cell_rows: BTreeMap<String, (usize, Vec<Row>)>,
    boundary: Vec<(usize, usize, bool)>,
    clamp: Lift,
    z_planes: Option<[f64; 2]>,
}
impl ParentField {
    pub fn solve(&self, job: &Json) -> Result<Json> {
        let domain = IntegrationDomain {
            scalar_dofs: self.ndof,
            cell_count: self.elements.len(),
            z_bounds: self
                .z_planes
                .ok_or("unsupported clamp: no physical z end planes")?,
            boundary: &self.boundary,
            essential: &self.clamp,
        };
        solve_sampled(&domain, job, |i, p| self.elements[i].sample(p))
    }
}
pub fn build(model: &Json) -> Result<ParentField> {
    let support = recover_complete_parent(model)?;
    let ndof = support.homogeneous_controls.len();
    let [_, nv, nw] = support.control_shape;
    let blocks = decode_blocks(model)?;
    let faces = enumerate_faces(&blocks)?;
    let blocks: BTreeMap<_, _> = blocks.iter().map(|b| (b.id.as_str(), b)).collect();
    let mut elements = Vec::new();
    let mut records = Vec::new();
    let mut cell_rows = BTreeMap::new();
    let mut indices = BTreeMap::new();
    for span in &support.source_spans {
        let r = &span.restriction;
        if r.axes != [0, 1, 2] {
            return Err("unsupported parent field span orientation");
        }
        let mut factors = Vec::new();
        for a in 0..3 {
            let knot = support.knots[a]
                .iter()
                .map(|&x| Q::new(i128::from(x), 1))
                .collect::<Result<Vec<_>>>()?;
            let [lo, hi] = r.intervals[a];
            let mut restriction = axis_rows(
                &knot,
                support.degrees[a],
                [Q::new(i128::from(lo), 1)?, Q::new(i128::from(hi), 1)?],
            )?;
            if r.reversals[a] {
                restriction.reverse();
            }
            factors.push(restriction);
        }
        let axes = [&factors[0], &factors[1], &factors[2]];
        let [(u, au), (v, av), (w, aw)] = axes.map(|rs| numeric_axis(rs));
        let mut dofs = Vec::new();
        let mut controls = Vec::new();
        for &i in &u {
            for &j in &v {
                for &k in &w {
                    let id = (i * nv + j) * nw + k;
                    dofs.push(id);
                    controls.push(support.homogeneous_controls[id].map(|x| x as f64));
                }
            }
        }
        let e = ExtractedElement {
            axes: [au, av, aw],
            dofs,
            field_weights: controls.iter().map(|c| c[3]).collect(),
            geometry_controls: controls,
        };
        e.validate()?;
        let block = blocks[span.source_block_id.as_str()];
        let gauge = Q::new(span.gauge[0], span.gauge[1])?;
        let mut polynomial = Vec::new();
        let mut rational = Vec::new();
        let mut physical_rows = Vec::new();
        for a in axes[0] {
            for b in axes[1] {
                for c in axes[2] {
                    let mut row = Row::new();
                    let mut control = [Q::zero(); 4];
                    for (&i, &lu) in a {
                        for (&j, &lv) in b {
                            for (&k, &lw) in c {
                                let id = (i * nv + j) * nw + k;
                                let basis = lu.checked_mul(lv)?.checked_mul(lw)?;
                                row.insert(id, basis);
                                for (dst, h) in
                                    control.iter_mut().zip(support.homogeneous_controls[id])
                                {
                                    *dst = dst.checked_add(basis.checked_mul(Q::new(h, 1)?)?)?;
                                }
                            }
                        }
                    }
                    let expected = block.controls[polynomial.len()];
                    for (&h, x) in control.iter().zip(expected) {
                        if h.checked_mul(gauge)? != Q::new(i128::from(x), 1)? {
                            return Err("canonical field extraction disagrees with the source net");
                        }
                    }
                    let physical = row
                        .iter()
                        .map(|(&id, &v)| {
                            Ok((
                                id,
                                v.checked_mul(Q::new(support.homogeneous_controls[id][3], 1)?)?
                                    .checked_div(control[3])?,
                            ))
                        })
                        .collect::<Result<Row>>()?;
                    polynomial.push(row_json(&row));
                    rational.push(row_json(&physical));
                    physical_rows.push(physical);
                }
            }
        }
        indices.insert(span.source_block_id.as_str(), elements.len());
        records.push(object([
            ("source_block", Json::String(span.source_block_id.clone())),
            ("E", Json::Array(polynomial)),
            ("C", Json::Array(rational)),
        ]));
        cell_rows.insert(
            span.source_block_id.clone(),
            (elements.len(), physical_rows),
        );
        elements.push(e);
    }
    let free_roots: Vec<_> = (0..ndof).filter(|i| i % nw != 0).collect();
    let lines = (0..ndof)
        .map(|i| LiftLine {
            dof: i,
            terms: free_roots
                .binary_search(&i)
                .map(|j| vec![(j, Q::one())])
                .unwrap_or_default(),
            offset: Q::zero(),
        })
        .collect();
    let boundary = faces
        .exterior
        .iter()
        .map(|s| (indices[s.block.as_str()], s.axis, s.upper))
        .collect();
    let h = &support.homogeneous_controls;
    let planes = [0, nw - 1]
        .map(|k| Q::new(h[k][2], h[k][3]))
        .into_iter()
        .collect::<Result<Vec<_>>>()?;
    let regular = h.iter().enumerate().all(|(i, p)| {
        i % nw != 0 && i % nw != nw - 1
            || Q::new(p[2], p[3]) == Ok(planes[usize::from(i % nw != 0)])
    });
    let z_planes = (regular && planes[0].cmp_value(&planes[1]).is_lt()).then(|| {
        [
            crate::inherited_space::val(planes[0]),
            crate::inherited_space::val(planes[1]),
        ]
    });
    let extraction = object([
        ("schema", Json::String("canonical_parent_field_v1".into())),
        ("support", Json::String(support.id)),
        ("scalar_dofs", Json::Number(ndof.to_string())),
        ("elements", Json::Array(records)),
    ]);
    Ok(ParentField {
        elements,
        ndof,
        extraction,
        cell_rows,
        boundary,
        clamp: Lift { free_roots, lines },
        z_planes,
    })
}

/// Compare the canonical NURBS space to its exact embedding in the full-face
/// C0 coefficients. Both integrals use identical source geometry evaluations.
pub fn check_galerkin(model: &Json, order: usize) -> Result<Json> {
    use crate::extracted_space::{assemble_sampled, compare_assembled};
    let field = build(model)?;
    let c0 = crate::inherited_space::build_axial_space(model, 1, 0)?;
    let mut lift = vec![Row::new(); c0.ndof];
    let mut element_order = Vec::new();
    for cell in &c0.cells {
        let (element, rows) = field
            .cell_rows
            .get(&cell.id)
            .ok_or("source field lineage differs")?;
        element_order.push(*element);
        if rows.len() != cell.extraction.len() {
            return Err("incompatible C0 degree");
        }
        for (c, r) in cell.extraction.iter().zip(rows) {
            if c.len() != 1 {
                return Err("source C0 row is not nodal");
            }
            let (&i, &v) = c.first_key_value().unwrap();
            if v != Q::one() || (!lift[i].is_empty() && lift[i] != *r) {
                return Err("parent transfer changed on shared coefficient");
            }
            lift[i] = r.clone();
        }
    }
    let (_, kernel) = crate::spline_constraints::nullspace(field.ndof, &lift)?;
    if !kernel.is_empty() {
        return Err("parent coefficient transfer lacks full column rank");
    }
    let body = |x: [f64; 3]| [1. + x[0], x[2] - x[1], 2. * x[1] + x[0]];
    let large = assemble_sampled(
        c0.ndof,
        c0.cells.len(),
        1.5,
        1.,
        order,
        |i, p| c0.cells[i].sample(p),
        body,
    )?;
    let small = assemble_sampled(
        field.ndof,
        c0.cells.len(),
        1.5,
        1.,
        order,
        |i, p| {
            let (common, _) = c0.cells[i].sample(p)?;
            let (canonical, basis) = field.elements[element_order[i]].sample(p)?;
            if common
                .position
                .iter()
                .zip(canonical.position)
                .any(|(a, b)| (a - b).abs() > 1e-10 * (1. + a.abs()))
                || common
                    .jacobian_columns
                    .iter()
                    .flatten()
                    .zip(canonical.jacobian_columns.iter().flatten())
                    .any(|(a, b)| (a - b).abs() > 1e-10 * (1. + a.abs()))
            {
                return Err("parent and source quadrature geometry differ");
            }
            Ok((common, basis))
        },
        body,
    )?;
    let result = compare_assembled(&large, &small, c0.ndof, field.ndof, &lift)?;
    if result.stiffness_relative_frobenius_defect > 1e-10 || result.force_relative_l2_defect > 1e-10
    {
        return Err("parent Galerkin identity exceeds tolerance");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("c0_scalar_dofs", Json::Number(c0.ndof.to_string())),
        (
            "canonical_scalar_dofs",
            Json::Number(field.ndof.to_string()),
        ),
        ("exact_transfer_rank", Json::Number(field.ndof.to_string())),
        (
            "relative_stiffness_defect",
            Json::Number(result.stiffness_relative_frobenius_defect.to_string()),
        ),
        (
            "relative_load_defect",
            Json::Number(result.force_relative_l2_defect.to_string()),
        ),
    ]))
}

/// Small-system oracle interchange. The native assembly remains block CSR.
pub fn stiffness_package(model: &Json, order: usize) -> Result<Json> {
    use crate::extracted_space::assemble_elements;
    let parent = recover_complete_parent(model)?;
    let field = build(model)?;
    if field.ndof > 256 {
        return Err("indeterminate: small-system spectrum audit budget");
    }
    if field.z_planes.is_none() {
        return Err("unsupported spectrum clamp orientation");
    }
    let material = assemble_elements(&field.elements, field.ndof, 2.6, 0.3, order)?;
    let m = material.stiffness;
    let n = |v: f64| Json::Number(v.to_string());
    let array =
        |v: Vec<usize>| Json::Array(v.into_iter().map(|i| Json::Number(i.to_string())).collect());
    let values = m
        .values
        .into_iter()
        .map(|block| {
            Json::Array(
                block
                    .into_iter()
                    .map(|r| Json::Array(r.map(n).to_vec()))
                    .collect(),
            )
        })
        .collect();
    let points = parent
        .homogeneous_controls
        .into_iter()
        .map(|h| Json::Array((0..3).map(|j| n(h[j] as f64 / h[3] as f64)).collect()))
        .collect();
    Ok(object([
        ("schema", Json::String("canonical_elasticity_csr_v1".into())),
        ("offsets", array(m.offsets)),
        ("columns", array(m.columns)),
        ("values", Json::Array(values)),
        ("control_positions_mm", Json::Array(points)),
        ("free_scalar_roots", array(field.clamp.free_roots)),
        ("volume_mm3", n(m.integrated_volume)),
        ("young_mpa", n(2.6)),
        ("poisson", n(0.3)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn q(n: i128, d: i128) -> Q {
        Q::new(n, d).unwrap()
    }
    #[test]
    fn subinterval_restriction_uses_actual_source_knots() {
        let knot = [q(-2, 1), q(-2, 1), q(2, 1), q(5, 1), q(5, 1)];
        let map = axis_rows(&knot, 1, [q(-1, 1), q(0, 1)]).unwrap();
        assert_eq!(
            map,
            vec![
                Row::from([(0, q(3, 4)), (1, q(1, 4))]),
                Row::from([(0, q(1, 2)), (1, q(1, 2))])
            ]
        );
        assert!(axis_rows(&knot, 1, [q(-1, 1), q(3, 1)]).is_err());
        assert!(axis_rows(&knot, 1, [q(0, 1), q(0, 1)]).is_err());
    }
    #[test]
    fn canonical_parent_stiffness_and_load_pull_back_from_c0() {
        let recipe =
            crate::json::parse(include_str!("../tests/references/nurbs_strip.recipe.json"))
                .unwrap();
        let model = crate::nurbs_sketch_json::build(&recipe).unwrap();
        let report = check_galerkin(&model, 5).unwrap();
        assert_eq!(
            report
                .field("exact_transfer_rank")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            20
        );
        assert_eq!(
            report
                .field("c0_scalar_dofs")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            28
        );
    }
}
