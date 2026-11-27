//! Native multi-support field on audited complete-face extrusions.
//!
//! Source cells remain integration cells. Jet equations are assembled in
//! physical extrusion length; their exact null space is the field extraction.
//! Each original homogeneous section is reused without fitting or knot removal.
use crate::{
    analysis_space::{face_axes, UnionFind},
    bezier::basis,
    bracket::BracketBlock,
    bracket_json::decode_blocks,
    constraint_system::Ratio as Q,
    extracted_space::{ElementBasisSample, ElementGeometrySample},
    incidence::{enumerate_faces, BlockSide, Interface},
    interchange::audit,
    json::{object, Json},
    spline_constraints::{add, nullspace, Result, Row},
    RationalBezierMap,
};
use std::collections::{BTreeMap, BTreeSet};

type V = [f64; 3];
fn index(d: [usize; 3], at: [usize; 3]) -> usize {
    (at[0] * (d[1] + 1) + at[1]) * (d[2] + 1) + at[2]
}
pub(crate) fn val(r: Q) -> f64 {
    r.numerator() as f64 / r.denominator() as f64
}
pub(crate) fn ratio_div(a: Q, b: Q) -> Result<Q> {
    a.checked_mul(Q::new(b.denominator(), b.numerator())?)
}
fn n(i: impl ToString) -> Json {
    Json::Number(i.to_string())
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: V, b: V) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[derive(Debug)]
pub struct InheritedCell {
    pub id: String,
    pub z_bounds: [f64; 2],
    pub degrees: [usize; 3],
    pub active: Vec<usize>,
    /// Local rational Bernstein rows indexed by global free scalar coefficient.
    pub extraction: Vec<Row>,
    geometry: RationalBezierMap,
    weights: Vec<f64>,
    stencils: Vec<Vec<(usize, f64)>>,
}

impl InheritedCell {
    pub fn sample(&self, uvw: V) -> Result<(ElementGeometrySample, ElementBasisSample)> {
        let g = self
            .geometry
            .evaluate(uvw)
            .map_err(|_| "geometry evaluation unavailable")?;
        let j = g.jacobian_columns;
        let cof = [cross(j[1], j[2]), cross(j[2], j[0]), cross(j[0], j[1])];
        let det = dot(j[0], cof[0]);
        if !det.is_finite() || det <= 0. {
            return Err("nonpositive sampled geometry Jacobian");
        }
        let b: [_; 3] = std::array::from_fn(|a| basis(self.degrees[a], uvw[a]));
        let mut values = vec![0.; self.active.len()];
        let mut grads = vec![[0.; 3]; self.active.len()];
        let mut w = 0.;
        let mut dw = [0.; 3];
        for i in 0..=self.degrees[0] {
            for j in 0..=self.degrees[1] {
                for k in 0..=self.degrees[2] {
                    let l = index(self.degrees, [i, j, k]);
                    let weight = self.weights[l];
                    let p = weight * b[0].0[i] * b[1].0[j] * b[2].0[k];
                    let dp = [
                        weight * b[0].1[i] * b[1].0[j] * b[2].0[k],
                        weight * b[0].0[i] * b[1].1[j] * b[2].0[k],
                        weight * b[0].0[i] * b[1].0[j] * b[2].1[k],
                    ];
                    w += p;
                    for a in 0..3 {
                        dw[a] += dp[a];
                    }
                    for &(d, c) in &self.stencils[l] {
                        values[d] += p * c;
                        for a in 0..3 {
                            grads[d][a] += dp[a] * c;
                        }
                    }
                }
            }
        }
        if w <= 0. || !w.is_finite() {
            return Err("invalid rational denominator");
        }
        for (r, d) in values.iter_mut().zip(&mut grads) {
            *r /= w;
            let ref_grad: [f64; 3] = std::array::from_fn(|a| (d[a] - *r * dw[a]) / w);
            *d = std::array::from_fn(|c| {
                (ref_grad[0] * cof[0][c] + ref_grad[1] * cof[1][c] + ref_grad[2] * cof[2][c]) / det
            });
        }
        Ok((
            ElementGeometrySample {
                position: g.position,
                jacobian_columns: g.jacobian_columns,
                determinant: det,
            },
            ElementBasisSample {
                dofs: self.active.clone(),
                values,
                physical_gradients: grads,
            },
        ))
    }
}

pub struct InheritedSpace {
    pub cells: Vec<InheritedCell>,
    pub ndof: usize,
    pub control_xyz: Vec<V>,
    pub c0_ndof: usize,
    pub common_parent_faces: usize,
    pub retained_faces: usize,
    pub transverse_faces: usize,
    pub interface_evidence: Vec<Json>,
    /// Deterministic audit trail for automatically selected regularity.
    pub selection_audit: Option<Json>,
    /// Full C0 coefficient rows, in the same canonical ordering for any continuity.
    pub lift: Vec<Row>,
}

fn exact_json(r: Q) -> Json {
    Json::Array(
        [r.numerator(), r.denominator()]
            .into_iter()
            .map(|v| Json::String(v.to_string()))
            .collect(),
    )
}
impl InheritedSpace {
    pub fn document(&self) -> Json {
        let cells = self
            .cells
            .iter()
            .map(|cell| {
                object([
                    ("id", Json::String(cell.id.clone())),
                    ("degrees", Json::Array(cell.degrees.map(n).to_vec())),
                    (
                        "rational_extraction",
                        Json::Array(
                            cell.extraction
                                .iter()
                                .map(|r| {
                                    Json::Array(
                                        r.iter()
                                            .map(|(&i, &c)| Json::Array(vec![n(i), exact_json(c)]))
                                            .collect(),
                                    )
                                })
                                .collect(),
                        ),
                    ),
                ])
            })
            .collect();
        let mut result = object([
            (
                "schema",
                Json::String("inherited_extrusion_field_v1".into()),
            ),
            ("scalar_dofs", n(self.ndof)),
            ("c0_scalar_dofs", n(self.c0_ndof)),
            ("inherited_join_count", n(self.common_parent_faces)),
            ("transverse_c1_join_count", n(self.transverse_faces)),
            ("retained_c0_join_count", n(self.retained_faces)),
            ("interfaces", Json::Array(self.interface_evidence.clone())),
            ("cells", Json::Array(cells)),
        ]);
        if let Some(selection) = &self.selection_audit {
            if let Json::Object(fields) = &mut result {
                fields.insert("continuity_selection".into(), selection.clone());
            }
        }
        result
    }
}

/// A section in x,y held fixed over an affine physical z interval, including weights.
pub(crate) struct Section<'a> {
    pub block: &'a BracketBlock,
    pub d: [usize; 3],
    pub z: [Q; 2],
    pub shear: [Q; 2],
    pub offset: usize,
}
impl Section<'_> {
    fn len(&self) -> usize {
        self.d.iter().map(|v| v + 1).product()
    }
    fn point(&self, at: [usize; 3]) -> Result<[Q; 3]> {
        let h = self.block.controls[index(self.block.degrees, [at[0], at[1], 0])];
        let z = self.z[0].checked_add(
            self.z[1]
                .checked_sub(self.z[0])?
                .checked_mul(Q::new(at[2] as i128, self.d[2] as i128)?)?,
        )?;
        let dz = z.checked_sub(self.z[0])?;
        Ok([
            Q::new(h[0].into(), h[3].into())?.checked_add(self.shear[0].checked_mul(dz)?)?,
            Q::new(h[1].into(), h[3].into())?.checked_add(self.shear[1].checked_mul(dz)?)?,
            z,
        ])
    }
    pub fn flat(&self, at: [usize; 3]) -> usize {
        self.offset + index(self.d, at)
    }
    fn height(&self) -> Result<Q> {
        self.z[1].checked_sub(self.z[0])
    }
}

fn section(block: &BracketBlock, p: usize, offset: usize) -> Result<Section<'_>> {
    if block.degrees[2] != 1 {
        return Err("unsupported: source w is not a degree-one extrusion");
    }
    let mut bounds = None;
    let mut direction = None;
    for pair in block.controls.as_chunks::<2>().0.iter() {
        let (a, b) = (pair[0], pair[1]);
        if a[3] <= 0 || a[3] != b[3] {
            return Err("unsupported: varying homogeneous section along w");
        }
        let z = [
            Q::new(a[2].into(), a[3].into())?,
            Q::new(b[2].into(), b[3].into())?,
        ];
        if z[0].cmp_value(&z[1]).is_ge() || bounds.is_some_and(|v| v != z) {
            return Err("unsupported: extrusion requires common ascending physical z planes");
        }
        let dz = z[1].checked_sub(z[0])?;
        let mut shear = [Q::zero(); 2];
        for c in 0..2 {
            shear[c] = ratio_div(
                Q::new(i128::from(b[c]) - i128::from(a[c]), i128::from(a[3]))?,
                dz,
            )?;
        }
        if direction.is_some_and(|old| old != shear) {
            return Err("unsupported: non-translational ruled section");
        }
        direction = Some(shear);
        bounds = Some(z);
    }
    Ok(Section {
        block,
        d: [block.degrees[0], block.degrees[1], p],
        z: bounds.ok_or("empty extrusion")?,
        shear: direction.ok_or("empty extrusion")?,
        offset,
    })
}

#[derive(Clone, Copy)]
struct CompilePolicy {
    match_chart_jets: bool,
    record_decisions: bool,
}

pub(crate) fn face_coord(side: &BlockSide, i: usize, j: usize, d: [usize; 3]) -> [usize; 3] {
    let axes = face_axes(side.axis);
    let mut at = [0; 3];
    at[axes[0]] = i;
    at[axes[1]] = j;
    at[side.axis] = if side.upper { d[side.axis] } else { 0 };
    at
}

fn paired_coords(
    link: &Interface,
    a: &Section,
    b: &Section,
) -> Result<Vec<([usize; 3], [usize; 3])>> {
    let [s, t] = face_axes(link.first.axis).map(|axis| a.d[axis]);
    let dims = face_axes(link.second.axis).map(|axis| b.d[axis]);
    let o = link.orientation;
    if dims != if o.swap { [t, s] } else { [s, t] } {
        return Err("unsupported: different field trace degrees");
    }
    let mut pairs = Vec::new();
    for i in 0..=s {
        for j in 0..=t {
            let k = if o.flip_first { s - i } else { i };
            let l = if o.flip_second { t - j } else { j };
            let [x, y] = if o.swap { [l, k] } else { [k, l] };
            pairs.push((
                face_coord(&link.first, i, j, a.d),
                face_coord(&link.second, x, y, b.d),
            ));
        }
    }
    Ok(pairs)
}

fn choose(n: usize, k: usize) -> i128 {
    match (n, k) {
        (2, 1) => 2,
        (3, 1 | 2) => 3,
        _ => 1,
    }
}

/// Enforce homogeneous numerator derivatives divided by each trace's gauge.
/// A common positive denominator with the same jets implies the rational field identity.
fn jet_row(
    row: &mut Row,
    s: &Section,
    mut at: [usize; 3],
    side: &BlockSide,
    order: usize,
    interval: (Q, i128),
    ids: &[usize],
) -> Result<()> {
    let axis = side.axis;
    let p = s.d[axis];
    let face_weight = s.block.controls[index(s.block.degrees, [at[0], at[1], 0])][3];
    let (length, sign) = interval;
    let mut scale = Q::new(sign, 1)?;
    for _ in 0..order {
        scale = ratio_div(scale, length)?;
    }
    for k in 0..=order {
        at[axis] = if side.upper { p - order + k } else { k };
        let h = s.block.controls[index(s.block.degrees, [at[0], at[1], 0])];
        let coefficient = choose(order, k) * if (order - k).is_multiple_of(2) { 1 } else { -1 };
        let c =
            Q::new(coefficient * i128::from(h[3]), i128::from(face_weight))?.checked_mul(scale)?;
        add(row, ids[s.flat(at)], c)?;
    }
    Ok(())
}

/// Identify an affine change of the normal coordinate from the complete first
/// homogeneous jet. Face coefficients already share a constant gauge by
/// incidence. The ratio compares outward normals with opposite signs.
fn normal_scale(link: &Interface, a: &Section, b: &Section) -> Result<Option<Q>> {
    let outward = |section: &Section, side: &BlockSide, mut at: [usize; 3]| -> Result<[Q; 4]> {
        at[2] = usize::from(at[2] == section.d[2]);
        let block = section.block;
        let face = block.controls[index(block.degrees, at)];
        if side.upper {
            at[side.axis] -= 1;
        } else {
            at[side.axis] += 1;
        }
        let inner = block.controls[index(block.degrees, at)];
        let mut result = [Q::zero(); 4];
        for (k, value) in result.iter_mut().enumerate() {
            *value = Q::new(
                (i128::from(face[k]) - i128::from(inner[k])) * block.degrees[side.axis] as i128,
                i128::from(face[3]),
            )?;
        }
        Ok(result)
    };
    let mut scale = None;
    for (i, j) in paired_coords(link, a, b)? {
        if i[2] != 0 && i[2] != a.d[2] {
            continue;
        }
        let first = outward(a, &link.first, i)?;
        let second = outward(b, &link.second, j)?;
        for (x, y) in first.into_iter().zip(second) {
            if y == Q::zero() {
                if x != Q::zero() {
                    return Ok(None);
                }
            } else {
                let r = ratio_div(Q::zero().checked_sub(x)?, y)?;
                if r.numerator() <= 0 || scale.is_some_and(|v| v != r) {
                    return Ok(None);
                }
                scale = Some(r);
            }
        }
    }
    Ok(scale)
}

pub fn build_axial_space(model: &Json, degree: usize, continuity: usize) -> Result<InheritedSpace> {
    build_space(model, degree, continuity, false, &[])
}

/// Each selected block pair names an existing lateral interface to enforce physically C1.
pub fn build_space(
    model: &Json,
    degree: usize,
    continuity: usize,
    restore_strip: bool,
    transverse: &[(String, String)],
) -> Result<InheritedSpace> {
    build_plane_space(model, degree, continuity, restore_strip, transverse, 0)
}

/// Raise the two section degrees to at least `plane_degree` before solving for
/// the smooth field. The original source geometry and source topology remain audited.
pub fn build_plane_space(
    model: &Json,
    degree: usize,
    continuity: usize,
    restore_strip: bool,
    transverse: &[(String, String)],
    plane_degree: usize,
) -> Result<InheritedSpace> {
    compile_extrusion(
        model,
        degree,
        continuity,
        restore_strip,
        transverse,
        plane_degree,
        CompilePolicy {
            match_chart_jets: false,
            record_decisions: false,
        },
    )
}

fn compile_extrusion(
    model: &Json,
    degree: usize,
    continuity: usize,
    restore_strip: bool,
    transverse: &[(String, String)],
    plane_degree: usize,
    policy: CompilePolicy,
) -> Result<InheritedSpace> {
    let CompilePolicy {
        match_chart_jets,
        record_decisions,
    } = policy;
    if plane_degree > 3 {
        return Err("unsupported transverse degree");
    }
    if !(1..=3).contains(&degree) || continuity > degree {
        return Err("invalid field degree or continuity");
    }
    if audit(model)?.field("status")?.text()? != "valid" {
        return Err("source geometry is not certified");
    }
    let mut blocks = decode_blocks(model)?;
    for block in &mut blocks {
        let target = [
            block.degrees[0].max(plane_degree),
            block.degrees[1].max(plane_degree),
            block.degrees[2],
        ];
        if target != block.degrees {
            let net = crate::nurbs_space::elevate_to(block, target);
            let gcd = |mut a: i128, mut b: i128| {
                while b != 0 {
                    (a, b) = (b, a % b);
                }
                a.abs()
            };
            let divisor = net.iter().flatten().copied().fold(0, gcd);
            if divisor == 0 {
                return Err("zero elevated net");
            }
            let mut controls = Vec::with_capacity(net.len());
            for h in net {
                let mut point = [0i64; 4];
                for k in 0..4 {
                    point[k] = i64::try_from(h[k] / divisor)
                        .map_err(|_| "indeterminate: elevated net exceeds i64")?;
                }
                controls.push(point);
            }
            block.degrees = target;
            block.controls = controls;
        }
    }
    if blocks.is_empty() || blocks.len() > 2048 {
        return Err("indeterminate: extrusion block budget");
    }
    let links = enumerate_faces(&blocks)?;
    let mut offset = 0;
    let mut sections = Vec::new();
    let mut by_id = BTreeMap::new();
    for (i, b) in blocks.iter().enumerate() {
        let s = section(b, degree, offset)?;
        offset += s.len();
        sections.push(s);
        by_id.insert(b.id.as_str(), i);
    }
    if offset > 120000 {
        return Err("indeterminate: raw scalar field budget");
    }
    let mut uf = UnionFind::new(offset);
    let mut shared_plane = BTreeSet::new();
    for link in &links.interfaces {
        let a = &sections[by_id[link.first.block.as_str()]];
        let b = &sections[by_id[link.second.block.as_str()]];
        for (i, j) in paired_coords(link, a, b)? {
            uf.join(a.flat(i), b.flat(j));
            if link.first.axis == 2 {
                shared_plane.insert(a.flat(i));
            }
        }
    }
    let mut c0 = BTreeMap::new();
    let mut ids = Vec::new();
    for i in 0..offset {
        let next = c0.len();
        ids.push(*c0.entry(uf.find(i)).or_insert(next));
    }
    let n0 = c0.len();
    let mut xyz = vec![None; n0];
    let mut priority = vec![1; n0];
    for s in &sections {
        for i in 0..=s.d[0] {
            for j in 0..=s.d[1] {
                for k in 0..=degree {
                    let at = [i, j, k];
                    let id = ids[s.flat(at)];
                    let pt = s.point(at)?;
                    if xyz[id].is_some_and(|v| v != pt) {
                        return Err("field coefficient trace positions differ");
                    }
                    xyz[id] = Some(pt);
                    if k == 0 || k == degree {
                        priority[id] = 2;
                    }
                }
            }
        }
    }
    for i in shared_plane {
        priority[ids[i]] = 0;
    }
    // Pivot on internal planes first, retaining bottom/top trace roots whenever possible.
    let mut order: Vec<_> = (0..n0).collect();
    order.sort_by_key(|&i| (priority[i], i));
    let mut priority_id = vec![0; n0];
    for (j, &i) in order.iter().enumerate() {
        priority_id[i] = j;
    }
    let pids: Vec<_> = ids.iter().map(|&i| priority_id[i]).collect();
    let mut eq = Vec::new();
    let mut lateral = Vec::new();
    let mut evidence = Vec::new();
    let mut inherited = 0;
    let mut coupled = 0;
    let wanted: BTreeSet<_> = transverse
        .iter()
        .map(|(a, b)| if a < b { (a, b) } else { (b, a) })
        .collect();
    if transverse.len() > 256 {
        return Err("indeterminate: transverse interface budget");
    }
    if wanted.len() != transverse.len() {
        return Err("duplicate transverse join request");
    }
    let strip_parent = if restore_strip {
        let parent = crate::spline_support::recover_nurbs_strip_parent(model)?;
        if plane_degree > parent.degrees[0] {
            return Err("unsupported elevated source-curve continuity");
        }
        Some(parent)
    } else {
        None
    };
    for link in &links.interfaces {
        let ai = by_id[link.first.block.as_str()];
        let bi = by_id[link.second.block.as_str()];
        let (a, b) = (&sections[ai], &sections[bi]);
        let key = if link.first.block < link.second.block {
            (&link.first.block, &link.second.block)
        } else {
            (&link.second.block, &link.first.block)
        };
        let mut record = object([
            ("first", link_side(&link.first)),
            ("second", link_side(&link.second)),
            ("role", Json::String("whole_face_join".into())),
            ("field_order", n(0)),
            (
                "basis_evidence",
                Json::String("oriented_proportional_homogeneous_face".into()),
            ),
        ]);
        if record_decisions {
            interface_reason(
                &mut record,
                if link.first.axis == 2 && link.second.axis == 2 {
                    "C0 retained: section maps or physical extrusion directions do not certify a common-parent join"
                } else {
                    "C0 retained: no qualified lateral smoothing rule applies to this interface"
                },
            );
        }
        let mut set_evidence = |role: &str, method: &str, r: usize| {
            if let Json::Object(fields) = &mut record {
                fields.insert("role".into(), Json::String(role.into()));
                fields.insert("basis_evidence".into(), Json::String(method.into()));
                fields.insert("field_order".into(), n(r));
            }
        };
        if match_chart_jets && link.first.axis < 2 && link.second.axis < 2 {
            if let Some(scale) = normal_scale(link, a, b)? {
                let sign_b = if link.first.upper == link.second.upper {
                    1
                } else {
                    -1
                };
                for (i, j) in paired_coords(link, a, b)? {
                    let mut row = Row::new();
                    jet_row(
                        &mut row,
                        a,
                        i,
                        &link.first,
                        1,
                        (scale, a.d[link.first.axis] as i128),
                        &pids,
                    )?;
                    jet_row(
                        &mut row,
                        b,
                        j,
                        &link.second,
                        1,
                        (Q::one(), sign_b * b.d[link.second.axis] as i128),
                        &pids,
                    )?;
                    lateral.push(row);
                }
                set_evidence(
                    "distinct_support_join",
                    "exact_proportional_complete_homogeneous_normal_jet",
                    1,
                );
                if record_decisions {
                    if let Json::Object(fields) = &mut record {
                        fields.insert(
                            "normal_interval_ratio".into(),
                            Json::Array(vec![n(scale.numerator()), n(scale.denominator())]),
                        );
                    }
                }
                coupled += 1;
                evidence.push(record);
                continue;
            }
        }
        if wanted.contains(&key) {
            if strip_parent.is_some() && link.first.axis == 0 && link.second.axis == 0 {
                return Err("duplicate source-parent and transverse join request");
            }
            if link.first.axis == 2 || link.second.axis == 2 {
                return Err("expected lateral transverse join");
            }
            crate::transverse_jets::append(link, a, b, &pids, &mut lateral)?;
            coupled += 1;
            set_evidence(
                "distinct_support_join",
                "exact_physical_gradient_unisolvence",
                1,
            );
            if record_decisions {
                interface_reason(
                    &mut record,
                    "qualified cubic radial face: exact physical-gradient constraints applied",
                );
            }
            evidence.push(record);
            continue;
        }
        let (r, ha, hb) = if link.first.axis == 2
            && link.second.axis == 2
            && a.shear == b.shear
            && continuity > 0
        {
            // Incidence proves the entire rational section matches up to a constant
            // homogeneous gauge and the section check proves extrusion in physical z.
            (continuity, a.height()?, b.height()?)
        } else if restore_strip && link.first.axis == 0 && link.second.axis == 0 {
            let parent = strip_parent.as_ref().ok_or("missing source parent")?;
            let restriction = |id: &str| {
                parent
                    .source_spans
                    .iter()
                    .find(|span| span.source_block_id == id)
                    .map(|span| &span.restriction)
                    .ok_or("missing source restriction")
            };
            let ar = restriction(&link.first.block)?;
            let br = restriction(&link.second.block)?;
            let [x, y] = ar.intervals[0];
            let [s, t] = br.intervals[0];
            let seam = ar.intervals[0][usize::from(link.first.upper ^ ar.reversals[0])];
            if seam != br.intervals[0][usize::from(link.second.upper ^ br.reversals[0])] {
                return Err("parent restrictions disagree at interface");
            }
            let knots = model.field("lower")?.field("knots")?.items()?;
            let m = knots
                .iter()
                .map(Json::number::<i64>)
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .filter(|&v| v == seam)
                .count();
            let smooth = a.d[0].checked_sub(m).ok_or("invalid knot multiplicity")?;
            if smooth == 0 {
                if record_decisions {
                    interface_reason(
                        &mut record,
                        "C0 retained: recovered source knot multiplicity permits no higher order",
                    );
                }
                evidence.push(record);
                continue;
            }
            (
                smooth,
                Q::new(i128::from(y) - i128::from(x), 1)?,
                Q::new(i128::from(t) - i128::from(s), 1)?,
            )
        } else {
            evidence.push(record);
            continue;
        };
        if link.first.upper == link.second.upper {
            return Err("unsupported reversed normal parent interval");
        }
        for (i, j) in paired_coords(link, a, b)? {
            for k in 1..=r {
                let mut row = Row::new();
                jet_row(&mut row, a, i, &link.first, k, (ha, 1), &pids)?;
                jet_row(&mut row, b, j, &link.second, k, (hb, -1), &pids)?;
                eq.push(row);
            }
        }
        set_evidence(
            "common_parent_subdivision",
            if link.first.axis == 2 {
                "same_homogeneous_section_affine_physical_z"
            } else {
                "recovered_source_knots_weighted_numerator_jets"
            },
            r,
        );
        if record_decisions {
            interface_reason(
                &mut record,
                "exact inherited-parent coefficient and mapping certificate",
            );
        }
        evidence.push(record);
        inherited += 1;
    }
    if !match_chart_jets && coupled != wanted.len() {
        return Err("requested transverse interface not found");
    }
    let (prior_lift, free) = crate::spline_constraints::joint_kernel(n0, &eq, &lateral)?;
    let mut free_original: Vec<_> = free.iter().map(|&i| order[i]).collect();
    free_original.sort_unstable();
    let free_indices: BTreeMap<_, _> = free_original
        .iter()
        .enumerate()
        .map(|(i, &r)| (r, i))
        .collect();
    let coords: Vec<_> = free_original.iter().map(|&i| xyz[i].unwrap()).collect();
    let lift: Vec<Row> = priority_id
        .iter()
        .map(|&i| {
            prior_lift[i]
                .iter()
                .map(|(&k, &v)| (free_indices[&order[free[k]]], v))
                .collect()
        })
        .collect();
    for (i, row) in lift.iter().enumerate() {
        let mut h = [Q::zero(); 4];
        for (&j, &k) in row {
            for (v, x) in h.iter_mut().zip(coords[j].into_iter().chain([Q::one()])) {
                *v = v.checked_add(k.checked_mul(x)?)?;
            }
        }
        let pos = xyz[i].unwrap();
        if h != [pos[0], pos[1], pos[2], Q::one()] {
            return Err("jet space lost an affine coordinate");
        }
    }
    let mut cells = Vec::new();
    for s in sections {
        let extraction: Vec<_> = ids[s.offset..s.offset + s.len()]
            .iter()
            .map(|&i| lift[i].clone())
            .collect();
        let active: Vec<_> = extraction
            .iter()
            .flat_map(|r| r.keys().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let local: BTreeMap<_, _> = active.iter().enumerate().map(|(i, &j)| (j, i)).collect();
        let stencils = extraction
            .iter()
            .map(|r| r.iter().map(|(i, &c)| (local[i], val(c))).collect())
            .collect();
        let weights = s
            .block
            .controls
            .as_chunks::<2>()
            .0
            .iter()
            .flat_map(|pair| std::iter::repeat_n(pair[0][3] as f64, degree + 1))
            .collect();
        let geometry = RationalBezierMap::new(
            s.block.degrees,
            s.block
                .controls
                .iter()
                .map(|h| h.map(|v| v as f64))
                .collect(),
        )
        .map_err(|_| "invalid source map")?;
        cells.push(InheritedCell {
            id: s.block.id.clone(),
            degrees: s.d,
            z_bounds: s.z.map(val),
            extraction,
            active,
            weights,
            stencils,
            geometry,
        });
    }
    Ok(InheritedSpace {
        cells,
        ndof: free.len(),
        control_xyz: coords.into_iter().map(|p| p.map(val)).collect(),
        c0_ndof: n0,
        common_parent_faces: inherited,
        retained_faces: links.interfaces.len() - inherited - coupled,
        transverse_faces: coupled,
        interface_evidence: evidence,
        selection_audit: None,
        lift,
    })
}

/// Select only continuity supported by exact inherited mappings and
/// certificates. Uncertified lateral interfaces remain C0 with an audit reason.
pub fn build_automatic_space(model: &Json, degree: usize) -> Result<InheritedSpace> {
    if !(1..=3).contains(&degree) {
        return Err("invalid field degree");
    }
    let qualified_pocket = crate::pocket_convergence::matches_qualified_layout(model)?;
    let radial_joins = if qualified_pocket {
        crate::pocket_convergence::radial_interfaces_after_qualification(model)?
            .into_iter()
            .map(|face| (face.first_block, face.second_block))
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    let selected_degree = degree;
    let selected_plane_degree = if qualified_pocket { 3 } else { 0 };
    let restore_strip = model.field("schema")?.text()? == "nurbs_sketch_strip_v1";
    let parent_decision = if restore_strip {
        "exact recovered source-parent restrictions available"
    } else {
        "no supported exact source-parent recovery; retain C0 on uncertified interfaces"
    };
    // Translation across a shared section proves an exact C1 match. Higher
    // axial order is inherited only when the source-parent knot certificate
    // establishes it.
    let axial_order = usize::from(selected_degree > 0);
    let mut space = compile_extrusion(
        model,
        selected_degree,
        axial_order,
        restore_strip,
        &radial_joins,
        selected_plane_degree,
        CompilePolicy {
            match_chart_jets: false,
            record_decisions: true,
        },
    )?;
    space.selection_audit = Some(object([
        (
            "policy",
            Json::String("conservative_exact_inheritance_v1".into()),
        ),
        (
            "outcome",
            Json::String(
                if space.common_parent_faces + space.transverse_faces > 0 {
                    "exact_certificates_applied"
                } else {
                    "retain_c0"
                }
                .into(),
            ),
        ),
        ("requested_degree", n(degree)),
        ("selected_degree", n(selected_degree)),
        ("axial_order_cap", n(selected_degree.min(1))),
        ("plane_degree", n(selected_plane_degree)),
        (
            "parent_mapping_decision",
            Json::String(parent_decision.into()),
        ),
        ("qualified_radial_join_count", n(radial_joins.len())),
        (
            "geometry_qualification",
            Json::String(
                if qualified_pocket {
                    "exact qualified cubic pocket layout"
                } else {
                    "outside the qualified pocket layout"
                }
                .into(),
            ),
        ),
        ("fallback_reason", Json::Null),
        (
            "contact_physics",
            Json::String("unsupported; no contact model is assembled".into()),
        ),
    ]));
    Ok(space)
}

fn link_side(s: &BlockSide) -> Json {
    Json::Array(vec![
        Json::String(s.block.clone()),
        Json::String(format!(
            "{}{}",
            ["u", "v", "w"][s.axis],
            usize::from(s.upper)
        )),
    ])
}

fn interface_reason(record: &mut Json, reason: &str) {
    if let Json::Object(fields) = record {
        fields.insert("reason".into(), Json::String(reason.into()));
    }
}

/// Compile an independently replayable field descriptor, retaining geometry, boundaries and job plan.
pub fn space_from_plan(model: &Json, plan: &Json) -> Result<InheritedSpace> {
    if plan.field("schema")?.text()? == "inherited_field_auto_plan_v1" {
        plan.keys(&["schema", "axial_degree"])?;
        return build_automatic_space(model, plan.field("axial_degree")?.number()?);
    }
    if plan.field("schema")?.text()? == "matched_chart_jets_v1" {
        plan.keys(&["schema", "axial_degree", "axial_continuity"])?;
        return compile_extrusion(
            model,
            plan.field("axial_degree")?.number()?,
            plan.field("axial_continuity")?.number()?,
            false,
            &[],
            0,
            CompilePolicy {
                match_chart_jets: true,
                record_decisions: false,
            },
        );
    }
    let plane = if let Ok(p) = plan.field("plane_degree") {
        plan.keys(&[
            "schema",
            "axial_degree",
            "axial_continuity",
            "restore_strip",
            "transverse_interfaces",
            "plane_degree",
        ])?;
        p.number()?
    } else {
        plan.keys(&[
            "schema",
            "axial_degree",
            "axial_continuity",
            "restore_strip",
            "transverse_interfaces",
        ])?;
        0
    };
    if plan.field("schema")?.text()? != "inherited_field_plan_v1" {
        return Err("unsupported inherited field plan");
    }
    let restore = match plan.field("restore_strip")? {
        Json::Bool(v) => *v,
        _ => return Err("restore_strip must be boolean"),
    };
    let pairs = plan
        .field("transverse_interfaces")?
        .items()?
        .iter()
        .map(|v| {
            let a = v.items()?;
            if a.len() != 2 {
                return Err("expected two named blocks per transverse join");
            }
            Ok((a[0].text()?.to_owned(), a[1].text()?.to_owned()))
        })
        .collect::<Result<Vec<_>>>()?;
    build_plane_space(
        model,
        plan.field("axial_degree")?.number()?,
        plan.field("axial_continuity")?.number()?,
        restore,
        &pairs,
        plane,
    )
}
pub fn automatic_plan(degree: usize) -> Json {
    object([
        (
            "schema",
            Json::String("inherited_field_auto_plan_v1".into()),
        ),
        ("axial_degree", n(degree)),
    ])
}
pub fn compile_automatic(model: &Json, degree: usize) -> Result<Json> {
    compile(model, &automatic_plan(degree))
}
pub fn compile(model: &Json, plan: &Json) -> Result<Json> {
    let space = space_from_plan(model, plan)?;
    Ok(object([
        (
            "schema",
            Json::String("solidcraft_inherited_analysis_v1".into()),
        ),
        ("geometry", model.clone()),
        ("plan", plan.clone()),
        ("field", space.document()),
    ]))
}
pub fn replay(package: &Json) -> Result<InheritedSpace> {
    package.keys(&["schema", "geometry", "plan", "field"])?;
    let space = space_from_plan(package.field("geometry")?, package.field("plan")?)?;
    if package.field("schema")?.text()? != "solidcraft_inherited_analysis_v1"
        || package.field("field")? != &space.document()
    {
        return Err("inherited field differs from canonical replay");
    }
    Ok(space)
}
pub fn audit_package(package: &Json) -> Result<Json> {
    let field = replay(package)?;
    Ok(object([
        ("status", Json::String("valid".into())),
        ("scalar_dofs", n(field.ndof)),
    ]))
}

impl InheritedSpace {
    /// Homogeneous essential data on selected full source faces. Every
    /// Bernstein trace row is constrained even when an eliminated coefficient
    /// depends on several roots whose own coordinates lie off the face.
    pub fn zero_face_lift(&self, faces: &[BlockSide]) -> Result<crate::constraint_system::Lift> {
        use crate::constraint_system::{Lift, LiftLine};
        let cells: BTreeMap<_, _> = self.cells.iter().map(|c| (c.id.as_str(), c)).collect();
        let mut eq = Vec::new();
        for side in faces {
            let c = cells
                .get(side.block.as_str())
                .ok_or("unknown essential face")?;
            if side.axis > 2 {
                return Err("invalid essential face axis");
            }
            let [s, t] = face_axes(side.axis).map(|d| c.degrees[d]);
            for i in 0..=s {
                for j in 0..=t {
                    eq.push(
                        c.extraction[index(c.degrees, face_coord(side, i, j, c.degrees))].clone(),
                    );
                }
            }
        }
        if eq.is_empty() {
            return Err("empty essential selection");
        }
        let (rows, free_roots) = nullspace(self.ndof, &eq)?;
        Ok(Lift {
            free_roots,
            lines: rows
                .into_iter()
                .enumerate()
                .map(|(dof, r)| LiftLine {
                    dof,
                    terms: r.into_iter().collect(),
                    offset: Q::zero(),
                })
                .collect(),
        })
    }
}
