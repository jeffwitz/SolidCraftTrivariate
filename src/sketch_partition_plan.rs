//! Exact planar events for two transverse sketch extrusions.
use crate::{
    constraint_system::Ratio,
    json::{object, Json},
    nurbs_curve::{extract_active_bezier, BezierSpan},
    sketch_junction::{plan_loop, Continuity},
};
use std::collections::{BTreeMap, BTreeSet};

type Error = &'static str;
type H = [Ratio; 3];

pub(crate) struct SourceCurve {
    pub degree: usize,
    pub knots: Vec<i64>,
    pub controls: Vec<[i64; 3]>,
}

impl SourceCurve {
    pub(crate) fn parse(curve: &Json) -> Result<Self, Error> {
        curve.keys(&["degree", "knots", "controls"])?;
        Ok(Self {
            degree: curve.field("degree")?.number()?,
            knots: curve
                .field("knots")?
                .items()?
                .iter()
                .map(Json::number)
                .collect::<Result<_, _>>()?,
            controls: curve
                .field("controls")?
                .items()?
                .iter()
                .map(triple)
                .collect::<Result<_, _>>()?,
        })
    }
}

fn triple(value: &Json) -> Result<[i64; 3], Error> {
    let [a, b, c] = value.items()? else {
        return Err("expected three coordinates");
    };
    Ok([a.number()?, b.number()?, c.number()?])
}
fn int(x: i64) -> Ratio {
    Ratio::new(x.into(), 1).expect("integer")
}
fn num(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn rational(x: Ratio) -> Json {
    Json::Array(vec![num(x.numerator()), num(x.denominator())])
}
fn affine_point(h: H) -> Result<[Ratio; 2], Error> {
    let inv = Ratio::new(h[2].denominator(), h[2].numerator())?;
    Ok([h[0].checked_mul(inv)?, h[1].checked_mul(inv)?])
}
#[derive(Clone, Copy)]
pub(crate) struct Crossing {
    pub span: usize,
    pub local: Ratio,
    pub position: [Ratio; 2],
}

struct RootSearch {
    level: Ratio,
    span: usize,
    found: BTreeMap<[Ratio; 2], Crossing>,
}

impl RootSearch {
    fn visit(&mut self, net: &[H], lo: Ratio, hi: Ratio, depth: u8) -> Result<(), Error> {
        let coefficients = net
            .iter()
            .map(|h| h[1].checked_sub(self.level.checked_mul(h[2])?))
            .collect::<Result<Vec<_>, Error>>()?;
        let signs: Vec<_> = coefficients
            .iter()
            .map(|r| r.numerator().signum())
            .collect();
        let last = net.len() - 1;
        for (i, near, t) in [(0, 1, lo), (last, last - 1, hi)] {
            if signs[i] == 0 {
                if signs[near] == 0 {
                    return Err("unresolved tangency or coincident sketch cut");
                }
                let position = affine_point(net[i])?;
                self.found.entry(position).or_insert(Crossing {
                    span: self.span,
                    local: t,
                    position,
                });
            }
        }
        if signs.iter().all(|&s| s >= 0) || signs.iter().all(|&s| s <= 0) {
            return Ok(());
        }
        if let Some(t) = crate::sketch_rational_event::simple_zero(&coefficients) {
            let complement = Ratio::one().checked_sub(t)?;
            let mut row = net.to_vec();
            for n in (1..=last).rev() {
                for i in 0..n {
                    let next = row[i + 1];
                    for (v, other) in row[i].iter_mut().zip(next) {
                        *v = v
                            .checked_mul(complement)?
                            .checked_add(other.checked_mul(t)?)?;
                    }
                }
            }
            let position = affine_point(row[0])?;
            let local = lo.checked_add(hi.checked_sub(lo)?.checked_mul(t)?)?;
            self.found.entry(position).or_insert(Crossing {
                span: self.span,
                local,
                position,
            });
            return Ok(());
        }
        if depth == 0 {
            return Err(
                "transverse intersection lacks a certified rational event in the search budget",
            );
        }
        let half = Ratio::new(1, 2)?;
        let mut row = net.to_vec();
        let mut left = vec![row[0]];
        let mut right = vec![row[last]];
        for n in (1..=last).rev() {
            for i in 0..n {
                let next = row[i + 1];
                for (v, other) in row[i].iter_mut().zip(next) {
                    *v = v.checked_add(other)?.checked_mul(half)?;
                }
            }
            left.push(row[0]);
            right.push(row[n - 1]);
        }
        right.reverse();
        let middle = lo.checked_add(hi)?.checked_mul(half)?;
        self.visit(&left, lo, middle, depth - 1)?;
        self.visit(&right, middle, hi, depth - 1)
    }
}
pub(crate) fn crossings(spans: &[BezierSpan], level: Ratio) -> Result<Vec<Crossing>, Error> {
    let mut search = RootSearch {
        level,
        span: 0,
        found: BTreeMap::new(),
    };
    for (index, span) in spans.iter().enumerate() {
        search.span = index;
        let net: Vec<_> = span.controls.iter().map(|p| p.map(int)).collect();
        search.visit(&net, Ratio::zero(), Ratio::one(), 12)?;
    }
    if search.found.len() != 2 {
        return Err("expected two distinct transverse intersections per sketch level");
    }
    Ok(search.found.into_values().collect())
}

fn rectangle_levels(spans: &[BezierSpan]) -> Result<[[Ratio; 2]; 2], Error> {
    let mut lines: [BTreeSet<Ratio>; 2] = Default::default();
    for span in spans {
        let pts = span
            .controls
            .iter()
            .map(|p| affine_point(p.map(int)))
            .collect::<Result<Vec<_>, _>>()?;
        let axis = (0..2)
            .find(|&d| pts.iter().all(|p| p[d] == pts[0][d]))
            .ok_or("transverse sketch has a curved or oblique side")?;
        let free = 1 - axis;
        let direction = pts[0][free].cmp_value(&pts[pts.len() - 1][free]);
        if direction.is_eq()
            || pts.windows(2).any(|w| {
                let step = w[0][free].cmp_value(&w[1][free]);
                !step.is_eq() && step != direction
            })
        {
            return Err("transverse sketch side folds or collapses");
        }
        lines[axis].insert(pts[0][axis]);
    }
    let levels = |s: BTreeSet<Ratio>| -> Result<[Ratio; 2], Error> {
        let mut ordered: Vec<_> = s.into_iter().collect();
        ordered.sort_by(Ratio::cmp_value);
        ordered
            .try_into()
            .map_err(|_| "transverse sketch needs exactly two distinct lines per axis")
    };
    let [y, z] = lines;
    Ok([levels(y)?, levels(z)?])
}

pub(crate) fn plan(recipe: &Json) -> Result<(Json, Json), Error> {
    recipe.keys(&["family", "bounds_xyz_mm", "sketches"])?;
    if recipe.field("family")?.text()? != "g2_cross_plate_v1" {
        return Err("unexpected transverse sketch family");
    }
    let bounds = recipe.field("bounds_xyz_mm")?;
    let [xb, _, zb] = bounds.items()? else {
        return Err("expected stock bounds on three axes");
    };
    let extent = |v: &Json| -> Result<i64, Error> {
        let [a, b] = v.items()? else {
            return Err("expected stock interval");
        };
        if a.number::<i64>()? != 0 {
            return Err("unsupported extrusion origin");
        }
        b.number()
    };
    let mut ids = BTreeSet::new();
    let mut profile = None;
    let mut rect = None;
    let mut source_events = Vec::new();
    for sketch in recipe.field("sketches")?.items()? {
        sketch.keys(&[
            "id",
            "operation",
            "frame",
            "direction",
            "length_mm",
            "loops",
        ])?;
        let id = sketch.field("id")?.text()?;
        if id.is_empty() || !ids.insert(id) || sketch.field("operation")?.text()? != "subtract" {
            return Err("invalid transverse sketch operation");
        }
        let frame = sketch.field("frame")?;
        frame.keys(&["origin_mm", "u", "v"])?;
        if triple(frame.field("origin_mm")?)? != [0; 3] {
            return Err("unsupported transverse sketch frame origin");
        }
        let (u, v, d) = (
            triple(frame.field("u")?)?,
            triple(frame.field("v")?)?,
            triple(sketch.field("direction")?)?,
        );
        let axial = match (u, v, d) {
            ([1, 0, 0], [0, 1, 0], [0, 0, 1]) => true,
            ([0, 1, 0], [0, 0, 1], [1, 0, 0]) => false,
            _ => return Err("unsupported transverse sketch plane or extrusion direction"),
        };
        if sketch.field("length_mm")?.number::<i64>()? != extent(if axial { zb } else { xb })? {
            return Err("extrusion must traverse the full stock");
        }
        let [lp] = sketch.field("loops")?.items()? else {
            return Err("one closed contour is required per transverse operation");
        };
        let loop_id = lp.field("id")?.text()?;
        if loop_id.is_empty() {
            return Err("missing transverse sketch loop identifier");
        }
        if let Ok(definition) = lp.field("primitive") {
            let stored_curve = lp.field("curve").ok();
            if stored_curve.is_some() {
                lp.keys(&["id", "primitive", "curve"])?;
            } else {
                lp.keys(&["id", "primitive"])?;
            }
            let primitive = crate::primitive_template_json::primitive(definition)?;
            let crate::primitive_template::Primitive::Rectangle { min, max } = primitive else {
                return Err("unsupported native cutter for transverse plate");
            };
            if axial || rect.is_some() {
                return Err("expected one transverse native rectangle");
            }
            let template = primitive.plan(crate::primitive_template::Operation::Extrude {
                length: extent(xb)?,
            })?;
            let generated = crate::primitive_template_json::contour(&template);
            if stored_curve.is_some_and(|curve| *curve != generated) {
                return Err("stored NURBS does not match native primitive provenance");
            }
            let junctions = template
                .junctions
                .iter()
                .enumerate()
                .map(|(i, kind)| {
                    object([
                        ("right_span", num(i)),
                        ("knot", num(template.contour.knots[i + 1])),
                        ("multiplicity", num(if i == 0 { 2 } else { 1 })),
                        ("closure", Json::Bool(i == 0)),
                        (
                            "continuity",
                            Json::String(crate::primitive_template_json::continuity(*kind).into()),
                        ),
                    ])
                })
                .collect();
            rect = Some(([[min[0], max[0]], [min[1], max[1]]], id, loop_id));
            source_events.push(object([
                ("sketch", Json::String(id.into())),
                ("loop", Json::String(loop_id.into())),
                ("curve", generated),
                (
                    "provenance",
                    Json::String("native_rectangle_template".into()),
                ),
                ("junctions", Json::Array(junctions)),
            ]));
            continue;
        }
        lp.keys(&["id", "curve"])?;
        let raw = lp.field("curve")?;
        let curve = SourceCurve::parse(raw)?;
        let spans = extract_active_bezier(curve.degree, &curve.knots, &curve.controls)?;
        let joins = plan_loop(curve.degree, &curve.knots, &curve.controls)?;
        let mut corners = 0;
        let mut junctions = Vec::new();
        for j in joins {
            let class = match j.continuity {
                Continuity::Curvature => "G2",
                Continuity::Tangent => "G1",
                Continuity::Corner => {
                    corners += 1;
                    "G0"
                }
                _ => return Err("disconnected or singular transverse sketch loop"),
            };
            if axial && j.continuity != Continuity::Curvature {
                return Err("axial free contour needs G2 at all junctions");
            }
            junctions.push(object([
                ("right_span", num(j.right_span)),
                ("knot", num(j.parameter)),
                ("multiplicity", num(j.multiplicity)),
                ("closure", Json::Bool(j.closure)),
                ("continuity", Json::String(class.into())),
            ]));
        }
        if axial {
            if curve.degree != 3 || profile.is_some() {
                return Err("expected one cubic free contour");
            }
            profile = Some((raw, spans, id, loop_id));
        } else {
            if corners != 4 || rect.is_some() {
                return Err("transverse contour needs four geometric corners");
            }
            rect = Some((rectangle_levels(&spans)?, id, loop_id));
        }
        source_events.push(object([
            ("sketch", Json::String(id.into())),
            ("loop", Json::String(loop_id.into())),
            ("provenance", Json::String("imported_nurbs_jets".into())),
            ("junctions", Json::Array(junctions)),
        ]));
    }
    let (raw, spans, profile_id, profile_loop) = profile.ok_or("missing axial free contour")?;
    let ([ys, zs], rect_id, rect_loop) = rect.ok_or("missing transverse rectangular contour")?;
    let mut induced = Vec::new();
    for y in ys {
        for root in crossings(&spans, y)? {
            let [a, b] = spans[root.span].interval;
            let t = int(a).checked_add(root.local.checked_mul(int(b).checked_sub(int(a))?)?)?;
            induced.push(object([
                ("profile", Json::String(profile_id.into())),
                ("profile_loop", Json::String(profile_loop.into())),
                ("transverse_sketch", Json::String(rect_id.into())),
                ("transverse_loop", Json::String(rect_loop.into())),
                ("y_mm", rational(y)),
                ("span", num(root.span)),
                ("parameter", rational(t)),
                ("xy_mm", Json::Array(root.position.map(rational).into())),
            ]));
        }
    }
    if zs.iter().any(|v| v.denominator() != 1) {
        return Err("nonintegral transverse axial layers outside this partition family");
    }
    let normalized = object([
        ("family", Json::String("g2_cross_plate_v1".into())),
        ("bounds_xyz_mm", bounds.clone()),
        ("profile", raw.clone()),
        ("rect_y_mm", Json::Array(ys.map(rational).into())),
        (
            "rect_z_mm",
            Json::Array(zs.map(|v| num(v.numerator())).into()),
        ),
    ]);
    let report = object([
        ("sketches", Json::Array(source_events)),
        ("crossings", Json::Array(induced)),
        ("z_cut_levels_mm", Json::Array(zs.map(rational).into())),
    ]);
    Ok((normalized, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::parse;

    #[test]
    fn four_interior_crossings_are_planned_from_the_two_input_curves() {
        let input = parse(include_str!(
            "../tests/references/g2_cross_sketches.recipe.json"
        ))
        .unwrap();
        let (recipe, report) = plan(&input).unwrap();
        assert_eq!(
            recipe.field("rect_y_mm").unwrap().to_string(),
            "[[719,36],[1081,36]]"
        );
        let parameters: BTreeSet<_> = report
            .field("crossings")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|e| e.field("parameter").unwrap().to_string())
            .collect();
        assert_eq!(
            parameters,
            ["[5,2]", "[7,2]", "[13,2]", "[15,2]"]
                .map(String::from)
                .into()
        );
        let native = parse(include_str!(
            "../tests/references/g2_cross_native.recipe.json"
        ))
        .unwrap();
        let (native_recipe, native_report) = plan(&native).unwrap();
        assert_eq!(native_recipe, recipe);
        assert_eq!(native_report.field("crossings"), report.field("crossings"));
        assert_eq!(
            native_report.field("sketches").unwrap().items().unwrap()[1]
                .field("provenance")
                .unwrap()
                .text()
                .unwrap(),
            "native_rectangle_template"
        );
    }

    #[test]
    fn ordinary_y_levels_return_an_unresolved_exact_event() {
        let mut input = parse(include_str!(
            "../tests/references/g2_cross_sketches.recipe.json"
        ))
        .unwrap();
        let Json::Object(root) = &mut input else {
            unreachable!()
        };
        let Json::Array(sketches) = root.get_mut("sketches").unwrap() else {
            unreachable!()
        };
        let Json::Object(rect) = &mut sketches[1] else {
            unreachable!()
        };
        rect.insert("loops".into(),parse(r#"[{"id":"box","curve":{"degree":1,"knots":[0,0,1,2,3,4,4],"controls":[[20,4,1],[30,4,1],[30,8,1],[20,8,1],[20,4,1]]}}]"#).unwrap());
        assert!(plan(&input)
            .unwrap_err()
            .contains("certified rational event"));
    }
}
