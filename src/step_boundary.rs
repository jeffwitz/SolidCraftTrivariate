//! AP214 exchange of a certified rational B-rep with optional interior voids.

use crate::{
    boundary_shell::{boundary_shell, BoundaryShell},
    json::Json,
};
use std::fmt::Write as _;

/// Proleptic Gregorian UTC timestamp for the Part 21 exchange header.
fn exchange_timestamp() -> Result<String, &'static str> {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|_| "STEP clock predates the Unix epoch")?
        .as_secs();
    let days = seconds / 86400;
    // Work in March-based 400-year cycles (146097 days).
    let day = days + 719468;
    let era = day / 146097;
    let doe = day % 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let month = (mp + 2) % 12 + 1;
    let year = yoe + era * 400 + u64::from(month <= 2);
    if year > 9999 {
        return Err("STEP date exceeds four-digit year");
    }
    let date = doy - (153 * mp + 2) / 5 + 1;
    Ok(format!(
        "{year:04}-{month:02}-{date:02}T{:02}:{:02}:{:02}Z",
        seconds / 3600 % 24,
        seconds / 60 % 60,
        seconds % 60
    ))
}

struct Writer {
    next: usize,
    data: String,
}
impl Writer {
    fn entity(&mut self, value: String) -> usize {
        let id = self.next;
        self.next += 1;
        writeln!(self.data, "#{id} = {value};").unwrap();
        id
    }
}

fn number(value: f64) -> String {
    if value == 0.0 {
        "0.".to_owned()
    } else {
        format!("{value:.16E}")
    }
}
fn point4(h: [i128; 4]) -> Result<[f64; 3], &'static str> {
    if h[3] <= 0 {
        return Err("invalid_boundary_weight");
    }
    let p = [
        h[0] as f64 / h[3] as f64,
        h[1] as f64 / h[3] as f64,
        h[2] as f64 / h[3] as f64,
    ];
    if p.iter().any(|x| !x.is_finite()) {
        return Err("boundary_coordinate_out_of_range");
    }
    Ok(p)
}
fn coord(p: [f64; 3]) -> String {
    format!("({},{},{})", number(p[0]), number(p[1]), number(p[2]))
}
fn step_text(s: &str) -> String {
    s.replace('\'', "''")
}
fn cartesian(w: &mut Writer, p: [f64; 3]) -> usize {
    w.entity(format!("CARTESIAN_POINT('',{})", coord(p)))
}
fn knots(degree: usize) -> String {
    format!("({},{}),(0.,1.)", degree + 1, degree + 1)
}
fn curve_entity(w: &mut Writer, controls: &[[i64; 4]]) -> Result<usize, &'static str> {
    if controls.len() < 2 || controls.len() > 4 || controls.iter().any(|p| p[3] <= 0) {
        return Err("unsupported_boundary_curve");
    }
    let max_weight = controls.iter().map(|p| p[3]).max().unwrap() as f64;
    let mut poles = Vec::new();
    let mut weights = Vec::new();
    for p in controls {
        poles.push(cartesian(
            w,
            [
                p[0] as f64 / p[3] as f64,
                p[1] as f64 / p[3] as f64,
                p[2] as f64 / p[3] as f64,
            ],
        ));
        weights.push(number(p[3] as f64 / max_weight));
    }
    let degree = controls.len() - 1;
    let pole_list = poles
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let entity = format!("(BOUNDED_CURVE() B_SPLINE_CURVE({degree},({pole_list}),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS({},.PIECEWISE_BEZIER_KNOTS.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({})) REPRESENTATION_ITEM(''))", knots(degree), weights.join(","));
    Ok(w.entity(entity))
}
fn surface_entity(
    w: &mut Writer,
    face: &crate::boundary_shell::ShellFace,
    compact: bool,
) -> Result<usize, &'static str> {
    if compact && face.feature.starts_with("bore/") && face.feature.ends_with("/wall") {
        if let Ok(Some(c)) = crate::step_cylinder::cylinder_surface(face) {
            let real =
                |v: crate::constraint_system::Ratio| v.numerator() as f64 / v.denominator() as f64;
            let origin = cartesian(w, c.center.map(real));
            let axis = w.entity(format!("DIRECTION('',{})", coord(c.axis.map(f64::from))));
            let reference = w.entity(format!("DIRECTION('',{})", coord(c.reference.map(real))));
            let placement = w.entity(format!(
                "AXIS2_PLACEMENT_3D('',#{origin},#{axis},#{reference})"
            ));
            return Ok(w.entity(format!(
                "CYLINDRICAL_SURFACE('',#{placement},{})",
                number(real(c.radius))
            )));
        }
    }
    let [du, dv] = face.degrees;
    if du > 3 || dv > 3 || face.controls.len() != (du + 1) * (dv + 1) {
        return Err("unsupported_boundary_surface");
    }
    let controls = face.controls.chunks(dv + 1).collect::<Vec<_>>();
    if face.controls.iter().any(|p| p[3] <= 0) {
        return Err("invalid_boundary_weight");
    }
    let max_weight = face.controls.iter().map(|p| p[3]).max().unwrap() as f64;
    let mut rows = Vec::new();
    let mut weight_rows = Vec::new();
    for row in controls {
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for p in row {
            poles.push(format!(
                "#{}",
                cartesian(
                    w,
                    [
                        p[0] as f64 / p[3] as f64,
                        p[1] as f64 / p[3] as f64,
                        p[2] as f64 / p[3] as f64
                    ]
                )
            ));
            weights.push(number(p[3] as f64 / max_weight));
        }
        rows.push(format!("({})", poles.join(",")));
        weight_rows.push(format!("({})", weights.join(",")));
    }
    let entity = format!("(BOUNDED_SURFACE() B_SPLINE_SURFACE({du},{dv},({}),.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS(({},{}) ,({},{}) ,(0.,1.),(0.,1.),.PIECEWISE_BEZIER_KNOTS.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({})) REPRESENTATION_ITEM('') SURFACE())",
        rows.join(","), du+1,du+1,dv+1,dv+1,weight_rows.join(","));
    Ok(w.entity(entity))
}

/// Serialize the audited boundary of `model` as a STEP AP214 manifold solid.
pub fn export_step(model: &Json) -> Result<String, &'static str> {
    let shell = boundary_shell(model)?;
    let source = model.field("source_model").unwrap_or(model);
    let compact = source
        .field("partition_plan")
        .and_then(|p| p.field("method"))
        .and_then(Json::text)
        == Ok(crate::orthogonal_sketch_json::COMPACT_METHOD);
    write_shell(shell, compact)
}

fn write_shell(shell: BoundaryShell, compact: bool) -> Result<String, &'static str> {
    let mut w = Writer {
        next: 11,
        data: String::new(),
    };
    let mut vertex_entities = Vec::new();
    for &h in &shell.vertices {
        let p = cartesian(&mut w, point4(h)?);
        vertex_entities.push(w.entity(format!("VERTEX_POINT('',#{p})")));
    }
    let mut edge_refs = Vec::new();
    for edge in &shell.edges {
        let curve = curve_entity(&mut w, &edge.controls)?;
        let edge_curve = w.entity(format!(
            "EDGE_CURVE('',#{},#{},#{curve},.T.)",
            vertex_entities[edge.vertices[0]], vertex_entities[edge.vertices[1]]
        ));
        edge_refs.push(edge_curve);
    }
    let mut face_refs = Vec::new();
    let mut is_void = vec![false; shell.faces.len()];
    for component in shell.components.iter().skip(1) {
        for &face in component {
            is_void[face] = true;
        }
    }
    for (face_index, face) in shell.faces.iter().enumerate() {
        // The oriented void wrapper reverses an independently outward-facing
        // closed shell. Source loops are already outward from the material.
        let reverse = is_void[face_index];
        let same_sense = face.same_sense ^ reverse;
        if face.loop_edges.len() < 3 {
            return Err("invalid_boundary_loop");
        }
        let mut oriented = Vec::new();
        for i in 0..face.loop_edges.len() {
            let (edge, forward) = if reverse {
                let (e, f) = face.loop_edges[face.loop_edges.len() - 1 - i];
                (e, !f)
            } else {
                face.loop_edges[i]
            };
            oriented.push(w.entity(format!(
                "ORIENTED_EDGE('',*,*,#{},.{}.)",
                edge_refs[edge],
                if forward { "T" } else { "F" }
            )));
        }
        let loop_id = w.entity(format!(
            "EDGE_LOOP('',({}))",
            oriented
                .iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        ));
        let bound = w.entity(format!(
            "FACE_BOUND('',#{loop_id},.{}.)",
            if same_sense { "T" } else { "F" }
        ));
        let surface = surface_entity(&mut w, face, compact)?;
        let label = step_text(&format!("{}:{}:{}", face.feature, face.block, face.side));
        face_refs.push(w.entity(format!(
            "ADVANCED_FACE('{label}',(#{bound}),#{surface},.{}.)",
            if same_sense { "T" } else { "F" }
        )));
    }
    let closed = shell
        .components
        .iter()
        .map(|component| {
            let members = component
                .iter()
                .map(|&i| format!("#{}", face_refs[i]))
                .collect::<Vec<_>>()
                .join(",");
            w.entity(format!("CLOSED_SHELL('',({members}))"))
        })
        .collect::<Vec<_>>();
    let outer = closed[0];
    let solid = if closed.len() == 1 {
        w.entity(format!(
            "MANIFOLD_SOLID_BREP('SolidCraft boundary',#{outer})"
        ))
    } else {
        let voids = closed[1..]
            .iter()
            .map(|c| {
                let id = w.entity(format!("ORIENTED_CLOSED_SHELL('',*,#{c},.F.)"));
                format!("#{id}")
            })
            .collect::<Vec<_>>()
            .join(",");
        w.entity(format!(
            "BREP_WITH_VOIDS('SolidCraft boundary',#{outer},({voids}))"
        ))
    };
    let origin = cartesian(&mut w, [0.0; 3]);
    let z_axis = w.entity("DIRECTION('',(0.,0.,1.))".into());
    let x_axis = w.entity("DIRECTION('',(1.,0.,0.))".into());
    let placement = w.entity(format!(
        "AXIS2_PLACEMENT_3D('',#{origin},#{z_axis},#{x_axis})"
    ));
    let context_id = w.next;
    let units_id = context_id + 3;
    let uncertainty_id = units_id + 3;
    writeln!(w.data, "#{context_id} = ( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{uncertainty_id})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{units_id},#{},#{})) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') );", units_id+1,units_id+2).unwrap();
    writeln!(
        w.data,
        "#{units_id} = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );"
    )
    .unwrap();
    writeln!(
        w.data,
        "#{} = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );",
        units_id + 1
    )
    .unwrap();
    writeln!(
        w.data,
        "#{} = ( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() );",
        units_id + 2
    )
    .unwrap();
    writeln!(w.data, "#{uncertainty_id} = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{units_id},'distance_accuracy_value','confusion accuracy');").unwrap();
    let product_category = uncertainty_id + 1;
    writeln!(
        w.data,
        "#{product_category} = PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#7));"
    )
    .unwrap();
    let header = "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('SolidCraft rational boundary'),'2;1');\nFILE_NAME('SolidCraft boundary','@UTC@',('SolidCraft'),('SolidCraft'),'SolidCraftTrivariate','SolidCraft','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n#1 = APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#2);\n#2 = APPLICATION_CONTEXT('core data for automotive mechanical design processes');\n#3 = SHAPE_DEFINITION_REPRESENTATION(#4,#10);\n#4 = PRODUCT_DEFINITION_SHAPE('','',#5);\n#5 = PRODUCT_DEFINITION('design','',#6,#9);\n#6 = PRODUCT_DEFINITION_FORMATION('','',#7);\n#7 = PRODUCT('SolidCraft rational boundary','SolidCraft rational boundary','',(#8));\n#8 = PRODUCT_CONTEXT('',#2,'mechanical');\n#9 = PRODUCT_DEFINITION_CONTEXT('part definition',#2,'design');\n";
    let mut output = header.replace("@UTC@", &exchange_timestamp()?);
    writeln!(
        output,
        "#10 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{placement},#{solid}),#{context_id});"
    )
    .unwrap();
    output.push_str(&w.data);
    output.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    Ok(output)
}
