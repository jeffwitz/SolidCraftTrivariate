//! Thin filesystem adapter for the versioned native construction service.

use solidcraft_trivariate::interchange::{audit, build, SECTION_FAMILIES};
use solidcraft_trivariate::json::{object, parse, Json};
use solidcraft_trivariate::RationalBezierMap;
use std::env;
use std::fs::{self, OpenOptions};
use std::io::{ErrorKind, Write};
use std::path::Path;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};

fn read(path: &str) -> Result<Json, String> {
    parse(&fs::read_to_string(path).map_err(|e| e.to_string())?).map_err(str::to_owned)
}

fn publish(target: &Path, contents: &str) -> std::io::Result<()> {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = target
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(dir)?;
    loop {
        let name = format!(
            ".trivariate-{}-{}.tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        );
        let temp = dir.join(name);
        if temp == target {
            continue;
        }
        let mut file = match OpenOptions::new().write(true).create_new(true).open(&temp) {
            Err(e) if e.kind() == ErrorKind::AlreadyExists => continue,
            other => other?,
        };
        let result = (|| {
            file.write_all(contents.as_bytes())?;
            if !contents.ends_with('\n') {
                file.write_all(b"\n")?;
            }
            file.sync_all()?;
            drop(file);
            fs::rename(&temp, target)
        })();
        if result.is_err() {
            let _ = fs::remove_file(temp);
        }
        return result;
    }
}

const COMMANDS_AND_EVALUATORS: &str = r#"{"commands":["build","audit","assemble-patches","reduce-axial-knots","plan-sketch","crossed-boolean","evaluate","capabilities","export-step","certify-map","export-analysis","export-sketch-analysis","export-nurbs-space","export-adaptive","inspect-travel","solve-job","differentiate-job","import-solidcraft","import-solidcraft-sketch","export-field-vtu","export-geometry-vtu","edit-bore-intervals","edit-s3-stock-rim","edit-fieldcraft-bore-rims"],"evaluators":{"rational_bezier_tensor":{"degrees":[[1,1,1],[1,1,2],[1,1,3],[1,2,1],[1,2,2],[1,2,3],[1,3,1],[1,3,2],[1,3,3],[2,1,1],[2,1,2],[2,1,3],[2,2,1],[2,2,2],[2,2,3],[2,3,1],[2,3,2],[2,3,3],[3,1,1],[3,1,2],[3,1,3],[3,2,1],[3,2,2],[3,2,3],[3,3,1],[3,3,2],[3,3,3]],"domain":[[0,1],[0,1],[0,1]]}}}"#;

fn capabilities() -> Result<Json, String> {
    let mut result = parse(COMMANDS_AND_EVALUATORS)?;
    let Json::Object(fields) = &mut result else {
        return Err("invalid capabilities".into());
    };
    let families = [
        "affine_pair_v1",
        "annular_pair_v1",
        "filleted_bracket_v1",
        "movable_bore_v1",
        "pocket_depth_v1",
        "exact_bore_v1",
        "row_bores_v1",
        "overlap_bore_v1",
        "axial_boolean_v1",
        "crescent_pocket_v1",
        "parameterized_bracket_v1",
        "travelling_bore_v1",
        "travelling_bore_bounded_v1",
        "ruled_section_sweep_v1",
        "coaxial_feature_program_v1",
        "circular_bend_v1",
        "orthogonal_feature_program_v1",
        "polygon_prism_v1",
        "polygonal_perforated_prism_v1",
        "polygon_circular_bores_v1",
        "polygon_interval_bores_v1",
        "polygon_counterbores_v1",
        "polygon_countersinks_v1",
        "nurbs_multiloop_v1",
        "nurbs_edge_loops_v1",
        "g2_cross_plate_v1",
        "s3_face_pocket_v1",
        "primitive_sweep_v1",
        "nurbs_sketch_strip_v1",
        "sketch_layer_partition_v1",
        "sketch_boolean_layers_v1",
        "orthogonal_sketch_features_v1",
        "nurbs_cross_bores_v1",
        "nurbs_cross_bores_bounded_v1",
        "nurbs_cross_oval_v1",
        "rational_sector_prism_v1",
        "affine_image_v1",
    ]
    .into_iter()
    .chain(SECTION_FAMILIES)
    .map(|s| Json::String(s.to_owned()))
    .collect();
    fields.insert("families".into(), Json::Array(families));
    Ok(result)
}

fn diagnostic(status: &str, reason: &str) -> Json {
    object([
        ("status", Json::String(status.into())),
        ("reason", Json::String(reason.into())),
    ])
}

fn map_data<T: std::str::FromStr>(map: &Json) -> Result<([usize; 3], Vec<[T; 4]>), &'static str> {
    map.keys(&["degrees", "controls"])?;
    let degrees = map
        .field("degrees")?
        .items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "expected three degrees")?;
    let mut controls = Vec::new();
    for point in map.field("controls")?.items()? {
        controls.push(
            point
                .items()?
                .iter()
                .map(Json::number)
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .map_err(|_| "expected four homogeneous coordinates")?,
        );
    }
    Ok((degrees, controls))
}

fn certify_map(map: &Json, leaves: &str, depth: &str) -> Result<Json, &'static str> {
    use solidcraft_trivariate::spatial_certificate::{
        certify_spatial, CertificateBudget, CertificateStatus,
    };
    let (degrees, controls) = map_data::<i64>(map)?;
    let budget = CertificateBudget {
        max_leaves: leaves.parse().map_err(|_| "invalid leaf budget")?,
        max_depth: depth.parse().map_err(|_| "invalid depth budget")?,
    };
    let result = certify_spatial(degrees, &controls, budget);
    let numeric = |x: f64| Json::Number(x.to_string());
    let integer = |x: usize| Json::Number(x.to_string());
    let witness = result.witness.map(|w| {
        object([
            (
                "axes",
                Json::Array(
                    w.domain
                        .axes
                        .map(|interval| {
                            object([
                                ("numerator", integer(interval.numerator as usize)),
                                ("depth", integer(interval.depth.into())),
                            ])
                        })
                        .into(),
                ),
            ),
            (
                "upper_corner",
                Json::Array(w.upper_corner.map(Json::Bool).into()),
            ),
            (
                "determinant_upper_bound_mm3",
                numeric(w.determinant_upper_bound_mm3),
            ),
        ])
    });
    let status = match result.status {
        CertificateStatus::Valid => "valid",
        CertificateStatus::Invalid => "invalid",
        CertificateStatus::Indeterminate => "indeterminate",
    };
    Ok(object([
        ("status", Json::String(status.into())),
        (
            "quantity",
            Json::String("local Jacobian determinant on the closed unit cube".into()),
        ),
        (
            "method",
            Json::String("outward_interval_tensor_Bernstein".into()),
        ),
        (
            "determinant_lower_bound_mm3",
            result
                .determinant_lower_bound_mm3
                .map_or(Json::Null, numeric),
        ),
        ("witness", witness.unwrap_or(Json::Null)),
        ("leaves_used", integer(result.leaves_used)),
        (
            "max_depth_reached",
            integer(result.max_depth_reached.into()),
        ),
        (
            "reason",
            result.reason.map_or(Json::Null, |s| Json::String(s.into())),
        ),
    ]))
}

fn inspect_travel(center: &str) -> Result<Json, &'static str> {
    use solidcraft_trivariate::travelling_bore_events::analyze_travelling_bore;
    let r = analyze_travelling_bore(
        center
            .parse()
            .map_err(|_| "invalid integer center in micrometres")?,
    )?;
    let a = r.areas;
    let n = |x: f64| Json::Number(x.to_string());
    Ok(object([
        ("status", Json::String("evaluated".into())),
        ("center_x_um", Json::Number(r.center_um.to_string())),
        ("ideal_volume_mm3", n(r.volume_mm3)),
        ("dvolume_dcenter_mm2", n(r.dvolume_dcenter_mm2)),
        (
            "domain_signature",
            object([
                (
                    "bore_plate",
                    Json::String(format!("{:?}", r.physical.bore_plate)),
                ),
                (
                    "bore_pocket",
                    Json::String(format!("{:?}", r.physical.bore_pocket)),
                ),
            ]),
        ),
        (
            "partition_locally_fixed",
            Json::Bool(r.representation.locally_fixed_partition),
        ),
        (
            "center_sides_in_x_order",
            Json::Array(
                r.representation
                    .center_sides
                    .map(|s| Json::Number(s.to_string()))
                    .into(),
            ),
        ),
        (
            "integer_cut_ordinates",
            Json::Bool(r.representation.integer_cut_ordinates),
        ),
        (
            "geometry_derivative_status",
            Json::String("unavailable".into()),
        ),
        (
            "tangencies",
            Json::Array(
                r.events
                    .into_iter()
                    .flatten()
                    .map(|e| Json::String(e.feature.into()))
                    .collect(),
            ),
        ),
        (
            "ideal_boundary_areas_mm2",
            object([
                ("lower_bore", n(a.bore_wall_below_pocket_mm2)),
                ("upper_bore", n(a.bore_wall_above_pocket_mm2)),
                ("pocket_walls", n(a.pocket_vertical_walls_mm2)),
                ("pocket_floor", n(a.pocket_floor_mm2)),
                ("top", n(a.plate_top_mm2)),
                ("bottom", n(a.plate_bottom_mm2)),
                ("outer_vertical", n(a.plate_outer_vertical_walls_mm2)),
            ]),
        ),
    ]))
}

fn run(args: &[String]) -> Result<Json, String> {
    match args {
        [cmd] if cmd=="capabilities"=>capabilities(),
        [cmd, path, state, target] if matches!(cmd.as_str(),"edit-fieldcraft-bore-rims"|"edit-s3-stock-rim") => {
            use solidcraft_trivariate::{fieldcraft_rims, s3_stock_rim};
            let source = read(path)?;
            let path = Path::new(target);
            let step = match path.extension().and_then(|s|s.to_str()) {
                Some("json") => false, Some("step") => true,
                _ => return Err("expected .json or .step output".into()),
            };
            let contents = if cmd == "edit-fieldcraft-bore-rims" {
                if step { fieldcraft_rims::export_step(&source,state)? }
                else { fieldcraft_rims::construct(&source,state)?.to_string() }
            } else {
                let state = match state.as_str() {
                    "sharp" => s3_stock_rim::RimState::Sharp,
                    "chamfer" => s3_stock_rim::RimState::Chamfer,
                    "fillet" => s3_stock_rim::RimState::Fillet,
                    _ => return Err("expected sharp, chamfer or fillet".into()),
                };
                if step { s3_stock_rim::export_step(&source,state)? }
                else { s3_stock_rim::construct(&source,state)?.to_string() }
            };
            publish(path,&contents).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("rim_state",Json::String(state.to_owned()))]))
        }
        [cmd,path,op,target] if cmd=="crossed-boolean"=>{
            let model=solidcraft_trivariate::g2_cross_plate::boolean::report(&read(path)?,op)?;
            publish(Path::new(target),&model.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("blocks",Json::Number(model.field("blocks")?.items()?.len().to_string()))]))
        }
        [cmd, model, layout, target] if cmd == "assemble-patches" => {
            let model = read(model)?;
            audit(&model)?;
            let cover = solidcraft_trivariate::nurbs_tensor_assembly::assemble_cover(
                &model, &read(layout)?,
            )?;
            publish(Path::new(target), &cover.to_string()).map_err(|e| e.to_string())?;
            Ok(object([("patch_count", Json::Number(cover.field("patches")?.items()?.len().to_string()))]))
        }
        [cmd, path, target] if cmd == "reduce-axial-knots" => {
            let compact = solidcraft_trivariate::nurbs_knot_reduction::reduce_patch_document(
                &read(path)?,
            )?;
            publish(Path::new(target), &compact.to_string()).map_err(|e| e.to_string())?;
            Ok(object([("patch_count", Json::Number(compact.field("patches")?.items()?.len().to_string()))]))
        }
        [cmd,path] if cmd=="plan-sketch"=>Ok(solidcraft_trivariate::interchange::plan_sketch(&read(path)?)?),
        [cmd,center] if cmd=="inspect-travel"=>Ok(inspect_travel(center)?),
        [cmd,path,target] if cmd=="export-nurbs-space"=>{
            let package=solidcraft_trivariate::nurbs_space::export(&read(path)?)?;
            publish(Path::new(target),&package.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("scalar_ndof",package.field("scalar_ndof")?.clone())]))
        }
        [cmd,path] if cmd=="audit"=>Ok(audit(&read(path)?).unwrap_or_else(|reason| diagnostic(if reason.ends_with("_indeterminate_bernstein") {"indeterminate"} else {"invalid"},reason))),
        [cmd,recipe,target] if cmd=="build"=>{
            let model=build(&read(recipe)?)?;
            let report=audit(&model)?;
            if report.field("status")?.text()?=="valid" {publish(Path::new(target),&model.to_string()).map_err(|e|e.to_string())?;}
            Ok(report)
        }
        [cmd,source,target] if cmd=="import-solidcraft" || cmd=="import-solidcraft-sketch"=>{
            let input=read(source)?;
            let model=if cmd=="import-solidcraft-sketch" {
                solidcraft_trivariate::solidcraft_sketch_features::import_sketch_timeline(&input)?
            } else {
                solidcraft_trivariate::solidcraft_bridge::import_script(&input)?
            };
            let report=audit(&model)?;
            publish(Path::new(target),&model.to_string()).map_err(|e|e.to_string())?;
            Ok(report)
        }
        [cmd,model,job,direction,target] if cmd=="differentiate-job"=>{
            let result=solidcraft_trivariate::mechanical_job::differentiate_job(&read(model)?,&read(job)?,&read(direction)?)?;
            publish(Path::new(target),&result.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("evaluated".into())),("schema",result.field("schema")?.clone())]))
        }
        [cmd,model,edits,target] if cmd=="edit-bore-intervals"=>{
            let (candidate, report) = solidcraft_trivariate::interval_bores_edit::edit(&read(model)?,&read(edits)?)?;
            publish(Path::new(target),&candidate.to_string()).map_err(|e|e.to_string())?;
            Ok(report)
        }
        [cmd,model,job,target] if cmd=="solve-job"=>{
            let solution=solidcraft_trivariate::mechanical_job::solve_job(&read(model)?,&read(job)?)?;
            publish(Path::new(target),&solution.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("solved".into())),("schema",solution.field("schema")?.clone()),
                ("relative_residual",solution.field("relative_residual")?.clone())]))
        }
        [cmd,path,samples,target] if cmd=="export-geometry-vtu"=>{
            let n=samples.parse().map_err(|_| "invalid geometry subdivision count")?;
            let vtu=solidcraft_trivariate::geometry_vtu::export_geometry_vtu(&read(path)?,n)?;
            publish(Path::new(target),&vtu).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("format",Json::String("VTK XML UnstructuredGrid".into()))]))
        }
        [cmd,path,samples,target] if cmd=="export-field-vtu"=>{
            let n=samples.parse().map_err(|_| "invalid field subdivision count")?;
            let vtu=solidcraft_trivariate::field_vtu::export_field_vtu(&read(path)?,n)?;
            publish(Path::new(target),&vtu).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("format",Json::String("VTK XML UnstructuredGrid".into()))]))
        }
        [cmd,path,order,target] if cmd=="export-analysis"=>{
            let package=solidcraft_trivariate::solver_package::export_analysis(&read(path)?,order.parse().map_err(|_| "invalid quadrature order")?)?;
            publish(Path::new(target),&package.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("schema",package.field("schema")?.clone())]))
        }
        [cmd,path,order,target] if cmd=="export-sketch-analysis"=>{
            let package=solidcraft_trivariate::solver_package::export_sketch_analysis(&read(path)?,order.parse().map_err(|_| "invalid quadrature order")?)?;
            publish(Path::new(target),&package.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("schema",package.field("schema")?.clone())]))
        }
        [cmd,model,plan,order,target] if cmd=="export-adaptive"=>{
            let package=solidcraft_trivariate::adaptive_package::export_adaptive_package(&read(model)?,&read(plan)?,order.parse().map_err(|_|"invalid quadrature order")?)?;
            publish(Path::new(target),&package.to_string()).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("schema",package.field("schema")?.clone()),("package_extension",package.field("package_extension")?.clone())]))
        }
        [cmd,path,target] if cmd=="export-step"=>{
            let model=read(path)?;
            let step=solidcraft_trivariate::step_boundary::export_step(&model)?;
            publish(Path::new(target), &step).map_err(|e|e.to_string())?;
            Ok(object([("status",Json::String("exported".into())),("format",Json::String("STEP AP214 boundary".into()))]))
        }
        [cmd,path,leaves,depth] if cmd=="certify-map"=>Ok(certify_map(&read(path)?,leaves,depth)?),
        [cmd,path,u,v,w] if cmd=="evaluate"=>{
            let (degrees,controls)=map_data::<f64>(&read(path)?)?;
            let map=match RationalBezierMap::new(degrees,controls) {
                Ok(map)=>map,
                Err(solidcraft_trivariate::MapError::UnsupportedDegrees)=>return Ok(diagnostic("unsupported","unsupported degree family")),
                Err(error)=>return Err(error.to_string()),
            };
            let parameter=|s:&str|s.parse::<f64>().map_err(|_|"invalid parameter");
            let res=map.evaluate([parameter(u)?,parameter(v)?,parameter(w)?]).map_err(|e|e.to_string())?;
            let vector=|a:[f64;3]|Json::Array(a.into_iter().map(|x|Json::Number(x.to_string())).collect());
            Ok(object([("status",Json::String("evaluated".into())),("position",vector(res.position)),
                ("jacobian_columns",Json::Array(res.jacobian_columns.into_iter().map(vector).collect()))]))
        }
        _=>Err("usage: solidcraft-trivariate assemble-patches MODEL LAYOUT OUTPUT | reduce-axial-knots PATCHES OUTPUT | capabilities | crossed-boolean RECIPE OP OUTPUT | plan-sketch RECIPE | build RECIPE OUTPUT | audit MODEL | evaluate MAP U V W | export-step MODEL OUTPUT | certify-map MAP MAX_LEAVES MAX_DEPTH | export-analysis MODEL GAUSS_ORDER OUTPUT | export-sketch-analysis MODEL GAUSS_ORDER OUTPUT | export-nurbs-space MODEL OUTPUT | export-adaptive MODEL PLAN GAUSS_ORDER OUTPUT | import-solidcraft SCRIPT OUTPUT | import-solidcraft-sketch TIMELINE OUTPUT | solve-job MODEL JOB OUTPUT | differentiate-job MODEL JOB DIRECTION OUTPUT | export-field-vtu RESULT SUBDIVISIONS OUTPUT | export-geometry-vtu MODEL SUBDIVISIONS OUTPUT | edit-bore-intervals MODEL EDITS OUTPUT | edit-s3-stock-rim MODEL sharp|chamfer|fillet OUTPUT.json|OUTPUT.step | edit-fieldcraft-bore-rims MODEL STATE OUTPUT.json|OUTPUT.step | inspect-travel CENTER_X_UM".into()),
    }
}

fn main() -> ExitCode {
    let args: Vec<_> = env::args().skip(1).collect();
    let report = run(&args).unwrap_or_else(|e| {
        object([
            (
                "status",
                Json::String(
                    if e.ends_with("_indeterminate_bernstein") {
                        "indeterminate"
                    } else {
                        "error"
                    }
                    .into(),
                ),
            ),
            ("reason", Json::String(e)),
        ])
    });
    println!("{report}");
    match report.field("status").and_then(Json::text) {
        Ok("valid" | "evaluated" | "exported" | "solved" | "updated") | Err(_) => ExitCode::SUCCESS,
        _ => ExitCode::FAILURE,
    }
}
