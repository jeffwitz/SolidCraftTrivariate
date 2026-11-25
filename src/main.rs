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

fn publish(target: &Path, model: &Json) -> std::io::Result<()> {
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
            writeln!(file, "{model}")?;
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

const COMMANDS_AND_EVALUATORS: &str = r#"{"commands":["build","audit","evaluate","capabilities"],"evaluators":{"rational_bezier_tensor":{"degrees":[[1,1,1],[2,1,1]],"domain":[[0,1],[0,1],[0,1]]}}}"#;

fn capabilities() -> Result<Json,String> {
    let mut result = parse(COMMANDS_AND_EVALUATORS)?;
    let Json::Object(fields) = &mut result else {return Err("invalid capabilities".into())};
    let families = ["affine_pair_v1", "annular_pair_v1", "filleted_bracket_v1", "movable_bore_v1", "pocket_depth_v1"]
        .into_iter().chain(SECTION_FAMILIES).map(|s| Json::String(s.to_owned())).collect();
    fields.insert("families".into(), Json::Array(families));
    Ok(result)
}

fn diagnostic(status: &str, reason: &str) -> Json {
    object([
        ("status", Json::String(status.into())),
        ("reason", Json::String(reason.into())),
    ])
}

fn run(args: &[String]) -> Result<Json, String> {
    match args {
        [cmd] if cmd=="capabilities"=>capabilities(),
        [cmd,path] if cmd=="audit"=>Ok(audit(&read(path)?).unwrap_or_else(|reason| diagnostic("invalid",reason))),
        [cmd,recipe,target] if cmd=="build"=>{
            let model=build(&read(recipe)?)?;
            let report=audit(&model)?;
            if report.field("status")?.text()?=="valid" {publish(Path::new(target),&model).map_err(|e|e.to_string())?;}
            Ok(report)
        }
        [cmd,path,u,v,w] if cmd=="evaluate"=>{
            let m=read(path)?;m.keys(&["degrees","controls"])?;
            let degrees=m.field("degrees")?.items()?.iter().map(Json::number::<usize>).collect::<Result<Vec<_>,_>>()?
                .try_into().map_err(|_|"expected three degrees")?;
            let mut controls=vec![];
            for p in m.field("controls")?.items()? {
                controls.push(p.items()?.iter().map(Json::number::<f64>).collect::<Result<Vec<_>,_>>()?
                    .try_into().map_err(|_|"expected four homogeneous coordinates")?);
            }
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
        _=>Err("usage: solidcraft-trivariate capabilities | build RECIPE OUTPUT | audit MODEL | evaluate MAP U V W".into()),
    }
}

fn main() -> ExitCode {
    let args: Vec<_> = env::args().skip(1).collect();
    let report = run(&args).unwrap_or_else(|e| {
        object([
            ("status", Json::String("error".into())),
            ("reason", Json::String(e)),
        ])
    });
    println!("{report}");
    match report.field("status").and_then(Json::text) {
        Ok("valid" | "evaluated") | Err(_) => ExitCode::SUCCESS,
        _ => ExitCode::FAILURE,
    }
}
