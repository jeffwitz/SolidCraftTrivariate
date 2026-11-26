//! Rebuild the published FieldCraft ten-chart bracket and emit its fixed-
//! partition, quarter-millimetre two-bore chamfer as JSON under `build/`.
use solidcraft_trivariate::{
    bracket::BracketBlock,
    bracket_edge_edit::{edit_compact_bracket, BoreEdge, EdgeTreatment},
    interchange,
    json::{object, parse, Json},
};
use std::{fs, path::PathBuf};

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn s(value: impl Into<String>) -> Json {
    Json::String(value.into())
}
fn arr<T>(values: impl IntoIterator<Item = T>, f: impl FnMut(T) -> Json) -> Json {
    Json::Array(values.into_iter().map(f).collect())
}
fn ints<const N: usize>(value: &Json) -> [i64; N] {
    value
        .items()
        .unwrap()
        .iter()
        .map(|x| x.number().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}
fn source_blocks(model: &Json) -> Vec<BracketBlock> {
    model
        .field("blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|b| BracketBlock {
            id: b.field("id").unwrap().text().unwrap().to_owned(),
            degrees: b
                .field("degrees")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .map(|x| x.number().unwrap())
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
            controls: b
                .field("controls")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .map(ints)
                .collect(),
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let source_path = root.join("docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json");
    let source = parse(&fs::read_to_string(&source_path)?)?;
    let audit = interchange::audit(&source)?;
    let parents = source_blocks(&source);
    let bores = [
        BoreEdge {
            chart_prefix: "wing_east_",
            axis: 1,
            center_quarters: [160, 16, 80],
            radius_quarters: 16,
        },
        BoreEdge {
            chart_prefix: "wing_north_",
            axis: 0,
            center_quarters: [16, 128, 80],
            radius_quarters: 16,
        },
    ];
    let state = std::env::args().nth(1).unwrap_or_else(|| "chamfer".into());
    let treatment = match state.as_str() {
        "sharp" => EdgeTreatment::Sharp,
        "chamfer" => EdgeTreatment::Chamfer { width_quarters: 1 },
        "fillet" => EdgeTreatment::Fillet { radius_quarters: 1 },
        _ => return Err("expected sharp, chamfer, or fillet".into()),
    };
    let edited = edit_compact_bracket(&parents, bores, treatment)?;
    let output = object([
        ("schema", s("fixed_partition_bore_edge_v1")),
        ("state", s(state.clone())),
        (
            "source_model",
            s("docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"),
        ),
        ("source_audit_status", s(audit.field("status")?.text()?)),
        ("units", s("mm; recipe parameters use quarter-mm ticks")),
        ("logical_blocks", arr(edited.logical_blocks, s)),
        ("interfaces", n(edited.incidence.interfaces.len())),
        (
            "determinant_lower_bound_mm3",
            n(edited.determinant_lower_bound_mm3),
        ),
        (
            "singular_boundary_qualified",
            Json::Bool(edited.singular_boundary_qualified),
        ),
        (
            "jacobian_scope",
            s(if edited.singular_boundary_qualified {
                "zero at certified bore/cap contact rails; strictly positive away from those rails"
            } else {
                "strictly positive on every closed span"
            }),
        ),
        (
            "bores",
            arr(edited.bores, |b| {
                object([
                    ("chart_prefix", s(b.chart_prefix)),
                    ("axis", n(b.axis)),
                    ("center_quarters", arr(b.center_quarters, n)),
                    ("radius_quarters", n(b.radius_quarters)),
                    ("both_axial_ends", Json::Bool(true)),
                ])
            }),
        ),
        (
            "edge_parameter_quarters",
            n(match treatment {
                EdgeTreatment::Sharp => 0,
                EdgeTreatment::Chamfer { width_quarters } => width_quarters,
                EdgeTreatment::Fillet { radius_quarters } => radius_quarters,
            }),
        ),
        (
            "elements",
            arr(edited.elements, |e| {
                object([
                    ("logical_id", s(e.logical_id)),
                    ("span_index", n(e.span_index)),
                    ("station_quarters", arr(e.station_quarters, n)),
                    (
                        "parameter_interval_32nds",
                        arr(e.parameter_interval_32nds, n),
                    ),
                    ("degrees", arr(e.degrees, n)),
                    ("controls_xyzw", arr(e.controls, |p| arr(p, n))),
                ])
            }),
        ),
        (
            "fillet_status",
            s(solidcraft_trivariate::bracket_edge_edit::round_fillet_fixed_partition_status()),
        ),
    ]);
    let out = root.join(format!("build/bracket-edge-edit-{state}.json"));
    fs::create_dir_all(out.parent().unwrap())?;
    fs::write(&out, format!("{}\n", output))?;
    println!("{}", out.display());
    Ok(())
}
