//! Independently replay an X2 model and refinement plan, solve it natively,
//! and save a field separate from the exchange package.
use solidcraft_trivariate::{
    adaptive_package::parse_refinement_plan,
    adaptive_partition::AdaptivePartition,
    adaptive_space::build_adaptive_space,
    analysis_space::{build_analysis_space, parameterized_bracket_space},
    boundary_traction::{integrate_traction, LoadPatch},
    bracket::BracketBlock,
    bracket_audit::audit_bracket,
    bracket_parameters_audit::audit_bracket_parameters,
    bracket_solve::bracket_tip_displacement,
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::assemble,
    interchange,
    json::{object, parse, Json},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    path::Path,
};

fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn decode_blocks(model: &Json) -> Result<Vec<BracketBlock>, &'static str> {
    let mut result = Vec::new();
    for item in model.field("blocks")?.items()? {
        item.keys(&["id", "degrees", "controls"])?;
        let degrees = item
            .field("degrees")?
            .items()?
            .iter()
            .map(Json::number)
            .collect::<Result<Vec<usize>, _>>()?
            .try_into()
            .map_err(|_| "invalid bracket degrees")?;
        let mut controls = Vec::new();
        for p in item.field("controls")?.items()? {
            controls.push(
                p.items()?
                    .iter()
                    .map(Json::number::<i64>)
                    .collect::<Result<Vec<_>, _>>()?
                    .try_into()
                    .map_err(|_| "invalid bracket control")?,
            );
        }
        result.push(BracketBlock {
            id: item.field("id")?.text()?.to_owned(),
            degrees,
            controls,
        });
    }
    result.sort_by(|a, b| a.id.cmp(&b.id));
    Ok(result)
}
fn run(
    model_path: &str,
    plan_path: &str,
    order: usize,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let model = parse(&fs::read_to_string(model_path)?)?;
    let plan_json = parse(&fs::read_to_string(plan_path)?)?;
    let plan = parse_refinement_plan(&plan_json)?;
    let audit = interchange::audit(&model)?;
    if audit.field("status")?.text()? != "valid" {
        return Err("native adaptive solve requires an audited model".into());
    }
    let schema = model.field("schema")?.text()?;
    if schema != "filleted_bracket_v1" && schema != "parameterized_bracket_v1" {
        return Err("unsupported native adaptive model".into());
    }
    let blocks = decode_blocks(&model)?;
    let (certificate, parent) = if schema == "filleted_bracket_v1" {
        let c = audit_bracket(&blocks)?;
        let p = build_analysis_space(&blocks)?;
        (c, p)
    } else {
        let [rx, ry] = model.field("radii_sixteenths_mm")?.items()? else {
            return Err("invalid radius ticks".into());
        };
        let rx = rx.number()?;
        let ry = ry.number()?;
        let c = audit_bracket_parameters(rx, ry, &blocks)?.certificate;
        let p = parameterized_bracket_space(rx, ry, &blocks)?;
        (c, p)
    };
    let root_ids = parent
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut partition = AdaptivePartition::new(
        &root_ids,
        &certificate.incidence,
        plan.max_depth,
        plan.max_cells,
    )?;
    for step in &plan.steps {
        partition.refine(step)?;
    }
    let adaptive = build_adaptive_space(&parent, &partition)?;
    let space = &adaptive.space;
    let tags = certificate
        .design_faces
        .iter()
        .cloned()
        .collect::<BTreeMap<_, _>>();
    let exterior = partition
        .exterior_patches()
        .iter()
        .map(|patch| solidcraft_trivariate::incidence::BlockSide {
            block: adaptive.leaf_blocks[&patch.side.leaf].clone(),
            axis: patch.side.root_side.axis,
            upper: patch.side.root_side.upper,
        })
        .collect::<BTreeSet<_>>();
    let load_patches = partition
        .exterior_patches()
        .iter()
        .filter(|patch| {
            tags.get(&patch.side.root_side)
                .is_some_and(|feature| feature == "x60")
        })
        .map(|patch| {
            LoadPatch::whole(solidcraft_trivariate::incidence::BlockSide {
                block: adaptive.leaf_blocks[&patch.side.leaf].clone(),
                axis: patch.side.root_side.axis,
                upper: patch.side.root_side.upper,
            })
        })
        .collect::<Vec<_>>();
    let load = integrate_traction(space, &exterior, &load_patches, order, |_, _| {
        [0., -1000. / 320., 0.]
    })?;
    if (load.area_mm2 - 320.).abs() > 1e-8 {
        return Err("native loaded area mismatch".into());
    }
    let stiffness = assemble(space, 210000., 0.3, order)?;
    let prescribed = space.boundary_dofs["y50"]
        .iter()
        .map(|&d| (d, Ratio::zero()))
        .collect::<Vec<_>>();
    let lift = resolve(space.ndof, &adaptive.constraints, &prescribed)?;
    let solved = solve_constrained(
        &stiffness,
        [&lift, &lift, &lift],
        &load.nodal_force_n,
        1e-10,
        20000,
    )?;
    let tip = bracket_tip_displacement(space, &solved.displacement_mm)?;
    let result = object([
        (
            "schema",
            native_text("solidcraft_adaptive_native_solution_v1"),
        ),
        ("source_schema", native_text(schema)),
        ("source_model", native_text(model_path)),
        ("refinement_plan", plan_json),
        ("quadrature_order", n(order)),
        ("scalar_ndof", n(space.ndof)),
        ("full_vector_dofs", n(solved.displacement_mm.len())),
        (
            "independent_scalar_lift_masters",
            n(adaptive.lift.free_roots.len()),
        ),
        ("clamped_scalar_masters", n(lift.free_roots.len())),
        ("volume_mm3", n(stiffness.integrated_volume)),
        ("x60_area_mm2", n(load.area_mm2)),
        ("resultant_n", Json::Array(load.resultant_n.map(n).into())),
        ("strain_energy_n_mm", n(solved.strain_energy_n_mm)),
        ("true_residual_relative", n(solved.relative_residual)),
        ("iterations", n(solved.iterations)),
        ("tip_displacement_mm", Json::Array(tip.map(n).into())),
        (
            "full_interleaved_displacement_mm",
            Json::Array(solved.displacement_mm.iter().copied().map(n).collect()),
        ),
    ]);
    let target = Path::new(output);
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent)?;
    }
    let tmp = target.with_extension("tmp");
    fs::write(&tmp, format!("{result}\n"))?;
    fs::rename(tmp, target)?;
    Ok(())
}
fn native_text(s: &str) -> Json {
    Json::String(s.into())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    let [model, plan, order, output] = args.as_slice() else {
        return Err("usage: cargo run --release --example solve_adaptive_plan -- MODEL PLAN GAUSS_ORDER OUTPUT_JSON".into());
    };
    run(model, plan, order.parse()?, output)
}
