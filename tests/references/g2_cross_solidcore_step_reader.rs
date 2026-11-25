use solidcore::Kernel;
use solidcore::mass::integrator::properties::brep_volume_default;
use solidcore::query::brep::BrepClassifier;
use solidcore::math::Point3;
use solidcore::intersect::point_class::PointClassType;
use solidcore::error::EntityRef;
use solidcore_step::{p21, topology};
use solidcore_convert::brep::to_kernel;
use std::{env, fs, process};

fn run() -> Result<(), String> {
    let path = env::args().nth(1).ok_or("usage: g2-trivariate-solidcore-reader STEP")?;
    let text = fs::read_to_string(&path).map_err(|e| format!("read STEP: {e}"))?;
    let exchange = p21::parse(&text).map_err(|e| format!("SolidCore Part 21 parse: {e}"))?;
    let solids = topology::read_solids(&exchange);
    println!("exchange_entities={}", exchange.entities.len());
    println!("read_solids={}", solids.solids.len());
    println!("read_failures={}", solids.failures.len());
    for failure in &solids.failures {
        println!("read_failure={} {}", failure.item, failure.message);
    }
    if solids.solids.len() != 1 || !solids.failures.is_empty() {
        return Err("expected exactly one fully read solid".into());
    }
    let imported = &solids.solids[0];
    let record = &imported.brep;
    println!("name={}", imported.name);
    println!("sheet={}", imported.sheet);
    println!("record_regions={}", record.regions.len());
    println!("record_shells={}", record.shells.len());
    println!("record_faces={}", record.faces.len());
    println!("record_edges={}", record.edges.len());
    println!("record_vertices={}", record.vertices.len());
    println!("record_tolerance_mm={:.17}", record.tolerance);
    println!("reader_notes={}", imported.notes.join(" | "));

    let mut kernel = Kernel::default();
    let (mut brep, _) = to_kernel(&mut kernel, record, false)
        .map_err(|e| format!("SolidCore record-to-kernel conversion: {e}"))?;
    println!("kernel_faces={}", brep.face_list().len());
    println!("kernel_edges={}", brep.edge_list().len());
    println!("kernel_vertices={}", brep.vertex_list().len());
    let volume = brep_volume_default(&kernel, &mut brep, 1e-8)
        .map_err(|e| format!("SolidCore imported STEP volume: {e:?}"))?;
    println!("volume_mm3={volume:.17}");

    let probes = [
        ("material", [10.0, 10.0, 10.0]),
        ("periodic_bore_only", [75.0, 25.0, 2.0]),
        ("overlap_both_bores", [75.0, 25.0, 6.0]),
        ("rectangular_bore", [50.0, 25.0, 6.0]),
        ("stock_boundary", [0.0, 40.0, 10.0]),
    ];
    let mut classifier = BrepClassifier::of_brep(&brep);
    for (name, point) in probes {
        let class = classifier.point_3d_classify(
            &kernel, Point3::new(point[0], point[1], point[2]), 2e-7, true,
        ).map_err(|e| format!("SolidCore classify {name}: {e:?}"))?;
        let label = match (class.point_class(), class.object()) {
            (PointClassType::Region, Some(EntityRef::Region(r))) => {
                if brep.region(r).is_ok_and(|x| x.is_void == 0) { "Inside" } else { "Outside" }
            }
            (PointClassType::Region, _) => "RegionUnknown",
            (PointClassType::Face | PointClassType::Edge | PointClassType::Vertex, _) => "OnBoundary",
            (PointClassType::Nothing, _) => "Outside",
            _ => "Unclassified",
        };
        println!("classification_{name}={label}");
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        process::exit(1);
    }
}
