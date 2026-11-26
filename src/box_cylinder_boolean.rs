//! Exact fixed-chart subtraction of a centered vertical cylinder from a box.
//!
//! The public operands are distinct. A validated operand pair is converted to
//! an orthogonal sketch only inside this constructor so the familiar four
//! circular sectors and their exact wall seams are reused.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::enumerate_faces,
    orthogonal_sketch::{self, Bore, OrthogonalSketch},
    sketch_layers::SketchLoop,
    spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus},
};
use std::collections::BTreeSet;

pub const SUITE_FAMILY: &str = "box_cylinder_boolean_suite_v1";
pub const SUITE_OPERATIONS: [&str; 7] = [
    "box",
    "cylinder",
    "union",
    "intersection",
    "a_minus_b",
    "b_minus_a",
    "symmetric_difference",
];

#[derive(Debug, Clone, PartialEq)]
pub struct BooleanState {
    pub operation: &'static str,
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub volume_constant_mm3: i64,
    pub volume_pi_coefficient_mm3: i64,
    pub connected_components: usize,
}

fn box_block(stock: BoxOperand) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for x in [stock.min_mm[0], stock.max_mm[0]] {
        for y in [stock.min_mm[1], stock.max_mm[1]] {
            for z in [stock.min_mm[2], stock.max_mm[2]] {
                controls.push([i64::from(x), i64::from(y), i64::from(z), 1]);
            }
        }
    }
    BracketBlock {
        id: "box_stock".into(),
        degrees: [1, 1, 1],
        controls,
    }
}

/// Build one square-core disk layer. Its four circular outer rows are copied
/// from the matching bore charts so the artificial extrusion has exact traces.
fn disk_layer(
    layer: &str,
    z: [i32; 2],
    stock: BoxOperand,
    radius: i32,
    ring: &[BracketBlock],
) -> Result<Vec<BracketBlock>, &'static str> {
    let ring_by_id: std::collections::BTreeMap<_, _> = ring
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let sector_names = ["east", "north", "west", "south"];
    if ring.len() != 4 || ring.iter().any(|b| b.degrees != [2, 1, 1]) {
        return Err("four_sector_ring_layout_required");
    }
    let center = [
        (i64::from(stock.min_mm[0]) + i64::from(stock.max_mm[0])) / 2,
        (i64::from(stock.min_mm[1]) + i64::from(stock.max_mm[1])) / 2,
    ];
    let half_num = i64::from(radius); // half-width is radius/4; weight is four.
    let mut result = Vec::with_capacity(5);

    // One biquadratic square core: linear geometry elevated exactly to match
    // the quadratic inner sides of the four sectors.
    let mut core_controls = Vec::with_capacity(18);
    for i in 0..3 {
        for j in 0..3 {
            let x = 4_i64 * center[0] + (i as i64 - 1) * half_num;
            let y = 4_i64 * center[1] + (j as i64 - 1) * half_num;
            for level in z {
                core_controls.push([x, y, i64::from(level) * 4, 4]);
            }
        }
    }
    result.push(BracketBlock {
        id: format!("disk_{layer}_core"),
        degrees: [2, 2, 1],
        controls: core_controls,
    });

    for (sector, name) in sector_names.iter().enumerate() {
        let source = ring_by_id
            .get(format!("cylinder_{name}").as_str())
            .ok_or("missing_matching_ring_sector")?;
        // Construct each angular row as circle-to-core, with z fastest.
        let mut controls = Vec::with_capacity(12);
        for angle in 0..3 {
            let source_index = angle * 4 + 2;
            let [x, y, _, weight] = source.controls[source_index];
            let (inner_x, inner_y) = match sector {
                0 => (
                    4 * center[0] + half_num,
                    4 * center[1] + (angle as i64 - 1) * half_num,
                ),
                1 => (
                    4 * center[0] + (1 - angle as i64) * half_num,
                    4 * center[1] + half_num,
                ),
                2 => (
                    4 * center[0] - half_num,
                    4 * center[1] + (1 - angle as i64) * half_num,
                ),
                _ => (
                    4 * center[0] + (angle as i64 - 1) * half_num,
                    4 * center[1] - half_num,
                ),
            };
            for (px, py, w) in [(x, y, weight), (inner_x, inner_y, 4)] {
                for level in z {
                    controls.push([
                        px,
                        py,
                        i64::from(level)
                            .checked_mul(w)
                            .ok_or("world_z_translation_overflow")?,
                        w,
                    ]);
                }
            }
        }
        result.push(BracketBlock {
            id: format!("disk_{layer}_sector_{name}"),
            degrees: [2, 1, 1],
            controls,
        });
    }
    Ok(result)
}

fn joined_components(blocks: &[BracketBlock], incidence: &crate::incidence::Incidence) -> usize {
    let mut parent: Vec<usize> = (0..blocks.len()).collect();
    let index: std::collections::BTreeMap<_, _> = blocks
        .iter()
        .enumerate()
        .map(|(i, block)| (block.id.as_str(), i))
        .collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for interface in &incidence.interfaces {
        let Some(&a) = index.get(interface.first.block.as_str()) else {
            continue;
        };
        let Some(&b) = index.get(interface.second.block.as_str()) else {
            continue;
        };
        let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
        parent[ra] = rb;
    }
    (0..blocks.len())
        .map(|i| root(&mut parent, i))
        .collect::<BTreeSet<_>>()
        .len()
}

fn certify_boolean_blocks(
    blocks: &[BracketBlock],
) -> Result<(BracketCertificate, usize), &'static str> {
    let incidence = enumerate_faces(blocks)?;
    let mut lower = f64::INFINITY;
    for block in blocks {
        let certificate =
            certify_spatial(block.degrees, &block.controls, CertificateBudget::default());
        if certificate.status != CertificateStatus::Valid {
            return Err("boolean_chart_jacobian_not_certified");
        }
        lower = lower.min(
            certificate
                .determinant_lower_bound_mm3
                .ok_or("missing_jacobian_bound")?,
        );
    }
    if !lower.is_finite() || lower <= 0.0 {
        return Err("boolean_chart_jacobian_not_certified");
    }
    let design_faces = incidence
        .exterior
        .iter()
        .map(|side| {
            (
                side.clone(),
                format!("boundary_{}_{}_{}", side.block, side.axis, side.upper),
            )
        })
        .collect();
    let components = joined_components(blocks, &incidence);
    Ok((
        BracketCertificate {
            incidence,
            determinant_lower_bound_mm3: lower,
            radial_error_upper_mm: 0.0,
            design_faces,
        },
        components,
    ))
}

/// Build a regularized Boolean state for a centered cylinder extending through both stock caps.
pub fn build_boolean_state(
    operation: &str,
    stock: BoxOperand,
    tool: ThroughCylinderOperand,
) -> Result<BooleanState, &'static str> {
    let operation = SUITE_OPERATIONS
        .iter()
        .copied()
        .find(|op| *op == operation)
        .ok_or("unsupported_box_cylinder_boolean_operation")?;
    let bore = subtract_through_cylinder(stock, tool)?;
    let stock_volume = bore.volume_constant_mm3;
    let overlap = -bore.volume_pi_coefficient_mm3;
    let tool_height = i64::from(tool.interval_mm[1]) - i64::from(tool.interval_mm[0]);
    let cylinder = (i64::from(tool.radius_mm).pow(2))
        .checked_mul(tool_height)
        .ok_or("cylinder_volume_overflow")?;
    let outside = cylinder - overlap;
    let ring = bore.blocks;
    let mut blocks = match operation {
        "box" => vec![box_block(stock)],
        "union" | "a_minus_b" | "symmetric_difference" => ring.clone(),
        _ => Vec::new(),
    };
    let lower = [tool.interval_mm[0], stock.min_mm[2]];
    let middle = [stock.min_mm[2], stock.max_mm[2]];
    let upper = [stock.max_mm[2], tool.interval_mm[1]];
    let caps = matches!(
        operation,
        "union" | "cylinder" | "b_minus_a" | "symmetric_difference"
    );
    let core = matches!(operation, "union" | "cylinder" | "intersection");
    for (enabled, name, interval) in [
        (caps, "lower", lower),
        (core, "middle", middle),
        (caps, "upper", upper),
    ] {
        if enabled {
            blocks.extend(disk_layer(name, interval, stock, tool.radius_mm, &ring)?);
        }
    }
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    if blocks
        .iter()
        .map(|b| b.id.as_str())
        .collect::<BTreeSet<_>>()
        .len()
        != blocks.len()
    {
        return Err("duplicate_boolean_block_identity");
    }
    let (certificate, connected_components) = certify_boolean_blocks(&blocks)?;
    let (volume_constant_mm3, volume_pi_coefficient_mm3) = match operation {
        "box" => (stock_volume, 0),
        "cylinder" => (0, cylinder),
        "union" => (stock_volume, outside),
        "intersection" => (0, overlap),
        "a_minus_b" => (stock_volume, -overlap),
        "b_minus_a" => (0, outside),
        "symmetric_difference" => (stock_volume, outside - overlap),
        _ => unreachable!(),
    };
    Ok(BooleanState {
        operation,
        blocks,
        certificate,
        volume_constant_mm3,
        volume_pi_coefficient_mm3,
        connected_components,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoxOperand {
    pub min_mm: [i32; 3],
    pub max_mm: [i32; 3],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThroughCylinderOperand {
    pub center_mm: [i32; 2],
    pub radius_mm: i32,
    pub interval_mm: [i32; 2],
}

#[derive(Debug, Clone, PartialEq)]
pub struct BoxCylinderBoolean {
    pub stock: BoxOperand,
    pub tool: ThroughCylinderOperand,
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub volume_constant_mm3: i64,
    pub volume_pi_coefficient_mm3: i64,
    pub planner_method: &'static str,
}

/// Construct the regularized difference `stock \ vertical_through_cylinder`.
pub fn subtract_through_cylinder(
    stock: BoxOperand,
    tool: ThroughCylinderOperand,
) -> Result<BoxCylinderBoolean, &'static str> {
    let lo = stock.min_mm;
    let hi = stock.max_mm;
    let span = [
        hi[0].checked_sub(lo[0]).ok_or("stock_extent_overflow")?,
        hi[1].checked_sub(lo[1]).ok_or("stock_extent_overflow")?,
        hi[2].checked_sub(lo[2]).ok_or("stock_extent_overflow")?,
    ];
    if span.iter().any(|value| *value <= 0) {
        return Err("stock_extents_must_be_positive");
    }
    if span[2] < 2 {
        return Err("stock_height_must_allow_an_interior_seed");
    }
    if tool.interval_mm[0] >= lo[2] || tool.interval_mm[1] <= hi[2] {
        return Err("cylinder_must_strictly_cross_both_stock_caps");
    }
    if tool.interval_mm[0] >= tool.interval_mm[1] {
        return Err("cylinder_interval_must_be_positive");
    }
    if tool.radius_mm <= 0
        || i64::from(tool.radius_mm) * 2 >= i64::from(span[0])
        || i64::from(tool.radius_mm) * 2 >= i64::from(span[1])
    {
        return Err("cylinder_must_fit_strictly_inside_stock_section");
    }
    if i64::from(tool.center_mm[0]) * 2 != i64::from(lo[0]) + i64::from(hi[0])
        || i64::from(tool.center_mm[1]) * 2 != i64::from(lo[1]) + i64::from(hi[1])
    {
        return Err("cylinder_axis_must_be_centered_in_stock");
    }

    let outer = SketchLoop {
        id: "stock".into(),
        vertices: vec![
            [lo[0], lo[1]],
            [hi[0], lo[1]],
            [hi[0], hi[1]],
            [lo[0], hi[1]],
        ],
    };
    let sketch = OrthogonalSketch {
        outer,
        height: span[2],
        fillets: Vec::new(),
        bores: vec![Bore {
            id: "cylinder".into(),
            axis: 2,
            center: [tool.center_mm[0], tool.center_mm[1], span[2] / 2],
            radius: tool.radius_mm,
        }],
    };
    let mut result = orthogonal_sketch::build(&sketch)?;
    // The sketch constructor uses z=[0,h]. The public box operand is expressed
    // in world coordinates, so translate its homogeneous controls exactly.
    if lo[2] != 0 {
        for block in &mut result.blocks {
            for point in &mut block.controls {
                let offset = i64::from(lo[2])
                    .checked_mul(point[3])
                    .ok_or("world_z_translation_overflow")?;
                point[2] = point[2]
                    .checked_add(offset)
                    .ok_or("world_z_translation_overflow")?;
            }
        }
    }
    let area = i64::from(span[0])
        .checked_mul(i64::from(span[1]))
        .ok_or("stock_measure_overflow")?;
    let volume_constant_mm3 = area
        .checked_mul(i64::from(span[2]))
        .ok_or("stock_measure_overflow")?;
    let volume_pi_coefficient_mm3 = -i64::from(tool.radius_mm)
        .checked_mul(i64::from(tool.radius_mm))
        .and_then(|r2| r2.checked_mul(i64::from(span[2])))
        .ok_or("cylinder_measure_overflow")?;
    let planner_method = if result.blocks.len() == 4 {
        "four_sector_local_sweeps"
    } else {
        "orthogonal_sketch_auto"
    };
    Ok(BoxCylinderBoolean {
        stock,
        tool,
        blocks: result.blocks,
        certificate: result.certificate,
        volume_constant_mm3,
        volume_pi_coefficient_mm3,
        planner_method,
    })
}
