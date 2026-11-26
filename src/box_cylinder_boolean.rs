//! Exact fixed-chart subtraction of a centered vertical cylinder from a box.
//!
//! The public operands are distinct. A validated operand pair is converted to
//! an orthogonal sketch only inside this constructor so the familiar four
//! circular sectors and their exact wall seams are reused.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    orthogonal_sketch::{self, Bore, OrthogonalSketch},
    sketch_layers::SketchLoop,
};

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
