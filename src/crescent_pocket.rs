//! Layered bore/pocket intersection model built from the certified G7 ring.
//!
//! All geometry is homogeneous integer data. The frozen rational value of
//! sqrt(3) and the complete G7 sector nets are reused through `overlap_blocks`.

use crate::{bracket::BracketBlock, overlap_bore::overlap_blocks};
use std::convert::TryFrom;

const MAX_DEPTH_TICKS: i32 = 16;

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("crescent coefficient overflow")
}
fn sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("crescent coefficient overflow")
}
fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("crescent coefficient overflow")
}
fn out(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "crescent control exceeds i64")
}
#[cfg(test)]
fn index(i: usize, j: usize, k: usize) -> usize {
    (i * 3 + j) * 2 + k
}

#[cfg(test)]
fn certify_patch_planarity(controls: &[[i64; 4]]) -> Result<(), &'static str> {
    if controls.len() != 18 {
        return Err("invalid crescent patch control count");
    }
    let (pairs, remainder) = controls.as_chunks::<2>();
    if !remainder.is_empty() {
        return Err("invalid crescent patch control count");
    }
    let planar = pairs
        .iter()
        .map(|pair| [pair[0][0], pair[0][1], pair[0][3]])
        .collect::<Vec<_>>();
    crate::planar::certify_planar([2, 2], &planar)
        .map(|_| ())
        .map_err(|error| match error {
            crate::planar::PlanarError::Invalid => "crescent patch has nonpositive Jacobian",
            crate::planar::PlanarError::Indeterminate => "crescent patch regularity indeterminate",
        })
}

fn source_arc(blocks: &[BracketBlock], id: &str) -> Result<[[i128; 4]; 3], &'static str> {
    let block = blocks
        .iter()
        .find(|b| b.id == id)
        .ok_or("missing G7 source sector")?;
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("invalid G7 source sector");
    }
    let mut arc = [[0_i128; 4]; 3];
    for (i, control) in arc.iter_mut().enumerate() {
        let p = block.controls[(i * 2 + 1) * 2];
        if p[3] <= 0 {
            return Err("nonpositive G7 source weight");
        }
        *control = p.map(i128::from);
    }
    Ok(arc)
}

fn reversed(mut arc: [[i128; 4]; 3]) -> [[i128; 4]; 3] {
    arc.reverse();
    arc
}

fn right_circle_to_bc(mut arc: [[i128; 4]; 3]) -> Result<[[i128; 4]; 3], &'static str> {
    // The source arc runs from (3,0) to (2,s). Reflect about x=1 to get B->C.
    for p in &mut arc {
        p[0] = sub(mul(2, p[3])?, p[0])?;
        for component in p.iter_mut() {
            *component = mul(3, *component)?;
        }
    }
    Ok(arc)
}

type Curve = [[i128; 4]; 3];
type Vertex = [i128; 4];

fn scale_curve(curve: Curve, scale: i128) -> Result<Curve, &'static str> {
    let mut result = curve;
    for p in &mut result {
        for component in p.iter_mut() {
            *component = mul(*component, scale)?;
        }
    }
    Ok(result)
}

fn split_half(curve: Curve, scale: i128) -> Result<(Curve, Curve), &'static str> {
    // The de Casteljau midpoint has denominator four for quadratic Bernstein data.
    // `scale` must be divisible by four; both child nets use one common multiplier.
    if scale % 4 != 0 {
        return Err("split scale must clear the quadratic midpoint denominator");
    }
    let factor = scale / 4;
    let mut left = [[0_i128; 4]; 3];
    let mut right = [[0_i128; 4]; 3];
    for c in 0..4 {
        left[0][c] = mul(scale, curve[0][c])?;
        left[1][c] = mul(mul(2, factor)?, add(curve[0][c], curve[1][c])?)?;
        let midpoint = add(add(curve[0][c], mul(2, curve[1][c])?)?, curve[2][c])?;
        left[2][c] = mul(factor, midpoint)?;
        right[0][c] = left[2][c];
        right[1][c] = mul(mul(2, factor)?, add(curve[1][c], curve[2][c])?)?;
        right[2][c] = mul(scale, curve[2][c])?;
    }
    Ok((left, right))
}

fn line_curve(start: Vertex, end: Vertex) -> Result<Curve, &'static str> {
    let mut result = [[0_i128; 4]; 3];
    result[0] = start;
    result[2] = end;
    for c in 0..4 {
        let sum = add(start[c], end[c])?;
        if sum % 2 != 0 {
            return Err("line elevation is not integral at common scale");
        }
        result[1][c] = sum / 2;
    }
    Ok(result)
}

fn patch_controls(
    bottom: Curve,
    right: Curve,
    top: Curve,
    left: Curve,
) -> Result<Vec<[i64; 4]>, &'static str> {
    if bottom[0] != left[0] || bottom[2] != right[0] || right[2] != top[2] || top[0] != left[2] {
        return Err("crescent patch corner traces disagree");
    }
    let mut net = [[[0_i128; 4]; 3]; 3];
    for i in 0..3 {
        net[i][0] = bottom[i];
        net[i][2] = top[i];
    }
    net[0].copy_from_slice(&left);
    net[2].copy_from_slice(&right);
    // Homogeneous Coons center for a biquadratic patch with four quadratic edges.
    let mut center = [0_i128; 4];
    for (c, component) in center.iter_mut().enumerate() {
        let edge_sum = add(
            add(net[1][0][c], net[1][2][c])?,
            add(net[0][1][c], net[2][1][c])?,
        )?;
        let corner_sum = add(
            add(net[0][0][c], net[0][2][c])?,
            add(net[2][0][c], net[2][2][c])?,
        )?;
        *component = sub(mul(2, edge_sum)?, corner_sum)?;
    }
    net[1][1] = center;
    if net.iter().flatten().any(|p| p[3] <= 0) {
        return Err("nonpositive crescent patch weight");
    }
    let mut controls = Vec::with_capacity(18);
    for row in net {
        for point in row {
            for z in [0_i128, 8] {
                let mut p = point;
                p[2] = mul(z, p[3])?;
                controls.push(
                    p.map(out)
                        .into_iter()
                        .collect::<Result<Vec<_>, _>>()?
                        .try_into()
                        .map_err(|_| "invalid crescent point")?,
                );
            }
        }
    }
    Ok(controls)
}

fn upper_crescent_patches(g7: &[BracketBlock]) -> Result<[BracketBlock; 3], &'static str> {
    // Pull the two exposed left-circle arcs and transform the right-circle arc.
    // The boundary curves receive one common scale per source edge; all seam
    // traces below are degree-elevated from the same homogeneous vertices.
    let arc_ad = scale_curve(
        reversed(source_arc(g7, "overlap_bore_left_sector_01_z0_8")?),
        8,
    )?; // A->D
    let arc_dc = scale_curve(
        reversed(source_arc(g7, "overlap_bore_left_sector_00_z0_8")?),
        8,
    )?; // D->C
    let arc_bc = right_circle_to_bc(source_arc(g7, "overlap_bore_right_sector_02_z0_8")?)?; // B->C
    let (arc_bf, arc_fc) = split_half(arc_bc, 8)?;
    let [a, d] = [arc_ad[0], arc_ad[2]];
    let [d2, c] = [arc_dc[0], arc_dc[2]];
    let [b, f] = [arc_bf[0], arc_bf[2]];
    let [f2, c2] = [arc_fc[0], arc_fc[2]];
    if d != d2 || f != f2 || c != c2 {
        return Err("crescent arc endpoints disagree");
    }
    let base = a[3] / 24;
    if base <= 0 || mul(24, base)? != a[3] {
        return Err("invalid G7 scale for crescent vertices");
    }
    if b[3] != mul(72, base)? || c[3] != mul(96, base)? || f[3] != mul(78, base)? {
        return Err("crescent source arc weights changed");
    }
    let e = [mul(-2, mul(48, base)?)?, 0, 0, mul(48, base)?];
    let o_weight = mul(64, base)?;
    let ox = mul(-7, o_weight)?;
    let oy = mul(7, o_weight)?;
    if ox % 4 != 0 || oy % 8 != 0 {
        return Err("interior crescent vertex is not integral");
    }
    let o = [ox / 4, oy / 8, 0, o_weight];
    let line = |p, q| line_curve(p, q);
    let upper = [
        BracketBlock {
            id: "g10_crescent_upper_a".into(),
            degrees: [2, 2, 1],
            controls: patch_controls(line(a, e)?, line(e, o)?, line(d, o)?, arc_ad)?,
        },
        BracketBlock {
            id: "g10_crescent_upper_b".into(),
            degrees: [2, 2, 1],
            controls: patch_controls(line(e, b)?, arc_bf, line(o, f)?, line(e, o)?)?,
        },
        BracketBlock {
            id: "g10_crescent_upper_c".into(),
            degrees: [2, 2, 1],
            controls: patch_controls(line(o, f)?, arc_fc, arc_dc, line(o, d)?)?,
        },
    ];
    for patch in &upper {
        if patch.controls.iter().any(|p| p[3] <= 0) {
            return Err("nonpositive crescent patch weight");
        }
    }
    Ok(upper)
}

fn lower_crescent_patches(upper: &[BracketBlock; 3]) -> Result<[BracketBlock; 3], &'static str> {
    let mut lower: [BracketBlock; 3] = std::array::from_fn(|i| BracketBlock {
        id: format!("g10_crescent_lower_{}", ["a", "b", "c"][i]),
        degrees: [2, 2, 1],
        controls: vec![[0; 4]; 18],
    });
    // Reflection y -> -y followed by reversing u preserves positive orientation.
    for (source, target) in upper.iter().zip(lower.iter_mut()) {
        let (source_rows, source_tail) = source.controls.as_chunks::<6>();
        let (target_rows, target_tail) = target.controls.as_chunks_mut::<6>();
        if !source_tail.is_empty() || !target_tail.is_empty() {
            return Err("invalid crescent patch controls");
        }
        for (source_row, target_row) in source_rows.iter().rev().zip(target_rows) {
            for (src, dst) in source_row.iter().zip(target_row.iter_mut()) {
                *dst = *src;
                dst[1] = dst[1]
                    .checked_neg()
                    .ok_or("crescent coefficient overflow")?;
            }
        }
    }
    Ok(lower)
}

fn extrude_layer(
    block: &BracketBlock,
    bottom_tick: i32,
    top_tick: i32,
    id: String,
) -> Result<BracketBlock, &'static str> {
    if bottom_tick >= top_tick {
        return Err("zero height crescent layer");
    }
    let mut result = block.clone();
    result.id = id;
    for (control_index, p) in result.controls.iter_mut().enumerate() {
        let w = i128::from(p[3]);
        let upper = control_index % 2 == 1;
        let z_tick = if upper { top_tick } else { bottom_tick };
        let n = mul(i128::from(z_tick), w)?;
        if n % 2 != 0 {
            return Err("half-mm layer is not exactly representable");
        }
        p[2] = out(n / 2)?;
    }
    Ok(result)
}

fn layer_id(base: &str, suffix: &str) -> String {
    format!("g10_ring_{base}_{suffix}")
}

fn source_sector_id(id: &str) -> Result<(&'static str, &str), &'static str> {
    let side = if id.contains("_right_sector_") {
        "right"
    } else if id.contains("_left_sector_") {
        "left"
    } else {
        return Err("invalid G7 sector id");
    };
    let (_, tail) = id.split_once("_sector_").ok_or("invalid G7 sector id")?;
    let (sector, _) = tail.split_once("_z0_8").ok_or("invalid G7 sector id")?;
    Ok((side, sector))
}

/// Build the plate with a through bore and an intersecting blind pocket.
/// `depth_half_mm` is a half-millimetre integer in the closed range 0..=16.
/// The returned collection contains only positive-height layers.
pub fn crescent_pocket_blocks(depth_half_mm: i32) -> Result<Vec<BracketBlock>, &'static str> {
    if !(0..=MAX_DEPTH_TICKS).contains(&depth_half_mm) {
        return Err("invalid crescent pocket depth");
    }
    let g7 = overlap_blocks()?;
    if g7.len() != 8 {
        return Err("invalid G7 sector count");
    }
    let split_tick = MAX_DEPTH_TICKS - depth_half_mm;
    let crescent_patches = if split_tick > 0 {
        let upper = upper_crescent_patches(&g7)?;
        let lower = lower_crescent_patches(&upper)?;
        Some((upper, lower))
    } else {
        None
    };
    let mut result = Vec::with_capacity(22);
    if split_tick > 0 {
        for source in &g7 {
            let mut renamed = source.clone();
            let (side, sector) = source_sector_id(&source.id)?;
            renamed.id = layer_id(&format!("lower_{side}_sector_{sector}"), "z0_split");
            result.push(extrude_layer(&renamed, 0, split_tick, renamed.id.clone())?);
        }
        let (upper, lower) = crescent_patches.unwrap();
        for patch in upper.into_iter().chain(lower) {
            result.push(extrude_layer(&patch, 0, split_tick, patch.id.clone())?);
        }
    }
    if split_tick < MAX_DEPTH_TICKS {
        for source in &g7 {
            let mut renamed = source.clone();
            let (side, sector) = source_sector_id(&source.id)?;
            renamed.id = layer_id(&format!("upper_{side}_sector_{sector}"), "split_8");
            result.push(extrude_layer(
                &renamed,
                split_tick,
                MAX_DEPTH_TICKS,
                renamed.id.clone(),
            )?);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g10_depth_regimes_omit_zero_height_layers() {
        let zero = crescent_pocket_blocks(0).unwrap();
        let full = crescent_pocket_blocks(16).unwrap();
        let blind = crescent_pocket_blocks(10).unwrap();
        assert_eq!(zero.len(), 14);
        assert_eq!(full.len(), 8);
        assert_eq!(blind.len(), 22);
        assert!(zero.iter().all(|b| !b.id.contains("upper_right")));
        assert!(full.iter().all(|b| !b.id.contains("lower")));
        assert!(blind
            .iter()
            .all(|b| b.controls.len() == b.degrees.iter().map(|d| d + 1).product::<usize>()));
    }

    #[test]
    fn six_crescent_patches_have_positive_weights_and_exact_shared_seams() {
        let g7 = overlap_blocks().unwrap();
        let patches = upper_crescent_patches(&g7).unwrap();
        assert_eq!(
            patches.iter().map(|b| b.id.as_str()).collect::<Vec<_>>(),
            vec![
                "g10_crescent_upper_a",
                "g10_crescent_upper_b",
                "g10_crescent_upper_c"
            ]
        );
        assert!(patches.iter().all(|p| p.degrees == [2, 2, 1]
            && p.controls.len() == 18
            && p.controls.iter().all(|h| h[3] > 0)));
        for patch in &patches {
            certify_patch_planarity(&patch.controls).unwrap();
        }
        // The baseline split uses the corner ratios A:B:C:D=1:3:4:4/3.
        let (a, b, c) = (&patches[0], &patches[1], &patches[2]);
        let wa = a.controls[index(0, 0, 0)][3];
        let wb = b.controls[index(2, 0, 0)][3];
        let wc = c.controls[index(2, 2, 0)][3];
        let wd = a.controls[index(0, 2, 0)][3];
        assert_eq!(i128::from(wb), i128::from(wa) * 3);
        assert_eq!(i128::from(wc), i128::from(wa) * 4);
        assert_eq!(i128::from(wd) * 3, i128::from(wa) * 4);
        // A-E-O and O-F seams are exactly the same elevated homogeneous lines.
        for q in 0..3 {
            assert_eq!(a.controls[index(2, q, 0)], b.controls[index(0, q, 0)]);
            assert_eq!(b.controls[index(q, 2, 0)], c.controls[index(q, 0, 0)]);
        }
        let lower = lower_crescent_patches(&patches).unwrap();
        for (up, down) in patches.iter().zip(&lower) {
            assert_eq!(down.degrees, up.degrees);
            certify_patch_planarity(&down.controls).unwrap();
            let (lower_rows, _) = down.controls.as_chunks::<6>();
            let (upper_rows, _) = up.controls.as_chunks::<6>();
            for (lower_row, upper_row) in lower_rows.iter().zip(upper_rows.iter().rev()) {
                let (lower_edges, _) = lower_row.as_chunks::<2>();
                let (upper_edges, _) = upper_row.as_chunks::<2>();
                for (u, v) in lower_edges.iter().zip(upper_edges.iter()) {
                    assert_eq!(u[0][0], v[0][0]);
                    assert_eq!(u[0][1], -v[0][1]);
                    assert_eq!(u[0][3], v[0][3]);
                }
            }
        }
    }
}
