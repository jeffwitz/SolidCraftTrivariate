//! Surface loads on whole or partial parameter faces of a certified mesh.
use crate::{
    analysis_space::AnalysisSpace,
    elasticity::{basis_geometry, gauss, Vector},
    incidence::BlockSide,
};
use std::collections::{BTreeMap, BTreeSet};

/// Coordinates use increasing local parameter axes after omitting the normal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LoadPatch {
    pub side: BlockSide,
    pub depth: u8,
    pub intervals: [[u32; 2]; 2],
}

impl LoadPatch {
    pub fn whole(side: BlockSide) -> Self {
        Self {
            side,
            depth: 0,
            intervals: [[0, 1]; 2],
        }
    }

    fn bounds(&self) -> Result<[[f64; 2]; 2], &'static str> {
        if self.depth > 30 {
            return Err("surface patch depth exceeds the dyadic coordinate range");
        }
        let scale = 1_u32 << self.depth;
        if self
            .intervals
            .iter()
            .any(|&[lo, hi]| lo >= hi || hi > scale)
        {
            return Err("empty or invalid surface patch rectangle");
        }
        Ok(self
            .intervals
            .map(|r| r.map(|n| f64::from(n) / f64::from(scale))))
    }
}

pub struct SurfaceLoad {
    pub nodal_force_n: Vec<f64>,
    pub area_mm2: f64,
    pub resultant_n: Vector,
    /// Moment about the physical coordinate origin.
    pub moment_n_mm: Vector,
}

fn cross(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}

/// Integrate a traction (N/mm²), evaluated at the physical position and
/// outward unit normal, into the interleaved full nodal force. The caller
/// supplies the certified exterior sides of these analysis blocks.
pub fn integrate_traction(
    space: &AnalysisSpace,
    exterior: &BTreeSet<BlockSide>,
    patches: &[LoadPatch],
    order: usize,
    traction: impl Fn(Vector, Vector) -> Vector,
) -> Result<SurfaceLoad, &'static str> {
    let quadrature = gauss(order)?;
    if patches.is_empty() {
        return Err("empty traction selection");
    }
    let mut known = BTreeMap::new();
    for b in &space.blocks {
        if known.insert(&b.id, b).is_some() {
            return Err("duplicate analysis block");
        }
    }
    let mut used = BTreeMap::<&BlockSide, Vec<[[f64; 2]; 2]>>::new();
    let len = space
        .ndof
        .checked_mul(3)
        .ok_or("force vector dimension overflow")?;
    let mut result = SurfaceLoad {
        nodal_force_n: vec![0.; len],
        area_mm2: 0.,
        resultant_n: [0.; 3],
        moment_n_mm: [0.; 3],
    };
    for patch in patches {
        let face = &patch.side;
        if face.axis >= 3 || !exterior.contains(face) {
            return Err("traction requires a certified exterior side");
        }
        let b = known.get(&face.block).ok_or("unknown traction block")?;
        if b.dofs.len() != 27 || b.dofs.iter().any(|&d| d >= space.ndof) {
            return Err("invalid traction block DOFs");
        }
        let rect = patch.bounds()?;
        let previous = used.entry(face).or_default();
        if previous
            .iter()
            .any(|r| (0..2).all(|i| rect[i][0] < r[i][1] && r[i][0] < rect[i][1]))
        {
            return Err("overlapping traction patches");
        }
        previous.push(rect);
        let axes = match face.axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let sign = if face.upper ^ (face.axis == 1) {
            1.
        } else {
            -1.
        };
        let width = rect.map(|[lo, hi]| hi - lo);
        for &(u, wu) in &quadrature {
            for &(v, wv) in &quadrature {
                let mut t = [0.; 3];
                t[face.axis] = f64::from(u8::from(face.upper));
                t[axes[0]] = rect[0][0] + width[0] * u;
                t[axes[1]] = rect[1][0] + width[1] * v;
                let sample = basis_geometry(b, t)?;
                let area = cross(
                    sample.jacobian_columns[axes[0]],
                    sample.jacobian_columns[axes[1]],
                );
                let jac = area[0].hypot(area[1]).hypot(area[2]);
                if jac <= 0. || !jac.is_finite() {
                    return Err("invalid surface Jacobian");
                }
                let measure = jac * wu * wv * width[0] * width[1];
                let normal = area.map(|s| sign * s / jac);
                let df = traction(sample.position, normal).map(|s| s * measure);
                if df.iter().any(|x| !x.is_finite()) {
                    return Err("nonfinite integrated traction");
                }
                let moment = cross(sample.position, df);
                result.area_mm2 += measure;
                for d in 0..3 {
                    result.resultant_n[d] += df[d];
                    result.moment_n_mm[d] += moment[d];
                }
                for (&id, value) in b.dofs.iter().zip(sample.values) {
                    for (d, f) in df.iter().enumerate() {
                        result.nodal_force_n[3 * id + d] += value * f;
                    }
                }
            }
        }
    }
    if !result.area_mm2.is_finite()
        || result
            .nodal_force_n
            .iter()
            .chain(&result.resultant_n)
            .chain(&result.moment_n_mm)
            .any(|v| !v.is_finite())
    {
        return Err("nonfinite assembled surface load");
    }
    Ok(result)
}
