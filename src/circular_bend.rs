//! Circular bends of certified planar sections, including a periodic fourth quadrant.
//! With section coordinates (ξ,η), the physical point at angle θ is
//! ((R+ξ)cosθ, (R+ξ)sinθ, −η). The two-dimensional orientation and the
//! positive angular direction yield positive three-dimensional orientation.
use crate::{bracket::BracketBlock, sections::section_blocks};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BendSection {
    Rectangle,
    Disk,
    Annulus,
    PiercedSquare,
}

impl BendSection {
    pub fn family(self) -> &'static str {
        match self {
            Self::Rectangle => "reference_rectangle_v1",
            Self::Disk => "core_disk_v1",
            Self::Annulus => "full_annulus_v1",
            Self::PiercedSquare => "pierced_square_v1",
        }
    }
    pub fn parse(name: &str) -> Result<Self, &'static str> {
        [
            Self::Rectangle,
            Self::Disk,
            Self::Annulus,
            Self::PiercedSquare,
        ]
        .into_iter()
        .find(|v| v.family() == name)
        .ok_or("unsupported circular bend section")
    }
    pub fn extent_mm(self) -> i32 {
        match self {
            Self::Rectangle | Self::PiercedSquare => 10,
            Self::Disk | Self::Annulus => 4,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BendRecipe {
    pub section: BendSection,
    pub major_radius_mm: i32,
    pub quadrants: u8,
}
impl BendRecipe {
    pub fn validate(self) -> Result<(), &'static str> {
        if !(1..=4).contains(&self.quadrants) || !(1..=64).contains(&self.major_radius_mm) {
            return Err("circular_bend_parameter_range");
        }
        if self.major_radius_mm <= self.section.extent_mm() {
            return Err("circular_bend_axis_intersection");
        }
        Ok(())
    }

    /// Inverse to the generating section and polar angle, measured in radians.
    /// The fourth-quadrant seam uses angle zero; section membership is separate.
    pub fn section_coordinates(self, point: [f64; 3]) -> Result<[f64; 3], &'static str> {
        self.validate()?;
        let [x, y, z] = point;
        if point.iter().any(|q| !q.is_finite()) || (x == 0. && y == 0.) {
            return Err("invalid circular bend query");
        }
        let theta = y.atan2(x).rem_euclid(std::f64::consts::TAU);
        if theta > f64::from(self.quadrants) * std::f64::consts::FRAC_PI_2 {
            return Err("point outside bend angle");
        }
        Ok([x.hypot(y) - f64::from(self.major_radius_mm), -z, theta])
    }
}

pub fn build_circular_bend(recipe: BendRecipe) -> Result<Vec<BracketBlock>, &'static str> {
    recipe.validate()?;
    let source = section_blocks(recipe.section.family())?;
    let mut result = Vec::with_capacity(source.len() * usize::from(recipe.quadrants));
    for quadrant in 0..recipe.quadrants {
        for block in &source {
            let mut controls = Vec::new();
            for pair in block.controls.as_chunks::<2>().0 {
                let [x, y, _, w] = pair[0].map(i128::from);
                let r = i128::from(recipe.major_radius_mm) * w + x;
                // Homogeneous Bernstein coefficients for (1-t²,2t,1+t²).
                for (c, s, d) in [(1, 0, 1), (1, 1, 1), (0, 2, 2)] {
                    let (mut a, mut b) = (r * c, r * s);
                    for _ in 0..quadrant {
                        (a, b) = (-b, a);
                    }
                    let mut control = [0; 4];
                    for (out, v) in control.iter_mut().zip([a, b, -y * d, w * d]) {
                        *out = i64::try_from(v).map_err(|_| "circular_bend_control_range")?;
                    }
                    controls.push(control);
                }
            }
            result.push(BracketBlock {
                id: format!("bend{quadrant}_{}", block.id),
                degrees: [block.degrees[0], block.degrees[1], 2],
                controls,
            });
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RationalBezierMap;
    fn map(b: &BracketBlock) -> RationalBezierMap {
        RationalBezierMap::new(
            b.degrees,
            b.controls.iter().map(|x| x.map(|v| v as f64)).collect(),
        )
        .unwrap()
    }
    #[test]
    fn physical_points_follow_the_rotated_section_and_seam() {
        for section in [
            BendSection::Rectangle,
            BendSection::Disk,
            BendSection::Annulus,
            BendSection::PiercedSquare,
        ] {
            for r in [section.extent_mm() + 1, 64] {
                let recipe = BendRecipe {
                    section,
                    major_radius_mm: r,
                    quadrants: 4,
                };
                let source = section_blocks(section.family()).unwrap();
                let blocks = build_circular_bend(recipe).unwrap();
                for (index, b) in blocks.iter().enumerate() {
                    let base = map(&source[index % source.len()]);
                    let curved = map(b);
                    let quadrant = index / source.len();
                    for t in [
                        [0., 0., 0.],
                        [1., 1., 1.],
                        [0.3, 0.7, 0.4],
                        [0.125, 0.625, 0.875],
                    ] {
                        let p = curved.evaluate(t).unwrap().position;
                        let uv = base.evaluate(t).unwrap().position;
                        let [xi, eta, angle] = recipe.section_coordinates(p).unwrap();
                        let expected =
                            quadrant as f64 * std::f64::consts::FRAC_PI_2 + 2. * t[2].atan();
                        assert!((xi - uv[0]).abs() < 1e-12 && (eta - uv[1]).abs() < 1e-12);
                        assert!(
                            (angle.sin() - expected.sin()).abs() < 1e-14
                                && (angle.cos() - expected.cos()).abs() < 1e-14
                        );
                    }
                }
            }
        }
    }
}
