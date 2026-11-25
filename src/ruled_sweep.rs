//! Ruled extrusion between planar sections related by a complex similarity.
//! The intermediate complex multiplier is affine in normalized height.
use crate::{bracket::BracketBlock, sections::section_blocks};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SweepSection {
    Rectangle,
    PiercedSquare,
}
impl SweepSection {
    pub fn family(self) -> &'static str {
        match self {
            Self::Rectangle => "reference_rectangle_v1",
            Self::PiercedSquare => "pierced_square_v1",
        }
    }
    pub fn parse(name: &str) -> Result<Self, &'static str> {
        match name {
            "reference_rectangle_v1" => Ok(Self::Rectangle),
            "pierced_square_v1" => Ok(Self::PiercedSquare),
            _ => Err("unsupported ruled sweep section"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SweepRecipe {
    pub section: SweepSection,
    /// Real and imaginary components of the final planar multiplier, in eighths.
    pub top_multiplier_eighths: [i32; 2],
    pub top_translation_eighths_mm: [i32; 2],
    pub height_mm: i32,
}
impl SweepRecipe {
    pub fn validate(&self) -> Result<(), &'static str> {
        if !(1..=32).contains(&self.height_mm)
            || self
                .top_multiplier_eighths
                .iter()
                .any(|v| !(-16..=16).contains(v))
            || self
                .top_translation_eighths_mm
                .iter()
                .any(|v| !(-32..=32).contains(v))
        {
            return Err("ruled_sweep_parameter_range");
        }
        let [a, b] = self.top_multiplier_eighths;
        if b == 0 && a <= 0 {
            return Err("ruled_sweep_collapsed_section");
        }
        Ok(())
    }

    /// Physical inverse to the section coordinates and the fractional height.
    /// Membership in the chosen planar section is checked separately.
    pub fn inverse_coordinates(&self, point: [f64; 3]) -> Result<[f64; 3], &'static str> {
        self.validate()?;
        if point.iter().any(|c| !c.is_finite()) {
            return Err("nonfinite sweep point");
        }
        let w = point[2] / f64::from(self.height_mm);
        if !(0.0..=1.).contains(&w) {
            return Err("sweep point outside axial domain");
        }
        let [a, b] = self.top_multiplier_eighths.map(|v| f64::from(v) / 8.);
        let [tx, ty] = self.top_translation_eighths_mm.map(|v| f64::from(v) / 8.);
        let (re, im) = (1. + (a - 1.) * w, b * w);
        let (x, y) = (point[0] - tx * w, point[1] - ty * w);
        let den = re * re + im * im;
        let output = [(re * x + im * y) / den, (re * y - im * x) / den, w];
        if output.iter().all(|v| v.is_finite()) {
            Ok(output)
        } else {
            Err("ruled_sweep_inverse_range")
        }
    }
}

/// Both slices have identical weights. The stored linear axial interpolation
/// therefore represents the intended family throughout each block.
pub fn build_ruled_sweep(recipe: SweepRecipe) -> Result<Vec<BracketBlock>, &'static str> {
    recipe.validate()?;
    let [a, b] = recipe.top_multiplier_eighths.map(i128::from);
    let [tx, ty] = recipe.top_translation_eighths_mm.map(i128::from);
    let h = i128::from(recipe.height_mm);
    let mut blocks = section_blocks(recipe.section.family())?;
    for block in &mut blocks {
        let mut controls = Vec::with_capacity(block.controls.len());
        for pair in block.controls.as_chunks::<2>().0 {
            let [x, y, _, w] = pair[0].map(i128::from);
            controls.extend([
                [8 * x, 8 * y, 0, 8 * w],
                [
                    a * x - b * y + tx * w,
                    b * x + a * y + ty * w,
                    8 * h * w,
                    8 * w,
                ],
            ]);
        }
        // All operands above are bounded i64 coordinates times small recipe
        // integers. One common gcd can be removed without changing any map.
        let common = controls.iter().flatten().copied().fold(0_i128, |mut a, v| {
            let mut b = v.abs();
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a
        });
        for (out, v) in block
            .controls
            .iter_mut()
            .flatten()
            .zip(controls.into_iter().flatten())
        {
            *out = i64::try_from(v / common).map_err(|_| "ruled_sweep_coefficient_range")?;
        }
    }
    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RationalBezierMap;
    #[test]
    fn varying_height_sections_match_their_analytic_similarities() {
        for section in [SweepSection::Rectangle, SweepSection::PiercedSquare] {
            for (a, b, h, tx, ty) in [
                (8, 0, 8, 0, 0),
                (0, 8, 12, 4, -8),
                (-8, 1, 7, -16, 24),
                (-8, -8, 16, 8, 4),
                (16, 0, 1, -32, 32),
                (4, -12, 32, 32, -32),
            ] {
                let r = SweepRecipe {
                    section,
                    height_mm: h,
                    top_multiplier_eighths: [a, b],
                    top_translation_eighths_mm: [tx, ty],
                };
                let source = section_blocks(section.family()).unwrap();
                let blocks = build_ruled_sweep(r).unwrap();
                for (base, volume) in source.into_iter().zip(blocks) {
                    let map = |b: BracketBlock| {
                        RationalBezierMap::new(
                            b.degrees,
                            b.controls.iter().map(|p| p.map(|n| n as f64)).collect(),
                        )
                        .unwrap()
                    };
                    let (src, dst) = (map(base), map(volume));
                    for j in 0..100 {
                        let t = [
                            f64::from(j % 11) / 10.,
                            f64::from(j % 7) / 6.,
                            f64::from(j % 17) / 16.,
                        ];
                        let x = src.evaluate(t).unwrap().position;
                        let y = dst.evaluate(t).unwrap().position;
                        assert!((y[2] - f64::from(h) * t[2]).abs() < 1e-13);
                        // The physical query uses the known exact cap coordinate.
                        // Direct binary64 evaluation can round it outside a closed cap.
                        let z = if t[2] == 1. { f64::from(h) } else { y[2] };
                        let inv = r.inverse_coordinates([y[0], y[1], z]).unwrap();
                        assert!(
                            (x[0] - inv[0]).abs() < 1e-12
                                && (x[1] - inv[1]).abs() < 1e-12
                                && (inv[2] - t[2]).abs() < 1e-14
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn collapsed_and_out_of_range_sweeps_are_refused() {
        let r = SweepRecipe {
            section: SweepSection::Rectangle,
            height_mm: 8,
            top_multiplier_eighths: [0, 0],
            top_translation_eighths_mm: [0; 2],
        };
        for a in -16..=0 {
            assert_eq!(
                build_ruled_sweep(SweepRecipe {
                    top_multiplier_eighths: [a, 0],
                    ..r
                }),
                Err("ruled_sweep_collapsed_section")
            );
        }
        assert_eq!(
            build_ruled_sweep(SweepRecipe { height_mm: 33, ..r }),
            Err("ruled_sweep_parameter_range")
        );
    }
}
