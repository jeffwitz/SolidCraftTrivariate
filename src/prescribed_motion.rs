//! Exact essential values of affine physical motions on rational control points.
use crate::{constraint_system::Ratio, json::Json};

#[derive(Clone, Copy)]
pub(crate) struct Motion {
    coefficients: [[Ratio; 4]; 3],
}

impl Motion {
    pub(crate) fn zero() -> Self {
        Self {
            coefficients: [[Ratio::zero(); 4]; 3],
        }
    }

    pub(crate) fn parse(input: &Json) -> Result<Self, &'static str> {
        input.keys(&["denominator", "numerators"])?;
        let den = i128::from(input.field("denominator")?.number::<i64>()?);
        if den <= 0 {
            return Err("nonpositive affine motion denominator");
        }
        let rows = input.field("numerators")?.items()?;
        if rows.len() != 3 {
            return Err("expected three affine motion components");
        }
        let mut motion = Self::zero();
        for (out, row) in motion.coefficients.iter_mut().zip(rows) {
            let terms = row.items()?;
            if terms.len() != 4 {
                return Err("expected affine offset and three slopes");
            }
            for (a, v) in out.iter_mut().zip(terms) {
                *a = Ratio::new(i128::from(v.number::<i64>()?), den)?;
            }
        }
        Ok(motion)
    }

    pub(crate) fn set_component(&mut self, c: usize, from: Self) {
        self.coefficients[c] = from.coefficients[c];
    }

    pub(crate) fn as_f64(self) -> [[f64; 4]; 3] {
        self.coefficients
            .map(|row| row.map(|v| v.numerator() as f64 / v.denominator() as f64))
    }

    /// Material derivative of this spatial support law along an affine velocity:
    /// the constant term drops out and the gradient multiplies all four velocity coefficients.
    pub(crate) fn along(self, velocity: Self) -> Result<Self, &'static str> {
        let mut out = Self::zero();
        for (component, row) in out.coefficients.iter_mut().enumerate() {
            for (j, slot) in row.iter_mut().enumerate() {
                for k in 0..3 {
                    *slot = slot.checked_add(
                        self.coefficients[component][k + 1]
                            .checked_mul(velocity.coefficients[k][j])?,
                    )?;
                }
            }
        }
        Ok(out)
    }

    /// Extract translation and infinitesimal rotation, which have zero linear strain.
    pub(crate) fn rigid_part(self) -> Result<Self, &'static str> {
        let mut r = Self::zero();
        let minus_half = Ratio::new(-1, 2)?;
        let half = Ratio::new(1, 2)?;
        for i in 0..3 {
            r.coefficients[i][0] = self.coefficients[i][0];
            for j in 0..3 {
                r.coefficients[i][1 + j] = self.coefficients[i][1 + j]
                    .checked_mul(half)?
                    .checked_add(self.coefficients[j][1 + i].checked_mul(minus_half)?)?;
            }
        }
        Ok(r)
    }

    pub(crate) fn value(&self, c: usize, h: [i128; 4]) -> Result<Ratio, &'static str> {
        let a = self.coefficients[c];
        let mut sum = a[0];
        for j in 0..3 {
            if a[1 + j] != Ratio::zero() {
                sum = sum.checked_add(a[1 + j].checked_mul(Ratio::new(h[j], h[3])?)?)?;
            }
        }
        Ok(sum)
    }
}
